#!/usr/bin/env python3
"""
clean_spectra.py
----------------
Apply matchms spectral cleaning to MassBank3 .txt files before DB import.

Cleaning steps applied (in order):
  1. normalize_intensities      — scale peaks to 0–1
  2. select_by_mz               — keep peaks in [mz_min, 2000]
  3. select_by_relative_intensity — keep peaks >= 0.5% of base peak
  4. remove_peaks_around_precursor_mz — remove ±17 Da window (precursor noise)
  5. require_minimum_number_of_peaks  — drop spectra with < 3 peaks after cleaning

Intensities are rescaled to MassBank3 convention (rel.int. 0–999) in the output.

Usage:
    conda run -n matchms python3 clean_spectra.py \
        --indir  /tmp/library_out_final \
        --outdir /lustre/.../mycomsdata/CBIO \
        [--min-peaks 3] [--min-rel-intensity 0.005] [--mz-min 50]
"""

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
from matchms import Spectrum
from matchms import filtering as msfilters


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir',  required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--min-peaks',         type=int,   default=1)
    ap.add_argument('--min-rel-intensity', type=float, default=0.005)
    ap.add_argument('--mz-min',            type=float, default=50.0)
    ap.add_argument('--remove-precursor-window', type=float, default=0.0,
                    help='Remove peaks within ±N Da of precursor m/z (0 = disabled; timsTOF does not need this)')
    return ap.parse_args()


def parse_massbank_peaks(text: str):
    """Parse m/z and intensity columns from a MassBank3 .txt string."""
    mzs, ints = [], []
    in_peaks = False
    for line in text.splitlines():
        if line.startswith('PK$PEAK:'):
            in_peaks = True
            continue
        if in_peaks:
            if line.startswith('//') or line.startswith('PK$'):
                break
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mzs.append(float(parts[0]))
                    ints.append(float(parts[1]))
                except ValueError:
                    pass
    return np.array(mzs), np.array(ints)


def parse_precursor_mz(text: str) -> Optional[float]:
    for line in text.splitlines():
        if line.startswith('MS$FOCUSED_ION: PRECURSOR_M/Z'):
            try:
                return float(line.split('PRECURSOR_M/Z')[1].strip())
            except ValueError:
                pass
    return None


def clean_spectrum(mzs: np.ndarray, ints: np.ndarray,
                   precursor_mz: Optional[float], args) -> Optional[tuple]:
    """Apply cleaning filters; return (mzs, ints) or None if spectrum is rejected."""
    if len(mzs) == 0:
        return None

    # Build a matchms Spectrum object (metadata optional — we just need peak filtering)
    metadata = {}
    if precursor_mz is not None:
        metadata['precursor_mz'] = precursor_mz

    spec = Spectrum(mz=mzs, intensities=ints, metadata=metadata)

    spec = msfilters.normalize_intensities(spec)
    if spec is None: return None

    spec = msfilters.select_by_mz(spec, mz_from=args.mz_min, mz_to=2000.0)
    if spec is None: return None

    spec = msfilters.select_by_relative_intensity(
        spec, intensity_from=args.min_rel_intensity, intensity_to=1.0)
    if spec is None: return None

    if args.remove_precursor_window > 0 and precursor_mz is not None:
        spec = msfilters.remove_peaks_around_precursor_mz(
            spec, mz_tolerance=args.remove_precursor_window)
    if spec is None: return None

    spec = msfilters.require_minimum_number_of_peaks(spec, n_required=args.min_peaks)
    if spec is None: return None

    return spec.peaks.mz, spec.peaks.intensities


def rewrite_massbank_file(inpath: Path, outpath: Path,
                          clean_mzs: np.ndarray, clean_ints: np.ndarray):
    """Preserve all metadata from original file; replace only the peak block."""
    lines = inpath.read_text().splitlines()

    # Renormalise intensities to 0–1 after filtering, then to 0–999
    max_i = clean_ints.max() if len(clean_ints) else 1.0
    ints_norm = clean_ints / max_i
    peak_lines = [
        f"  {mz:.5f} {inten:.6g} {int(round(inten * 999))}"
        for mz, inten in zip(clean_mzs, ints_norm)
    ]

    new_lines = []
    in_peaks  = False
    for line in lines:
        if line.startswith('PK$NUM_PEAK:'):
            new_lines.append(f'PK$NUM_PEAK: {len(clean_mzs)}')
            continue
        if line.startswith('PK$PEAK:'):
            new_lines.append('PK$PEAK: m/z int. rel.int.')
            new_lines.extend(peak_lines)
            in_peaks = True
            continue
        if in_peaks:
            if line.startswith('//'):
                new_lines.append('//')
                in_peaks = False
            continue
        new_lines.append(line)

    outpath.write_text('\n'.join(new_lines))


def main():
    args   = parse_args()
    indir  = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = sorted(indir.glob('*.txt'))
    print(f"Input files:  {len(files)}")
    print(f"Settings:")
    print(f"  mz_min={args.mz_min}, min_rel_intensity={args.min_rel_intensity}")
    print(f"  precursor_window=±{args.remove_precursor_window} Da")
    print(f"  min_peaks={args.min_peaks}")
    print()

    kept = dropped = 0
    peak_before = peak_after = 0
    dropped_names = []

    for fpath in files:
        text = fpath.read_text()
        mzs, ints = parse_massbank_peaks(text)
        precursor_mz = parse_precursor_mz(text)

        n_before = len(mzs)
        peak_before += n_before

        result = clean_spectrum(mzs, ints, precursor_mz, args)

        if result is None:
            name = next((l.split(':',1)[1].strip() for l in text.splitlines()
                         if l.startswith('CH$NAME:')), fpath.name)
            dropped_names.append(name)
            dropped += 1
            continue

        clean_mzs, clean_ints = result
        peak_after += len(clean_mzs)

        rewrite_massbank_file(fpath, outdir / fpath.name, clean_mzs, clean_ints)
        kept += 1

    pct = 100 * peak_after // peak_before if peak_before else 0
    print(f"Results:")
    print(f"  Kept:    {kept}")
    print(f"  Dropped: {dropped}  (< {args.min_peaks} peaks after cleaning)")
    print(f"  Peaks:   {peak_before} → {peak_after} ({pct}% retained)")

    if dropped_names:
        print(f"\nDropped spectra:")
        for n in dropped_names:
            print(f"  {n}")

    print(f"\nOutput written to: {outdir}")


if __name__ == '__main__':
    main()
