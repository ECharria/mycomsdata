#!/usr/bin/env python3
"""
convert_msp_to_massbank.py
--------------------------
Convert an MSP‑like library file into MassBank3‑compatible text and insert a
PK$SPLASH line for every record.


"""

import argparse
from pathlib import Path
from typing import List, Tuple

# splash import comes from the spectra‑hash package
from splash import Splash
from splash import Spectrum
from splash import SpectrumType
from typing import Callable, Any  # new import at the top

from collections import Counter
import re

# Monoisotopic masses (Dalton)
MONOISOTOPIC = dict(H=1.007825032, C=12.0, N=14.003074,  O=15.99491462,
                    P=30.973762,   S=31.972071,  Cl=34.968853, Br=78.918338,
                    F=18.9984032,  I=126.90447)

formula_re = re.compile(r'InChI=1S?/([^/]+)')

def formula_and_exactmass_from_inchi(inchi: str):
    m = formula_re.match(inchi)
    if not m:
        return None, None          # couldn’t parse — leave placeholders
    formula = m.group(1)           # e.g.  C22H22O6
    # split elements and counts
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    counts = Counter({el: int(n) if n else 1 for el, n in tokens})

    exact_mass = sum(MONOISOTOPIC[el] * n for el, n in counts.items())
    return formula, f"{exact_mass:.4f}"

def parse_custom_msp(lines: List[str]) -> List[dict]:
    """Parse the custom MSP‑like input format into a list of dictionaries."""
    records, current, peaks = [], {}, []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        def take(key: str, transform: Callable[[str], Any] = str) -> None:
            current[key] = transform(line.split(":", 1)[1].strip())

        # ─── field handlers ────────────────────────────────────────────────
        if line.startswith("CHARGE:"):
            # ### FLUSH the previous record (if complete) ###
            if current and peaks:
                current["PEAKS"] = peaks
                records.append(current)
            elif current and len(peaks) == 0 :
                print(f"Record starting with {current.get('COMPOUND_NAME',f'index {len(records) +1}')  } has no peaks.")
        
            # start a fresh one
            current, peaks = {}, []
            take("CHARGE", int)
            

        elif line.startswith("IONMODE:"):
            take("IONMODE", str.upper)
        elif line.startswith("SMILES:"):
            take("SMILES")
        elif line.startswith("INCHI:"):
            take("INCHI")
        elif line.startswith("PUBMED:"):
            take("PUBMED")
        elif line.startswith("CCS:"):
            take("CCS", float)
        elif line.startswith("COL_ENERGY1:"):
            take("COL_ENERGY1", float)
        elif line.startswith("COL_ENERGY2:"):
            take("COL_ENERGY2", float)
        elif line.startswith("MS_LEVEL:"):
            take("MS_LEVEL", int)
        elif line.startswith("INSTRUMENT_TYPE:"):
            take("INSTRUMENT_TYPE")
        elif line.startswith("COMPOUND_NAME:"):
            take("COMPOUND_NAME")
        elif line.startswith("PRECURSOR_MZ:"):
            take("PRECURSOR_MZ", float)
        elif line.startswith("INCHIKEY:"):
            take("INCHIKEY")
        elif line.startswith("NUM PEAKS:"):
            take("NUM_PEAKS", int)
        else:
            # assumed to be a peak row: "<mz> <intensity>"
            try:
                mz, intensity = map(float, line.split())
                peaks.append((mz, intensity))
            except ValueError:
                pass
        
    if current and peaks:
        current["PEAKS"] = peaks
        records.append(current)
    
    return records

def splash_for_peaks(peaks: List[Tuple[float, float]]) -> str:
    """Generate a SPLASH identifier from a list of (mz, intensity)."""
    # Spectra‑hash expects intensities scaled to arbitrary units – raw is fine.
    if len(peaks) == 0 :
        raise ValueError(f"Record starting with has no peaks.")
        print(f"Record starting with {current.get('COMPOUND_NAME','<unknown>')} has no peaks.")
    spectrum = Spectrum(peaks, SpectrumType.MS)
    return Splash().splash(spectrum)


def record_to_massbank(record: dict, accession_prefix: str, index: int) -> Tuple[str, str]:
    """Convert one parsed record into a MassBank3‑compatible text block."""
    accession = f"MSBNK-{accession_prefix}{index:06d}"
    title = (
        f"{record.get('COMPOUND_NAME')}; "
        f"{record.get('INSTRUMENT_TYPE')}; "
        f"MS{record.get('MS_LEVEL')}; "
        f"CE:{record.get('COL_ENERGY1')}; "
        f"[M+H]+"
    )
        
    formula, exact = formula_and_exactmass_from_inchi(record.get("INCHI", ""))

    mb = [
        f"ACCESSION: {accession}",
        f"RECORD_TITLE: {title}",
        f"DATE: 2025.08.05",
        "AUTHORS: Unknown",
        "LICENSE: CC BY",
        "COPYRIGHT: Unknown",
        "COMMENT: Converted from MSP",
        f"CH$NAME: {record.get('COMPOUND_NAME')}",
        f"CH$FORMULA: {formula or 'CH4'}",
        f"CH$EXACT_MASS: {exact or '16'}",
        f"CH$SMILES: {record.get('SMILES')or 'C'}",
        f"CH$IUPAC: {record.get('INCHI')or'InChI=1S/CH4/h1H4'}",
        f"CH$LINK: INCHIKEY {record.get('INCHIKEY')or 'VNWKTOKETHGBQD-UHFFFAOYSA-N'}",
        f"AC$INSTRUMENT: {record.get('INSTRUMENT_TYPE', 'Unknown')}",
        f"AC$INSTRUMENT_TYPE: {record.get('INSTRUMENT_TYPE')}",
        f"AC$MASS_SPECTROMETRY: MS_TYPE MS{record.get('MS_LEVEL')}",
        f"AC$MASS_SPECTROMETRY: ION_MODE {record.get('IONMODE')}",
        f"AC$MASS_SPECTROMETRY: COLLISION_ENERGY {record.get('COL_ENERGY1')}",
        f"MS$FOCUSED_ION: PRECURSOR_M/Z {record.get('PRECURSOR_MZ')}",
        "MS$DATA_PROCESSING: Converted with convert_msp_to_massbank.py",
        # --- SPLASH ---------------------------------------------------------
        f"PK$SPLASH: {splash_for_peaks(record['PEAKS'])}",
        # --------------------------------------------------------------------
        f"PK$NUM_PEAK: {record.get('NUM_PEAKS')}",
        "PK$PEAK: m/z int. rel.int.",
    ]

    max_intensity = max(i for _, i in record["PEAKS"]) or 1.0
    for mz, intensity in record["PEAKS"]:
        rel_int = int(round(intensity / max_intensity * 999))
        mb.append(f"  {mz:.5f} {intensity:.6g} {rel_int}")

    mb.append("//")
    mb_text= "\n".join(mb)
    return accession, mb_text

# ---- main() rewritten ------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Convert MSP to individual MassBank files")
    ap.add_argument("--input",  default="input.msp", help="Input MSP file")
    ap.add_argument("--outdir", default="massbank_records", help="Directory to write *.txt files")
    ap.add_argument("--accession-prefix", default="HZI-CBIO-AA-", help="ACCESSION prefix")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse the MSP once
    with open(args.input) as f:
        records = parse_custom_msp(f.readlines())

    # Convert and save each record
    for idx, rec in enumerate(records, start=1):
        accession, mb_text = record_to_massbank(rec, args.accession_prefix, idx)
        (outdir / f"{accession}.txt").write_text(mb_text)

    print(f"Wrote {len(records)} MassBank files to {outdir.resolve()}")



if __name__ == "__main__":
    main()