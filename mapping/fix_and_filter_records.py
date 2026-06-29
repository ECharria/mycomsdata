#!/usr/bin/env python3
"""
fix_and_filter_records.py
--------------------------
Post-processing pass on MassBank3 .txt files before DB import.

Fixes applied:
  1. Deduplication   — same SPLASH → keep lowest accession, drop extras
  2. Adduct fix      — derive correct PRECURSOR_TYPE from precursor m/z vs exact mass
  3. RECORD_TITLE    — regenerate to match corrected adduct
  4. Unknowns flagged in a report (diff > tolerance and no standard adduct found)

Usage:
    python3 fix_and_filter_records.py \
        --indir  /tmp/library_cleaned_v3 \
        --outdir /tmp/library_fixed \
        --filter species,publication        # comma-sep fields that must be present
"""

import argparse
import re
import shutil
from collections import defaultdict
from pathlib import Path

# Adduct m/z offsets from neutral mass (positive mode only for now)
ADDUCTS = {
    '[M+H]+':       1.007276,
    '[M+Na]+':      22.989218,
    '[M+K]+':       38.963158,
    '[M+NH4]+':     18.034164,
    '[M+H-H2O]+':   1.007276 - 18.010565,
    '[M+H-2H2O]+':  1.007276 - 36.021130,
    '[M+2H]2+':     (2 * 1.007276) / 2,
    '[M+2Na-H]+':   2 * 22.989218 - 1.007276,
    '[M-H]-':      -1.007276,
    '[M+Cl]-':      34.969402,
    '[M+HCOO]-':    44.998201,
}
# Dimer adducts — checked separately (prec vs 2*mass + offset)
DIMER_ADDUCTS = {
    '[2M+H]+':  1.007276,
    '[2M+Na]+': 22.989218,
}
TOL = 0.02  # Da


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir',  required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--filter', default='SP$SCIENTIFIC_NAME,PUBLICATION',
                    help='Comma-separated field prefixes; records missing any are excluded')
    ap.add_argument('--dry-run', action='store_true',
                    help='Report what would be done without writing output')
    return ap.parse_args()


def read_fields(text: str) -> dict:
    """Extract relevant scalar fields from a MassBank3 text."""
    d = {}
    for line in text.splitlines():
        for key, attr in [
            ('ACCESSION:', 'accession'),
            ('RECORD_TITLE:', 'title'),
            ('CH$NAME:', 'name'),
            ('CH$EXACT_MASS:', 'exact_mass'),
            ('MS$FOCUSED_ION: PRECURSOR_M/Z', 'prec_mz'),
            ('MS$FOCUSED_ION: PRECURSOR_TYPE', 'adduct'),
            ('AC$MASS_SPECTROMETRY: ION_MODE', 'ion_mode'),
            ('PK$SPLASH:', 'splash'),
        ]:
            if line.startswith(key) and attr not in d:
                d[attr] = line[len(key):].strip()
    return d


def derive_adduct(prec_mz: float, exact_mass: float):
    """Return (best_adduct, diff_da) for the closest matching adduct."""
    best_adduct = None
    best_diff = 999.0

    for adduct, offset in ADDUCTS.items():
        diff = abs(prec_mz - (exact_mass + offset))
        if diff < best_diff:
            best_diff = diff
            best_adduct = adduct

    # Check dimers
    for adduct, offset in DIMER_ADDUCTS.items():
        diff = abs(prec_mz - (2 * exact_mass + offset))
        if diff < best_diff:
            best_diff = diff
            best_adduct = adduct

    return best_adduct, best_diff


def fix_record(text: str, correct_adduct: str, original_adduct: str) -> str:
    """Replace PRECURSOR_TYPE and regenerate RECORD_TITLE with new adduct."""
    lines = text.splitlines()
    new_lines = []
    for line in lines:
        # Fix PRECURSOR_TYPE
        if line.startswith('MS$FOCUSED_ION: PRECURSOR_TYPE'):
            new_lines.append(f'MS$FOCUSED_ION: PRECURSOR_TYPE {correct_adduct}')
            continue
        # Fix RECORD_TITLE: replace adduct token at end
        if line.startswith('RECORD_TITLE:') and original_adduct in line:
            new_lines.append(line.replace(original_adduct, correct_adduct))
            continue
        new_lines.append(line)
    return '\n'.join(new_lines)


def main():
    args = parse_args()
    indir  = Path(args.indir)
    outdir = Path(args.outdir)
    if not args.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)

    filter_fields = [f.strip() for f in args.filter.split(',') if f.strip()]

    files = sorted(indir.glob('*.txt'))
    print(f"Input files:  {len(files)}")
    print(f"Required fields: {filter_fields}")
    print()

    # ---- Pass 1: filter by required fields and collect field data ----
    retained = []
    filtered_out = 0
    for fpath in files:
        text = fpath.read_text()
        lines = text.splitlines()
        missing = [f for f in filter_fields if not any(l.startswith(f) for l in lines)]
        if missing:
            filtered_out += 1
            continue
        retained.append((fpath, text))

    print(f"After field filter:  {len(retained)} kept, {filtered_out} excluded")

    # ---- Pass 2: deduplicate by SPLASH (keep lowest accession) ----
    splash_groups = defaultdict(list)
    for fpath, text in retained:
        d = read_fields(text)
        splash_groups[d.get('splash', '')].append((d.get('accession', ''), fpath, text))

    kept_set = {}
    dup_count = 0
    for splash, items in splash_groups.items():
        items_sorted = sorted(items, key=lambda x: x[0])
        acc, fpath, text = items_sorted[0]
        kept_set[acc] = (fpath, text)
        dup_count += len(items_sorted) - 1

    print(f"After deduplication: {len(kept_set)} unique, {dup_count} duplicates removed")

    # ---- Pass 3: adduct correction ----
    adduct_changes = []
    adduct_unknown = []

    final_records = {}
    for acc, (fpath, text) in kept_set.items():
        d = read_fields(text)
        try:
            prec  = float(d.get('prec_mz', ''))
            mass  = float(d.get('exact_mass', ''))
        except (ValueError, TypeError):
            final_records[acc] = (fpath, text)
            continue

        current_adduct = d.get('adduct', '[M+H]+')
        best_adduct, best_diff = derive_adduct(prec, mass)

        if best_diff > TOL:
            adduct_unknown.append({
                'accession': acc,
                'name': d.get('name', ''),
                'prec_mz': prec,
                'exact_mass': mass,
                'diff': prec - mass,
                'best_guess': best_adduct,
                'best_diff': best_diff,
            })
            final_records[acc] = (fpath, text)  # keep as-is, flag for review
        elif best_adduct != current_adduct:
            new_text = fix_record(text, best_adduct, current_adduct)
            adduct_changes.append({
                'accession': acc,
                'name': d.get('name', ''),
                'from': current_adduct,
                'to': best_adduct,
                'diff_da': best_diff,
            })
            final_records[acc] = (fpath, new_text)
        else:
            final_records[acc] = (fpath, text)

    print(f"Adduct corrections:  {len(adduct_changes)}")
    print(f"Unknown adducts:     {len(adduct_unknown)} (kept with flag)")
    print()

    # ---- Write output ----
    if not args.dry_run:
        for acc, (fpath, text) in final_records.items():
            dest = outdir / fpath.name
            dest.write_text(text)
        print(f"Wrote {len(final_records)} records to {outdir}")
    else:
        print(f"[dry-run] Would write {len(final_records)} records")

    # ---- Reports ----
    print("\n=== Adduct corrections ===")
    from_counts = {}
    for c in adduct_changes:
        k = (c['from'], c['to'])
        from_counts[k] = from_counts.get(k, 0) + 1
    for (f, t), n in sorted(from_counts.items(), key=lambda x: -x[1]):
        print(f"  {n:4d}x  {f} → {t}")

    if adduct_unknown:
        print(f"\n=== Unknown adducts (diff > {TOL} Da) — needs review ===")
        print(f"{'Accession':<30} {'Name':<40} {'prec':>10} {'mass':>10} {'diff':>10} {'best':>15} {'Δ':>8}")
        for u in adduct_unknown:
            print(f"  {u['accession']:<28} {u['name'][:38]:<40} "
                  f"{u['prec_mz']:10.4f} {u['exact_mass']:10.4f} "
                  f"{u['diff']:10.4f} {u['best_guess']:>15} {u['best_diff']:8.4f}")


if __name__ == '__main__':
    main()
