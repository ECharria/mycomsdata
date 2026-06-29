#!/usr/bin/env python3
"""
library_to_massbank.py
----------------------
Convert a Bruker MetaboScape .library file into individual MassBank3-compatible
.txt files ready for import into MycoMSBase.

Structure enrichment: for entries that lack SMILES/InChI in the .library file,
a companion Excel spreadsheet (--excel) is used to look up compound info by the
numeric STO ID embedded in the library Name prefix (e.g. "4204_Flavipucine").

Usage:
    python3 library_to_massbank.py \
        --input   /path/to/file.library \
        --excel   /path/to/compounds.xlsx \
        --outdir  /path/to/output_dir \
        --start   130 \
        --accession-prefix HZI-CBIO-AA-
"""

import argparse
import hashlib
import re
from pathlib import Path
from collections import Counter
from typing import List, Tuple, Optional, Dict

from splash import Splash, Spectrum, SpectrumType

# ---------------------------------------------------------------------------
# Monoisotopic masses
# ---------------------------------------------------------------------------
MONOISOTOPIC = {
    'H': 1.007825032, 'C': 12.0,        'N': 14.003074,
    'O': 15.99491462, 'P': 30.973762,   'S': 31.972071,
    'Cl': 34.968853,  'Br': 78.918338,  'F': 18.9984032,
    'I': 126.90447,   'Si': 27.9769265,
}

def exact_mass_from_formula(formula: str) -> Optional[str]:
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    counts = Counter({el: int(n) if n else 1 for el, n in tokens if el in MONOISOTOPIC})
    if not counts:
        return None
    return f"{sum(MONOISOTOPIC[el] * n for el, n in counts.items()):.4f}"

def inchikey_from_inchi(inchi: str, _cache: dict = {}) -> Optional[str]:
    """Derive InChIKey from InChI. Tries rdkit, then PubChem REST API (network), then None."""
    if inchi in _cache:
        return _cache[inchi]

    # 1. rdkit (local, no network)
    try:
        from rdkit.Chem.inchi import InchiToInchiKey
        key = InchiToInchiKey(inchi)
        _cache[inchi] = key
        return key
    except Exception:
        pass

    # 2. pubchempy (network)
    try:
        import pubchempy as pcp
        results = pcp.get_compounds(inchi, 'inchi')
        if results:
            key = results[0].inchikey
            _cache[inchi] = key
            return key
    except Exception:
        pass

    # 3. PubChem REST API directly (network, no extra dependencies)
    try:
        import urllib.request, urllib.parse, json, time
        encoded = urllib.parse.quote(inchi, safe='')
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/property/InChIKey/JSON'
        req = urllib.request.Request(
            url,
            data=f'inchi={encoded}'.encode(),
            headers={'Content-Type': 'application/x-www-form-urlencoded'},
        )
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read())
        key = data['PropertyTable']['Properties'][0]['InChIKey']
        _cache[inchi] = key
        time.sleep(0.12)  # respect PubChem rate limit
        return key
    except Exception:
        pass

    _cache[inchi] = None
    return None

# ---------------------------------------------------------------------------
# Strip leading "NNN_" prefix from MetaboScape Name field
# ---------------------------------------------------------------------------
_prefix_re = re.compile(r'^(\d+)_')

def parse_name(raw: str) -> Tuple[Optional[int], str]:
    """Return (sto_id_or_None, clean_name)."""
    m = _prefix_re.match(raw.strip())
    if m:
        return int(m.group(1)), raw[m.end():].strip()
    return None, raw.strip()

# ---------------------------------------------------------------------------
# Load Excel compound table
# ---------------------------------------------------------------------------
def load_excel(path: Path) -> Dict[int, dict]:
    """Return dict keyed by STO Myxobase key compound (int) -> row dict."""
    try:
        import openpyxl
    except ImportError:
        raise SystemExit("openpyxl required: pip install openpyxl")

    wb = openpyxl.load_workbook(str(path))
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    headers = rows[0]
    result = {}
    for row in rows[1:]:
        d = dict(zip(headers, row))
        sto_id = d.get('STO Myxobase key compound')
        if sto_id is not None:
            try:
                result[int(sto_id)] = d
            except (ValueError, TypeError):
                pass
    return result

# ---------------------------------------------------------------------------
# Parser for the Bruker .library format
# ---------------------------------------------------------------------------
def parse_library(path: Path) -> List[dict]:
    text = path.read_text(encoding='utf-8', errors='replace')
    lines = [l.rstrip('\r\n') for l in text.splitlines()]

    records: List[dict] = []
    current: dict = {}
    peaks: List[Tuple[float, float]] = []
    in_peaks = False

    def flush():
        if current and peaks:
            current['PEAKS'] = peaks[:]
            records.append(dict(current))
        elif current:
            print(f"  WARNING: skipping '{current.get('NAME', '?')}' — no peaks")

    for line in lines:
        stripped = line.strip()
        if not stripped:
            in_peaks = False
            continue

        if stripped.startswith('Name:'):
            flush()
            raw_name = stripped.split(':', 1)[1].strip()
            sto_id, clean = parse_name(raw_name)
            current = {'NAME': clean, 'STO_ID': sto_id}
            peaks = []
            in_peaks = False

        elif ':' in stripped and not stripped[0].isdigit():
            in_peaks = False
            key, _, val = stripped.partition(':')
            key = key.strip()
            val = val.strip()
            if key == 'Formula':       current['FORMULA'] = val
            elif key == 'Smiles':      current['SMILES'] = val
            elif key == 'InChI':
                # InChI field may contain a second colon — re-split
                current['INCHI'] = stripped.split(':', 1)[1].strip()
            elif key == 'InChIKey':    current['INCHIKEY'] = val
            elif key == 'InstType':    current['INST_TYPE'] = val
            elif key == 'IoniMethod':  current['IONI_METHOD'] = val
            elif key == 'IonPolarity': current['ION_POLARITY'] = val.upper()
            elif key == 'MSMS':        current['MS_LEVEL'] = int(val)
            elif key == 'PreIon':      current['PRECURSOR_MZ'] = float(val)
            elif key == 'CCS':         current['CCS'] = float(val)
            elif key == 'ColEnergy':   current['COL_ENERGY'] = float(val)
            elif key == 'Num Peaks':
                current['NUM_PEAKS'] = int(val)
                in_peaks = True
            # Structure (hex-encoded molfile) and CommentSpec ignored

        elif in_peaks:
            parts = stripped.split()
            try:
                for i in range(0, len(parts) - 1, 2):
                    peaks.append((float(parts[i]), float(parts[i + 1])))
            except ValueError:
                pass

    flush()
    return records

# ---------------------------------------------------------------------------
# Load compound metadata from the live MycoMSBase PostgreSQL database
# ---------------------------------------------------------------------------
def load_db_metadata(tsv_path: str) -> Dict[str, dict]:
    """Load curated species/publication from a pipe-delimited TSV: name|species|publication."""
    result = {}
    for line in Path(tsv_path).read_text().splitlines():
        parts = line.split('|')
        if len(parts) >= 3:
            name, species, pub = parts[0].strip(), parts[1].strip(), parts[2].strip()
            if name and (species or pub):
                result[name] = {'species': species, 'doi': pub}
    return result

def load_db_classes(tsv_path: str) -> Dict[str, str]:
    """Load curated biosynthetic classes from a pipe-delimited TSV: name|class."""
    result = {}
    for line in Path(tsv_path).read_text().splitlines():
        parts = line.split('|')
        if len(parts) >= 2:
            name, cls = parts[0].strip(), parts[1].strip()
            if name and cls:
                result[name] = cls
    return result

def load_db_compounds(tsv_path: Optional[str] = None) -> Dict[str, dict]:
    """Return dict keyed by lowercase compound name -> {smiles, inchi, inchikey, formula}.

    If tsv_path is given, load from a pre-exported TSV (pipe-delimited:
    name|smiles|inchi|inchikey|formula).  Otherwise try psycopg2 direct connection.
    """
    if tsv_path and Path(tsv_path).exists():
        result = {}
        for line in Path(tsv_path).read_text().splitlines():
            parts = line.split('|')
            if len(parts) < 5:
                continue
            name, smiles, inchi, inchikey, formula = parts[:5]
            if name and name not in result:
                result[name] = {'smiles': smiles, 'inchi': inchi,
                                'inchikey': inchikey, 'formula': formula}
        return result

    try:
        import psycopg2
        conn = psycopg2.connect('dbname=mycomsbase user=mycomsbase host=localhost')
        cur = conn.cursor()
        cur.execute("""
            SELECT LOWER(cn.name), c.smiles, c.inchi, b.inchikey, c.formula
            FROM compound_name cn
            JOIN compound c ON c.id = cn.compound_id
            JOIN browse_options b ON b.massbank_id = cn.massbank_id
            WHERE c.smiles != 'C'
            GROUP BY cn.name, c.smiles, c.inchi, b.inchikey, c.formula
        """)
        result = {}
        for name, smiles, inchi, inchikey, formula in cur.fetchall():
            if name not in result:
                result[name] = {'smiles': smiles, 'inchi': inchi,
                                'inchikey': inchikey, 'formula': formula}
        conn.close()
        return result
    except Exception as e:
        print(f"  DB enrichment unavailable: {e}")
        return {}

# ---------------------------------------------------------------------------
# NP Atlas metadata lookup (producer organism + original reference)
# ---------------------------------------------------------------------------
def lookup_npatlas_metadata(name: str, _cache: dict = {}) -> Optional[dict]:
    """Return {'species': str, 'doi': str} from NP Atlas for a compound name, or None."""
    if name in _cache:
        return _cache[name]
    import urllib.request, urllib.parse, json, time
    headers = {'User-Agent': 'Mozilla/5.0 (compatible; research-script/1.0)'}

    def _fetch_detail(npaid: str) -> Optional[dict]:
        url = f'https://www.npatlas.org/api/v1/compound/{npaid}'
        try:
            req = urllib.request.Request(url, headers=headers)
            data = json.loads(urllib.request.urlopen(req, timeout=10).read())
            time.sleep(0.1)
            org = data.get('origin_organism') or {}
            ref = data.get('origin_reference') or {}
            species = f"{org.get('genus','')} {org.get('species','')}".strip() if org else ''
            doi = ref.get('doi', '') if ref else ''
            if not (species or doi):
                return None
            return {'species': species, 'doi': f'doi:{doi}' if doi else ''}
        except Exception:
            return None

    # 1. Exact name match
    enc = urllib.parse.quote(name)
    try:
        req = urllib.request.Request(
            f'https://www.npatlas.org/api/v1/compounds?name={enc}&exact=true', headers=headers)
        hits = json.loads(urllib.request.urlopen(req, timeout=10).read())
        time.sleep(0.1)
        if hits:
            r = _fetch_detail(hits[0]['npaid'])
            _cache[name] = r
            return r
    except Exception:
        pass

    # 2. basicSearch POST (fuzzy)
    try:
        req = urllib.request.Request(
            f'https://www.npatlas.org/api/v1/compounds/basicSearch?name={enc}&limit=5',
            method='POST', headers=headers)
        hits = json.loads(urllib.request.urlopen(req, timeout=10).read())
        time.sleep(0.1)
        name_l = name.lower()
        match = next((h for h in hits if h.get('original_name','').lower() == name_l), None)
        if not match:
            match = next((h for h in hits if name_l in h.get('original_name','').lower()), None)
        if match:
            r = _fetch_detail(match['npaid'])
            _cache[name] = r
            return r
    except Exception:
        pass

    _cache[name] = None
    return None

# ---------------------------------------------------------------------------
# Name-based lookup via PubChem REST API
# ---------------------------------------------------------------------------
def classify_by_npclassifier(smiles: str, _cache: dict = {}) -> Optional[str]:
    """Return NPClassifier pathway string(s) for a SMILES, e.g. 'Polyketides' or None."""
    if smiles in _cache:
        return _cache[smiles]
    import urllib.request, urllib.parse, json, time
    enc = urllib.parse.quote(smiles)
    url = f'https://npclassifier.gnps2.org/classify?smiles={enc}'
    try:
        headers = {'User-Agent': 'Mozilla/5.0 (compatible; research-script/1.0)'}
        req = urllib.request.Request(url, headers=headers)
        data = json.loads(urllib.request.urlopen(req, timeout=15).read())
        pathway = data.get('pathway_results', [])
        result = ' / '.join(pathway) if pathway else None
        _cache[smiles] = result
        time.sleep(0.15)
        return result
    except Exception:
        _cache[smiles] = None
        return None

def lookup_by_name_pubchem(name: str, _cache: dict = {}) -> Optional[dict]:
    """Query PubChem for SMILES/InChI/InChIKey/formula by compound name."""
    if name in _cache:
        return _cache[name]
    import urllib.request, urllib.parse, json, time
    enc = urllib.parse.quote(name)
    url = (f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc}'
           f'/property/InChIKey,InChI,IsomericSMILES,MolecularFormula/JSON')
    try:
        resp = urllib.request.urlopen(url, timeout=10)
        p = json.loads(resp.read())['PropertyTable']['Properties'][0]
        result = {
            'smiles':   p.get('IsomericSMILES', ''),
            'inchi':    p.get('InChI', ''),
            'inchikey': p.get('InChIKey', ''),
            'formula':  p.get('MolecularFormula', ''),
        }
        _cache[name] = result
        time.sleep(0.12)
        return result
    except Exception:
        _cache[name] = None
        return None

# ---------------------------------------------------------------------------
# Enrich records: .library → DB by name → Excel by STO ID → PubChem by name
# ---------------------------------------------------------------------------
def enrich_records(records: List[dict], excel: Dict[int, dict],
                   db_compounds: Dict[str, dict]) -> None:
    from_lib = from_db = from_excel = from_pubchem = still_unknown = 0

    # Collect unique unknown names first so we can batch-report progress
    for rec in records:
        # Already has real structure from the .library file itself
        if rec.get('SMILES') and rec['SMILES'] != 'C':
            from_lib += 1
            continue

        # Try DB lookup by compound name
        name_key = rec['NAME'].lower()
        if name_key in db_compounds:
            row = db_compounds[name_key]
            rec['SMILES']   = row['smiles']
            rec['INCHI']    = row['inchi']
            rec['INCHIKEY'] = row['inchikey']
            rec.setdefault('FORMULA', row['formula'])
            from_db += 1
            continue

        # Try Excel lookup by STO ID
        sto_id = rec.get('STO_ID')
        if sto_id and sto_id in excel:
            row = excel[sto_id]
            inchi   = row.get('InChI') or ''
            formula = (row.get('STO Formula') or '').strip()
            if inchi:
                rec['INCHI'] = inchi
                rec.setdefault('FORMULA', formula)
                ikey = inchikey_from_inchi(inchi)
                if ikey:
                    rec['INCHIKEY'] = ikey
                from_excel += 1
                continue

        # Fall back to PubChem name lookup (network)
        hit = lookup_by_name_pubchem(rec['NAME'])
        if hit and (hit.get('smiles') or hit.get('inchi')):
            if hit.get('smiles'):
                rec['SMILES'] = hit['smiles']
            if hit.get('inchi'):
                rec['INCHI'] = hit['inchi']
            if hit.get('inchikey'):
                rec['INCHIKEY'] = hit['inchikey']
            if hit.get('formula'):
                rec.setdefault('FORMULA', hit['formula'])
            from_pubchem += 1
            continue

        still_unknown += 1

    print(f"  {from_lib} from .library, {from_db} from DB, "
          f"{from_excel} from Excel, {from_pubchem} from PubChem name lookup, "
          f"{still_unknown} still unknown")

    # Metadata pass: producer organism + publication — curated DB first, then NP Atlas
    print("Enriching organism/publication metadata...")
    meta_db = meta_npa = 0
    for rec in records:
        if rec.get('SPECIES') or rec.get('PUBLICATION'):
            continue
        npa = lookup_npatlas_metadata(rec['NAME'])
        if npa:
            if npa.get('species'):
                rec['SPECIES'] = npa['species']
            if npa.get('doi'):
                rec['PUBLICATION'] = npa['doi']
            meta_npa += 1
    print(f"  {meta_npa} records enriched with NP Atlas organism/publication")

    # Biosynthetic class pass: NPClassifier by SMILES for records that have a real structure
    print("Deriving biosynthetic class via NPClassifier...")
    class_hits = 0
    seen_smiles: dict = {}  # smiles -> class (cache across records with same compound)
    for rec in records:
        if rec.get('COMPOUND_CLASS'):
            continue
        smiles = rec.get('SMILES', '')
        if not smiles or smiles == 'C':
            continue
        if smiles not in seen_smiles:
            seen_smiles[smiles] = classify_by_npclassifier(smiles)
        cls = seen_smiles[smiles]
        if cls:
            rec['COMPOUND_CLASS'] = cls
            class_hits += 1
    print(f"  {class_hits} records classified via NPClassifier")

# ---------------------------------------------------------------------------
# SPLASH
# ---------------------------------------------------------------------------
def make_splash(peaks: List[Tuple[float, float]]) -> str:
    return Splash().splash(Spectrum(peaks, SpectrumType.MS))

# ---------------------------------------------------------------------------
# MassBank record builder
# ---------------------------------------------------------------------------
def to_massbank(record: dict, accession: str) -> str:
    name     = record['NAME']
    formula  = (record.get('FORMULA') or 'CH4').strip()
    # Normalize formula: remove count=1 from single-letter elements like N1→N
    formula  = re.sub(r'([A-Z][a-z]?)1(?!\d)', r'\1', formula)
    exact    = exact_mass_from_formula(formula) or '16.0313'
    smiles   = record.get('SMILES') or 'C'
    inchi    = record.get('INCHI') or 'InChI=1S/CH4/h1H4'
    inchikey = record.get('INCHIKEY') or 'N/A'
    ms_level = record.get('MS_LEVEL', 2)
    ce       = record.get('COL_ENERGY', 0)
    polarity = 'POSITIVE' if record.get('ION_POLARITY', 'POS').startswith('POS') else 'NEGATIVE'
    precursor = record.get('PRECURSOR_MZ', 0.0)
    peaks    = record['PEAKS']

    instrument      = 'timsTOF Pro, Bruker Daltonics'
    instrument_type = 'timsTOF'

    title = f"{name}; {instrument_type}; MS{ms_level}; CE:{ce}; [M+H]+"

    max_int = max(i for _, i in peaks) or 1.0
    peak_lines = [
        f"  {mz:.5f} {intensity:.6g} {int(round(intensity / max_int * 999))}"
        for mz, intensity in peaks
    ]

    ccs_line = (
        [f"AC$MASS_SPECTROMETRY: CCS {record['CCS']}"]
        if 'CCS' in record else []
    )

    species          = (record.get('SPECIES') or '').strip()
    publication      = (record.get('PUBLICATION') or '').strip()
    compound_class   = (record.get('COMPOUND_CLASS') or '').strip()
    species_line     = [f"SP$SCIENTIFIC_NAME: {species}"]         if species        else []
    pub_line         = [f"PUBLICATION: {publication}"]             if publication    else []
    class_line       = [f"CH$COMPOUND_CLASS: {compound_class}"]   if compound_class else []

    lines = [
        f"ACCESSION: {accession}",
        f"RECORD_TITLE: {title}",
        "DATE: 2026.06.26",
        "AUTHORS: Charrià-Girón E, Surup F",
        "LICENSE: CC BY",
        "COPYRIGHT: Helmholtz Centre for Infection Research (HZI)",
        *pub_line,
        *species_line,
        "COMMENT: Converted from Bruker MetaboScape .library export",
        f"CH$NAME: {name}",
        *class_line,
        f"CH$FORMULA: {formula}",
        f"CH$EXACT_MASS: {exact}",
        f"CH$SMILES: {smiles}",
        f"CH$IUPAC: {inchi}",
        f"CH$LINK: INCHIKEY {inchikey}",
        f"AC$INSTRUMENT: {instrument}",
        f"AC$INSTRUMENT_TYPE: {instrument_type}",
        f"AC$MASS_SPECTROMETRY: MS_TYPE MS{ms_level}",
        f"AC$MASS_SPECTROMETRY: ION_MODE {polarity}",
        f"AC$MASS_SPECTROMETRY: COLLISION_ENERGY {ce}",
        *ccs_line,
        f"MS$FOCUSED_ION: BASE_PEAK {precursor:.5f}",
        f"MS$FOCUSED_ION: PRECURSOR_M/Z {precursor:.5f}",
        "MS$FOCUSED_ION: PRECURSOR_TYPE [M+H]+",
        "MS$DATA_PROCESSING: Converted with library_to_massbank.py",
        f"PK$SPLASH: {make_splash(peaks)}",
        f"PK$NUM_PEAK: {len(peaks)}",
        "PK$PEAK: m/z int. rel.int.",
        *peak_lines,
        "//",
    ]
    return "\n".join(lines)

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input',  required=True, help='Input .library file')
    ap.add_argument('--excel',  default=None,  help='Companion Excel with compound info')
    ap.add_argument('--outdir', required=True, help='Output directory for .txt files')
    ap.add_argument('--start',  type=int, default=1, help='First accession number')
    ap.add_argument('--accession-prefix', default='HZI-CBIO-AA-')
    ap.add_argument('--db-tsv', default=None,
                    help='Pre-exported TSV from DB (name|smiles|inchi|inchikey|formula)')
    ap.add_argument('--db-meta-tsv', default=None,
                    help='Pre-exported TSV with curated metadata (name|species|publication)')
    ap.add_argument('--db-class-tsv', default=None,
                    help='Pre-exported TSV with curated biosynthetic classes (name|class)')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = parse_library(Path(args.input))
    print(f"Parsed {len(records)} records from library")

    excel = {}
    if args.excel:
        excel = load_excel(Path(args.excel))
        print(f"Loaded {len(excel)} compounds from Excel")

    # Pre-load curated biosynthetic classes from DB export
    if args.db_class_tsv and Path(args.db_class_tsv).exists():
        db_classes = load_db_classes(args.db_class_tsv)
        print(f"Loaded {len(db_classes)} curated class entries from DB")
        for rec in records:
            cls = db_classes.get(rec['NAME'].lower())
            if cls:
                rec['COMPOUND_CLASS'] = cls

    # Pre-load curated species/publication metadata from DB export
    if args.db_meta_tsv and Path(args.db_meta_tsv).exists():
        db_meta = load_db_metadata(args.db_meta_tsv)
        print(f"Loaded {len(db_meta)} curated metadata entries from DB")
        for rec in records:
            m = db_meta.get(rec['NAME'].lower())
            if m:
                if m.get('species'):
                    rec['SPECIES'] = m['species']
                if m.get('doi'):
                    rec['PUBLICATION'] = m['doi']

    print("Enriching structures...")
    db_compounds = load_db_compounds(args.db_tsv)
    print(f"  Loaded {len(db_compounds)} named compounds from DB")
    enrich_records(records, excel, db_compounds)

    written = 0
    for i, rec in enumerate(records, start=args.start):
        accession = f"MSBNK-{args.accession_prefix}{i:06d}"
        mb = to_massbank(rec, accession)
        (outdir / f"{accession}.txt").write_text(mb)
        written += 1

    print(f"Wrote {written} MassBank files to {outdir.resolve()}")
    print(f"Accessions: MSBNK-{args.accession_prefix}{args.start:06d} "
          f"→ MSBNK-{args.accession_prefix}{args.start + written - 1:06d}")

if __name__ == '__main__':
    main()
