"""Pin exact inherited Git blobs and LF artifact bytes; checking never writes."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
REV='b1c353a48c7062c6900a622c761c4700de74706b'
INCOMING_REV='bdbe68566e8c56d08cdc9032229f6e50cccb47fc'
SOURCES=[
    '05-knowledge/results/session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md',
    '05-knowledge/results/session114_continuation_symmetry_20261030/lrc/RESONANT_TAIL.md',
    '05-knowledge/results/session114_continuation_symmetry_20261030/lrc_audit/GENERAL_AUDIT.md',
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/regular_parent_symmetry_20261029/PROOFS.md',
    '05-knowledge/results/session37_statistics_parity_cuts_20261023/arithmetic/ARITHMETIC_CUT_COMPLETION.md',
    '01-canon/theorems/THM-468-tournament-determinant-floor.md',
    '01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md',
    '01-canon/theorems/THM-080-transfer-matrix-odd-walsh.md',
    '01-canon/theorems/THM-3043-lrc-tight-instances-are-not-only-APs-and-what-FC2-does-not-transfer.md',
    '01-canon/theorems/THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md',
    '01-canon/theorems/THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary.md',
    '05-knowledge/results/switching_branch_integrality_20261020/PAIRING_BOUNDARY.md',
    '05-knowledge/results/diagonal_palette_holonomy_20261011g.md',
    '05-knowledge/results/odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md',
    '07-reflections/lrc-dyadic-fibre-grid-ternary-carry-and-boolean-quartet-codex-20260815.md',
]
INCOMING=[
    '05-knowledge/results/session115_odd_eighth_classifications_20261030/lrc/HARMONIC_COMPLETION.md',
    '05-knowledge/results/session115_odd_eighth_classifications_20261030/quadratic/ODD_EIGHTH_PHASES.md',
    '05-knowledge/results/session115_odd_eighth_classifications_20261030/switching/ODD_EIGHTH_GROUPS.md',
    '05-knowledge/results/odd_eighths_uniform_families_20261030/LRC.md',
    '05-knowledge/results/odd_eighths_uniform_families_20261030/PHASE_BRIDGE.md',
    '05-knowledge/results/session115_odd_eighth_classifications_20261030/SYNTHESIS.md',
    '05-knowledge/results/entry_cosets_resonant_echoes_20261030/LRC.md',
]


def sha(data): return hashlib.sha256(data).hexdigest()


def blob(ref,path):
    return subprocess.check_output(['git','show',f'{ref}:{path}'],cwd=ROOT)


def longpath(path):
    value=str(path.resolve())
    return Path('\\\\?\\'+value) if os.name=='nt' and not value.startswith('\\\\?\\') else Path(value)


def build():
    artifacts={}
    for file in sorted(HERE.rglob('*')):
        if not file.is_file() or file.name=='manifest.json' or '__pycache__' in file.parts: continue
        data=file.read_bytes()
        if b'\r' in data: raise RuntimeError(f'Non-LF artifact: {file.name}')
        artifacts[file.relative_to(HERE).as_posix()]=sha(data)
    return {'source_revision':REV,'incoming_revision':INCOMING_REV,
            'hash_basis':'Exact Git source blobs; exact LF bytes for owned artifacts.',
            'sources':{p:sha(blob(REV,p)) for p in sorted(SOURCES)},
            'incoming_sources':{p:sha(blob(INCOMING_REV,p)) for p in sorted(INCOMING)},
            'artifacts':artifacts,
            'scope':'Same-agent proof audit and separate exact native algorithms; no formal proof build.',
            'runtime':'Python standard library; regular_parent_symmetry/verify.py is the sole inherited runtime import.'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true'); parser.add_argument('--staged',action='store_true')
    args=parser.parse_args(); result=build(); target=HERE/'manifest.json'
    if args.write: target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    elif json.loads(target.read_text(encoding='utf-8'))!=result:
        raise RuntimeError('Manifest mismatch')
    for section in ['sources','incoming_sources']:
        for name,digest in result[section].items():
            if sha(longpath(ROOT/name).read_bytes().replace(b'\r\n',b'\n'))!=digest:
                raise RuntimeError(f'Current source differs from pin: {name}')
    if args.staged:
        for name,digest in result['artifacts'].items():
            path=(HERE/name).relative_to(ROOT).as_posix()
            if sha(blob('',path))!=digest: raise RuntimeError(f'Staged artifact differs: {path}')
        if blob('',target.relative_to(ROOT).as_posix())!=target.read_bytes():
            raise RuntimeError('Staged manifest differs')
    print(json.dumps({'base_sources':len(SOURCES),'incoming_sources':len(INCOMING),
                      'artifacts':len(result['artifacts']),'write':args.write,'staged':args.staged}))


if __name__=='__main__': main()
