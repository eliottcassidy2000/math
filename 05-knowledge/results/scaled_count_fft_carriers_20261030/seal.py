"""Pin inherited Git blobs and exact packet bytes; checking never mutates."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
REV='b51ef5cd739375e6475ceafd42e36e86858c70ca'
INCOMING_REV='d32897c77717c6d8a5f4fa332082feb26a256e0b'
INCOMING=[
    '05-knowledge/results/session112_counting_interfaces_20261030/arithmetic/MODULAR_OBSERVERS.md',
    '05-knowledge/results/session112_counting_interfaces_20261030/norms/NORM_SIGNS.md',
    '05-knowledge/results/session112_counting_interfaces_20261030/norms/NORM_DEFECT.md',
    '05-knowledge/results/session112_counting_interfaces_20261030/arithmetic/POWER_OBSERVERS.md',
    '05-knowledge/results/pinned_ramified_counts_20261030/README.md',
    '05-knowledge/results/pinned_ramified_counts_20261030/ARITHMETIC_AND_SOURCE.md',
]
SOURCES=[
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
    '05-knowledge/results/regular_parent_symmetry_20261029/PROOFS.md',
    '05-knowledge/results/prime_power_plus_one_20261029/README.md',
    '05-knowledge/results/prime_power_plus_one_20261029/ARITHMETIC.md',
    '05-knowledge/results/prime_power_plus_one_20261029/GROWTH.md',
    '05-knowledge/results/session111_scale_connections_20261030/README.md',
    '05-knowledge/results/anchor_content_flat_quotients_20261030/README.md',
    '05-knowledge/results/anchor_content_flat_quotients_20261030/TOURNAMENTS.md',
    '05-knowledge/results/anchor_content_flat_quotients_20261030/COLLATZ.md',
    '05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md',
    '05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md',
    '05-knowledge/results/arithmetic_braids2_20260917_fano_code.md',
    '05-knowledge/results/four_vertex_glue_fermat_20261029/README.md',
    '05-knowledge/results/metagraph_scheduled_recursion_20261029/README.md',
    '05-knowledge/results/cursor_cm_phase_20261028/TRACE_CARRIER.md',
    '05-knowledge/results/cursor_cm_phase_20261028/RAMANUJAN.md',
    '05-knowledge/results/pinned_discrepancy_reset_tails_20261029/FINITE_DEPTH.md',
]


def sha(data): return hashlib.sha256(data).hexdigest()


def blob(ref,path):
    return subprocess.check_output(['git','show',f'{ref}:{path}'],cwd=ROOT)


def longpath(path):
    resolved=str(path.resolve())
    return Path('\\\\?\\'+resolved) if os.name=='nt' and not resolved.startswith('\\\\?\\') else Path(resolved)


def build():
    artifacts={}
    for file in sorted(HERE.rglob('*')):
        if not file.is_file() or file.name=='manifest.json' or '__pycache__' in file.parts: continue
        data=file.read_bytes()
        if file.suffix!='.pdf' and b'\r' in data:
            raise RuntimeError(f'Non-LF text artifact: {file.name}')
        artifacts[file.relative_to(HERE).as_posix()]=sha(data)
    return {'source_revision':REV,
            'hash_basis':'Inherited sources: exact Git blobs. Packet text: exact LF bytes. PDF: exact binary bytes.',
            'sources':{p:sha(blob(REV,p)) for p in sorted(SOURCES)},
            'incoming_revision':INCOMING_REV,
            'incoming_nonruntime_sources':{p:sha(blob(INCOMING_REV,p)) for p in sorted(INCOMING)},
            'artifacts':artifacts,
            'scope':'Same-agent proof audit and alternative native calculations; separate pinned external Lean build.',
            'runtime':'Python standard library; native regular_parent_symmetry verifier is the sole inherited runtime import.'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--staged',action='store_true')
    args=parser.parse_args()
    result=build(); target=HERE/'manifest.json'
    if args.write:
        target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    elif json.loads(target.read_text(encoding='utf-8'))!=result:
        raise RuntimeError('Manifest does not match current source or artifact bytes')
    for path,expected in result['sources'].items():
        if sha(longpath(ROOT/path).read_bytes().replace(b'\r\n',b'\n'))!=expected:
            raise RuntimeError(f'Inherited source changed: {path}')
    if args.staged:
        for name,expected in result['artifacts'].items():
            path=(HERE/name).relative_to(ROOT).as_posix()
            if sha(blob('',path))!=expected: raise RuntimeError(f'Staged bytes differ: {path}')
        path=target.relative_to(ROOT).as_posix()
        if blob('',path)!=target.read_bytes(): raise RuntimeError('Staged manifest differs')
    print(json.dumps({'sources':len(result['sources']),'artifacts':len(result['artifacts']),
                      'write':args.write,'staged':args.staged}))


if __name__=='__main__': main()
