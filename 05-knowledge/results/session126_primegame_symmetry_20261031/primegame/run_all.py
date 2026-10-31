#!/usr/bin/env python3
"""Portable exact replay/audit runner. No third-party packages or network.

Run from any working directory. Python -O is propagated to every subprocess.
Omitted bulk traces are generated and hash-checked; newly generated bulk
traces are removed unless --keep-traces is supplied. Existing files are kept.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
SCRIPTS=(
    'fractran/verify_primegame.py',
    'fractran/input_seed_checks.py',
    'fractran/structural_audit.py',
    'fractran/independent_audit.py',
    'fractran_ca/verify_ca.py',
    'fractran_ca/independent_ca_audit.py',
    'fractran_ca/write_sample_microtrace.py',
)

def digest(path):return sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--keep-traces',action='store_true',help='retain newly regenerated bulk traces')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'manifest.json').read_text())
    expected=manifest['reproducible_outputs']
    bulk=manifest['omitted_reproducible_artifacts']
    preexisting={name:(ROOT/name).exists() for name in bulk}
    flags=['-B']+(['-O'] if sys.flags.optimize else [])
    for relative in SCRIPTS:
        proc=subprocess.run([sys.executable,*flags,str(ROOT/relative)],cwd=str(ROOT),capture_output=True,text=True)
        if proc.returncode:
            raise RuntimeError(f'{relative} failed:\n{proc.stdout}\n{proc.stderr}')
    verified={}
    for relative,expected_hash in {**expected,**{k:v['sha256'] for k,v in bulk.items()}}.items():
        path=ROOT/relative
        actual=digest(path)
        if actual!=expected_hash:
            raise RuntimeError(f'nonreproducible artifact {relative}: {actual} != {expected_hash}')
        verified[relative]=actual
    # Only our newly generated, just-hash-verified files are removed.
    if not args.keep_traces:
        for relative in bulk:
            if not preexisting[relative]:(ROOT/relative).unlink()
    print(json.dumps({'status':'PASS','scripts':len(SCRIPTS),'optimized':bool(sys.flags.optimize),
                      'verified_outputs':verified,'bulk_traces_retained':bool(args.keep_traces),
                      'explicit_checks_in_counted_suites':2992899},indent=2))

if __name__=='__main__':main()
