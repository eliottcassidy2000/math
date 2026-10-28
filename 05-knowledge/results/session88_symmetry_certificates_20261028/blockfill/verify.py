#!/usr/bin/env python3
"""Replay all finite checks; no pre-existing receipt files are required."""
from hashlib import sha256
from pathlib import Path
import json
import sys
import check_blockfill as blockfill
import check_adaptive as adaptive
import check_named_dictionary as named
import check_transport_replay as transport
import check_factored_operator as factored

HERE=Path(__file__).resolve().parent
# Hashes of normalized semantic results. Timing/benchmark measurements are
# deliberately excluded. These pin full output hashes, census/audit results,
# operation counts, scope metadata and hostile witnesses, not just class totals.
EXPECTED={
    'full_blockfill': 'c4651aa7b2e154f060a778bcf8a61267ed5835e6c490f93286ed5b41c3124cd5',
    'adaptive': 'b34d7cd870d9b18f648fbc978739571b319507ae9e40e7da29af5aff9be319aa',
    'named_dictionary': '99fc23e108eefaa105c601b493191910cd9535cfc84319445d1a924c49a0e54d',
    'transport_replay': '0640c8972adf891d54905054c391ef20de17ea92f47fc0d060a8625e07e4aab9',
    'factored_operator': '3f883910577595db2ef4fc7b9f2a579c8c22a2bd935052630e163ddeb6b698ab',
}


def deterministic(value):
    if isinstance(value,dict):
        return {str(k):deterministic(v) for k,v in value.items()
                if k not in {'seconds','preparation_seconds','benchmarks'}}
    if isinstance(value,list): return [deterministic(v) for v in value]
    if isinstance(value,tuple): return [deterministic(v) for v in value]
    return value


def semantic_hash(value):
    serialized=json.dumps(deterministic(value),sort_keys=True,separators=(',',':'))
    return sha256(serialized.encode()).hexdigest()


def main():
    # The baseline loader has already consumed its optional staging flag.
    blockfill.need(len(sys.argv)==1,'verify.py accepts only optional --baseline-root PATH')
    results=[]
    for name,run in [('full_blockfill',blockfill.main),('adaptive',adaptive.main),
                     ('named_dictionary',named.main),('transport_replay',transport.main),
                     ('factored_operator',factored.main)]:
        actual=run(); digest=semantic_hash(actual)
        blockfill.need(digest==EXPECTED[name],'semantic replay mismatch: '+name)
        results.append({'check':name,'semantic_sha256':digest,'matches_expected':True})
    sources={p.name:sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))}
    return {'status':'PASS','python_optimization_level':sys.flags.optimize,
            'checks':results,'requires_preexisting_receipts':False,'source_sha256':sources}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
