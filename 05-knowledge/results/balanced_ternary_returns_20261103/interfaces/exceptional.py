#!/usr/bin/env python3
"""Tiny exact controls for the proved five-family intersection lemma."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    families = ((3,-1),(6,1),(5,-1),(13,3),(14,3))
    hits = []
    comparisons = 0
    for m in range(2,65):
        target = 15*2**m-1
        for k in range(1,67):
            for index,(a,b) in enumerate(families,1):
                comparisons += 1
                if target == a*2**k+b:
                    hits.append({'m':m,'k':k,'family':index,'H':target})
    need(hits == [{'m':2,'k':2,'family':5,'H':59}], 'only bounded intersection is59')
    need(15*2-1 == 13*2+3 == 29, 'excluded numerical boundary m1')
    sources = [
        '05-knowledge/results/two_entry_funnels_20261103/READY.md',
        '05-knowledge/results/ready_return_pentagon_20261103/OBSTRUCTION.md',
        '05-knowledge/results/orbit_spectrum_session_20261103/coverage/NOTE.md']
    print(json.dumps({'status':'PASS','scope':'FINITE-EXACT controls m2..64,k1..66; all-exponent proof is EXCEPTIONAL.md.',
                      'comparisons':comparisons,'hits':hits,
                      'excluded_numeric_boundary':{'m':1,'k':1,'family':4,'H':29},
                      'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
