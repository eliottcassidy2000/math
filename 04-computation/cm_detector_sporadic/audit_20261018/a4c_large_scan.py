"""Audit extension of check 4/8: (B1') and (B1'_ord) survivors over 0 < |m| <= M for large M (chunked),
residue tables mod p^2 for 7 <= p <= 199 as in a4_detector_scan.py.  usage: python a4c_large_scan.py M [families]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import seq, FAMILIES
from a4_detector_scan import tables, PRIMES

M = int(sys.argv[1])
fams = sys.argv[2:] or list(FAMILIES)
CH = 2_000_000
res = {}
for name in fams:
    t0 = time.time()
    u = seq(name, 420)
    tabs = [(p, *tables(u, p)) for p in PRIMES]
    sB, sO = [], []
    for lo in range(-M, M + 1, CH):
        ms = np.arange(lo, min(lo + CH, M + 1), dtype=np.int64)
        ms = ms[ms != 0]
        aB = np.ones(ms.size, dtype=bool)
        aO = np.ones(ms.size, dtype=bool)
        for p, pB, pO in tabs:
            idx = ms % (p * p)
            aO &= pO[idx]
            aB &= pB[idx]
        sB += ms[aB].tolist()
        sO += ms[aO].tolist()
    res[name] = dict(B1=sB, B1ord=sO)
    print(f'{name:8s} M={M}: (B1\') {sB}; (B1\'_ord) {sO}  [{time.time() - t0:.0f}s]', flush=True)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'a4c_large_scan_M{M}.json'), 'w') as f:
    json.dump(res, f)
