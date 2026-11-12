#!/usr/bin/env python3
"""S93 lane seeds: literal check (all 8! vertex orders, all cut points) of the gadgets of order 8 with kappa = 0 that
violate HYP-9662 (D): every one with |P| = 2 and r_X != -H(X), and 20 with |P| = 3."""
import numpy as np
from seedlib import out_from_colcode
from a_verify import literal
from c_neutral import scan, NCH, MAXREC

n = 8
reps = np.atleast_1d(np.loadtxt(f"cls_{n}.txt", dtype=np.int64))
rec = np.zeros((NCH, MAXREC, 4), np.int64); nrec = np.zeros(NCH, np.int64)
scan(n, reps, rec, nrec, np.zeros((NCH, 3), np.int64), np.zeros((NCH, 2), np.float64))
viol2, viol3 = [], []
for ch in range(NCH):
    for c in range(int(nrec[ch])):
        code, P, rX, H = (int(x) for x in rec[ch, c])
        k = bin(P).count("1")
        if k == 2 and rX != -H:
            viol2.append((code, P, rX, H))
        if k == 3:
            viol3.append((code, P, rX, H))
viol3.sort()
todo = viol2 + viol3[::max(1, len(viol3) // 20)][:20]
ok = True
for code, P, rX, H in todo:
    out = np.zeros(n, np.int64); out_from_colcode(n, code, out)
    Q2 = np.zeros((n, n), np.int64); U2 = np.zeros(n, np.int64); au = np.zeros((n, n), np.int64); H3 = np.zeros((n, n), np.int64)
    Hl, pi2 = literal(n, out, Q2, U2, au, H3)
    vs = [v for v in range(n) if (P >> v) & 1]
    w = sum(int(Q2[e, f]) for e in vs for f in vs if e != f)
    w3 = sum(int(H3[e, f]) for e in vs for f in vs if e != f)
    good = (Hl == H and Hl - w == rX and pi2 - 2 * w - 2 * w3 == 0)
    ok &= good
print(f"order 8: gadgets with kappa = 0: |P| = 2 and r_X != -H(X): {len(viol2)}; |P| = 3: {len(viol3)}; |P| >= 4: {sum(1 for ch in range(NCH) for c in range(int(nrec[ch])) if bin(int(rec[ch, c, 1])).count('1') >= 4)}")
print(f"literal enumeration confirms kappa = 0, H and r_X on all {len(viol2)} with |P| = 2 and on {len(todo) - len(viol2)} with |P| = 3: {ok}")
print("the |P| = 2 violations (colcode, P mask, r_X, H(X)):", viol2[:12])
