#!/usr/bin/env python3
"""Verifier: is 'no two-module family is additively separable' (summary) true beyond transitive/nested modules?
Search for a 2x2 zero mixed difference  sum_{j,k} c_jk dP_j dP'_k = 0  with dP = P(C1) - P(C1') != 0 and
dP' = P(C2) - P(C2') != 0, C1, C1' (resp. C2, C2') tournaments of the SAME order (so not nested), in a strong context.
Then H(T0[C1,C2]) - H(T0[C1',C2]) - H(T0[C1,C2']) + H(T0[C1',C2']) = 0: verified by direct counting."""
import itertools, sys
import numpy as np
from vlib import *
from v4_balanced_modules import walks, ordered_path_partitions

def profiles(k):
    seen = {}
    for m in range(1 << (k * (k - 1) // 2)):
        A = adjm_from_mask(k, m)
        P = tuple(ordered_path_partitions(k, A)[1:])
        seen.setdefault(P, m)
    return seen

prof = {k: profiles(k) for k in (3, 4, 5)}
for k in prof:
    print(f"order {k}: {len(prof[k])} distinct ordered-path-partition profiles: {sorted(prof[k])}")
found = 0
for n0 in (3, 4, 5):
    pi, pj = pair_arrays(n0)
    seenc = set()
    for m in range(1 << (n0 * (n0 - 1) // 2)):
        adj = np.zeros(n0, np.int64); radj = np.zeros(n0, np.int64)
        adj_from_mask(n0, m, pi, pj, adj, radj)
        if not is_strong(n0, adj, radj):
            continue
        for u in range(n0):
            for v in range(n0):
                if u == v: continue
                out = np.zeros((6, 6), np.int64)
                walks(n0, adj, u, v, 5, 5, out)
                key = tuple(out[1:, 1:].flatten().tolist())
                if key in seenc: continue
                seenc.add(key)
                c = out[1:, 1:]
                for k1 in (4, 5):
                    for k2 in (4, 5):
                        L1 = sorted(prof[k1]); L2 = sorted(prof[k2])
                        for a, b in itertools.combinations(L1, 2):
                            d1 = np.array(a) - np.array(b)
                            for a2, b2 in itertools.combinations(L2, 2):
                                d2 = np.array(a2) - np.array(b2)
                                val = int(d1 @ c[:k1, :k2] @ d2)
                                if val == 0:
                                    found += 1
                                    if found <= 3:
                                        print(f"ZERO mixed difference: context order {n0} mask {m} ports ({u},{v}); modules of order "
                                              f"{k1}: masks {prof[k1][a]}, {prof[k1][b]} (profiles {a}, {b}); order {k2}: masks "
                                              f"{prof[k2][a2]}, {prof[k2][b2]} (profiles {a2}, {b2})")
                                        # direct check
                                        def build(mc1, mc2):
                                            C1 = adjm_from_mask(k1, mc1); C2 = adjm_from_mask(k2, mc2)
                                            A0 = adjm_from_mask(n0, m)
                                            verts = []
                                            for x in range(n0):
                                                if x == u: verts += [('u', i) for i in range(k1)]
                                                elif x == v: verts += [('v', i) for i in range(k2)]
                                                else: verts.append(('o', x))
                                            N = len(verts)
                                            ad = np.zeros(N, np.int64)
                                            for p in range(N):
                                                for q in range(N):
                                                    if p == q: continue
                                                    (tp, ip), (tq, iq) = verts[p], verts[q]
                                                    if tp == 'u' and tq == 'u': ar = C1[ip][iq]
                                                    elif tp == 'v' and tq == 'v': ar = C2[ip][iq]
                                                    else:
                                                        xp = u if tp == 'u' else (v if tp == 'v' else ip)
                                                        xq = u if tq == 'u' else (v if tq == 'v' else iq)
                                                        ar = A0[xp][xq]
                                                    if ar: ad[p] |= 1 << q
                                            return int(hp(N, ad, np.zeros((1 << N, N), np.int64)))
                                        h = [build(prof[k1][a], prof[k2][a2]), build(prof[k1][b], prof[k2][a2]),
                                             build(prof[k1][a], prof[k2][b2]), build(prof[k1][b], prof[k2][b2])]
                                        print(f"   direct counts H = {h}; mixed difference {h[0]-h[1]-h[2]+h[3]}")
    print(f"context order {n0}: distinct c-matrices {len(seenc)}; zero mixed differences found so far: {found}")
    sys.stdout.flush()
