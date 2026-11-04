#!/usr/bin/env python3
"""cpq_defect_k7.py -- the C15 (3,5,1^7) row at n = 15 (HYP-9463 (2) pre-registered test: predicted defect 8652),
memory-lean version of cpq_defect_census.py: coloured types are streamed per base type, full objects are reduced to
128-bit BLAKE2b digests of their nauty certificates stored in numpy arrays, and the graph side records the sign-even
flag of every object (recomputed per object; constancy per type is asserted at the end).
Usage: python3 cpq_defect_k7.py [p q k]   (default 3 5 7)
"""
import sys, itertools, time, hashlib
import numpy as np
import pynauty
from cpq_defect_census import circulant_tournament_types, circulant_graph_types, build, reversal_parity, all_types

def digest(cert):
    d = hashlib.blake2b(cert, digest_size=16).digest()
    return int.from_bytes(d[:8], 'little'), int.from_bytes(d[8:], 'little')

def coloured_stream(k, directed):
    """yield (adj, col, even) for every 4-coloured type on k points (even = colour-preserving automorphisms all
    reversal-even; graphs only), streaming per base type."""
    for adj in all_types(k, directed):
        edges = sorted((u, v) for u in adj for v in adj[u] if u < v) if not directed else None
        seen = set()
        for col in itertools.product(range(4), repeat=k):
            cells = [set(v for v in range(k) if col[v] == c) for c in range(4)]
            present = tuple(c for c in range(4) if cells[c])
            g = pynauty.Graph(k, directed=directed, adjacency_dict=adj, vertex_coloring=[cells[c] for c in present])
            key = (tuple((c, len(cells[c])) for c in present), pynauty.certificate(g))
            if key in seen: continue
            seen.add(key)
            if directed: yield adj, col, True
            else:
                gens = pynauty.autgrp(g)[0]
                yield adj, col, all(reversal_parity(gen, edges) == 0 for gen in gens)

def main():
    p, q, k = (int(a) for a in sys.argv[1:4]) if len(sys.argv) >= 4 else (3, 5, 7)
    n = p + q + k; t0 = time.time()
    tT, tQ = circulant_tournament_types(p), circulant_tournament_types(q)
    gP = [C for C, ev in circulant_graph_types(p) if ev]; gQ = [C for C, ev in circulant_graph_types(q) if ev]
    print("(%d,%d,1^%d) n=%d: circulant classes t=(%d,%d) e=(%d,%d)" % (p, q, k, n, len(tT), len(tQ), len(gP), len(gQ)), flush=True)
    # ---- tournaments
    hashes = []; cnt = 0
    for adjF, col, _ in coloured_stream(k, True):
        for Sp in tT:
            for Sq in tQ:
                for cross in (0, 1):
                    adj = build(p, q, k, Sp, Sq, cross, adjF, col, True)
                    hashes.append(digest(pynauty.certificate(pynauty.Graph(n, directed=True, adjacency_dict=adj))))
                    cnt += 1
                    if cnt % 1000000 == 0: print("  tournaments: %d objects, %.0fs" % (cnt, time.time() - t0), flush=True)
    H = np.array(hashes, dtype=np.uint64); del hashes
    R = len(H); tau_T = len(np.unique(H, axis=0)); E_T = R - tau_T
    print("tournaments: R=%d tau_T=%d E_T=%d (%.0fs)" % (R, tau_T, E_T, time.time() - t0), flush=True); del H
    # ---- graphs (N-even orbits only): hash + sign flag per object
    hashes = []; flags = []; cnt = 0
    for adjF, col, ev in coloured_stream(k, False):
        if not ev: continue
        for Cp in gP:
            for Cq in gQ:
                for cross in (0, 1):
                    adj = build(p, q, k, Cp, Cq, cross, adjF, col, False)
                    g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
                    hashes.append(digest(pynauty.certificate(g)))
                    gens = pynauty.autgrp(g)[0]
                    edges = [(u, v) for u in adj for v in adj[u] if u < v]
                    flags.append(all(reversal_parity(gen, edges) == 0 for gen in gens))
                    cnt += 1
                    if cnt % 1000000 == 0: print("  graphs: %d objects, %.0fs" % (cnt, time.time() - t0), flush=True)
    H = np.array(hashes, dtype=np.uint64); Fl = np.array(flags, dtype=np.uint8); del hashes, flags
    RG = len(H)
    uniq, inv, counts = np.unique(H, axis=0, return_inverse=True, return_counts=True)
    inv = inv.reshape(-1)
    # flag constancy per type
    fmin = np.full(len(uniq), 2, dtype=np.uint8); fmax = np.zeros(len(uniq), dtype=np.uint8)
    np.minimum.at(fmin, inv, Fl); np.maximum.at(fmax, inv, Fl)
    assert np.array_equal(fmin, fmax), "sign flag not constant on a type: hash collision or bug"
    even = fmin.astype(bool)
    tau_E = int(even.sum()); E_even = int((counts[even] - 1).sum()); W = int((~even).sum()); E_wit = int((counts[~even] - 1).sum())
    defect = tau_T - tau_E
    print("graphs: RG=%d tau_E=%d E_even=%d W=%d E_wit=%d (%.0fs)" % (RG, tau_E, E_even, W, E_wit, time.time() - t0), flush=True)
    print("RESULT (%d,%d,1^%d) n=%d: R=%d RG=%d tau_T=%d tau_E=%d defect=%d W=%d E_T=%d E_even=%d E_wit=%d identity_ok=%s theoremR_ok=%s" % (
        p, q, k, n, R, RG, tau_T, tau_E, defect, W, E_T, E_even, E_wit, defect == W + E_even + E_wit - E_T, R == RG), flush=True)

if __name__ == "__main__":
    main()
