#!/usr/bin/env python3
"""increment_ladder.py -- Claim A at the new vertex: for X of order n and a cut (A, B) defining v,
H(X + v) - H(X) = 2 * sum over odd cycles C through v of H((X + v) - V(C)).  Tabulates, for every (X, cut) with
increment k = (H(X+v) - h)/2 in {1, 2, 3, 4}: the number of odd cycles through v and the multiset of complement counts;
checks the identity; and lists the minimal positive increment per h."""
import sys, itertools
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from ext_spectrum_true import tournaments, ext_spectrum
from spine_numeration import ham_paths

def odd_cycles_through(n, adj, v):
    """list of vertex sets of odd directed cycles through v in the tournament (adj out-masks, n vertices)."""
    out = []
    def rec(x, S, length):
        nxt = adj[x]
        if (nxt >> v) & 1 and length >= 3 and length % 2 == 1: out.append(S)
        m = nxt & ~S & ~(1 << v)
        while m:
            y = (m & -m).bit_length() - 1; m &= m - 1
            rec(y, S | (1 << y), length + 1)
    rec(v, 1 << v, 1)
    return out

def sub_count(n, adj, S):
    """Hamiltonian paths of the sub-tournament on vertex set S (bitmask); 1 if empty."""
    verts = [x for x in range(n) if (S >> x) & 1]
    if not verts: return 1
    idx = {x: i for i, x in enumerate(verts)}
    a = [sum(1 << idx[y] for y in verts if (adj[x] >> y) & 1) for x in verts]
    return ham_paths(len(verts), a)

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(3, NMAX + 1):
        table = defaultdict(Counter); checked = 0; bad = 0; minpos = {}
        for adj in tournaments(n):
            Hs = ext_spectrum(n, adj); h = int(Hs[0])
            for A in range(1 << n):
                k2 = int(Hs[A]) - h
                if k2 == 0: continue
                k = k2 // 2
                if k < minpos.get(h, 10 ** 9): minpos[h] = k
                if k <= 4:
                    adj2 = [adj[u] | (0 if (A >> u) & 1 else (1 << n)) for u in range(n)] + [A]
                    cyc = odd_cycles_through(n + 1, adj2, n)
                    full = (1 << (n + 1)) - 1
                    comps = sorted(sub_count(n + 1, adj2, full & ~S) for S in cyc)
                    checked += 1
                    if 2 * sum(comps) != k2: bad += 1
                    table[k][(h, len(cyc), tuple(comps))] += 1
        print("n=%d: Claim A at v checked on %d cuts with increment <= 8: mismatches %d" % (n, checked, bad))
        for k in (1, 2, 3, 4):
            keys = sorted(table[k])
            print("   increment 2k=%d: (h, #odd cycles through v, complement counts) -> cases: %s" % (2 * k, [(key, table[k][key]) for key in keys][:10]))
        print("   minimal positive increment k per h (h <= 40):", {h: minpos[h] for h in sorted(minpos) if h <= 40})
