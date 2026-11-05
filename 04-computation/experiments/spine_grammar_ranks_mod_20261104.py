#!/usr/bin/env python3
"""spine_grammar_ranks_mod_20261104.py -- opus-2026-11-04-S79

Modular (GF(p), p = 2^31 - 1 and p = 2^31 - 19) version of spine_grammar_ranks_20261104.py for large automata: the
reachable space (row vectors e_0 M_u), the observable space (column vectors M_r f) and the rank of their pairing are
computed with dense echelon bases over GF(p).  A rank over GF(p) is at most the rational rank; two primes agreeing is
the reported value (the exact script is used where it runs).
Usage: python3 spine_grammar_ranks_mod_20261104.py dump.txt
"""
import sys
import numpy as np
from collections import defaultdict

def parse(path):
    states = {}; trans = defaultdict(list); apex = defaultdict(list); letters = {}; nstates = 0
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts: continue
            if parts[0] == 'STATES': nstates = int(parts[1])
            elif parts[0] == 'S':
                sid = int(parts[1]); c = int(parts[2][2:]); comp = int(parts[3].split('=')[1]); states[sid] = (c, comp)
            elif parts[0] == 'T':
                sid = int(parts[1]); li = int(parts[2]); T, b, A = int(parts[3]), int(parts[4]), int(parts[5])
                lst = [(int(x.split('*')[0]), int(x.split('*')[1])) for x in parts[7:]]
                trans[(sid, li)] = lst; letters[(states[sid][0], li)] = (T, b, A)
            elif parts[0] == 'A':
                sid = int(parts[1]); a = int(parts[2]); lst = [(int(x.split('*')[0]), int(x.split('*')[1])) for x in parts[4:]]
                apex[(sid, a)] = lst
    return nstates, states, trans, apex, letters

class Echelon:
    """dense row-echelon basis over GF(p)"""
    def __init__(self, n, p): self.n = n; self.p = p; self.rows = []; self.pivots = []
    def add(self, v):
        v = v % self.p
        for row, piv in zip(self.rows, self.pivots):
            if v[piv]:
                v = (v - v[piv] * row) % self.p
        nz = np.nonzero(v)[0]
        if len(nz) == 0: return False
        piv = int(nz[0]); inv = pow(int(v[piv]), self.p - 2, self.p); v = (v * inv) % self.p
        self.rows.append(v); self.pivots.append(piv); return True

def main():
    path = sys.argv[1]
    nstates, states, trans, apex, letters = parse(path)
    K = 1 + max([a for (_, a) in apex] + [-1]); n = max(states) + 1
    letters_by_c = defaultdict(list)
    for (c, li) in letters: letters_by_c[c].append(li)
    # forward (row) transitions by state and letter; column (transpose) by target
    out_rows = defaultdict(list)   # (li) -> list of (s, t, m)
    for (s, li), lst in trans.items():
        for (t, m) in lst: out_rows[li].append((s, t, m))
    results = {}
    for p in (2147483647, 2147483629):
        # reachable space
        reach = Echelon(n, p); e0 = np.zeros(n, dtype=np.int64); e0[0] = 1; reach.add(e0.copy())
        frontier = [e0]; reach_by_c = defaultdict(list); reach_by_c[0].append(e0)
        while frontier:
            nxt = []
            for v in frontier:
                nzs = np.nonzero(v)[0]; c = states[int(nzs[0])][0]
                for li in letters_by_c[c]:
                    w = np.zeros(n, dtype=np.int64)
                    for s in nzs:
                        for (t, m) in trans.get((int(s), li), []): w[t] = (w[t] + int(v[s]) * m) % p
                    if not w.any(): continue
                    if reach.add(w.copy()):
                        nxt.append(w); reach_by_c[states[int(np.nonzero(w)[0][0])][0]].append(w)
            frontier = nxt
        # observable space: functionals; f = apex placement then complete indicator
        f = np.zeros(n, dtype=np.int64)
        for s, (c, cp) in states.items():
            if cp: f[s] = 1
        if '--empty1' in sys.argv: f[0] = 1    # H(empty word) = 1 (the empty tournament has one path): the empty cut is accepting
        for a in reversed(range(K)):
            g = np.zeros(n, dtype=np.int64)
            for s in states:
                tot = 0
                for (t, m) in apex.get((s, a), []): tot += m * int(f[t])
                g[s] = tot % p
            f = g
        obs = Echelon(n, p); obs.add(f.copy()); frontier = [f]
        all_letters = sorted(set(li for (c, li) in letters))
        while frontier:
            nxt = []
            for func in frontier:
                for li in all_letters:
                    g = np.zeros(n, dtype=np.int64)
                    for (s, t, m) in out_rows[li]:
                        if func[t]: g[s] = (g[s] + m * int(func[t])) % p
                    if not g.any(): continue
                    if obs.add(g.copy()): nxt.append(g)
            frontier = nxt
        R = np.array(reach.rows, dtype=np.int64); O = np.array(obs.rows, dtype=np.int64)
        def pairing_rank(R, O):
            # R O^T mod p without int64 overflow: split entries (< 2^31) into 16-bit halves
            R = R % p; O = O % p
            R1, R0 = R >> 16, R & 0xFFFF; O1, O0 = O >> 16, O & 0xFFFF
            A = (R1 @ O1.T) % p; B = ((R1 @ O0.T) + (R0 @ O1.T)) % p; C = (R0 @ O0.T) % p
            M = ((A * pow(2, 32, p)) % p + (B * pow(2, 16, p)) % p + C) % p
            E = Echelon(O.shape[0], p); r = 0
            for i in range(M.shape[0]):
                if E.add(M[i].copy()): r += 1
            return r
        hankel = pairing_rank(R, O)
        seams = {}
        cs = sorted(set(c for (c, cp) in states.values()))
        for c in cs:
            Ec = Echelon(n, p); rows = []
            for v in reach_by_c[c]:
                if Ec.add(v.copy()): rows.append(v)
            states_c = np.array([s for s, (cc, cp) in states.items() if cc == c])
            mask = np.zeros(n, dtype=np.int64); mask[states_c] = 1
            Oc = (O * mask) % p; Oc = Oc[[i for i in range(Oc.shape[0]) if Oc[i].any()]]
            seams[c] = (len(states_c), len(rows), pairing_rank(np.array(rows, dtype=np.int64), Oc) if rows and len(Oc) else 0)
        results[p] = (len(reach.rows), len(obs.rows), hankel, seams)
    (r1, o1, h1, s1), (r2, o2, h2, s2) = results[2147483647], results[2147483629]
    agree = (r1, o1, h1, s1) == (r2, o2, h2, s2)
    print(f"states {nstates}, K = {K}; reachable dim {r1}, observable dim {o1}, Hankel rank = {h1} (two primes agree: {agree})")
    for c, (nc, rd, sr) in sorted(s1.items()):
        print(f"  seam with c = {c} pending slots: {nc} cut states, reachable dim {rd}, seam rank = {sr}")

main()
