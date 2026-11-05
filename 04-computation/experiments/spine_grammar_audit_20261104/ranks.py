#!/usr/bin/env python3
"""ranks.py -- auditor's own exact linear algebra on the automaton dumps (format: S/T/A lines, state 0 = empty cut).
Computes, with Fractions: reachable space (row vectors e_0 M_u), observable space (column functionals M_r P f),
the Hankel rank of the pairing, for TWO alphabets:
  index alphabet : a letter is the index li, acting on every state s through lettersFor(c(s))[li] (the engine's encoding)
  triple alphabet: a letter is the triple (T,b,A), acting on the states s where it is a letter of lettersFor(c(s))
and the seam ranks (pairing restricted to reachable vectors supported on cut states with c pending slots).
Usage: python3 -I ranks.py dump.txt
"""
import sys
from fractions import Fraction
from collections import defaultdict

def parse(path):
    states = {}; trans = {}; triple = {}; apex = {}
    for line in open(path):
        p = line.split()
        if not p: continue
        if p[0] == 'S':
            sid = int(p[1]); c = int(p[2][2:]); comp = int(p[3].split('=')[1]); states[sid] = (c, comp)
        elif p[0] == 'T':
            sid = int(p[1]); li = int(p[2]); T, b, A = int(p[3]), int(p[4]), int(p[5])
            lst = [(int(x.split('*')[0]), int(x.split('*')[1])) for x in p[7:]]
            trans[(sid, li)] = lst; triple[(sid, li)] = (T, b, A)
        elif p[0] == 'A':
            sid = int(p[1]); a = int(p[2]); lst = [(int(x.split('*')[0]), int(x.split('*')[1])) for x in p[4:]]
            apex[(sid, a)] = lst
    return states, trans, triple, apex

MOD = None   # None = exact rationals; otherwise a prime p and all arithmetic is in GF(p) (rank mod p <= rational rank)

class Echelon:
    def __init__(self): self.rows = []   # list of (pivot, dict) ; pivots distinct
    def reduce(self, v):
        if MOD is None:
            v = {k: Fraction(x) for k, x in v.items() if x != 0}
            for piv, row in self.rows:
                if piv in v:
                    f = v[piv] / row[piv]
                    for k2, val in row.items():
                        nv = v.get(k2, Fraction(0)) - f * val
                        if nv == 0: v.pop(k2, None)
                        else: v[k2] = nv
            return v
        p = MOD
        v = {k: x % p for k, x in v.items() if x % p}
        for piv, row in self.rows:
            if piv in v:
                f = (v[piv] * pow(row[piv], p - 2, p)) % p
                for k2, val in row.items():
                    nv = (v.get(k2, 0) - f * val) % p
                    if nv == 0: v.pop(k2, None)
                    else: v[k2] = nv
        return v
    def add(self, v):
        v = self.reduce(v)
        if not v: return False
        piv = min(v); self.rows.append((piv, v)); return True
    def dim(self): return len(self.rows)

def rank_of_pairing(rows, cols):
    E = Echelon()
    for r in rows:
        vec = {}
        for j, cvec in enumerate(cols):
            s = sum(r[k] * cvec.get(k, 0) for k in r)
            if s: vec[j] = s
        E.add(vec)
    return E.dim()

EMPTY1 = False

def main():
    global MOD, EMPTY1
    args = [a for a in sys.argv[1:] if a != '--empty1']; EMPTY1 = ('--empty1' in sys.argv); sys.argv = [sys.argv[0]] + args
    if len(sys.argv) > 2: MOD = int(sys.argv[2])
    states, trans, triple, apex = parse(sys.argv[1])
    K = 1 + max([a for (_, a) in apex] + [-1])
    cs = sorted(set(c for c, _ in states.values()))
    letters_by_c = defaultdict(set); triples_by_c = defaultdict(set)
    for (s, li), tr in triple.items():
        letters_by_c[states[s][0]].add(li); triples_by_c[states[s][0]].add(tr)
    idx_of = {}   # (state, triple) -> li
    for (s, li), tr in triple.items(): idx_of[(s, tr)] = li
    all_li = sorted(set(li for (_, li) in trans)); all_tr = sorted(set(triple.values()))
    # --- reachable space (same for both alphabets): BFS from e_0, letters valid at the vector's c
    def row_apply(vec, li):
        out = defaultdict(int)
        for s, val in vec.items():
            for (t, m) in trans.get((s, li), []): out[t] += val * m
        return {k: v for k, v in out.items() if v}
    reach = Echelon(); e0 = {0: 1}; reach.add(e0); frontier = [e0]; reach_vecs = [e0]
    while frontier:
        nxt = []
        for vec in frontier:
            c = states[next(iter(vec))][0]
            assert all(states[s][0] == c for s in vec)
            for li in sorted(letters_by_c[c]):
                w = row_apply(vec, li)
                if w and reach.add(w): nxt.append(w); reach_vecs.append(w)
        frontier = nxt
    # --- final functional f: complete indicator after apex steps
    f = {s: 1 for s, (c, comp) in states.items() if comp}
    if EMPTY1 and K == 0: f[0] = 1   # convention H(empty word) = 1: the empty cut (no vertices) is accepting
    for a in reversed(range(K)):
        g = {}
        for s in states:
            tot = sum(m * f.get(t, 0) for (t, m) in apex.get((s, a), []))
            if tot: g[s] = tot
        f = g
    def col_apply_index(func, li):
        out = {}
        for s in states:
            tot = sum(m * func.get(t, 0) for (t, m) in trans.get((s, li), []))
            if tot: out[s] = tot
        return out
    def col_apply_triple(func, tr):
        out = {}
        for s in states:
            li = idx_of.get((s, tr))
            if li is None: continue
            tot = sum(m * func.get(t, 0) for (t, m) in trans.get((s, li), []))
            if tot: out[s] = tot
        return out
    def observable(apply, alphabet):
        obs = Echelon(); obs.add(f); frontier = [f]; vecs = [f]
        while frontier:
            nxt = []
            for func in frontier:
                for a in alphabet:
                    g = apply(func, a)
                    if g and obs.add(g): nxt.append(g); vecs.append(g)
            frontier = nxt
        return obs, vecs
    obs_i, O_i = observable(col_apply_index, all_li)
    obs_t, O_t = observable(col_apply_triple, all_tr)
    R = [dict(r) for _, r in reach.rows]
    Oi = [dict(r) for _, r in obs_i.rows]; Ot = [dict(r) for _, r in obs_t.rows]
    hi = rank_of_pairing(R, Oi); ht = rank_of_pairing(R, Ot)
    print(f"{sys.argv[1].split('/')[-1]}: {len(states)} states, K={K}; reachable dim {reach.dim()}; observable dim index-alphabet {obs_i.dim()}, triple-alphabet {obs_t.dim()}; Hankel rank index-alphabet {hi}, triple-alphabet {ht}")
    for c in cs:
        rows_c = [r for r in R if states[next(iter(r))][0] == c]
        nstates = sum(1 for s in states if states[s][0] == c)
        ri = rank_of_pairing(rows_c, Oi); rt = rank_of_pairing(rows_c, Ot)
        print(f"  seam c={c}: {nstates} cut states, reachable dim {len(rows_c)}, seam rank index-alphabet {ri}, triple-alphabet {rt}")

main()
