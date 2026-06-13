#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-469 lab: the sum-free grading dichotomy (why the Erdos-592 seam is 2-adic).
kind-pasteur-2026-06-11-S1.  Answers THM-464 D's sharp open note / THM-465 C caveat.

Parts:
  A. Range-verification of the hand lemma: v_p level sets are sum-free iff p=2
     (parity closure); leading-digit refinements p^v(u + pZ) are sum-free for ALL p.
  B. Diagonal census in the (sign, v_p) triangle quotient at n=3: count realizable
     mono triples (f,f,f) and total realizable triples for p=2 vs p=3.
     Prediction: p=2 mono count = 0 at every t (parity closure), p=3 positive from t=3.
  C. UNSAT-core extraction for the (sign,v3) feature instance at (3,4): WHICH forced
     classes clash through WHICH compositions (the concrete algebraic mechanism).
  D. HYP-2390 branching experiment: b=3 subgrid game at n=3 with full feature
     algebras F2=(sign,v2) vs F3v=(sign,v3).  Schur-arity explanation predicts the
     v2/v3 asymmetry PERSISTS at b=3 (it follows triangle 2-sums, not branching).
  E. HYP-2391 leading-digit rescue: F3LD=(sign, v3, leading 3-adic digit) at the
     binary game (3,4),(3,5): prediction SAT where bare F3v is feature-UNSAT.
     Plus (sign,v5) control at (3,6) (odd p, same escape: predicted feature-UNSAT).

Mistake guards: MISTAKE-067 (complete verifier only), feature-UNSAT semantics
(algebra too coarse, NOT problem-UNSAT), MISTAKE-062 (nothing hardcoded).
"""

import itertools, time, sys
from pysat.solvers import Glucose3

from erdos592_satverifier_frontier_macmini_s2 import leaves, SubgridVerifier, seed_subgrids
from erdos592_bidyadic_rule_macmini_s3 import v2, F2, independent_verify
from erdos592_ternary_seam_macmini_s3 import BVerifier, v3


# ---------------------------------------------------------------- feature maps

def vp(g, p):
    k = 0
    while g % p == 0:
        g //= p; k += 1
    return k


def sgn_vp(d, p):
    if d == 0:
        return ('=',)
    return ('+' if d > 0 else '-', vp(abs(d), p))


def F3v(x, y):
    return tuple(sgn_vp(y[i] - x[i], 3) for i in range(len(x)))


def F5v(x, y):
    return tuple(sgn_vp(y[i] - x[i], 5) for i in range(len(x)))


def sgn_vp_lead(d, p):
    if d == 0:
        return ('=',)
    g = abs(d)
    v = vp(g, p)
    lead = (g // p ** v) % p
    return ('+' if d > 0 else '-', v, lead)


def F3LD(x, y):
    return tuple(sgn_vp_lead(y[i] - x[i], 3) for i in range(len(x)))


# ------------------------------------------------- A. the lemma, range-checked

def part_A(limit=4096):
    print("=== A. sum-freeness of valuation classes, range [1,%d] ===" % limit, flush=True)
    for p in (2, 3, 5):
        viol = None
        for x in range(1, limit):
            if viol: break
            for y in range(x, limit - x):
                if vp(x, p) == vp(y, p) == vp(x + y, p):
                    viol = (x, y, x + y); break
        print(f"   bare v_{p} level sets sum-free: "
              + ("YES (no violation)" if viol is None else f"NO  e.g. {viol}"), flush=True)
    for p in (3, 5):
        viol = None
        for x in range(1, limit):
            if viol: break
            vx, lx = vp(x, p), 0
            lx = (x // p ** vx) % p
            for y in range(x, limit - x):
                vy = vp(y, p)
                if vy != vx: continue
                ly = (y // p ** vy) % p
                if ly != lx: continue
                s = x + y
                vs = vp(s, p)
                if vs == vx and (s // p ** vs) % p == lx:
                    viol = (x, y, s); break
        print(f"   v_{p}+leading-digit classes sum-free: "
              + ("YES (no violation)" if viol is None else f"NO  e.g. {viol}"), flush=True)
    print("   (hand proof: L_v sum-free <=> u1+u2 = 0 mod p for ALL units <=> |units|=1 <=> p=2;"
          "\n    refined classes: u+u = 2u != u mod p whenever u != 0 — sum-free for every p)", flush=True)


# ---------------------------------------- B. diagonal census in the quotient

def part_B(ts=(3, 4, 5), n=3):
    print("\n=== B. realizable triangle feature-triples, (sign,v_p) quotient, n=%d ===" % n, flush=True)
    for p in (2, 3):
        Fp = (lambda x, y: tuple(sgn_vp(y[i] - x[i], p) for i in range(len(x))))
        for t in ts:
            L = leaves(n, t); N = len(L)
            triples = set()
            mono = set()
            for i, j, k in itertools.combinations(range(N), 3):
                f1, f2, f3 = Fp(L[i], L[j]), Fp(L[j], L[k]), Fp(L[i], L[k])
                key = tuple(sorted((f1, f2, f3)))
                triples.add(key)
                if f1 == f2 == f3:
                    mono.add(f1)
            print(f"   p={p} t={t}: {len(triples):5d} realizable triples, "
                  f"{len(mono):3d} mono classes (f,f,f)"
                  + ("" if not mono else "   e.g. " + str(sorted(mono)[0])), flush=True)


# --------------------------------------------------- C. v3 UNSAT core at (3,4)

def build_feature_instance(n, t, F):
    """Return (clauses, qvar, kinds): triangle + seed clauses over feature vars."""
    L = leaves(n, t); N = len(L)
    qv = {}; cnt = [0]

    def q(x, y):
        if x > y:
            x, y = y, x
        key = F(x, y)
        if key not in qv:
            cnt[0] += 1; qv[key] = cnt[0]
        return qv[key]

    clauses = []; kinds = []
    seen = set()
    for i, j, k in itertools.combinations(range(N), 3):
        c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
        if c not in seen:
            seen.add(c)
            clauses.append(list(c)); kinds.append(('tri', (L[i], L[j], L[k])))
    for g in seed_subgrids(n, t):
        c = tuple(sorted(set(q(a, b) for a, b in itertools.combinations(g, 2))))
        if c not in seen:
            seen.add(c)
            clauses.append(list(c)); kinds.append(('hit', g))
    return clauses, qv, kinds


def part_C(n=3, t=4):
    print(f"\n=== C. UNSAT core of the (sign,v3) instance at ({n},{t}) ===", flush=True)
    clauses, qv, kinds = build_feature_instance(n, t, F3v)
    nv = max(qv.values())
    sol = Glucose3()
    sel = []
    for idx, c in enumerate(clauses):
        s = nv + 1 + idx
        sel.append(s)
        sol.add_clause(c + [-s])
    assert not sol.solve(assumptions=sel), "expected UNSAT"
    core = set(sol.get_core())
    active = [i for i, s in enumerate(sel) if s in core]
    # greedy minimization
    changed = True
    while changed:
        changed = False
        for i in list(active):
            trial = [sel[j] for j in active if j != i]
            if not sol.solve(assumptions=trial):
                active.remove(i); changed = True
    inv = {v: k for k, v in qv.items()}
    print(f"   minimized core: {len(active)} clauses "
          f"(of {len(clauses)} total)", flush=True)
    forced_on = []
    for i in active:
        kind, obj = kinds[i]
        lits = clauses[i]
        if kind == 'hit':
            feats = sorted(set(inv[abs(l)] for l in lits))
            print(f"   [hit ] subgrid root {obj[0]}..{obj[-1]} -> one of {len(feats)} classes:", flush=True)
            for f in feats:
                print(f"          {f}", flush=True)
        else:
            feats = [inv[abs(l)] for l in lits]
            print(f"   [tri ] {obj} forbids together:", flush=True)
            for f in feats:
                print(f"          {f}", flush=True)
    return active


# ----------------------------------- shared featured CEGAR solver (b-ary game)

def solve_featured_b(n, t, b, F, name, tlimit=900, verbose=True):
    """THM-465 feature-quotient CEGAR for the b-ary subgrid game (BVerifier).
    Returns (True, (edges, L, model, qvar)) | (False, None) feature-UNSAT | (None, None)."""
    L = leaves(n, t); N = len(L)
    qv = {}; cnt = [0]

    def q(x, y):
        if x > y:
            x, y = y, x
        key = F(x, y)
        if key not in qv:
            cnt[0] += 1; qv[key] = cnt[0]
        return qv[key]

    sol = Glucose3()
    seen = set()
    for i, j, k in itertools.combinations(range(N), 3):
        c = tuple(sorted(set((-q(L[i], L[j]), -q(L[i], L[k]), -q(L[j], L[k])))))
        if c not in seen:
            seen.add(c); sol.add_clause(list(c))
    ver = BVerifier(n, t, b)
    t0 = time.time(); added = 0
    while True:
        if time.time() - t0 > tlimit:
            if verbose: print(f"   TIMEOUT {name} b={b} ({n},{t}) lazy={added}", flush=True)
            return None, None
        if not sol.solve():
            if verbose: print(f"   feature-UNSAT {name} b={b} ({n},{t}) "
                              f"({cnt[0]} features, lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return False, None
        model = set(l for l in sol.get_model() if l > 0)
        edges = set()
        for i in range(N):
            for j in range(i + 1, N):
                if q(L[i], L[j]) in model:
                    edges.add((i, j))
        bad = ver.find(edges)
        if bad is None:
            if verbose: print(f"   SAT {name} b={b} ({n},{t}) ({len(edges)} edges, "
                              f"{cnt[0]} features, lazy={added}, {time.time()-t0:.1f}s)", flush=True)
            return True, (edges, L, model, qv)
        cl = sorted(set(q(a, c) for a, c in itertools.combinations(bad, 2)))
        sol.add_clause(cl)
        added += 1


def verify_b(n, t, b, edges):
    """Independent re-verification for the b-ary game (fresh verifier + triangle scan)."""
    L = leaves(n, t); N = len(L)
    adj = [[False] * N for _ in range(N)]
    for i, j in edges:
        adj[i][j] = adj[j][i] = True
    for i, j, k in itertools.combinations(range(N), 3):
        assert not (adj[i][j] and adj[i][k] and adj[j][k]), "triangle!"
    return BVerifier(n, t, b).find(set(edges)) is None


# ------------------------------------------------- D. the branching experiment

def part_D(ts=(4, 5, 6), tlim=1200):
    print("\n=== D. HYP-2390: b=3 game at n=3 — does the v2/v3 asymmetry persist? ===", flush=True)
    print("   (Schur-arity prediction: YES — triangles are 2-sums regardless of b)", flush=True)
    results = {}
    for name, F in (("F2 ", F2), ("F3v", F3v)):
        for t in ts:
            r, w = solve_featured_b(3, t, 3, F, name, tlimit=tlim)
            results[(name.strip(), t)] = ('SAT' if r else ('timeout' if r is None else 'feature-UNSAT'))
            if r:
                ok = verify_b(3, t, 3, w[0])
                print(f"      independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)
            if r is False:
                break  # antitone in t within this algebra: coarser data at larger t? NOT
                       # guaranteed for feature algebras -- but an UNSAT here already
                       # answers the asymmetry question at this t; stop to save time.
    print("   D summary:", results, flush=True)
    return results


# --------------------------------------------- E. the leading-digit rescue

def part_E(tlim=1200):
    print("\n=== E. HYP-2391: leading-digit rescue (binary game, n=3) ===", flush=True)
    from erdos592_bidyadic_rule_macmini_s3 import solve_featured
    results = {}
    print("   bare (sign,v3) controls (expected feature-UNSAT, re-run of THM-465 C):", flush=True)
    for t in (4,):
        r, w = solve_featured(3, t, F3v, "F3v")
        results[('F3v', t)] = 'SAT' if r else ('timeout' if r is None else 'feature-UNSAT')
    print("   (sign,v3,leading-digit) — sum-free classes for p=3:", flush=True)
    for t in (4, 5):
        r, w = solve_featured(3, t, F3LD, "F3LD", tlimit=tlim)
        results[('F3LD', t)] = 'SAT' if r else ('timeout' if r is None else 'feature-UNSAT')
        if r:
            edges, L, model, qv = w
            ok = independent_verify(3, t, edges)
            print(f"      independent re-verification: {'PASS' if ok else 'FAIL'}", flush=True)
    print("   (sign,v5) control at (3,6) (odd p, bare valuation — expected feature-UNSAT):", flush=True)
    r, w = solve_featured(3, 6, F5v, "F5v", tlimit=tlim)
    results[('F5v', 6)] = 'SAT' if r else ('timeout' if r is None else 'feature-UNSAT')
    print("   E summary:", results, flush=True)
    return results


def main():
    t0 = time.time()
    part_A()
    part_B()
    part_C()
    d = part_D()
    e = part_E()
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)
    print("FINAL: D(branching)=", d, flush=True)
    print("FINAL: E(rescue)  =", e, flush=True)


if __name__ == "__main__":
    main()
