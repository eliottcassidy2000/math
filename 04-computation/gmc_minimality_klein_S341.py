#!/usr/bin/env python3
"""
klein-2026-07-20-S341 -- IS THE GMC COUNTEREXAMPLE MINIMAL?  Bounded exhaustive searches
in (a) 2 real Gaussians, (b) fewer monomials, (c) lower degree -- each with a POSITIVE
CONTROL, because a null with no positive control is worthless (MISTAKE-162).
"""
from itertools import combinations, product as iproduct

def fact(n):
    r = 1
    for i in range(2, n + 1): r *= i
    return r
def pmul(A, B):
    C = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = tuple(x + y for x, y in zip(ma, mb)); C[m] = C.get(m, 0) + ca * cb
    return {m: c for m, c in C.items() if c}
def ppow(A, k, K):
    R = {(0,) * (2 * K): 1}
    for _ in range(k): R = pmul(R, A)
    return R
def E(A):
    tot = 0
    for m, c in A.items():
        ok, val = True, 1
        for k in range(0, len(m), 2):
            a, b = m[k], m[k + 1]
            if a != b: ok = False; break
            val *= fact(a)
        if ok: tot += c * val
    return tot

def monos(K, D):
    """all exponent tuples over K complex Gaussians with total degree <= D"""
    out = []
    def go(pos, cur, deg):
        if pos == 2 * K: out.append(tuple(cur)); return
        for e in range(D - deg + 1):
            cur.append(e); go(pos + 1, cur, deg + e); cur.pop()
    go(0, [], 0)
    return out

def is_counterexample(P, K, Mmax=8, Qdeg=2):
    """E[P^m]=0 for m=1..Mmax, and some small Q has E[Q P^m] != 0 at several large m"""
    pw = {}
    cur = {(0,) * (2 * K): 1}
    for m in range(1, Mmax + 1):
        cur = pmul(cur, P); pw[m] = cur
        if E(cur) != 0: return None
    for q in monos(K, Qdeg):
        if sum(q) == 0: continue
        Q = {q: 1}
        nz = [m for m in range(1, Mmax + 1) if E(pmul(Q, pw[m])) != 0]
        # demand nonvanishing at the top few m, not just small m
        if len([m for m in nz if m >= Mmax - 2]) == 3: return (q, nz)
    return None

# ------------------------------------------------------------------ the known one
K = 2
A_ = (1, 0, 0, 0); Ab = (0, 1, 0, 0); B_ = (0, 0, 1, 0); Bb = (0, 0, 0, 1)
def add(*ts):
    d = {}
    for t, c in ts: d[t] = d.get(t, 0) + c
    return {m: c for m, c in d.items() if c}
KNOWN = add((B_, 1), ((1, 1, 0, 0), -1), ((0, 0, 1, 1), 1), ((1, 1, 0, 1), -1))  # (1+Bb)(B-|A|^2)
print("=" * 80)
print("POSITIVE CONTROL -- the search predicate must FIRE on the known 4-term example")
print("=" * 80)
r = is_counterexample(KNOWN, 2)
print(f"  P' = (1+Bb)(B - |A|^2), 4 monomials, cubic -> predicate fires: {r is not None}")
if r: print(f"    witness Q exponent {r[0]} (= Bb), nonvanishing at m = {r[1]}")
assert r is not None, "POSITIVE CONTROL FAILED -- the search would be meaningless"

# ------------------------------------------------------------------ (a) n = 2 real
print("\n" + "=" * 80)
print("(a) 2 REAL GAUSSIANS (one complex A): exhaustive over small P")
print("=" * 80)
M1 = [m for m in monos(1, 4)]
print(f"  monomial pool: {len(M1)} (degree <= 4 in A, Ab);  P with <= 4 terms, coeffs in +-1,+-2")
hits = []
for size in (2, 3, 4):
    for combo in combinations(M1, size):
        for coeffs in iproduct((1, -1, 2, -2), repeat=size):
            P = {}
            for mm, c in zip(combo, coeffs): P[mm] = P.get(mm, 0) + c
            P = {m: c for m, c in P.items() if c}
            if not P: continue
            res = is_counterexample(P, 1)
            if res: hits.append((P, res))
    print(f"   size {size}: cumulative hits = {len(hits)}")
print(f"  RESULT: {len(hits)} counterexamples found in 2 real Gaussians"
      f"{'' if hits else '  -- NONE'}")
if hits: print(f"   example: {hits[0]}")

# ------------------------------------------------------------------ (b) fewer terms
print("\n" + "=" * 80)
print("(b) 4 REAL GAUSSIANS, FEWER THAN 4 MONOMIALS")
print("=" * 80)
M2 = [m for m in monos(2, 3)]
print(f"  monomial pool: {len(M2)} (degree <= 3 in A,Ab,B,Bb)")
for size in (2, 3):
    found = []
    for combo in combinations(M2, size):
        for coeffs in iproduct((1, -1, 2, -2), repeat=size):
            P = {}
            for mm, c in zip(combo, coeffs): P[mm] = P.get(mm, 0) + c
            P = {m: c for m, c in P.items() if c}
            if not P: continue
            res = is_counterexample(P, 2)
            if res: found.append((P, res)); break
        if found: break
    print(f"   {size} monomials: {'FOUND ' + str(found[0][0]) if found else 'none'}")

# ------------------------------------------------------------------ (c) lower degree
print("\n" + "=" * 80)
print("(c) 4 REAL GAUSSIANS, DEGREE <= 2 (is cubic necessary?)")
print("=" * 80)
M2q = [m for m in monos(2, 2)]
found = []
for size in (2, 3, 4):
    for combo in combinations(M2q, size):
        for coeffs in iproduct((1, -1, 2, -2), repeat=size):
            P = {}
            for mm, c in zip(combo, coeffs): P[mm] = P.get(mm, 0) + c
            P = {m: c for m, c in P.items() if c}
            if not P: continue
            res = is_counterexample(P, 2)
            if res: found.append((P, res)); break
        if found: break
    print(f"   degree <= 2, {size} monomials: {'FOUND ' + str(found[0][0]) if found else 'none'}")
    if found: break
print("\nDONE.")
