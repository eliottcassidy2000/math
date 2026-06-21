#!/usr/bin/env python3
# Wide rigorous verification of D_{p,q} <= 14/p over ALL coprime p/q in (1, 43/20], q<=Q.
# Also confirm sup(D*p)=20/7 and sup(D*q)=12/7, and the apex law D=0 iff 7|pq.
from fractions import Fraction as Fr
from math import gcd

def sec(y):
    f = y % 1
    s = int(7*f)
    return s if s < 7 else 6

def mu_row_counts(p, q):
    # exact cell counts on the common 7pq grid; c[i][j] integer, sum=7pq? Actually use
    # the integer-grid form: on [0,1), partition by both sec(qv) and sec(pv). Use breakpoints.
    bps = set()
    for c in (p, q):
        for t in range(0, 7*c+1):
            bps.add(Fr(t, 7*c))
    bps = sorted(bps)
    from collections import defaultdict
    cell = defaultdict(lambda: Fr(0))
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        mid = (a+b)/2
        cell[(sec(q*mid), sec(p*mid))] += b-a
    return cell

def D(p, q):
    cell = mu_row_counts(p, q)
    tot = Fr(0)
    for i in range(7):
        for j in range(7):
            tot += abs(cell.get((i,j), Fr(0)) - Fr(1,49))
    return tot

Q = 80
viol = 0
checked = 0
supDp = Fr(0); supDp_at = None
supDq = Fr(0); supDq_at = None
apex_fail = 0
for q in range(1, Q+1):
    for p in range(q+1, (Fr(43,20)*q).__floor__()+1):
        if gcd(p, q) != 1: continue
        if not (Fr(1,1) < Fr(p,q) <= Fr(43,20)): continue
        checked += 1
        d = D(p, q)
        if d > Fr(14, p):
            viol += 1
            if viol <= 5: print(f"  VIOLATION p/q={p}/{q}: D={float(d):.5f} > 14/p={float(Fr(14,p)):.5f}")
        if d*p > supDp: supDp = d*p; supDp_at = (p,q)
        if d*q > supDq: supDq = d*q; supDq_at = (p,q)
        # apex law: D==0 iff 7|pq
        is_apex = (d == 0)
        should = (7 % 1 == 0) and (p*q) % 7 == 0
        if is_apex != should:
            apex_fail += 1
            if apex_fail <= 5: print(f"  APEX-LAW FAIL p/q={p}/{q}: D==0? {is_apex} but 7|pq? {should}")

print(f"checked {checked} coprime p/q in (1,43/20], q<={Q}")
print(f"D<=14/p violations: {viol}  ({'PROOF HOLDS' if viol==0 else 'PROOF BROKEN'})")
print(f"sup D*p = {supDp} = {float(supDp):.5f} at {supDp_at}  (proof bound 14; tight sup claimed 20/7={float(Fr(20,7)):.5f})")
print(f"sup D*q = {supDq} = {float(supDq):.5f} at {supDq_at}  (claimed 12/7={float(Fr(12,7)):.5f})")
print(f"apex law (D=0 iff 7|pq) failures: {apex_fail}")
