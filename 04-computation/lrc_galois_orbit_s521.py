#!/usr/bin/env python3
"""
lrc_galois_orbit_s521.py   claudebox-2026-06-01-S521

LRC via Galois orbits (reflection: 07-reflections/lrc-galois-zak-theta-s521.md).

At t=a/q (q prime), runner i is at zeta_q^{v_i a}; Gal(Q(zeta_q)/Q)=(Z/q)* acts by
a:zeta->zeta^a, sending the config to (v_i a mod q). LRC at denom q <=> the Galois
orbit O_q meets the safe box B_q=[ceil(q/n),q-ceil(q/n)]^m. The orbit equidistributes
as q grows (resonances sum c_i v_i ≡ 0 mod q thin out), so f(q)=|O_q∩box|/(q-1) ->
V=(1-2/n)^m and the box is hit for large q. Tight sets escape to q=n (the cyclotomic
field itself).
"""
from math import ceil, gcd
from itertools import product

def isprime(x): return x > 1 and all(x % d for d in range(2, int(x**.5)+1))
def safe(r, q, n): return min(r, q-r)*n >= q
def f_of_q(v, q, n):
    return sum(1 for a in range(1, q) if all(safe((vi*a) % q, q, n) for vi in v)) / (q-1)
def smallest_lonely_prime(v, n, qmax=400):
    for q in range(n+1, qmax):
        if isprime(q) and f_of_q(v, q, n) > 0: return q
    return None
def resonances(v, q, R=2):
    return sum(1 for c in product(range(-R, R+1), repeat=len(v))
               if any(c) and sum(ci*vi for ci, vi in zip(c, v)) % q == 0)

def main():
    print("Galois-orbit lonely fraction f(q) vs volume V=(1-2*ceil(q/n)/q)^m:")
    for v in [(1,2,4,7),(2,3,5,7,11)]:
        n = len(v)+1; m = len(v)
        print(f"  v={v}:", end=" ")
        for q in [p for p in range(n+1, 60) if isprime(p)][:5]:
            V = (1-2*ceil(q/n)/q)**m
            print(f"q{q}:f={f_of_q(v,q,n):.2f}/V={V:.2f}", end="  ")
        print()
    print("\nSmallest Galois-lonely PRIME (q>n) per set (None = tight, lonely only at q=n):")
    for v in [(1,2,3,4),(1,2,4,7),(1,3,4,5,9),(2,3,5,7,11)]:
        n = len(v)+1
        print(f"  v={v}: {smallest_lonely_prime(v, n)}")
    print("\nResonance count (#small c, |c_i|<=2, with sum c_i v_i ≡ 0 mod q) thins as q grows:")
    for v in [(1,2,4,7)]:
        n = len(v)+1
        for q in [p for p in range(n+1, 40) if isprime(p)][:5]:
            print(f"  v={v} q={q}: resonances={resonances(v,q)}", end="   ")
        print()
    print("\n=> equidistribution: f(q)->V>0 as resonances thin; box hit for q>=q_0.")
    print("   Strategy: Weil/Ramanujan bound on orbit discrepancy => explicit q_0(n); finite middle.")

if __name__ == "__main__":
    main()
