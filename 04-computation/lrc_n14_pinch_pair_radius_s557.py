#!/usr/bin/env python3
"""
Exact verification of the Pinch Lemma / radius=r/s structure for LRC.
opus-2026-06-02-S557 (remote-control). Deductive backbone, not a search.

Claim (proved in the reflection): for M(S)=max_t min_i||v_i t|| < 1/2, the max is
attained at t*=m/(v_a+v_b) where a,b straddle the observer (frac(v_a t*)=M,
frac(v_b t*)=1-M); and M(S)=r/s with s=(v_a+v_b)/gcd(v_a,v_b), r>=1.
Consequences: counterexample => binding reduced-sum s>=15; tight (M=1/14) =>
s ≡ 0 (mod 14).
"""
from fractions import Fraction
from math import gcd
import random


def M_and_argmax(V):
    c = set()
    for v in V:
        for k in range(2 * v):
            c.add(Fraction(2 * k + 1, 2 * v) % 1)
    for i in range(len(V)):
        for j in range(len(V)):
            for s in (1, -1):
                d = V[i] + s * V[j]
                if d:
                    for k in range(abs(d) + 1):
                        c.add(Fraction(k, d) % 1)
    best, bt = Fraction(-1), None
    for t in c:
        mn = min(min((v * t) % 1, 1 - (v * t) % 1) for v in V)
        if mn > best:
            best, bt = mn, t
    return best, bt


def binding_pair(V, M, t):
    plus = [v for v in V if (v * t) % 1 == M]
    minus = [v for v in V if (v * t) % 1 == 1 - M]
    return plus, minus


def check(n_random=200, hi=40, seed=11):
    rng = random.Random(seed)
    ok = 0
    tot = 0
    for _ in range(n_random):
        V = tuple(sorted(rng.sample(range(1, hi + 1), 13)))
        M, t = M_and_argmax(V)
        plus, minus = binding_pair(V, M, t)
        tot += 1
        if M == Fraction(1, 2):
            ok += 1
            continue
        if not (plus and minus):
            print("  NO STRADDLE:", V, M, t)
            continue
        a, b = plus[0], minus[0]
        s = (a + b) // gcd(a, b)
        r = M * s
        cond = (r.denominator == 1) and ((a + b) * t == int((a + b) * t)) \
            and ((a - b) * t % 1 in (2 * M, 1 - 2 * M))
        if cond:
            ok += 1
        else:
            print("  FAIL:", V, M, t, (a, b), "s=", s, "r=", r)
    print(f"Pinch/r-over-s verified: {ok}/{tot} random 13-sets (speeds<= {hi})")


def critical():
    print("Critical configs — straddling pair and reduced sum s:")
    cases = [("AP@14", tuple(range(1, 14)), 14),
             ("V*@14", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24), 14),
             ("sporadic@5", (1, 3, 4, 7), 5),
             ("sporadic@6", (1, 3, 4, 5, 9), 6),
             ("sporadic@8", (1, 2, 3, 4, 5, 7, 12), 8)]
    for name, V, n in cases:
        thr = Fraction(1, n)
        # tight witnesses
        T = set()
        for v in V:
            for k in range(v + 1):
                for sgn in (-1, 1):
                    T.add(Fraction(n * k + sgn, n * v) % 1)
        wit = [t for t in T if min(min((v * t) % 1, 1 - (v * t) % 1) for v in V) == thr]
        for t in sorted(wit):
            plus, minus = binding_pair(V, thr, t)
            if plus and minus:
                a, b = plus[0], minus[0]
                s = (a + b) // gcd(a, b)
                print(f"  {name} t*={t}: straddle (+{a},-{b}) sum={a+b} reduced s={s} "
                      f"(s mod {n}={s % n}); (a-b)t*={(a-b)*t%1}=2/n? "
                      f"{(a-b)*t%1 in (2*thr,1-2*thr)}")


if __name__ == "__main__":
    check()
    print()
    critical()
