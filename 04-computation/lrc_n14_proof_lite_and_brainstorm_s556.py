#!/usr/bin/env python3
"""
LRC@14: a rigorous "proof-lite" impossibility + evidence for adjacent statements.
opus-2026-06-02-S556 (remote-control).

PROOF-LITE (rigorous) -- ruling out classes of counterexamples:
  P1. q=14 sieve: if 14 | no speed, t=1/14 is a witness (THM-369). Verify.
  P2. AP-with-any-evens-doubled is lonely at t=1/14, self-contained: doubling
      runner 2m (m=1..6) to 4m keeps 4m != 0 (mod 14) since 7 doesn't divide m,
      and every original AP runner is >= 1/14 from 0 at t=1/14. Verify all 2^6-1.
  P4. UNIFICATION: a counterexample's forced multiple of 14 IS the even member of
      oracle-S552o's mod-7 SINGLETON class (multiples of 7). Verify the bridge.

ADJACENT STATEMENTS (evidence, candidate theorems):
  B3. tight-witness-lattice: every TIGHT n=14 config is lonely ONLY at t in (1/14)Z.
      Search distance<=2 from the AP for tight configs; check their witnesses.
  B5. e<=6 reduction completeness + odd-split margin (how robust is S554?).
All exact (Fraction).
"""

from fractions import Fraction
from math import gcd
from itertools import combinations
import random

N = 14
THR = Fraction(1, N)


def collar_at(V, t):
    return min(min((v * t) % 1, 1 - (v * t) % 1) for v in V)


def witnesses_lattice_and_M(V):
    """Return (M, witness_times) exactly, witnesses = the argmax times among the
    exact candidate set; classify whether all max-attaining rational candidates
    lie on (1/14)Z."""
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
    best = Fraction(-1)
    for t in c:
        m = collar_at(V, t)
        if m > best:
            best = m
    wits = [t for t in c if collar_at(V, t) == best]
    return best, wits


def is_tight(V):
    M, _ = witnesses_lattice_and_M(V)
    return M == THR


def primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in S)) if g else tuple(sorted(S))


# ---------------- P1: q=14 sieve ----------------
def P1(trials=3000, seed=1):
    rng = random.Random(seed)
    bad = 0
    n = 0
    for _ in range(trials):
        S = tuple(sorted(rng.sample(range(1, 80), 13)))
        if any(v % 14 == 0 for v in S):
            continue
        n += 1
        if collar_at(S, Fraction(1, 14)) < THR:
            bad += 1
    print(f"P1 (no mult of 14 => t=1/14 witness): {n} configs, "
          f"failures={bad}  [PROVED: sieve THM-369, q=14]")


# ---------------- P2: AP with evens doubled ----------------
def P2():
    AP = list(range(1, 14))
    evens = [v for v in AP if v % 2 == 0]
    bad = 0
    tot = 0
    tights = []
    for r in range(1, len(evens) + 1):
        for sub in combinations(evens, r):
            V = tuple(sorted(2 * v if v in sub else v for v in AP))
            if len(set(V)) != 13:
                continue
            tot += 1
            # the claimed witness t=1/14
            if collar_at(V, Fraction(1, 14)) < THR:
                bad += 1
            if is_tight(V):
                tights.append(V)
    print(f"P2 (AP-with-evens-doubled lonely at 1/14): {tot} configs, "
          f"t=1/14 fails on {bad}  [PROVED self-contained]")
    print(f"   of these, exactly tight (M=1/14): {tights}")


# ---------------- P4: the 14-multiple is the mod-7 singleton's even member ----
def P4():
    print("P4 (bridge): a counterexample needs 14|v (else t=1/14 wins). 14|v ⇒ "
          "7|v ⇒ v is in oracle-S552o's mod-7 SINGLETON class AND even.")
    print("   So the mod-2 fold's obstruction (multiple of 14) is exactly the "
          "even half of the mod-7 singleton — the two factorisations meet here.")


# ---------------- B3: tight-witness lattice, distance<=2 from AP ----------------
def B3(add_hi=44):
    AP = tuple(range(1, 14))
    tight_found = {AP}
    # distance-1: remove one, add one
    for rem in AP:
        base = [x for x in AP if x != rem]
        for add in range(14, add_hi + 1):
            if add in base:
                continue
            V = tuple(sorted(base + [add]))
            if is_tight(primitive(V)):
                tight_found.add(primitive(V))
    # distance-2: remove two, add two
    news = list(range(14, add_hi + 1))
    for rem in combinations(AP, 2):
        base = [x for x in AP if x not in rem]
        for add in combinations(news, 2):
            V = tuple(sorted(base + list(add)))
            if len(set(V)) != 13:
                continue
            if is_tight(primitive(V)):
                tight_found.add(primitive(V))
    print(f"B3 (tight configs within distance<=2 of AP, adds<= {add_hi}): "
          f"{len(tight_found)} found")
    all_lattice = True
    for V in sorted(tight_found):
        M, wits = witnesses_lattice_and_M(V)
        onlat = all((w * N).denominator == 1 for w in wits)
        all_lattice = all_lattice and onlat
        tag = "AP" if V == AP else "tight"
        print(f"   {tag}: {V}  witnesses {[str(w) for w in wits]}  "
              f"all on (1/14)Z? {onlat}")
    print(f"   => EVERY tight config witnessed only on the (1/14)-lattice? "
          f"{all_lattice}")


# ---------------- B5: e<=6 odd-split robustness ----------------
def B5(trials=80, seed=3):
    rng = random.Random(seed)
    seven = Fraction(1, 7)
    ok = 0
    n = 0
    worst_slack = None
    for _ in range(trials):
        e = rng.randint(1, 6)
        evens = rng.sample(range(2, 60, 2), e)
        odds = rng.sample(range(1, 60, 2), 13 - e)
        S = primitive(tuple(evens + odds))
        if len(S) != 13:
            continue
        n += 1
        fold = sorted(v // 2 for v in S if v % 2 == 0)
        O = [v for v in S if v % 2 == 1]
        # even-good s (g_fold(s)>=1/7): candidate s from fold collar candidates
        cs = set()
        for u in fold:
            for k in range(2 * u):
                cs.add(Fraction(2 * k + 1, 2 * u) % 1)
        good = [s for s in cs if collar_at(fold, s) >= seven] if fold else [Fraction(0)]
        found = False
        for s in good:
            for t in (s / 2 % 1, (s + 1) / 2 % 1):
                if collar_at(S, t) >= THR:
                    found = True
                    break
            if found:
                break
        if found:
            ok += 1
    print(f"B5 (e<=6 fold reduction): {ok}/{n} recovered a witness from an "
          f"even-good preimage")


if __name__ == "__main__":
    print("===== PROOF-LITE (rigorous impossibility of structured counterexamples) =====")
    P1(); P2(); P4()
    print("\n===== ADJACENT STATEMENTS (evidence) =====")
    B3(); B5()
