#!/usr/bin/env python3
"""
lrc14_discrete_quintic_bonferroni_klein_S210.py

HYP-5758: the DISCRETE QUINTIC BONFERRONI CERTIFICATE (B5-cert) — the port of
THM-604's depth-5 truncation (and THM-660/661's moment-LP role) to Z_q kill
sets: the crux item of HYP-5732 (the aggregated modular supply route).

SETUP. S a covering 13-set, q > Vmax a modulus. Kill sets
    B_l = {p in Z_q : v_l·p mod q is OUTSIDE the middle band [c, q-c]},
    c = ceil(q/14)   (residue r safe iff 14r >= q and 14(q-r) >= q).
Coverage count C(p) = #{l-classes killing p} (classes: v_l ≡ ±v_k mod q merged,
their kill sets coincide). LM(q) = #{p != 0 : C(p) = 0} (live multipliers; a
live (q,p) is a DIRECT witness: t = p/q, M(S) >= 1/14; Lean consumer
LRCPairSumDispatch.mreach_ge_of_pairsum_band).

THE CERTIFICATE (THM-604's truncation identity, pointwise combinatorial):
    for integers c >= 0 and odd D:   1_{c=0} >= sum_{d<=D} (-1)^d C(c,d)
so, summing over p != 0,
    LM(q) >= B_D(S,q) := sum_{p!=0} f_D(C(p)),   f_D(c) = sum_{d<=D} (-1)^d C(c,d).
B_D(S,q) is an EXACT integer computed in one O(q) histogram pass; B_D > 0 is a
decidable loneliness certificate. iid heuristic: with 13 classes each of
density ~1/7, E[f5] -> B5(13) = 2052/16807 = 0.1221 > 0 while the cubic
B3(13) = -0.0991 < 0 — exactly why depth-1/2 ledgers (C1, C4-Hunter) broke
adversarially at klein-S209 while the truth stayed fat.

TESTS:
 [A] The S209 battlefield: does B5 certify where C1∪C4 died (the adversarial
     0-cert instances), on the @91 hard cluster, and on the mid-band census?
     Per instance: #q in (V,2V] with B5>0, the best B5, exact LM comparison,
     and the depth ladder B1/B3/B5/B7 at the best q.
 [B] SCALE: the B5-certified fraction of moduli vs V (V = 30..300); the
     empirical V0 where B5 starts certifying (the explicit-threshold shape:
     [B5-cert for V > V0] + [native_decide below V0, kps-S115 has [1,18]]).
 [C] STRUCTURE PENALTY: planted Schur triples / near-interval sets — how does
     max_q B5 degrade with E3(S)? (the LEM-015/rigidity dichotomy dial).
 [D] Histogram forensics: where B5 < 0 < LM, what does the coverage histogram
     look like vs Binomial(13, |K|/q) (which resonances inflate high coverage)?
"""

import random
from math import gcd, comb

random.seed(20260710)
QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


def coverage_hist_and_LM(S, q):
    """One pass: coverage counter over p in Z_q for the ±-merged classes.
    Returns (hist, LM, nclasses) with hist[c] = #{p != 0 : C(p) = c}."""
    def r_safe(r):
        return 14 * r >= q and 14 * (q - r) >= q
    cov = bytearray(q)
    seen = set()
    ncl = 0
    for v in S:
        r = v % q
        key = min(r, (q - r) % q)
        if key in seen:
            continue
        seen.add(key)
        ncl += 1
        if r == 0:
            for p in range(q):
                cov[p] += 1
            continue
        g = gcd(r, q)
        rr, qq = r // g, q // g
        inv = pow(rr, -1, qq)
        for m in range(qq):
            s = m * g
            if not r_safe(s):
                p0 = (m * inv) % qq
                for t in range(g):
                    cov[p0 + t * qq] += 1
    hist = [0] * (ncl + 1)
    for p in range(1, q):
        hist[cov[p]] += 1
    LM = hist[0]
    return hist, LM, ncl


def BD_from_hist(hist, D):
    """B_D = sum_p f_D(C(p)) from the coverage histogram (exact integer)."""
    tot = 0
    for c, n in enumerate(hist):
        if n == 0:
            continue
        fD = sum((-1) ** d * comb(c, d) for d in range(0, D + 1))
        tot += n * fD
    return tot


def gen_instance(V, style='slowheavy'):
    P = random.choice([(8, 9, 10, 12), (7, 9, 10, 11, 12), (11, 12, 13),
                       (10, 11, 12, 13), (9, 11, 13)])
    k = 13 - len(P)
    if k < 8:
        return None
    L = {V}
    missed = [q for q in QS if not any(p % q == 0 for p in P)]
    for q in missed:
        if any(u % q == 0 for u in L):
            continue
        lo, hi = -(-14 // q), V // q
        if lo > hi:
            return None
        L.add(q * random.randint(lo, hi))
    if style == 'slowheavy':
        for _ in range(3):
            if len(L) < k:
                L.add(random.randint(max(14, V // 14 + 1), max(16, 9 * V // 14 - 1)))
    while len(L) < k:
        L.add(random.randint(14, V))
    S = sorted(set(P) | L)
    if len(S) == 13 and is_covering(S):
        return S
    return None


def instance_report(S, name, qstep=1, verbose=True):
    V = max(S)
    qs = range(V + 1, 2 * V + 1, qstep)
    n_b5 = n_b3 = n_b1 = n_live = n_q = 0
    best = None
    for q in qs:
        hist, LM, ncl = coverage_hist_and_LM(S, q)
        b5 = BD_from_hist(hist, 5)
        b3 = BD_from_hist(hist, 3)
        b1 = BD_from_hist(hist, 1)
        n_q += 1
        n_live += 1 if LM > 0 else 0
        n_b5 += 1 if b5 > 0 else 0
        n_b3 += 1 if b3 > 0 else 0
        n_b1 += 1 if b1 > 0 else 0
        if best is None or b5 > best[1]:
            best = (q, b5, LM, ncl)
    if verbose:
        q, b5, LM, ncl = best
        hist, _, _ = coverage_hist_and_LM(S, q)
        b1, b3, b7 = BD_from_hist(hist, 1), BD_from_hist(hist, 3), BD_from_hist(hist, 7)
        print(f"  {name}: V={V}")
        print(f"    moduli q in (V,2V] (step {qstep}): {n_q}; live {n_live}; "
              f"B1-cert {n_b1}; B3-cert {n_b3}; B5-cert {n_b5}  "
              f"({100*n_b5/max(1,n_q):.0f}%)")
        print(f"    best q={q} (classes {ncl}): depth ladder B1={b1} B3={b3} "
              f"B5={b5} B7={b7} vs exact LM={LM}  [B5/q = {b5/q:.4f}]")
    return n_q, n_b5, n_live, best


def main():
    print("=" * 78)
    print("HYP-5758: the discrete quintic Bonferroni certificate (klein-S210)")
    print("=" * 78)
    print(f"\niid reference: B1(13)={1-13/7:.4f}  B3(13)={sum((-1)**d*comb(13,d)/7**d for d in range(4)):.4f}"
          f"  B5(13)={sum((-1)**d*comb(13,d)/7**d for d in range(6)):.4f}"
          f"  (6/7)^13={(6/7)**13:.4f}")

    # [A] the S209 battlefield
    print("\n[A] The S209 battlefield (C1∪C4 died here; exact truth was fat):")
    battlefield = [
        ([12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120], "adversarial-worst V=120 (0 certs)"),
        ([31, 33, 45, 48, 73, 76, 82, 86, 98, 102, 103, 104, 120], "adversarial V=120 (0 certs)"),
        ([62, 66, 69, 102, 109, 118, 120, 126, 130, 136, 159, 185, 200], "adversarial V=200 (0 certs)"),
        ([9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91], "7-structured @91"),
        ([9, 12, 64, 110, 123, 127, 146, 149, 155, 169, 182, 185, 200], "adversarial V=200 #2"),
        ([10, 46, 48, 59, 66, 148, 177, 181, 208, 213, 236, 261, 280], "adversarial V=280"),
    ]
    for S, name in battlefield:
        instance_report(S, name)

    # [B] scale: B5-certified fraction vs V
    print("\n[B] Scale sweep (random mid-band covering instances):")
    print(f"{'V':>5} {'inst':>5} {'avg %B5-cert':>13} {'min %B5':>8} {'min bestB5/q':>13} {'all inst >=1 B5?':>17}")
    for V in (30, 45, 60, 90, 120, 180, 260):
        rows = []
        tries = 0
        while len(rows) < 6 and tries < 60:
            tries += 1
            S = gen_instance(V)
            if S is None:
                continue
            step = 1 if V <= 120 else 2
            n_q, n_b5, n_live, best = instance_report(S, "", qstep=step, verbose=False)
            rows.append((n_q, n_b5, best[1] / best[0]))
        if not rows:
            continue
        fr = [b / a for a, b, _ in rows]
        allpos = all(b >= 1 for _, b, _ in rows)
        print(f"{V:>5} {len(rows):>5} {100*sum(fr)/len(fr):>12.1f}% {100*min(fr):>7.1f}% "
              f"{min(r[2] for r in rows):>13.4f} {str(allpos):>17}")

    # [C] structure penalty: E3 dial
    print("\n[C] Structure penalty (E3 = #Schur triples a+b=c in S):")
    def E3(S):
        st = set(S)
        return sum(1 for i, a in enumerate(S) for b in S[i:] if a + b in st)
    # base random covering vs planted-Schur variants at V=150
    tested = 0
    while tested < 6:
        S = gen_instance(150)
        if S is None:
            continue
        e3 = E3(S)
        n_q, n_b5, n_live, best = instance_report(S, "", verbose=False)
        print(f"   E3={e3:>3}  %B5-cert={100*n_b5/n_q:>5.1f}%  best B5/q={best[1]/best[0]:.4f}  S={S}")
        tested += 1
    # near-interval (high E3, covering): dilate {1..13} won't be covering; use {1..12, 14}-style
    for S in ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14],
              [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26],
              [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14]):
        if len(set(S)) == 13 and is_covering(S):
            e3 = E3(S)
            n_q, n_b5, n_live, best = instance_report(S, "", verbose=False)
            print(f"   E3={e3:>3}  %B5-cert={100*n_b5/n_q:>5.1f}%  best B5/q={best[1]/best[0]:.4f}  S={S} (structured)")

    # [D] forensics on one failure (if any): histogram vs binomial
    print("\n[D] Histogram forensics (worst instance/modulus with LM>0 but B5<=0, if any):")
    S = [9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91]
    V = max(S)
    shown = 0
    for q in range(V + 1, 2 * V + 1):
        hist, LM, ncl = coverage_hist_and_LM(S, q)
        b5 = BD_from_hist(hist, 5)
        if LM > 0 and b5 <= 0 and shown < 2:
            kdens = sum(c * n for c, n in enumerate(hist)) / ((q - 1) * ncl)
            from math import comb as CB
            bino = [(q - 1) * CB(ncl, c) * kdens ** c * (1 - kdens) ** (ncl - c) for c in range(ncl + 1)]
            print(f"   @91 q={q}: LM={LM} but B5={b5} (classes {ncl}, kill density {kdens:.3f})")
            print(f"      hist  = {hist}")
            print(f"      binom = {[round(x,1) for x in bino]}")
            shown += 1
    if shown == 0:
        print("   (none found on @91 — B5 tracks LM positivity everywhere there)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
