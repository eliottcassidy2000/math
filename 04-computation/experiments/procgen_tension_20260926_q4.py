#!/usr/bin/env python3
"""
procgen_tension_20260926_q4.py -- Q4: the negation duality nu on the strategy cube (q = 3).
u/d coordinates: U(sigma) = {r odd : sigma(r) != chi_{-4}(r)}; nu acts by U -> -U; self-dual <=> U = -U.
Classification of class (i) at levels 2..5 by nu-orbits, rho_max, components; self-dual census at level 6.
Run through the runner.
"""
from collections import Counter, deque
from fractions import Fraction

from procgen_tension_20260926_lib import (
    below_c, chi4_mask, check, claim, karp_batch, nu, rho_all, sigma_k_mask, bad_set, uset)


def nu_perm(k):
    """bit permutation of nu on masks: bit i (residue 2i+1) goes to the bit of -(2i+1); plus a global flip"""
    M = 1 << k
    return [((M - (2 * i + 1)) - 1) // 2 for i in range(1 << (k - 1))]


def nu_fast(mask, perm, full):
    out = 0
    for i, j in enumerate(perm):
        if (mask >> i) & 1:
            out |= 1 << j
    return out ^ full


def run():
    print("Q4. the negation duality on the strategy cube")
    summary = {}
    for k in (2, 3, 4, 5):
        M = 1 << k
        H = 1 << (k - 1)
        full = (1 << H) - 1
        perm = nu_perm(k)
        R = rho_all(k)
        NM = len(R)
        chi = chi4_mask(k)
        # nu via the fast bit permutation agrees with the definition, and preserves rho_max (Theorem E)
        for m in (range(NM) if k <= 4 else range(0, NM, 97)):
            check(nu_fast(m, perm, full) == nu(m, k), "nu_fast = nu")
        numap = [nu_fast(m, perm, full) for m in range(NM)]
        check(all(R[numap[m]] == R[m] for m in range(NM)), "rho_max is nu-invariant, k=%d" % k)
        # u/d lemma: U(nu sigma) = -U(sigma); self-dual <=> U symmetric
        for m in (range(NM) if k <= 4 else range(0, NM, 31)):
            U = set(uset(m, k))
            check(set(uset(numap[m], k)) == {M - r for r in U}, "U(nu sigma) = -U(sigma)")
            check((numap[m] == m) == (U == {M - r for r in U}), "self-dual <=> U = -U")
        C1 = [m for m in range(NM) if below_c(R[m])]
        C1s = set(C1)
        sd = [m for m in C1 if numap[m] == m]
        pairs = [m for m in C1 if numap[m] > m]
        check(len(sd) + 2 * len(pairs) == len(C1), "orbit count")
        # every self-dual strategy is at Hamming distance exactly 2^(k-2) from Collatz (mask 0) and 3n-1
        allsd = [m for m in range(NM) if numap[m] == m]
        check(len(allsd) == 2 ** (2 ** (k - 2)), "number of self-dual strategies")
        check(all(bin(m).count('1') == 2 ** (k - 2) and bin(m ^ full).count('1') == 2 ** (k - 2) for m in allsd),
              "self-dual strategies sit at distance 2^(k-2) from both constant strategies")
        # rho_max distributions
        dsd = Counter(R[m] for m in sd)
        dpr = Counter(R[m] for m in pairs)
        # components of class (i) under single sign flips
        comp = {}
        ncomp = 0
        for m in C1:
            if m in comp:
                continue
            q = deque([m])
            comp[m] = ncomp
            while q:
                x = q.popleft()
                for i in range(H):
                    y = x ^ (1 << i)
                    if y in C1s and y not in comp:
                        comp[y] = ncomp
                        q.append(y)
            ncomp += 1
        sizes = Counter(comp.values())
        comp_sd = Counter(comp[m] for m in sd)
        comps_with_sd = len(comp_sd)
        comps_nu_fixed = sum(1 for c in range(ncomp)
                             if comp[numap[next(m for m in C1 if comp[m] == c)]] == c)
        # nearest self-dual class-(i) strategy
        dist = Counter()
        for m in C1:
            dist[min(bin(m ^ t).count('1') for t in sd)] += 1
        # symmetric core: U cap -U (keep only the symmetric part of the u-set)
        core_ok = 0
        for m in C1:
            U = set(uset(m, k))
            core = {r for r in U if (M - r) in U}
            cm = chi
            for r in core:
                cm ^= 1 << ((r - 1) // 2)
            check(numap[cm] == cm, "the symmetric core is self-dual")
            core_ok += cm in C1s
        # is class (i) closed under turning one u into d?  (refuted: count the violating moves)
        moves = viol = 0
        for m in C1:
            for r in uset(m, k):
                moves += 1
                viol += (m ^ (1 << ((r - 1) // 2))) not in C1s
        summary[k] = (len(C1), len(sd), len(pairs), dsd, dpr, ncomp, sorted(sizes.values(), reverse=True)[:6],
                      comps_with_sd, comps_nu_fixed, dist, core_ok, moves, viol)
        print("    level %d: class (i) %d = %d self-dual + %d pairs; rho_max self-dual %s; paired %s"
              % (k, len(C1), len(sd), len(pairs), sorted(dsd.items()), sorted(dpr.items())))
        print("             components under single flips: %d (largest %s), %d contain a self-dual, %d are nu-invariant;"
              " distance to the nearest self-dual class-(i) strategy: %s; symmetric core class (i): %d/%d"
              % (ncomp, sorted(sizes.values(), reverse=True)[:6], comps_with_sd, comps_nu_fixed,
                 sorted(dist.items()), core_ok, len(C1)))
    check(summary[2][:3] == (1, 1, 0) and summary[3][:3] == (1, 1, 0) and summary[4][:3] == (16, 2, 7)
          and summary[5][:3] == (1052, 12, 520), "class-(i) nu-orbit census")
    claim(True, "class (i) = self-dual + pairs: 1 = 1 + 0 (k=2), 1 = 1 + 0 (k=3), 16 = 2 + 2*7 (k=4), "
          "1052 = 12 + 2*520 (k=5); rho_max is nu-invariant on all 65,812 strategies of levels 2..5")
    check(all(summary[k][5] == 1 for k in (2, 3, 4, 5)), "class (i) is flip-connected at levels 2..5")
    claim(True, "class (i) is a single component of the single-flip graph at each level 2..5 (so it contains chi_(-4) and "
          "every self-dual class-(i) strategy); at level 5 every class-(i) strategy is within %d flips of a self-dual "
          "class-(i) one" % max(summary[5][9]))
    check((summary[4][11], summary[4][12], summary[5][11], summary[5][12]) == (30, 4, 4794, 380)
          and (summary[5][10], summary[4][10]) == (972, 16), "refuted guesses: u->d closure, symmetric core")
    claim(True, "REFUTED guesses: turning one u into d leaves class (i) in 4 of 30 moves (k=4) and 380 of 4794 (k=5); "
          "the symmetric core U cap -U is class (i) for 16/16 (k=4) but only 972/1052 (k=5)")
    claim(True, "u/d coordinates: U(nu sigma) = -U(sigma) and sigma is self-dual iff its u-set is symmetric (checked on "
          "all strategies k <= 4, a 1/31 sample at k = 5); every one of the 2^(2^(k-2)) self-dual strategies is at "
          "Hamming distance exactly 2^(k-2) (Haar 1/2) from Collatz and from 3n-1")
    # sigma_k: its nearest self-dual class-(i) strategy is chi_{-4}, at distance 2^(k-2) - |Bad_k| (k = 2..5 exact,
    # the formula for all k is proved in the note)
    for k in (2, 3, 4, 5):
        H = 1 << (k - 1)
        R = rho_all(k)
        perm = nu_perm(k)
        full = (1 << H) - 1
        sd = [m for m in range(len(R)) if below_c(R[m]) and nu_fast(m, perm, full) == m]
        sk = sigma_k_mask(k)
        dmin = min(bin(sk ^ t).count('1') for t in sd)
        check(dmin == 2 ** (k - 2) - len(bad_set(k)) == bin(sk ^ chi4_mask(k)).count('1'),
              "nearest self-dual class-(i) to sigma_k is chi_{-4}")
    claim(True, "k = 2..5: the nearest self-dual class-(i) strategy to sigma_k is chi_{-4} (all-d), at distance "
          "2^(k-2) - |Bad_k|")
    # the three integer expanding cycles of 3n+1 on the negatives, and their nu-images (3n-1 on the positives):
    # every class-(i) strategy has sigma = - at an odd point of each negative cycle and sigma = + at an odd point of
    # each positive one (else the integer cycle's itinerary is an expanding closed walk of G_sigma)
    neg_cycles = [[-1], [-5, -7], [-17, -25, -37, -55, -41, -61, -91]]      # odd points
    for cyc in neg_cycles:
        pts = [x for x in cyc]
        y = pts[0]
        # confirm they are 3n+1 cycles: iterate T with sign + from the first odd point
        seen = [y]
        for _ in range(40):
            y = y // 2 if y % 2 == 0 else (3 * y + 1) // 2
            if y == pts[0]:
                break
            seen.append(y)
        check(y == pts[0] and sorted(x for x in seen if x % 2) == sorted(pts), "integer 3n+1 cycle")
    for k in (2, 3, 4, 5):
        M = 1 << k
        R = rho_all(k)
        for m in range(len(R)):
            if not below_c(R[m]):
                continue
            for cyc in neg_cycles:
                check(any((m >> (((x % M) - 1) // 2)) & 1 for x in cyc), "class (i) breaks a negative cycle")
                check(any(not ((m >> (((-x) % M - 1) // 2)) & 1) for x in cyc), "class (i) breaks a positive cycle")
    claim(True, "every class-(i) strategy of levels 2..5 puts a minus sign on an odd point of each of the 3n+1 cycles "
          "{-1}, {-5,-7,-10}, {-17,...,-34} and a plus sign on an odd point of each 3n-1 cycle {1}, {5,7,10}, {17,...,34}")
    # level 6: all 2^16 self-dual strategies, exact Karp
    k = 6
    M = 1 << k
    H = 1 << (k - 1)
    reps = [r for r in range(1, M, 4)]          # residues 1 mod 4; their negatives are 3 mod 4
    masks = []
    for bits in range(1 << len(reps)):
        m = 0
        for j, r in enumerate(reps):
            s = -1 if (bits >> j) & 1 else 1     # sigma(r)
            if s == -1:
                m |= 1 << ((r - 1) // 2)
            if -s == -1:                         # sigma(-r) = -sigma(r)
                m |= 1 << (((M - r) - 1) // 2)
        masks.append(m)
    perm = nu_perm(k)
    full = (1 << H) - 1
    for m in masks[:: 257]:
        check(nu(m, k) == m, "level-6 self-dual enumeration")
    cnt = Counter()
    ncls = 0
    for st in range(0, len(masks), 4096):
        rn, rd = karp_batch(k, masks[st:st + 4096])
        for a, b in zip(rn, rd):
            f = Fraction(int(a), int(b))
            if below_c(f):
                ncls += 1
                cnt[f] += 1
    print("    level 6: %d of the 65536 self-dual strategies are class (i); rho_max distribution %s"
          % (ncls, sorted(cnt.items())))
    claim(ncls > 0, "level 6 (exhaustive over the 2^16 self-dual strategies, exact Karp): %d self-dual class-(i) "
          "strategies; their largest rho_max is %s" % (ncls, max(cnt)))
    return summary


if __name__ == "__main__":
    run()
