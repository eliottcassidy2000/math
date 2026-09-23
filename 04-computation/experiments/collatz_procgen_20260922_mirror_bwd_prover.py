#!/usr/bin/env python3
"""collatz_procgen_20260922_mirror_bwd_prover.py -- exact hostility test for positive dyadic 3-adic points
of the backward E-game (Q2), used by the duality lane.

Reverse move x -> (2^k x - 1)/3, legal iff 2^k x = 1 mod 3 and the result is a 3-adic unit (for a
rational: numerator and denominator prime to 3; the value 0 is illegal).  A legal path of s moves with
K doublings has multiplier R = 2^K/3^s; y0 is hostile (in Bad_inf) iff R > 1 for every prefix of every
legal path.  For y0 = N/2^e write every path state as x_s = (2^K y0 - B_s)/3^s, so
      R_s = (x_s + beta'_s)/y0,  beta'_s = B_s/3^s >= (1 - 3^-s)/2.
FINITE CRITERION (PROVED in the note, sec. 5): let y0 < 3/2 and S* = max{S : 3^-S > 3 - 2 y0}.
  * integer states z >= 2 have R >= (2 + 1/3)/y0 > 1; the integer state 1 at depth S > S* has R > 1;
    positive integers are closed under legal moves (Lemma S_b);
  * non-integer states have K < e, so they form a finite tree (each k = 0 move divides R by 3);
  * an integer state z at depth tau < S* can reach 1 by depth S* only if z <= c_(S*-tau), c_r = (3^(r+1)-1)/2.
  So y0 < 3/2 is hostile iff an exhaustive search of (non-integer tree) + (integer states that can still
  reach 1 by depth S*) finds no prefix with R <= 1.  The search is complete and terminates.
SLOW-DESCENT SEARCH (for y0 > 3/2, where the criterion does not apply): first move k1 (all doublings
  first), then the greedy map G(m) = (2^k m - 1)/3 (k minimal legal) to 1; report a descent (R < 1 at some
  prefix) if one is found for k1 <= K1MAX.  A found descent is an explicit certificate (the class is not
  hostile); failure proves nothing.
Usage: python3 ..._mirror_bwd_prover.py  p/2^e ...  (or --census FILE)
"""
import sys
from fractions import Fraction as F

def v3(n):
    n = abs(n); c = 0
    while n and n % 3 == 0: n //= 3; c += 1
    return c

def is_unit3(x):
    return x != 0 and x.numerator % 3 != 0 and x.denominator % 3 != 0

def moves(x, kmax):
    """legal reverse moves (k, child) with k <= kmax"""
    out = []
    for k in range(0, kmax + 1):
        t = x * 2 ** k - 1
        if t.numerator % 3 != 0: continue         # 2^k x = 1 mod 3 (denominator is a power of 2)
        y = t / 3
        if is_unit3(y): out.append((k, y))
    return out

def c(r): return (3 ** (r + 1) - 1) // 2

def hostile_finite(y0, state_cap=2_500_000):
    """exact decision for a positive dyadic y0 < 3/2.  Depth-first over all legal paths, memoising every
    state with the least multiplier R = 2^K/3^s seen so far (a state reached again with R' >= R has a
    subtree whose multipliers are all >= those already checked, so it is not re-explored).
    Positive integer states are handled by the finite criterion (depth >= S* or value > c(S*-depth): safe).
    Returns ('hostile', #states) or ('descends', path) or ('undecided', reason)."""
    assert 0 < y0 < F(3, 2)
    Sstar = 0
    while F(1, 3 ** (Sstar + 1)) > 3 - 2 * y0: Sstar += 1
    best = {}                                   # state -> (K, s) with minimal 2^K/3^s
    stack = [(y0, 0, 0, ())]
    while stack:
        x, s, K, path = stack.pop()
        if s > 0 and 2 ** K <= 3 ** s: return ('descends', list(path))
        if x < 0:
            # negative: all descendants negative; follow least legal k until R <= 1
            xx, ss, KK, pp = x, s, K, list(path)
            for _ in range(20000):
                for k in (0, 1, 2, 3):
                    t = xx * 2 ** k - 1
                    if t.numerator % 3 == 0 and is_unit3(t / 3): break
                xx = (xx * 2 ** k - 1) / 3; ss += 1; KK += k; pp.append(k)
                if 2 ** KK <= 3 ** ss: return ('descends', pp)
            return ('undecided', 'negative branch without greedy descent')
        if x.denominator == 1:
            z = x.numerator
            if s >= Sstar or z > c(Sstar - s): continue   # safe (finite criterion)
            kmax = 1
            while (2 ** kmax * z - 1) // 3 <= c(max(Sstar - s - 1, 0)) + 1: kmax += 1
            kids = moves(x, kmax)
        else:
            d = x.denominator.bit_length() - 1
            # non-integer children: k < d; integer children: k >= d, only those that can still reach 1 by S*
            kmax = d
            while True:
                zc = (x * 2 ** kmax - 1) / 3
                if zc.denominator == 1 and (zc.numerator > c(max(Sstar - s - 1, 0)) + 1 or s + 1 >= Sstar): break
                kmax += 1
            kids = moves(x, kmax)
        for k, ch in kids:
            K2, s2 = K + k, s + 1
            if ch.denominator == 1 and ch > 0 and (s2 >= Sstar or ch.numerator > c(Sstar - s2)):
                if 2 ** K2 <= 3 ** s2: return ('descends', list(path) + [k])
                continue
            prev = best.get(ch)
            if prev is not None and 2 ** prev[0] * 3 ** s2 <= 2 ** K2 * 3 ** prev[1]: continue
            best[ch] = (K2, s2)
            if len(best) > state_cap: return ('undecided', f'state cap, S*={Sstar}')
            stack.append((ch, s2, K2, path + (k,)))
    return ('hostile', len(best))

def greedy_k(m):
    """minimal legal k at a positive integer m (3 not dividing m): 2^k m in {4,7} mod 9"""
    for k in range(0, 7):
        if (2 ** k * m) % 9 in (4, 7): return k
    raise ValueError

def slow_descent(y0, K1MAX=400, SMAX=20000):
    """first move k1 then greedy G to 1; returns the first (k1, s, K) with 2^K < 3^s, else None"""
    best = None
    for k1 in range(1, K1MAX + 1):
        t = y0 * 2 ** k1 - 1
        if t.numerator % 3: continue
        x = t / 3
        if not is_unit3(x) or x.denominator != 1: continue
        m = x.numerator; s = 1; K = k1
        if 2 ** K < 3 ** s: return (k1, s, K, [k1])
        ks = [k1]
        while m != 1 and s < SMAX:
            k = greedy_k(m); m = (2 ** k * m - 1) // 3; s += 1; K += k; ks.append(k)
            if 2 ** K < 3 ** s: return (k1, s, K, ks)
        # track the closest approach R_final
        r = F(2 ** K, 3 ** s)
        if best is None or r < best[0]: best = (r, k1, s)
    return ('none', float(best[0]) if best else None, best[1] if best else None)

def check_path(y0, ks):
    """independent replay: legality and multiplier"""
    x = y0; K = 0
    for s, k in enumerate(ks, 1):
        t = x * 2 ** k - 1
        if t.numerator % 3: return False
        x = t / 3; K += k
        if not is_unit3(x): return False
    return 2 ** K < 3 ** len(ks)

if __name__ == '__main__':
    args = sys.argv[1:]
    pts = []
    if args and args[0] == '--census':
        for line in open(args[1]):
            line = line.strip()
            if line: pts.append(F(line))
    else:
        pts = [F(a) for a in args]
    for y0 in pts:
        if y0 < F(3, 2):
            res = hostile_finite(y0)
            print(f"{str(y0):>28} = {float(y0):.6f}: {res[0]} ({res[1] if res[0] != 'descends' else 'path ' + str(res[1])})")
        else:
            res = slow_descent(y0)
            if res[0] == 'none':
                print(f"{str(y0):>28} = {float(y0):.6f}: >3/2, no greedy slow descent for k1<=400 (closest final R {res[1]:.4f} at k1={res[2]})")
            else:
                k1, s, K, ks = res
                print(f"{str(y0):>28} = {float(y0):.6f}: >3/2, DESCENDS: first move k1={k1}, then greedy; 2^{K} < 3^{s} "
                      f"(R={float(F(2**K,3**s)):.6f}); replay ok: {check_path(y0, ks)}")
