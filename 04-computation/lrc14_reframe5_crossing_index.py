#!/usr/bin/env python3
"""
REFRAME 5 part (2b): structural law for the crossing index j.

Frozen-core (C5) is REFUTED: the optimum genuinely MOVES, the dip is real.
So the binding sum-crossing (flank, w) at tau*=t/D, D=flank+w, gap = j/D really
is the operative quantity, and M(S)>=1/14 <=> j >= D/14 = (flank+w)/14.

NON-CIRCULAR target: a LOWER bound on j from (flank, w, 14, q) arithmetic alone.

Tests:
  (T1) j = M(S)*D : tabulate j, D, flank across all binding N=14 configs + resonant k.
  (T2) is j determined by  j = round-to-make-others-clear? Specifically, at the
       sum-crossing the gap is the distance of flank (= w) from integer, and that
       equals  ||flank * t/D||  with t chosen so that the SMALL runners {1..} clear.
  (T3) THE candidate non-circular law:  j >= ceil((flank+w)/14)  AND the gap j/D
       comes from  j = (number of the integer the pair crosses) so that
       j/D >= flank/(flank+w)*(1/?)... we just measure j vs flank, w, q.
  (T4) Does  j*14 >= D  reduce to a PROVABLE statement  j >= flank  ?  i.e. is
       j >= flank always, and is flank+w <= 14*flank (i.e. w <= 13*flank)?  If w<=13*flank
       AND j>=flank then j/D = j/(flank+w) >= flank/(14 flank) = 1/14.  TWO clean
       sub-inequalities, BOTH potentially non-circular.
"""
from fractions import Fraction as F
from math import gcd
import itertools


def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r


def g(S, t):
    return min(nrm(v * t) for v in S)


def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C


def Mall(S):
    b = F(0); ats = []
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; ats = [t]
        elif v == b:
            ats.append(t)
    return b, ats


def lcm(a, b):
    return a * b // gcd(a, b)


def binding_pair_with(S, tstar, w):
    """find a binding pair (flank,w) at tstar; return (flank, D=flank+w, kind, gapnum_over_D)."""
    val = g(S, tstar)
    mr = [v for v in S if nrm(v * tstar) == val]
    if w not in mr:
        return None
    for a in mr:
        if a == w:
            continue
        if nrm((a + w) * tstar) == 0:
            return (a, a + w, 'sum', val * (a + w))
        if nrm((w - a) * tstar) == 0 and w != a:
            return (a, w - a, 'diff', val * (w - a))
    return None


print("=" * 78)
print("REFRAME 5 (2b): crossing-index law  j = M(S)*D,  and the split j>=flank & w<=13 flank")
print("=" * 78)

# Collect all binding configs: single-drop + resonant k + two-drop, N=14 level.
rows = []
for q in range(2, 14):
    A = [v for v in range(1, 14) if v != q]
    L = lcm(q, 14)
    for k in range(1, 13):
        w = L * k
        if w in A:
            continue
        S = A + [w]
        MS, taus = Mall(S)
        info = None
        for t in taus:
            bp = binding_pair_with(S, t, w)
            if bp:
                info = (t, bp); break
        if info is None:
            continue  # w not binding -> dip 0, skip
        t, (flank, D, kind, jnum) = info
        # jnum = M(S)*D should be an integer (crossing index over D)
        jint = jnum.numerator if jnum.denominator == 1 else jnum
        rows.append(dict(q=q, k=k, w=w, MS=MS, flank=flank, D=D, kind=kind,
                         j=jint, jint=(jnum.denominator == 1), drop_set=A))

print(f"\nCollected {len(rows)} binding single-drop configs (w actually binds).")
print(f"\n{'q':>2} {'k':>2} {'w':>4} {'flank':>5} {'kind':>4} {'D':>5} {'j':>3} "
      f"{'M(S)':>8} {'j>=flank?':>9} {'w<=13*flank?':>12} {'j*14>=D?':>9}")
viol_jflank = 0; viol_w13 = 0; viol_main = 0
maxratio_w = None
for r in rows:
    jf = r['j'] >= r['flank']
    w13 = (r['kind'] == 'diff') or (r['w'] <= 13 * r['flank'])  # sum-crossing case
    main = F(r['j']) * 14 >= r['D']
    if not jf: viol_jflank += 1
    if not w13: viol_w13 += 1
    if not main: viol_main += 1
    rat = F(r['w'], r['flank'])
    if maxratio_w is None or rat > maxratio_w[0]:
        maxratio_w = (rat, r['q'], r['k'], r['flank'], r['w'])
    print(f"{r['q']:>2} {r['k']:>2} {r['w']:>4} {r['flank']:>5} {r['kind']:>4} {r['D']:>5} "
          f"{str(r['j']):>3} {str(r['MS']):>8} {str(jf):>9} {str(w13):>12} {str(main):>9}")

print(f"\n  violations  j>=flank: {viol_jflank}   w<=13*flank (sum only): {viol_w13}   "
      f"j*14>=D (THE goal): {viol_main}")
if maxratio_w:
    print(f"  max w/flank ratio: {float(maxratio_w[0]):.3f} "
          f"(q={maxratio_w[1]} k={maxratio_w[2]} flank={maxratio_w[3]} w={maxratio_w[4]})")

print("""
INTERPRETATION:
  If  j >= flank  AND  w <= 13*flank  both hold with NON-circular proofs, then
      j/D = j/(flank+w) >= flank/(flank + 13*flank) = flank/(14 flank) = 1/14,
  giving M(S) >= 1/14 WITHOUT assuming it.  So the question becomes whether
  these two sub-inequalities are themselves non-circular -- check below.
""")

# Are the two sub-inequalities themselves non-circular?
# j>=flank: j is the crossing index = round(D * tau*). It depends on tau* which depends
#   on the WHOLE set including the clearance of small runners. NOT obviously arithmetic.
# w<=13*flank: w = lcm(q,14)*k is arithmetic in q. flank is the binding small runner,
#   determined by the optimum. NOT obviously arithmetic either.
# So measure: is flank PREDICTABLE from (q,14) alone, independent of solving for M?
print("-" * 78)
print("Is 'flank' predictable from (q,14) arithmetic, or does it require solving M?")
print("-" * 78)
# tabulate flank as function of (q) for k=1 (the minimal, hardest case)
print(f"{'q':>3} {'w=lcm(q,14)':>11} {'flank(k=1)':>10} {'w/flank':>8} {'<=13?':>6}")
for r in rows:
    if r['k'] == 1:
        rat = F(r['w'], r['flank'])
        print(f"{r['q']:>3} {r['w']:>11} {r['flank']:>10} {float(rat):>8.3f} "
              f"{str(rat<=13):>6}")

print("\nDONE.")
