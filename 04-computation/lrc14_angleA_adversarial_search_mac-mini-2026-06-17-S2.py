#!/usr/bin/env python3
"""
ANGLE A (adversarial): hunt for the covering 13-set that MAXIMIZES DIP/SLACK,
to find whether the hand-in-hand reduction (DIP <= SLACK => M(S)>=1/14) can be
broken. Also hunt for any covering 13-set with M(S) < 1/14 (would refute LRC14).

We are honest: the reduction is a SUFFICIENT certificate. We test whether the
certificate ever FAILS (DIP>SLACK) and, separately, whether M(S) ever < 1/14.

Strategy:
  - All covering 13-sets of the form A u {w}, A a 12-subset of {1..13}, w=14m.
  - Slack is maximized by SMALLEST uncovered q; but ratio is stressed when the
    parked runner is the ONLY multiple of a LARGE q (small slack) yet still
    forced. We let the certificate pick the smallest uncovered q (best slack)
    AND also evaluate the worst-q variant (smallest slack) to see if a chooser
    could break it.
  - Exhaustive over A subset of {1..13} (size 12 -> 13 of them) x m in 1..N.
  - Plus brute random covering 13-sets with entries in 1..28 and a parked
    multiple of 14, to escape the {1..13}+14m scaffold.

stdlib only.
"""
from fractions import Fraction as F
from itertools import combinations
import random

def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at

def covers(S, q):
    return any(v % q == 0 for v in S)

def uncovered(S, qmax=14):
    return [q for q in range(2, qmax+1) if not covers(S, q)]

def smallest_uncovered(A, qmax=13):
    u = uncovered(A, qmax)
    return u[0] if u else None

def largest_uncovered(A, qmax=13):
    u = uncovered(A, qmax)
    return u[-1] if u else None

ONE14 = F(1, 14)

def analyze(S, A):
    """Return certificate stats for S=A u parked, choosing best-slack q."""
    MS, tS = M(S)
    MA, tA = M(A)
    dip = MA - MS
    qbest = smallest_uncovered(A, 13)   # largest slack
    qworst = largest_uncovered(A, 13)   # smallest slack
    out = dict(MS=MS, tS=tS, MA=MA, tA=tA, dip=dip,
               qbest=qbest, qworst=qworst)
    if qbest:
        sl = F(1,qbest)-ONE14
        out['slack_best'] = sl
        out['ratio_best'] = (dip/sl) if sl>0 else (F(0) if dip<=0 else None)
    if qworst:
        slw = F(1,qworst)-ONE14
        out['slack_worst'] = slw
        out['ratio_worst'] = (dip/slw) if slw>0 else (F(0) if dip<=0 else None)
    return out

def main():
    # ---- Exhaustive: A = 12-subset of {1..13}, w=14m ----
    base = list(range(1,14))
    worst_best = (F(0), None)     # max ratio under best-slack chooser
    worst_worst = (F(0), None)    # max ratio under worst-slack (adversarial q)
    minMS = (F(1), None)
    n = 0; broke = []
    examples = []
    for A in combinations(base, 12):
        A = list(A)
        for m in range(1, 25):
            w = 14*m
            S = sorted(set(A) | {w})
            if len(S) != 13:
                continue
            if uncovered(S, 14):   # must be a COVERING set
                continue
            if smallest_uncovered(A,13) is None:  # core must be easy
                continue
            st = analyze(S, A)
            n += 1
            if st['MS'] < minMS[0]:
                minMS = (st['MS'], (A, S, st))
            if st['MS'] < ONE14:
                broke.append((A,S,st))
            rb = st.get('ratio_best')
            if rb is not None and rb > worst_best[0]:
                worst_best = (rb, (A,S,st))
            rw = st.get('ratio_worst')
            if rw is not None and rw > worst_worst[0]:
                worst_worst = (rw, (A,S,st))
    print("="*78)
    print("EXHAUSTIVE: A=12-subset of {1..13}, parked w=14m, m=1..24, covering 13-sets")
    print("="*78)
    print(f"tested {n} covering 13-sets")
    print(f"GLOBAL MIN M(S) = {minMS[0]} ({float(minMS[0]):.6f}) >= 1/14 ? {minMS[0]>=ONE14}")
    if minMS[1]:
        A,S,st = minMS[1]
        print(f"   at S={S}, A={A}, tau*={st['tS']}")
    print(f"\nWORST ratio under BEST-slack chooser (smallest uncovered q): {worst_best[0]} "
          f"({float(worst_best[0]):.6f})")
    if worst_best[1]:
        A,S,st = worst_best[1]
        print(f"   S={S}")
        print(f"   A={A}, M(A)={st['MA']}, M(S)={st['MS']}, dip={st['dip']}")
        print(f"   qbest={st['qbest']}, slack={st['slack_best']} ({float(st['slack_best']):.6f})")
        print(f"   reduction certifies? dip<=slack: {st['dip']<=st['slack_best']}")
    print(f"\nWORST ratio under ADVERSARIAL (largest uncovered q, smallest slack): "
          f"{worst_worst[0]} ({float(worst_worst[0]):.6f})")
    if worst_worst[1]:
        A,S,st = worst_worst[1]
        print(f"   S={S}")
        print(f"   A={A}, M(A)={st['MA']}, M(S)={st['MS']}, dip={st['dip']}")
        print(f"   qworst={st['qworst']}, slack_worst={st['slack_worst']} ({float(st['slack_worst']):.6f})")
        print(f"   adversarial dip<=slack_worst? {st['dip']<=st['slack_worst']}  "
              f"(but real M(S)>=1/14? {st['MS']>=ONE14})")
    print(f"\ncovering 13-sets that break LRC14 (M(S)<1/14): {len(broke)}")

    # ---- Random covering 13-sets with a parked multiple of 14 ----
    print("\n" + "="*78)
    print("RANDOM: covering 13-sets, entries in 1..28, exactly one parked 14k")
    print("="*78)
    random.seed(12345)
    rn = 0; rbroke = []; rworst = (F(0), None); rminMS = (F(1), None)
    tries = 0
    while rn < 4000 and tries < 200000:
        tries += 1
        parked = 14*random.randint(1,2)   # 14 or 28
        pool = [x for x in range(1,29) if x % 14 != 0]
        rest = random.sample(pool, 12)
        S = sorted(set(rest) | {parked})
        if len(S) != 13:
            continue
        if uncovered(S,14):
            continue
        A = sorted(set(S) - {parked})
        if len(A) != 12:
            continue
        if smallest_uncovered(A,13) is None:
            continue
        st = analyze(S, A)
        rn += 1
        if st['MS'] < rminMS[0]:
            rminMS = (st['MS'], (A,S,st))
        if st['MS'] < ONE14:
            rbroke.append((A,S,st))
        rb = st.get('ratio_best')
        if rb is not None and rb > rworst[0]:
            rworst = (rb, (A,S,st))
    print(f"tested {rn} random covering 13-sets (parked=14 or 28)")
    print(f"GLOBAL MIN M(S) = {rminMS[0]} ({float(rminMS[0]):.6f}) >= 1/14 ? {rminMS[0]>=ONE14}")
    if rminMS[1]:
        A,S,st = rminMS[1]
        print(f"   at S={S}, tau*={st['tS']}")
    print(f"WORST best-slack ratio = {rworst[0]} ({float(rworst[0]):.6f})")
    if rworst[1]:
        A,S,st = rworst[1]
        print(f"   S={S}, A={A}")
        print(f"   M(A)={st['MA']}, M(S)={st['MS']}, dip={st['dip']}, qbest={st['qbest']}, "
              f"slack={st.get('slack_best')}")
        print(f"   dip<=slack? {st['dip']<=st['slack_best']}")
    print(f"covering 13-sets breaking LRC14 (M(S)<1/14): {len(rbroke)}")
    for c in rbroke[:5]:
        A,S,st = c
        print(f"   S={S} M(S)={st['MS']} ({float(st['MS']):.6f}) tau={st['tS']}")

if __name__ == "__main__":
    main()
