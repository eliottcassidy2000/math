"""
opus-2026-07-19-S407 (HYP-8005): three connection investigations.

(A) OBSERVER-TOURNAMENT REGULARITY BRIDGE: at each family's binding time t*,
    canon's half-turn comparator (THM-373/381) on the 14 positions (observer at 0
    + 13 runners) yields a tournament; c3 = C(14,3) - sum C(outdeg,2) counts its
    3-cycles (max 112 near-regular, 0 transitive).  Hypothesis: c3 anti-correlates
    with M; the extremals sit near-regular.  (The quantitative face of THM-640's
    'M-minimizer = H-maximizer' and of S399's 'tournament lens returns at the wall'.)
(B) SLACK-2 SINGLE-FAR EMPTINESS (my gate-table row): single-far families
    {1..N-2, N, x} never attain a first-gap-window value with slack != 1.
    Exhaustive N = 6..14, x <= 300 (gridmax-at-(N+1) prune, then exact scan).
(C) THE GHOST IS THE PENULTIMATE CONVERGENT: verify delta(ghost) = 1 at the
    maximizer for deep well (ghost 13 at 14/183), ladder3 (ghost 12 at 17/41),
    F4(31) (ghost 30 at 55/127) -- closing F1 + THM-1291 + THM-1292 into one loop.
"""
import random
from math import gcd
from fractions import Fraction
random.seed(8005)

def exact_max(V, qmax):
    bg, bq, wit = 0, 1, None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = q
            for v in V:
                r = (v*a) % q
                r = min(r, q-r)
                if r < m:
                    m = r
                    if m*bq < bg*q: break
            if m*bq > bg*q:
                bg, bq, wit = m, q, (a, q)
    return Fraction(bg, bq), wit

def c3_at(V, a, q):
    pos = [Fraction(0)] + [Fraction((v*a) % q, q) for v in V]
    n = len(pos)
    out = [0]*n
    for i in range(n):
        for j in range(n):
            if i == j: continue
            d = (pos[j] - pos[i]) % 1
            if d == 0 or d == Fraction(1,2):
                if i < j and d != 0: out[i] += 1   # tie-break half-turn by index
                elif d == 0 and i < j: out[i] += 1
            elif d < Fraction(1,2):
                out[i] += 1
    from math import comb
    return comb(n,3) - sum(comb(o,2) for o in out)

print("(A) tournament regularity at t* vs M:")
fams = [("AP13", list(range(1,14))), ("GW", list(range(1,12))+[13,24]),
        ("DW", list(range(1,13))+[182]), ("ladder3", list(range(1,12))+[13,36]),
        ("ladder5", list(range(1,12))+[13,60]), ("{1..12,14}", list(range(1,13))+[14]),
        ("{1..12,15}", list(range(1,13))+[15]), ("primes13", [2,3,5,7,11,13,17,19,23,29,31,37,41]),
        ("cluster N=5", list(range(5,18))), ("cluster N=20", list(range(20,33))),
        ("mixed", [3,7,10,14,22,25,31,38,44,52,57,63,70]),
        ("rand1", sorted(random.sample(range(1,150),13))),
        ("rand2", sorted(random.sample(range(1,150),13)))]
rows = []
for name, V in fams:
    M, (a, q) = exact_max(V, 300)
    c3 = c3_at(V, a, q)
    rows.append((float(M), c3, name, M))
    print(f"  {name}: M = {M} = {float(M):.4f}; c3(t*) = {c3} / 112")
rows.sort()
xs = [r[0] for r in rows]; ys = [r[1] for r in rows]
n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
den = (sum((x-mx)**2 for x in xs)*sum((y-my)**2 for y in ys))**0.5
print(f"  Pearson corr(M, c3) = {num/den:.3f}  (hypothesis: strongly negative)")

print("\n(B) slack-2 single-far emptiness, N = 6..14, x <= 300:")
hits = []
for N in range(6, 15):
    core = list(range(1, N-1)) + [N]
    for x in range(N+1, 301):
        if x == N: continue
        V = core + [x]
        # prune: gridmax at q = N+1 certifies M >= 1/(N+1) (not in window)
        ok = False
        for j in range(1, N+1):
            if all(min((v*j) % (N+1), (N+1)-(v*j) % (N+1)) >= 1 for v in V):
                ok = True; break
        if ok: continue
        M, (a, q) = exact_max(V, 2*x + 40)
        if M < Fraction(1, N+1) and M > Fraction(1, N+2):
            D = None
            for s in range(2, 3*(N+1)*8):
                if M*s == int(M*s): D = int(M*s); break
            s0 = int(D/M)
            slack = (N+1)*D - s0
            hits.append((N, x, M, D, s0, slack))
            print(f"  N={N} x={x}: M={M} in window; (D,s)=({D},{s0}), slack={slack}")
print(f"  window hits: {len(hits)}; slacks found: {sorted(set(h[5] for h in hits))} "
      f"(hypothesis: all = 1)")

print("\n(C) the ghost is the penultimate convergent (delta(ghost) = 1):")
for name, ghost, a, q in [("DW", 13, 14, 183), ("ladder3", 12, 17, 41),
                          ("F4(31)", 30, 55, 127)]:
    d = min((ghost*a) % q, q-(ghost*a) % q)
    print(f"  {name}: delta({ghost}) at {a}/{q} = {d}  (== 1: {d == 1})")
