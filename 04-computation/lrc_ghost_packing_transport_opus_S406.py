"""
opus-2026-07-19-S406: THM-1292 referee -- the F1 transport (tower seals -> covering-min).

(T1) GHOST-PACKING CEILING (L3-transport, new 5-line proof): if {1..m} subset V,
     m+1 not in V, and K*(m+1) in V, then M(V) <= K/(K(m+1)+1).
     Referee: exact M for structured families vs the ceiling; adversarial search for
     violations over random V with the hypothesis; tightness cases.
(T2) L1-TRANSPORT (witness four-liner at Q=183): at a = 14, danger classes are
     exactly the 13-multiples 13s with |s| <= 13; the deep well misses them all
     (182 = 13*14 sits at |s| = 14, exactly the floor) => M(DW) >= 14/183.
(T3) The composition check: ceiling value K/(K(m+1)+1) = M({1..m} u {K(m+1)})
     exactly (the ghost-subfamily attains its own ceiling) -- so the lemma =
     SubfamilyCap o generalized-THM-633, now with a uniform packing proof.
"""
import random
from math import gcd
from fractions import Fraction
random.seed(1292)

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

def ghost_params(V):
    m = 0
    while (m+1) in V: m += 1
    if m == 0: return None
    g = m + 1
    Ks = sorted(v // g for v in V if v % g == 0)
    return (m, g, Ks[0]) if Ks else None

print("(T1) ceiling referee:")
fams = [("deep well", list(range(1,13))+[182], 260),
        ("GW", list(range(1,12))+[13,24], 200),
        ("ladder3", list(range(1,12))+[13,36], 200),
        ("ladder4", list(range(1,12))+[13,48], 200),
        ("F4(31)", list(range(1,30))+[31,120], 200)]
for name, V, qm in fams:
    gp = ghost_params(V)
    if not gp: print(f"  {name}: no ghost structure"); continue
    m, g, K = gp
    ceil_ = Fraction(K, K*g+1)
    M, wit = exact_max(V, qm)
    print(f"  {name}: m={m} ghost={g} K={K}: M={M} <= ceiling {ceil_}: {M <= ceil_}"
          f"{'  TIGHT' if M == ceil_ else ''}")
viol = 0
for _ in range(3000):
    m = random.randint(3, 11)
    K = random.randint(1, 8)
    extra = random.sample([x for x in range(m+2, 160) if x % (m+1) or x == K*(m+1)], 3)
    V = list(range(1, m+1)) + [K*(m+1)] + extra
    if any(x % (m+1) == 0 and x // (m+1) < K for x in extra): continue
    gp = ghost_params(V)
    if not gp or gp[1] != m+1: continue
    K0 = gp[2]
    M, _ = exact_max(V, 90)
    if M > Fraction(K0, K0*(m+1)+1):
        viol += 1
        print("  VIOLATION:", V, M)
print(f"  adversarial: violations {viol}/~3000 random hypothesis-satisfying families")

print("\n(T2) L1-transport at Q=183, a=14:")
Q, a = 183, 14
danger = sorted(v for v in range(1, Q) if min((v*a) % Q, Q-(v*a) % Q) < 14)
mults13 = sorted(set((13*s) % Q for s in range(1, 14)) | set((-13*s) % Q for s in range(1, 14)))
print(f"  danger residues == 13-multiples with |s|<=13: {danger == mults13}")
DW = list(range(1,13)) + [182]
hits = [v for v in DW if (v % Q) in danger]
val = min(Fraction(min((v*a) % Q, Q-(v*a) % Q), Q) for v in DW)
print(f"  deep well hits danger: {hits} (expect []); value at 14/183 = {val} (expect 14/183)")

print("\n(T3) ghost-subfamily attains its own ceiling:")
for m, K in [(12,14),(11,3),(11,4),(9,5),(7,2)]:
    S = list(range(1, m+1)) + [K*(m+1)]
    M, wit = exact_max(S, 3*(K*(m+1)+1))
    print(f"  {{1..{m}}} u {{{K*(m+1)}}}: M = {M} vs K/(K(m+1)+1) = "
          f"{Fraction(K, K*(m+1)+1)}: {'EXACT' if M == Fraction(K, K*(m+1)+1) else 'DIFFERS'}")
