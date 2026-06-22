"""
LRC realizability: is the AP the UNIQUE realizable extremizer of the cover measure?
(kind-pasteur S31j)

Reframing: LRC(14) <=> the 13 danger zones D_i={t:||v_i t||<1/14} (each measure 1/7)
do NOT cover [0,1).  Equivalently the cover atom p0(E)=meas{x: all 6 inner sectors hit}
stays <= cap.  A counterexample over-covers (p0>cap).  REALIZABILITY claim: no integer
speed set beats the AP -- the AP/interval is the unique extremizer.

This test:
 - compute p0(E) (exact-ish, breakpoint sweep) and additive energy A(E)=#{a+b=c+d} for
   the AP/interval vs spread (Sidon-ish), random, and dilated sets at fixed k;
 - DETERMINE the direction: does the interval MAXIMIZE p0 (most equidistributed coverage)
   or MINIMIZE it?  Does p0 track additive energy or anti-track it?  This pins down which
   spectral extremality (Fejer concentration vs equidistribution) is the realizability.
"""
import random, itertools

def p0(E):
    E = sorted(set(e for e in E if e != 0))
    bset = {0.0, 1.0}
    for e in E:
        for j in range(8):
            b = j/7.0; m = 0
            while True:
                xv = (b+m)/e
                if xv >= 1: break
                if xv >= 0: bset.add(xv)
                m += 1
    B = sorted(bset); tot = 0.0
    for lo, hi in zip(B, B[1:]):
        if hi <= lo: continue
        mid = (lo+hi)/2
        hit = set(int((e*mid) % 1 * 7) for e in E) & set(range(1,7))
        if len(hit) == 6: tot += hi-lo
    return tot

def add_energy(E):
    E = list(E); from collections import Counter
    c = Counter(a+b for a in E for b in E)
    return sum(v*v for v in c.values())

def sidon_like(k, hi=60):
    # greedy Sidon set (all pairwise sums distinct) in [0,hi]
    S = [0]; sums = set()
    x = 1
    while len(S) < k and x <= hi:
        new = [x+s for s in S]
        if len(set(new)) == len(new) and not (set(new) & sums):
            sums |= set(new) | {2*x}
            # rebuild sums properly
            S.append(x)
            sums = set(a+b for a in S for b in S)
        x += 1
    return S[:k]

if __name__ == "__main__":
    k = 9
    interval = list(range(k))
    sidon = sidon_like(k)
    random.seed(1)
    rnd = [sorted(random.sample(range(1, 40), k-1)) for _ in range(3)]
    rnd = [[0]+r for r in rnd]
    cap9 = 1979/4004
    print(f"k={k}, cap_9={cap9:.4f}")
    print(f"{'set':28s} {'p0':>8s} {'A(E)':>8s}  note")
    def show(name, E):
        print(f"  {name:26s} {p0(E):8.4f} {add_energy(E):8d}")
    show("INTERVAL/AP "+str(interval), interval)
    show("Sidon "+str(sidon), sidon)
    for i,r in enumerate(rnd):
        show(f"random{i} "+str(r), r)
    # dilation invariance check
    show("2*interval", [2*x for x in interval])
    print("\n=> If INTERVAL has the LARGEST p0, the AP is the cover-extremizer (the realizable max);")
    print("   correlate p0 with A(E): tracks (Fejer concentration) or anti-tracks (equidistribution)?")
