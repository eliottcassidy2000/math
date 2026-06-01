from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce
from collections import Counter
import random
from lrc import frac_dist, D, witness_times, M_witness, n_of

def is_primitive(sp): return reduce(gcd,sp)==1

# ---------- PART 2: witness-time tournament ----------
def witness_stats(speeds):
    times = witness_times(speeds)
    Ts = list(times.keys())
    N = len(Ts)
    # D values
    dvals = {t: D(t,speeds) for t in Ts}
    Mt = max(Ts, key=lambda t: dvals[t])
    M = dvals[Mt]
    # is the source unique? (strictly beats all others)
    sources = [t for t in Ts if dvals[t]==M]
    # theoretical count estimate: sum over pairs of (v_i+v_j-1), before dedup
    raw = sum((speeds[i]+speeds[j]-1) for i,j in combinations(range(len(speeds)),2))
    return N, raw, len(sources), M, Mt, dvals

def part2(maxspeed=10, ms=(4,5,6), sample=None):
    for m in ms:
        sets=[list(c) for c in combinations(range(1,maxspeed+1),m) if is_primitive(list(c))]
        if sample and len(sets)>sample:
            random.seed(0); sets=random.sample(sets,sample)
        Ns=[]; raws=[]; nsrc=[]
        for sp in sets:
            N,raw,ns,M,Mt,_=witness_stats(sp)
            Ns.append(N); raws.append(raw); nsrc.append(ns)
        print(f"m={m} n={m+1}: {len(sets)} sets")
        print(f"  distinct witness times: mean {sum(Ns)/len(Ns):.1f}, min {min(Ns)}, max {max(Ns)}")
        print(f"  raw (pre-dedup) count:  mean {sum(raws)/len(raws):.1f}  (dedup ratio {sum(Ns)/sum(raws):.3f})")
        print(f"  C(m,2)*avg_sum approx: {len(list(combinations(range(m),2)))} pairs")
        print(f"  unique source (D=M attained once): {sum(1 for x in nsrc if x==1)}/{len(sets)}; multi-source max {max(nsrc)}")

# ---------- PART 4: meta-LRC on pairwise sums ----------
def meta_instance(speeds):
    """Build meta-runners = pairwise sums. Reduce to primitive distinct set, compute meta-M*n_meta."""
    sums = sorted(set(speeds[i]+speeds[j] for i,j in combinations(range(len(speeds)),2)))
    return sums

def part4(maxspeed=9, ms=(4,5), sample=None):
    for m in ms:
        sets=[list(c) for c in combinations(range(1,maxspeed+1),m) if is_primitive(list(c))]
        if sample and len(sets)>sample:
            random.seed(0); sets=random.sample(sets,sample)
        # Hypothesis A: does the original optimal time 1/n behave specially among meta-runners?
        # Hypothesis B: is M(meta-set) related to M(orig)?
        relate=Counter(); collisions=0
        for sp in sets:
            sums=meta_instance(sp)
            n=m+1
            Mo,to,_=M_witness(sp)
            # distinct pairwise sums count
            full = len(list(combinations(range(m),2)))
            if len(sums)<full: collisions+=1
            # meta as its own LRC instance (primitive reduce)
            g=reduce(gcd,sums); prim=[s//g for s in sums]
            if len(set(prim))==len(prim) and len(prim)>=1:
                Mm,tm,_=M_witness(prim) if len(prim)>=2 else (Fraction(1,2),Fraction(1,2),None)
                nm=len(prim)+1
                # compare M*n
                relate[(round(float(Mo*n),3))]+=0  # placeholder
        print(f"m={m} n={m+1}: {len(sets)} sets; sets with colliding pairwise sums: {collisions}")
        # Key structural test: is min pairwise sum vs n the discriminator?
        # Test the natural meta-claim: orig M*n=1 (tight) iff meta-set contains n as an element (a pair sums to n)
        contains_n=0; tight_and_contains=0; tight=0; contains_and_tight=0
        for sp in sets:
            sums=meta_instance(sp); n=m+1
            Mo,_,_=M_witness(sp)
            isT = (Mo*n==1)
            hasN = (n in sums)
            if hasN: contains_n+=1
            if isT: tight+=1
            if isT and hasN: tight_and_contains+=1
        print(f"  tight sets: {tight}; sets where some pair sums to exactly n: {contains_n}; tight AND pair-sums-to-n: {tight_and_contains}")
        print(f"  => tight => some pair sums to n: {tight==tight_and_contains}")

if __name__=="__main__":
    print("######## PART 2: WITNESS-TIME TOURNAMENT ########")
    part2(10,(4,5,6))
    print("\n######## PART 4: META-LRC ON PAIRWISE SUMS ########")
    part4(9,(4,5))
