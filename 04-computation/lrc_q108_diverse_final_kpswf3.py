#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_diverse_final_kpswf3.py   (kind-pasteur 2026-06-21, THREAD D final)

Pin down the 3 leads with real signal:

 (A') SCALING INVARIANCE on WIDE sets (k>=8 so p0>0): is measS7(c*E)==measS7(E) for
      gcd(c,7)=1?  And is it FALSE when 7|c? (structural lemma, search reduction)

 (B') GENUINE Bonferroni: the *truncated* IE sums P_t = sum_{j<=t}(-1)^j S_j.  Classical
      Bonferroni: even t -> UPPER bound on p0, odd t -> LOWER.  Check the bracketing
      actually holds (the earlier table's last entry was the FULL sum = p0, not a
      truncation -- circular).  Report the smallest even t whose P_t <= cap, i.e. how
      many miss-orders you must control for an elementary proof.

 (C') MULTIPLICATIVE-INDEPENDENCE anomaly: consecutive [1..k] is the WORST (max p0);
      multiplicatively-spread sets (primes, geometric) are FAR lower.  Conjecture:
      consecutive integers MAXIMIZE measS7 among k-subsets.  Probe exhaustively for
      small k whether [1..k] (or its dilates) is the unique maximizer, and quantify the
      gap to the runner-up.  If consecutive is the max, the LRC cover bound reduces to
      ONE family -> a clean finite check.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, comb

P = 7
def sector(yf): return int(P*yf)
def breakpoints(E):
    bp={Fr(0),Fr(1)}
    for e in E:
        if e==0: continue
        for t in range(0,P*e): bp.add(Fr(t,P*e))
    return sorted(bp)
def measS7(E):
    E=[int(e) for e in E if int(e)!=0]
    if len(set(e%P for e in E))<P and len(E)>=P:
        pass
    xs=breakpoints(E); tot=Fr(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        if len(set(sector((e*mid)%1) for e in E))==P: tot+=(b-a)
    return tot
def miss_prob(E,missing):
    E=[int(e) for e in E if int(e)!=0]
    xs=breakpoints(E); miss=set(missing); tot=Fr(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        if not (set(sector((e*mid)%1) for e in E)&miss): tot+=(b-a)
    return tot
def Svec(E):
    S=[Fr(0)]*(P+1); S[0]=Fr(1)
    for j in range(1,P+1):
        for T in itertools.combinations(range(P),j):
            S[j]+=miss_prob(E,set(T))
    return S

CAP={8:Fr(2243,5880),9:Fr(1979,4004),10:Fr(55,91)}

def main():
    print("#"*80)
    print("# THREAD D FINAL: scaling(wide), genuine Bonferroni, consecutive-maximizer")
    print("#"*80)

    # ---- (A') scaling on WIDE sets ----
    print("\n=== (A') SCALING INVARIANCE on wide sets (p0>0) ===")
    wide=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8,9],[1,2,4,5,8,10,11,13],
          [2,3,5,7,11,13,17,19,23]]
    allinv=True
    for E in wide:
        base=measS7(E); row=[]
        for c in [2,3,4,5,6,8,9,10,11,12,13]:
            if gcd(c,7)!=0 and c%7!=0:
                v=measS7([c*e for e in E]); row.append(v==base)
                if v!=base: allinv=False
        v7=measS7([7*e for e in E])
        print(f"  E(k={len(E)}) base={float(base):.5f}  inv(all c coprime7)={all(row)}"
              f"   measS7(7E)={float(v7):.5f} {'(DIFFERS-good)' if v7!=base else '(same)'}")
    print(f"  VERDICT scaling-invariance for gcd(c,7)=1: {'HOLDS (all wide tests)' if allinv else 'FAILS'}")

    # ---- (B') genuine Bonferroni truncations ----
    print("\n=== (B') GENUINE Bonferroni truncations P_t = sum_{j<=t}(-1)^j S_j ===")
    print("  classical: P_0>=P_2>=...>=p0>=...>=P_3>=P_1 ? (does bracketing hold here?)")
    for name,E in {"k8[1..8]":[1,2,3,4,5,6,7,8],"k9[1..9]":[1,2,3,4,5,6,7,8,9]}.items():
        S=Svec(E); part=[]; acc=Fr(0)
        for j in range(P+1):
            acc+=(-1)**j*S[j]; part.append(acc)
        p0=part[-1]; cap=CAP[len(E)]
        # check monotone bracketing
        evens=[(t,part[t]) for t in range(0,P,2)]   # EXCLUDE t=P (that's the full sum=p0)
        odds=[(t,part[t]) for t in range(1,P,2)]
        upper_ok=all(part[t]>=p0 for t in range(0,P,2))
        lower_ok=all(part[t]<=p0 for t in range(1,P,2))
        print(f"  {name}: p0={float(p0):.5f} cap={float(cap):.4f}")
        print(f"    even-trunc (upper?) {[(t,round(float(v),4)) for t,v in evens]}  bracket-upper={upper_ok}")
        print(f"    odd-trunc  (lower?) {[(t,round(float(v),4)) for t,v in odds]}  bracket-lower={lower_ok}")
        best_even=min(v for t,v in evens)
        print(f"    best genuine even-trunc UPPER (t<7) = {float(best_even):.4f}"
              f"  {'<=cap' if best_even<=cap else '> cap (Bonferroni too weak)'}")

    # ---- (C') consecutive-maximizer conjecture ----
    print("\n=== (C') CONSECUTIVE-MAXIMIZER: is [1..k] the max-p0 k-subset? ===")
    # exhaustive over k-subsets of [1..N] (dilation-normalize: fix min element=1, gcd=1)
    for k,N in [(8,12),(9,13)]:
        cons=list(range(1,k+1)); pc=measS7(cons)
        best=pc; bestE=cons; beat=0; total=0; ties=0
        # search subsets of {1..N} of size k containing 1, gcd 1, to bound work
        from math import gcd as g
        for combo in itertools.combinations(range(1,N+1),k):
            if combo[0]!=1: continue
            gg=0
            for e in combo: gg=g(gg,e)
            if gg!=1: continue
            total+=1
            v=measS7(combo)
            if v>best+Fr(1,10**9):
                best=v; bestE=combo; beat+=1
            elif abs(v-pc)<Fr(1,10**12):
                ties+=1
        print(f"  k={k},N={N}: consec p0={float(pc):.5f}  searched {total} subsets")
        print(f"     max found = {float(best):.5f} at {bestE}"
              f"   {'CONSEC IS MAX' if best<=pc+Fr(1,10**9) else 'CONSEC NOT MAX (beaten by '+str(beat)+')'}")
        # runner-up gap
        runner=Fr(0); runnerE=None
        for combo in itertools.combinations(range(1,N+1),k):
            if combo[0]!=1: continue
            gg=0
            for e in combo: gg=g(gg,e)
            if gg!=1: continue
            v=measS7(combo)
            if v<pc-Fr(1,10**12) and v>runner:
                runner=v; runnerE=combo
        print(f"     runner-up (strictly below consec) = {float(runner):.5f} at {runnerE}"
              f"   GAP={float(pc-runner):.5f}")

    print("\nDONE.")

if __name__=="__main__":
    main()
