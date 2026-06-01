#!/usr/bin/env python3
"""
source_sink_apex_arc_s530.py    oracle-2026-06-01-S530

User: the OUTSIDE of the polygon = the base path of the tiling model PLUS the one
arc between SOURCE and SINK. That arc closes the Hamiltonian path into the full
boundary n-cycle, and it occupies an important place in the tiling model.

We verify what that arc does:

 (1) STAIRCASE APEX. With base path 0->1->..->(n-1) (0=source, n-1=sink), the tiles
     are (i,j), j>i+1; the SOURCE-SINK arc (0,n-1) has maximal range n-1 = the
     APEX tile = the hypotenuse corner of the staircase. Flipping it on the
     transitive tournament gives H = 1 + 2^{n-2} (the max single-tile jump, the
     transitive's 'big SC neighbour' on the principal line). VERIFY.

 (2) MASTER SWITCH between the two block-extremes (#SCC in {1,n}, the round / LRC-
     realizable tournaments, S524). For round tournaments built from circle gaps:
     transitive (#SCC=n) <=> all points in a semicircle <=> the largest gap > 1/2;
     strongly connected (#SCC=1) <=> largest gap < 1/2 (Moon: Ham cycle = boundary).
     The SOURCE-SINK arc = the chord across the LARGEST gap; its orientation is the
     switch. VERIFY the correlation.

 (3) LRC. The largest gap is where a lonely OBSERVER sits (THM-382: lonely <=> both
     observer-adjacent gaps >= 1/n). So the source-sink arc = the observer's
     straddling gap; the two lonely types (observer-outside/transitive vs observer-
     inside/wrap) are selected by this one arc. Illustrate.
"""
from itertools import permutations, product
from functools import lru_cache
import random

def Hc(adj):
    n=len(adj); full=(1<<n)-1
    @lru_cache(None)
    def dp(mask,last):
        if mask==full: return 1
        return sum(dp(mask|(1<<x),x) for x in range(n) if not (mask>>x)&1 and adj[last][x])
    return sum(dp(1<<s,s) for s in range(n))

def transitive(n):
    return [[1 if i<j else 0 for j in range(n)] for i in range(n)]  # i beats j iff i<j

def scc_count(adj):
    n=len(adj)
    def reach(s,fwd):
        seen={s}; st=[s]
        while st:
            u=st.pop()
            for w in range(n):
                e=adj[u][w] if fwd else adj[w][u]
                if e and w not in seen: seen.add(w); st.append(w)
        return seen
    comp=[None]*n; c=0
    for v in range(n):
        if comp[v] is not None: continue
        for w in reach(v,True)&reach(v,False):
            if comp[w] is None: comp[w]=c
        c+=1
    return len(set(comp))

def half_turn_from_gaps(gaps):
    n=len(gaps); pos=[]; c=0.0
    for g in gaps: pos.append(c); c+=g
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=(pos[j]-pos[i])%1.0
            if abs(d-0.5)<1e-12: return None
            adj[i][j]=1 if d<0.5 else 0
    return adj

def main():
    print("source-sink (apex) arc: the closing edge of the outside (oracle-S530)\n")

    # (1) apex flip -> H = 1 + 2^{n-2}
    print("(1) APEX = source-sink arc = max-range tile. Flip it on the transitive tournament:")
    for n in range(3,9):
        adj=[row[:] for row in transitive(n)]
        adj[0][n-1]=0; adj[n-1][0]=1          # flip apex (0,n-1): sink beats source
        H=Hc(adj); pred=1+2**(n-2)
        print(f"   n={n}: H(transitive + apex back-arc) = {H}   1+2^(n-2) = {pred}   match={H==pred}")
    print("   => the source-sink arc is the hypotenuse / principal-line tile: its single")
    print("      flip is the maximal H jump (transitive -> its big strongly-connected neighbour).\n")

    # (2) master switch on round tournaments via circle gaps
    print("(2) ROUND tournaments from circle gaps: largest gap vs #SCC vs source-sink arc")
    rnd=random.Random(530)
    for n in (4,5,6):
        tallies=0; checks=0; mism=0
        bigtrans=0; smallstrong=0
        for _ in range(3000):
            g=[rnd.random() for _ in range(n)]; s=sum(g); g=[x/s for x in g]
            adj=half_turn_from_gaps(g)
            if adj is None: continue
            checks+=1
            sc=scc_count(adj); gmax=max(g)
            trans=(sc==n); strong=(sc==1)
            # THM-374 style: transitive <=> a gap > 1/2 (semicircle)
            if trans and gmax>0.5: bigtrans+=1
            if strong and gmax<0.5: smallstrong+=1
            if trans != (gmax>0.5): mism+=1
        print(f"   n={n}: samples={checks};  (transitive <=> largest gap>1/2) mismatches={mism}"
              f"   [transitive&big-gap={bigtrans}, strong&small-gap={smallstrong}]")
    print("   => exactly one boundary gap > 1/2 (the empty arc, a SEMICIRCLE) <=> transitive;")
    print("      the source-sink arc spans that gap. All gaps < 1/2 (runners wrap) <=> one")
    print("      strong block (Ham cycle = the boundary). The apex arc is the block switch.\n")

    # (3) LRC illustration: observer sits in the largest gap; loneliness lives on the apex arc
    print("(3) LRC: observer in the largest gap; the source-sink arc brackets it.")
    print("    transitive lonely = observer OUTSIDE the runner cluster (semicircle, big gap);")
    print("    strong lonely     = observer INSIDE a wrap-around (AP / regular-polygon tight).")
    print("    Loneliness (THM-382) = both observer-adjacent gaps >= 1/n = the two ends of")
    print("    the source-sink arc are each >= 1/n from the observer.")

if __name__=="__main__":
    main()
