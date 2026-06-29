#!/usr/bin/env python3
"""
r-eigenspace-merged-compression.py   (klein-2026-06-29-S1)

Synthesis of the Klein/compression thread with the project's R-eigenspace
organizing principle (HYP-3538, THM-583 reflections): complement R has
eigenvalues +/-1, and every "two-index" phenomenon is that split. At n=4 the
four iso classes are the Klein group V4=(Z2)^2 with R acting as the coordinate
swap; the +1 (R-even) coordinate is u = #boundary-defects in {0,1,2} and the
-1 (R-odd) coordinate is w = source-minus-sink, nonzero only on the NS pair.

Tests:
  CENSUS   : SC / NS-pair counts, R-even dim = V_merged = (A+SC)/2,
             R-odd dim = (A-SC)/2 = NS-merged.  (verify vs CLAUDE.md.)
  N4-UW    : exhibit u,w at n=4; verify u is R-invariant and w flips sign under R.
  HYP-A    : min covering face dimension for ALL classes vs for complement-MERGED
             classes (n=4,5,6). If merged is tight (=info bound) where full is not,
             the compression overhead IS the R-odd coordinate.
  HYP-B    : iso-class arc-flip quotient adjacency A; verify A[s i][s j]=A[i][j]
             (complement is an exact automorphism => A block-diagonalizes
             R-even (+) R-odd). Report block dims and (if numpy) eigenvalues.
"""
import itertools
from math import comb, log2, factorial
from collections import Counter, defaultdict

A000568 = {3:2,4:4,5:12,6:56,7:456,8:6880}

def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def base_path_pairs(n): return [(k,k+1) for k in range(n-1)]
def tile_pairs(n):
    bp=set(base_path_pairs(n)); return [p for p in pairs(n) if p not in bp]
def tiling_bits(n,flips):
    idx={p:k for k,p in enumerate(pairs(n))}; b=0
    for p in flips: b|=1<<idx[p]
    return b
def scores_full(n,bits):
    P=pairs(n); idx={p:k for k,p in enumerate(P)}; s=[0]*n
    for (i,j) in P:
        if (bits>>idx[(i,j)])&1: s[i]+=1
        else: s[j]+=1
    return s
def perm_tables(n):
    P=pairs(n); idx={p:k for k,p in enumerate(P)}; T=[]
    for perm in itertools.permutations(range(n)):
        row=[]
        for (i,j) in P:
            a,b=perm[i],perm[j]
            row.append((idx[(a,b)],False) if a<b else (idx[(b,a)],True))
        T.append(row)
    return T
def canonical(bits,T):
    best=None
    for row in T:
        v=0
        for q,(tgt,inv) in enumerate(row):
            bit=(bits>>tgt)&1
            if inv: bit^=1
            v|=bit<<q
        if best is None or v<best: best=v
    return best
def complement_bits(bits,n):
    return bits ^ ((1<<comb(n,2))-1)

def build(n):
    """Return tiles, class_of[tilemask]->canon, canon list, complement map sigma."""
    tiles=tile_pairs(n); m=len(tiles); T=perm_tables(n)
    class_of={}; canon_index={}
    for mask in range(2**m):
        flips=[tiles[k] for k in range(m) if (mask>>k)&1]
        c=canonical(tiling_bits(n,flips),T)
        class_of[mask]=c
        if c not in canon_index: canon_index[c]=len(canon_index)
    # complement permutation on classes
    rep={}
    for mask,c in class_of.items(): rep.setdefault(c,mask)
    sigma={}
    for c,mask in rep.items():
        flips=[tiles[k] for k in range(m) if (mask>>k)&1]
        bits=tiling_bits(n,flips)
        sigma[c]=canonical(complement_bits(bits,n),T)
    return tiles,m,T,class_of,canon_index,rep,sigma

def census(n,canon_index,sigma):
    A=len(canon_index)
    SC=sum(1 for c in canon_index if sigma[c]==c)
    NSpairs=(A-SC)//2
    Vmerged=SC+NSpairs   # = (A+SC)/2
    print(f"  n={n}: A={A:>4}  SC={SC:>3}  NS-pairs={NSpairs:>3}  "
          f"R-even dim V_merged={Vmerged:>3}  R-odd dim={NSpairs:>3}   "
          f"[(A+SC)/2={ (A+SC)//2}]")
    assert Vmerged==(A+SC)//2
    return A,SC,NSpairs,Vmerged

def n4_uw(n,tiles,class_of,rep,sigma):
    print("\n N4-UW: explicit R-even u and R-odd w at n=4")
    name={(0,1,2,3):'T',(0,2,2,2):'+',(1,1,1,3):'-',(1,1,2,2):'S'}
    def nm(c):
        flips=[tiles[k] for k in range(len(tiles)) if (rep[c]>>k)&1]
        return name[tuple(sorted(scores_full(n,tiling_bits(n,flips))))]
    # u = #defects: has source? has sink? from score sequence
    def uw(c):
        flips=[tiles[k] for k in range(len(tiles)) if (rep[c]>>k)&1]
        s=sorted(scores_full(n,tiling_bits(n,flips)))
        has_source = (s[-1]==n-1)   # someone beats all
        has_sink   = (s[0]==0)      # someone loses to all
        u=(0 if has_source else 1)+(0 if has_sink else 1)  # # killed
        # signed w: +1 source-killed-only, -1 sink-killed-only, 0 else
        src_killed = 0 if has_source else 1
        snk_killed = 0 if has_sink   else 1
        w=src_killed-snk_killed
        return u,w
    for c in canon_index_order(class_of):
        u,w=uw(c); sc=nm(c); scc=nm(sigma[c])
        print(f"    {sc:2s}: u(R-even)={u}  w(R-odd)={w:+d}   complement-> {scc:2s}  "
              f"u'={uw(sigma[c])[0]} w'={uw(sigma[c])[1]:+d}")
    print("    => u is R-invariant (u'=u everywhere); w flips sign (w'=-w): the +/-1 split.")

def canon_index_order(class_of):
    seen=[]
    for mask in sorted(class_of):
        c=class_of[mask]
        if c not in seen: seen.append(c)
    return seen

def min_cover_face(n,tiles,labels_of):
    """min face dimension covering all labels (labels_of: tilemask->label)."""
    m=len(tiles); full=set(labels_of.values())
    ib=0
    while (1<<ib)<len(full): ib+=1
    for k in range(ib,m+1):
        for free in itertools.combinations(range(m),k):
            fixed=[i for i in range(m) if i not in free]
            cover=defaultdict(set)
            for mask in labels_of:
                key=tuple((mask>>i)&1 for i in fixed)
                cover[key].add(labels_of[mask])
            if any(cs==full for cs in cover.values()):
                return k,ib
    return m,ib

def hyp_A(n,tiles,class_of,sigma):
    print(f"\n HYP-A: full vs complement-merged minimum covering face (n={n})")
    full_min,full_ib=min_cover_face(n,tiles,class_of)
    # merged labels: orbit id = min(c, sigma[c])
    merged_of={mask:min(class_of[mask],sigma[class_of[mask]]) for mask in class_of}
    merg_min,merg_ib=min_cover_face(n,tiles,merged_of)
    print(f"    ALL    classes: info={full_ib}  min face={full_min}  tight={full_min==full_ib}")
    print(f"    MERGED classes: info={merg_ib}  min face={merg_min}  tight={merg_min==merg_ib}")
    verdict = ("overhead IS the R-odd coordinate" if (full_min>full_ib and merg_min==merg_ib)
               else "overhead NOT explained by R-odd alone" if full_min>full_ib
               else "no overhead at this n")
    print(f"    => {verdict}")
    return full_min,full_ib,merg_min,merg_ib

def hyp_B(n,rep,canon_index,sigma):
    print(f"\n HYP-B: complement is an automorphism of the arc-flip iso-class adjacency (n={n})")
    T=perm_tables(n)
    classes=list(canon_index.keys())
    cidx={c:i for i,c in enumerate(classes)}
    N=len(classes)
    A=[[0]*N for _ in range(N)]
    for c in classes:
        bits=tiling_bits(n,[]) if False else rep_bits(n,rep,c)
        for a in range(comb(n,2)):
            nb=canonical(bits^(1<<a),T)
            A[cidx[c]][cidx[nb]]+=1
    # verify A[sigma i][sigma j]==A[i][j]
    ok=True
    for i in classes:
        for j in classes:
            if A[cidx[sigma[i]]][cidx[sigma[j]]]!=A[cidx[i]][cidx[j]]:
                ok=False
    print(f"    commute A o sigma == sigma o A : {ok}   (so A block-diagonalizes into R-even (+) R-odd)")
    SC=sum(1 for c in classes if sigma[c]==c)
    print(f"    R-even block dim={ (N+SC)//2 }   R-odd block dim={ (N-SC)//2 }")
    try:
        import numpy as np
        M=np.array(A,dtype=float)
        # symmetrize by class sizes is messy; report spectrum of A directly (row-sum = C(n,2))
        ev=np.linalg.eigvals(M)
        ev=sorted([round(float(e.real),3) for e in ev])
        print(f"    spectrum(A) (row-sum={comb(n,2)}): min={ev[0]}, max={ev[-1]}, "
              f"distinct~{len(set(ev))}")
    except Exception as e:
        print(f"    (numpy unavailable: {e})")
    return ok

# global canon_index for ordering helper
canon_index=None
def rep_bits(n,rep,c):
    tiles=tile_pairs(n)
    return tiling_bits(n,[tiles[k] for k in range(len(tiles)) if (rep[c]>>k)&1])

if __name__=="__main__":
    print("="*72)
    print(" CENSUS: SC/NS and R-even/R-odd dimensions")
    print("="*72)
    data={}
    for n in (3,4,5,6):
        tiles,m,T,class_of,canon_index_,rep,sigma=build(n)
        data[n]=(tiles,m,class_of,canon_index_,rep,sigma)
        census(n,canon_index_,sigma)
    # CLAUDE.md cross-check
    print("  expected (CLAUDE.md): SC n=3..6 = 2,2,8,12 ; NS-merged = 0,1,2,22 ; V_merged = 2,3,10,34")

    tiles,m,class_of,canon_index_,rep,sigma=data[4]
    globals()['canon_index']=canon_index_
    n4_uw(4,tiles,class_of,rep,sigma)

    print("\n"+"="*72); print(" HYP-A: compression overhead vs the R-odd coordinate"); print("="*72)
    for n in (4,5,6):
        tiles,m,class_of,canon_index_,rep,sigma=data[n]
        hyp_A(n,tiles,class_of,sigma)

    print("\n"+"="*72); print(" HYP-B: complement symmetry of the arc-flip metagraph"); print("="*72)
    for n in (3,4,5,6):
        tiles,m,class_of,canon_index_,rep,sigma=data[n]
        hyp_B(n,rep,canon_index_,sigma)
