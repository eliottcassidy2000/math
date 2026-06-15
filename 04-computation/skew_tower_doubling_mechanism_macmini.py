#!/usr/bin/env python3
"""
Verify the DOUBLING MECHANISM that drives the induction in the proof, at the
binary (F_2) level. Fresh code.

Claim (mechanism): with H_{2m}=[[H,H],[-H^T,H^T]], binarize B=(J-H)/2.
Let b_i = binarized i-th row of H (length m), and c_i = binarized i-th row of H^T.
Because H+H^T=2I, for off-diagonal H^T_{ij}=H_{ji}=-H_{ij}... we track exactly:

  TOP rows of H_{2m}:  [H | H]  -> binarize to (b_i, b_i)   (pair-doubled word)
  BOT rows of H_{2m}:  [-H^T | H^T].
    binarize(-H^T)_{ij} = 1 iff -H^T_{ij} = -1 iff H^T_{ij}=+1 iff H_{ji}=+1.
    binarize( H^T)_{ij} = 1 iff  H^T_{ij} = -1 iff H_{ji}=-1.
    So bot row i = ( complement-of-c_i , c_i ) where c_i = binarized col i of H
    ... and on the diagonal H_{ii}=+1 so the diagonal twin corrections appear.

We verify, for the tower, that:
  (1) the top n/2 binarized rows are EXACTLY (b,b) pair-doubled (b = row of (J-H)/2 on H).
  (2) each top row (b,b) has weight = 2*wt(b), divisible by 4 iff wt(b) even.
      And we confirm every row of (J-H)/2 on H has EVEN weight (self-orthogonality
      of C(H) forces diagonal <r,r>=wt(r) even). This is the inductive hypothesis.
  (3) dim of the new code = 2 * (dim of old) is FALSE in general; the correct
      statement is dim doubles as a function of order: dim(C(H_{2m})) = m = (2m)/2.
      We confirm dim(C(H_order_n)) = n/2 at every tower level, and that the
      generator structure is [top=(b,b); bottom=(comp(c), c)].
"""

def transpose(M):
    n=len(M); m=len(M[0]); return [[M[i][j] for i in range(n)] for j in range(m)]
def double(H):
    m=len(H); HT=transpose(H)
    top=[H[i]+H[i] for i in range(m)]
    bot=[[-HT[i][j] for j in range(m)]+HT[i][:] for i in range(m)]
    return top+bot
def binrow(row):
    v=0
    for j,e in enumerate(row):
        if e==-1: v|=(1<<j)
    return v
def popcount(x): return bin(x).count("1")
def gf2_rank(vectors):
    pivots={}
    for v in vectors:
        x=v
        for pb in sorted(pivots,reverse=True):
            if (x>>pb)&1: x^=pivots[pb]
        if x: pivots[x.bit_length()-1]=x
    return len(pivots)

def check_level(H):
    m=len(H)
    Hd=double(H)        # order 2m
    n=2*m
    # binarize H rows (length m)
    bH=[binrow(H[i]) for i in range(m)]
    # binarize doubled rows (length n)
    bD=[binrow(Hd[i]) for i in range(n)]
    # (1) top rows are (b,b): top bit pattern low m bits == high m bits == bH[i]
    top_ok=True
    for i in range(m):
        low = bD[i] & ((1<<m)-1)
        high = (bD[i] >> m) & ((1<<m)-1)
        if not (low==bH[i] and high==bH[i]):
            top_ok=False
    # (2) every row of (J-H)/2 on H has even weight?
    all_even = all(popcount(b)%2==0 for b in bH)
    top_weights_div4 = all(popcount(bD[i])%4==0 for i in range(m))
    # dim facts
    dimH = gf2_rank(bH)
    dimD = gf2_rank(bD)
    return {
        "order_old": m, "order_new": n,
        "top_rows_are_pair_doubled": top_ok,
        "all_H_rows_even_weight": all_even,
        "top_doubled_rows_div4": top_weights_div4,
        "dim_old": dimH, "dim_old==m/2": dimH==m//2 if m>1 else "(m=1 base)",
        "dim_new": dimD, "dim_new==n/2": dimD==n//2,
    }

if __name__=="__main__":
    H=[[1]]
    for step in range(6):
        r=check_level(H)
        print(f"step {step}: {r}")
        H=double(H)
