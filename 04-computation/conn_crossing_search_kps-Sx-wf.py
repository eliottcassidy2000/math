"""
Exploratory: which orientation-weighted CROSSING functional of the circular
tournament drawing equals a repo OCF datum (c3, alpha_2, H) exactly?

For each 4-subset a<b<c<d on the circle, the convex-position chords that
CROSS are (a,c) and (b,d). We tabulate the orientation of all 6 arcs among
{a,b,c,d} and look for a per-quad weight w(orientation) such that
   sum_quads w  == c3   (or another datum)  for ALL tilings.

We also test the rectilinear/2-page crossing reading and the 'directed crossing'
where a crossing counts only if the two chords are head-to-head etc.
"""
from itertools import combinations, product
import random

def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]
def adj(n,bits,T):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): A[k][k-1]=1
    for (a,b),bit in zip(T,bits):
        if bit==0: A[a][b]=1
        else: A[b][a]=1
    return A
def c3(A,n):
    t=0
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            for k in range(j+1,n+1):
                if (A[i][j]+A[i][k],A[j][i]+A[j][k],A[k][i]+A[k][j])==(1,1,1):t+=1
    return t
def scores(A,n): return [sum(A[v][w] for w in range(1,n+1)) for v in range(1,n+1)]

def quad_type(A,a,b,c,d):
    # encode orientation of the two crossing chords (a,c),(b,d) only -> 4 types
    s1 = A[a][c]      # 1 if a->c
    s2 = A[b][d]      # 1 if b->d
    return (s1,s2)

def main():
    out=[]
    def pr(*a):
        s=" ".join(str(x) for x in a); print(s); out.append(s)

    pr("=== crossing functional search ===")
    # Hypothesis families, each a function of the crossing-chord orientations (s1,s2)
    # over all C(n,4) crossings. We fit: does sum equal c3? alpha-stuff?
    # weight tables indexed by (s1,s2):
    families = {
        "count_s1ne s2 (alt)":     {(0,0):0,(0,1):1,(1,0):1,(1,1):0},
        "count_s1eq s2 (parallel)":{(0,0):1,(0,1):0,(1,0):0,(1,1):1},
        "count_(1,0)":             {(0,0):0,(0,1):0,(1,0):1,(1,1):0},
        "count_(0,1)":             {(0,0):0,(0,1):1,(1,0):0,(1,1):0},
        "count_(1,1)":             {(0,0):0,(0,1):0,(1,0):0,(1,1):1},
        "count_(0,0)":             {(0,0):1,(0,1):0,(1,0):0,(1,1):0},
    }
    for n in range(4,8):
        T=tiles(n); F=len(T)
        if F<=12:
            sample=list(product([0,1],repeat=F))
        else:
            random.seed(1); sample=[tuple(random.randint(0,1) for _ in range(F)) for _ in range(3000)]
        # check exact equality of each family with c3 across whole sample
        eqc3={k:True for k in families}
        # also: is each family a fixed AFFINE function of c3? track (val,c3) pairs span
        pairs={k:set() for k in families}
        for bits in sample:
            A=adj(n,bits,T)
            cc=c3(A,n)
            vals={k:0 for k in families}
            for (a,b,c,d) in combinations(range(1,n+1),4):
                t=quad_type(A,a,b,c,d)
                for k in families: vals[k]+=families[k][t]
            for k in families:
                if vals[k]!=cc: eqc3[k]=False
                pairs[k].add((vals[k],cc))
        pr(f"n={n} F={F} sample={len(sample)}")
        for k in families:
            # is val an exact affine function of c3? fit a*c3+b from 2 points
            P=sorted(pairs[k])
            affine="?"
            if len({c for v,c in P})>=2:
                (v0,c0),(v1,c1)=P[0],P[-1]
                if c1!=c0:
                    a_=(v1-v0)/(c1-c0)
                    b_=v0-a_*c0
                    ok=all(abs(v-(a_*c+b_))<1e-9 for v,c in P)
                    affine=f"val={a_:.3g}*c3+{b_:.3g}" if ok else "NOT affine in c3"
            pr(f"   {k:28s}: exact==c3? {eqc3[k]!s:5s}  {affine}")
    pr("")
    pr("Now test the DEGREE / quadratic-form structure of the best crossing count.")
    pr("c3 is a quadratic form in tile bits (THM-554/HYP-2707). We test whether")
    pr("the 'parallel crossing' count P (s1==s2) is ALSO a quadratic form, i.e.")
    pr("P(bits) - P(0) is degree<=2 in bits, by finite differences.")
    for n in [5,6]:
        T=tiles(n); F=len(T)
        def Pcount(bits):
            A=adj(n,bits,T); tot=0
            for (a,b,c,d) in combinations(range(1,n+1),4):
                s1=A[a][c]; s2=A[b][d]
                if s1==s2: tot+=1
            return tot
        # third finite difference should vanish for a quadratic
        base=[0]*F
        maxdeg=0
        random.seed(3)
        for _ in range(200):
            # pick up to 3 distinct coords, test 3rd mixed difference
            ids=random.sample(range(F),3)
            def ev(flips):
                b=base[:]
                for i in flips: b[i]^=1
                return Pcount(b)
            # mixed 3rd difference over coords ids[0],1,2
            s=0
            for mask in range(8):
                sub=[ids[j] for j in range(3) if mask&(1<<j)]
                sign=(-1)**bin(mask).count("1")
                s+=sign*ev(sub)
            if s!=0: maxdeg=max(maxdeg,3)
        pr(f"   n={n}: parallel-crossing 3rd-difference always 0? {'YES (deg<=2, QUADRATIC like c3)' if maxdeg<3 else 'NO (deg>=3)'}")
    with open("05-knowledge/results/conn_crossing_search_kps-Sx-wf.out","w") as f:
        f.write("\n".join(out)+"\n")

if __name__=="__main__":
    main()
