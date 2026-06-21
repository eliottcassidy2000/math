"""
Nail the closed form of the parallel-crossing count P of the circular drawing.
P = #{4-subsets a<b<c<d : crossing chords (a,c),(b,d) have SAME winding sign}.

Claims to test (exact, full enumeration n<=6, big sample n=7):
 (A) P is a quadratic form in the tile bits  [already VERIFIED -> recheck]
 (B) P is determined by the SCORE SEQUENCE alone?  (like c3 is)  TEST.
 (C) Closed form P = f(scores)?  Fit against C(s_v,2), inversions, etc.
 (D) Relation to c3: is P - c3 score-determined?
 (E) The 'directed crossing number' dcr = #crossings that are part of a directed
     4-cycle's two diagonals -> tie to alpha_2 / C4.
Also: report cr(K_n) convex = C(n,4) vs Guy rectilinear cr(K_n) (the MIN).
"""
from itertools import combinations, product
import random
from fractions import Fraction

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
def scores(A,n): return tuple(sum(A[v][w] for w in range(1,n+1)) for v in range(1,n+1))

def Pcount(A,n):
    tot=0
    for (a,b,c,d) in combinations(range(1,n+1),4):
        if A[a][c]==A[b][d]: tot+=1
    return tot

def all_crossings_full(A,n):
    """count crossing types using full convex order.
    For 4-subset w<x<y<z the crossing pair is (w,y),(x,z). Also count whether the
    4-set carries a directed 4-cycle."""
    Pc=0; cyc4=0
    for (w,x,y,z) in combinations(range(1,n+1),4):
        if A[w][y]==A[x][z]: Pc+=1
        sub=[w,x,y,z]
        for perm in [(0,1,2,3),(0,1,3,2),(0,2,1,3)]:
            p=[sub[i] for i in perm]
            if all(A[p[i]][p[(i+1)%4]] for i in range(4)) or all(A[p[(i+1)%4]][p[i]] for i in range(4)):
                cyc4+=1; break
    return Pc,cyc4

def main():
    out=[]
    def pr(*a):
        s=" ".join(str(x) for x in a); print(s); out.append(s)
    pr("=== closed form of circular parallel-crossing count P ===")
    # cr(K_n) Guy (rectilinear/optimal) known small values:
    guy={4:0,5:1,6:3,7:9,8:18}
    pr("convex crossing = C(n,4); Guy optimal cr(K_n) =", guy)
    for n in range(4,8):
        pr(f"--- n={n}, C(n,4)={len(list(combinations(range(n),4)))}, Guy cr={guy.get(n)} ---")
        T=tiles(n); F=len(T)
        if F<=12: sample=list(product([0,1],repeat=F))
        else:
            random.seed(2); sample=[tuple(random.randint(0,1) for _ in range(F)) for _ in range(5000)]
        # (B) is P score-determined? map score-seq(sorted? NO, ordered) -> set of P
        score_to_P={}
        c3_minus_P=set()
        P_minus_c3_scoredet=True
        score_P_ordered={}
        for bits in sample:
            A=adj(n,bits,T)
            P=Pcount(A,n); cc=c3(A,n); sc=scores(A,n)
            key=sc  # ordered score by vertex label
            score_to_P.setdefault(key,set()).add(P)
            score_P_ordered.setdefault(key,set()).add((P,cc))
        scoredet=all(len(v)==1 for v in score_to_P.values())
        pr(f"   (B) P determined by ORDERED score seq? {scoredet}")
        # P - c3 determined by score?
        diff_scoredet=all(len({P-cc for (P,cc) in v})==1 for v in score_P_ordered.values())
        pr(f"   (D) P - c3 determined by score seq? {diff_scoredet}")
        # (C) fit P = A0 + sum over vertices g(s_v) + winding term?
        # Try P = C(n,4)/2-ish + quadratic in scores. Collect (scores,P) and least-squares
        # over basis [1, sum C(s_v,2), sum s_v^2 with position weights...]. Keep simple:
        # we already know c3 = C(n,3)-sum C(s_v,2). Test if P is affine in (c3, n).
        # gather distinct (c3,P)
        cp=set()
        for bits in sample:
            A=adj(n,bits,T); cp.add((c3(A,n),Pcount(A,n)))
        # affine fit P=a*c3+b
        cps=sorted(cp)
        affine="(insufficient pts)"
        cs={c for c,p in cps}
        if len(cs)>=2:
            (c0,p0),(c1,p1)=cps[0],cps[-1]
            if c1!=c0:
                a_=Fraction(p1-p0,c1-c0); b_=Fraction(p0)-a_*c0
                ok=all(Fraction(p)==a_*c+b_ for c,p in cps)
                affine=f"P={a_}*c3+{b_}  exact:{ok}"
        pr(f"   (D') P affine in c3? {affine}")
        # show range
        Ps=sorted({p for c,p in cp})
        pr(f"   P range: [{Ps[0]},{Ps[-1]}], #distinct={len(Ps)}")
    pr("")
    pr("INTERPRETATION:")
    pr(" - raw circular crossing number = C(n,4) (orientation-free) -- TRIVIAL. PROVED.")
    pr(" - parallel-crossing P is a QUADRATIC form in tile bits (Clifford level), same")
    pr("   tier as c3 (THM-554 / HYP-2707 quadratic/Gauss-sum side of the wall).")
    pr(" - Whether P is score-determined is the decisive test above.")
    with open("05-knowledge/results/conn_crossing_closedform_kps-Sx-wf.out","w") as f:
        f.write("\n".join(out)+"\n")

if __name__=="__main__":
    main()
