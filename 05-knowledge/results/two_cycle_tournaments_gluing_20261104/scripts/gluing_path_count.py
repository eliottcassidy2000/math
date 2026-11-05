# path count of a vertex gluing Y = A (+)_s B from the pieces: five block patterns (see packet README / PROOF_CONVERSE §8).
#!/usr/bin/env python3
"""gluing_path_count.py (klein-2026-11-05): checks the five-term formula for the number of Hamiltonian paths of a vertex gluing
Y = A (+)_s B against every vertex-gluing two-cycle class of data/cls_{n}.txt (README section 10, Proposition).  Terms:
P1 = E_A s_{b1}(Z); P2 = h(Z) h(A) - e_{b1}(Z) (s_s(A) + s_{ak}(A)); P3 = Rend_A s_{b1}(Z); P4 = E_A R_{b1}(Z); P5 = the
Z1 X1 Z2 X2 two-path-cover sum; with X = V(A), Z = V(B) - s, E_A = e_s(A) + e_{ak}(A), Rend_A = sum over nonempty proper
X1 of (paths of A[X1] ending at s or a_k) h(A[X - X1]), R_{b1}(Z) = sum over nonempty Z1 not containing b1 of h(Z[Z1]) s_{b1}(Z[Z - Z1]).
Usage: gluing_path_count.py N1 N2."""
import sys, os, itertools
from cutvertex_lib import *
HERE = os.path.dirname(os.path.abspath(__file__))
def hp(rows, verts):
    """(h, dict start->count, dict end->count, dict (s,e)->count) for the sub-tournament on verts (verts nonempty)"""
    tot,SE=hpaths(rows,list(verts)); st={}; en={}
    for (s,e),c in SE.items(): st[s]=st.get(s,0)+c; en[e]=en.get(e,0)+c
    return tot,st,en,SE
def runs_of(rows,n):
    C=cycles(rows,n); A=[set((c[i],c[(i+1)%n]) for i in range(n)) for c in C]
    c=C[0]; diff=A[0]-A[1]; chg=sorted(i for i in range(n) if (c[i],c[(i+1)%n]) in diff)
    runs=[]
    for j in range(3):
        a=chg[j]+1; b=chg[(j+1)%3]; L=(b-a)%n+1; runs.append([c[(a+t)%n] for t in range(L)])
    return runs
D = os.path.join(HERE, '..', 'data') + '/'
for n in range(int(sys.argv[1]), int(sys.argv[2])+1):
    tested=bad=0
    for l in open(D+'cls_%d.txt'%n):
        h,d6=l.split()
        if h!='2': continue
        rows=rows_from_d6(d6); runs=runs_of(rows,n); sizes=[len(r) for r in runs]
        if 1 not in sizes: continue
        j=sizes.index(1); s=runs[j][0]; S2=runs[(j+1)%3]; PR=runs[(j+2)%3]
        X=[s]+S2; Z=PR; a_k=S2[-1]; b_1=PR[0]
        hY,_,_,_=hp(rows,range(n))
        hA,stA,enA,_=hp(rows,X); hZ,stZ,enZ,SEZ=hp(rows,Z)
        E_A=enA.get(s,0)+enA.get(a_k,0)
        def E(W):   # paths of A[W] ending at s or a_k
            if not W: return 0
            if len(W)==1: return 1 if W[0] in (s,a_k) else 0
            _,_,en,_=hp(rows,W); return en.get(s,0)+en.get(a_k,0)
        def H(W): return 1 if len(W)==1 else (hp(rows,W)[0] if W else 1)
        def S_excl(W):  # paths of A[W] starting at a vertex other than s, a_k (for a predecessor b_1)
            if len(W)==1: return 0 if W[0] in (s,a_k) else 1
            _,st,_,_=hp(rows,W); return sum(c for v,c in st.items() if v not in (s,a_k))
        P1=E_A*stZ.get(b_1,0)
        P2=hZ*hA-enZ.get(b_1,0)*(stA.get(s,0)+stA.get(a_k,0))
        # P3: X1 Z X2 with X1, X2 nonempty
        Rend=0
        for r in range(1,len(X)):
            for X1 in itertools.combinations(X,r):
                X2=[v for v in X if v not in X1]; Rend+=E(list(X1))*H(X2)
        P3=Rend*stZ.get(b_1,0)
        # P4: Z1 X Z2 (b_1 in Z2, Z1 nonempty): E_A * sum h(Z[Z1]) s_{b_1}(Z[Z2])
        R=0
        Zm=[v for v in Z if v!=b_1]
        for r in range(1,len(Zm)+1):
            for Z1 in itertools.combinations(Zm,r):
                Z2=[v for v in Z if v not in Z1]
                _,st2,_,_=hp(rows,Z2) if len(Z2)>1 else (1,{Z2[0]:1},{Z2[0]:1},None)
                R+=H(list(Z1))*st2.get(b_1,0)
        P4=E_A*R
        # P5: Z1 X1 Z2 X2: Z1 nonempty, b_1 in Z2, X1,X2 nonempty
        P5=0
        for r in range(1,len(Zm)+1):
            for Z1 in itertools.combinations(Zm,r):
                Z2=[v for v in Z if v not in Z1]; hZ1=H(list(Z1))
                if len(Z2)==1: ends={b_1:1}
                else:
                    _,_,_,SE2=hp(rows,Z2); ends={}
                    for (a,b),c in SE2.items():
                        if a==b_1: ends[b]=ends.get(b,0)+c
                inner=0
                for r2 in range(1,len(X)):
                    for X1 in itertools.combinations(X,r2):
                        X2=[v for v in X if v not in X1]; e1=E(list(X1))
                        if not e1: continue
                        for bt,c in ends.items():
                            inner+=e1*c*(H(X2) if bt!=b_1 else S_excl(X2))
                P5+=hZ1*inner
        tested+=1
        if hY!=P1+P2+P3+P4+P5: bad+=1; print("MISMATCH", n, d6, hY, (P1,P2,P3,P4,P5))
    print("order %d: vertex gluings tested %d, formula failures %d"%(n,tested,bad), flush=True)
