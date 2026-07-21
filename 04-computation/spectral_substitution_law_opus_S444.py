#!/usr/bin/env python3
"""
THE SPECTRAL SUBSTITUTION LAW (opus-S444). The skew matrix of the uniform substitution T[S^m]
(blow up each of T's k vertices by a copy of S, |S|=m) is  S_T (x) J_m + I_k (x) S_S  (Kronecker),
since between-block arcs follow T (all +1/-1 = sign times J_m) and within-block arcs are S_S.
J_m has eigenvalues m (once) and 0 (m-1 times), so the skew spectrum SPLITS:

   nonzero-spec(T[S^m]) = [ nonzero-spec(S), each with multiplicity k ]
                          U [ m * nonzero-spec(T) ]        (block-graph scaled by block size)

and the zero eigenvalues collapse to the parity-forced count ([n odd] => one 0). Verify on many
(T,S) pairs; this is the tournament/skew analogue of the lexicographic-product spectrum, and it
makes char_S of any substitution tree FACTOR over the seeds.
"""
import numpy as np, itertools

def skew(adj,n):
    return np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])
def substitute(T,k,S,m):
    N=k*m; big=[[0]*N for _ in range(N)]
    for i in range(k):
        for a in range(m):
            for c in range(m):
                if a!=c: big[i*m+a][i*m+c]=S[a][c]
        for j in range(k):
            if i!=j and T[i][j]:
                for a in range(m):
                    for c in range(m): big[i*m+a][j*m+c]=1
    return big,N
def nz_l2(M):
    n=len(M)
    S=np.array([[0.0 if i==j else (1.0 if M[i][j] else -1.0) for j in range(n)] for i in range(n)])
    ev=np.linalg.eigvals(S)
    return sorted(round(abs(e)**2,4) for e in ev if abs(e)>1e-7)

C3=[[0,1,0],[0,0,1],[1,0,0]]
T3=[[0,1,1],[0,0,1],[0,0,0]]                # transitive 3
C5=[[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1],[1,1,0,0,0]]  # rotational C5 (regular)
def rand_reg_or_any(seed_adj): return seed_adj

print("SPECTRAL SUBSTITUTION LAW: nz-spec(T[S^m]) = [nz-spec(S) x k] U [m * nz-spec(T)] ?")
print("="*70)
tests=[("C3","C3",C3,C3),("C3","T3",C3,T3),("T3","C3",T3,C3),("C3","C5",C3,C5),
       ("C5","C3",C5,C3),("T3","T3",T3,T3),("C5","C5",C5,C5)]
for Tn,Sn,T,S in tests:
    k=len(T); m=len(S)
    big,N=substitute(T,k,S,m)
    actual=nz_l2(big)
    # predicted:
    seed=nz_l2(S)*k                                        # nz-spec(S) each mult k (as lambda^2 list, repeated k times)
    block=[round(m*m*v,4) for v in nz_l2(T)]               # m*eigenvalue => lambda^2 scales by m^2
    pred=sorted(seed+block)
    ok = (actual==pred)
    print(f"  T={Tn}[{Sn}^{m}] (N={N}): law holds={ok}")
    if not ok:
        print(f"     actual={actual}\n     pred  ={pred}")

# consequence: char_S of a substitution tree factors over seeds -> closed-form tr(S^8) for octonion tests
print("\nCONSEQUENCE: C3[C3] (n=9) spectrum = seed sqrt3 (mult 3 pairs) + block 3*sqrt3 (=sqrt27) --")
big,N=substitute(C3,3,C3,3)
print(f"  nz-lambda^2(C3[C3]) = {nz_l2(big)}  = [3]*6 (seed, k=3) U [27]*2 (block m^2*3=9*3)")
print("  => tr(S^8) and every even moment of any substitution object is a CLOSED FORM in the seed")
print("     moments -- giving exact octonion-wall (degree-8) test objects at n=9,15,21,...")
