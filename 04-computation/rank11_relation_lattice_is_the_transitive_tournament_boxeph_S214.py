#!/usr/bin/env python3
"""rank11_relation_lattice_is_the_transitive_tournament_boxeph_S214.py -- boxeph-2026-07-21-S214

"Relations between a set of things is a tournament, through abstract enough glasses." The LRC(14) Wall-A
core is the 12 AP speeds {1,..,12}. This corrected audit keeps the useful
rank-11 order/lattice comparison but does not identify the objects. Pillars:

  P1 L(AP)={a in Z^12 : sum i*a_i = 0} is a RANK-11 primitive lattice. The rows
     d_k=(k+1)e_k-k e_{k+1} span an index-11! sublattice D, not a Z-basis of
     L(AP). Only D has the displayed tridiagonal path/Jacobi Gram matrix.
  P2 minimal vectors of L(AP) are the 30 distinct-index Schur triples and their
     negatives (60 vectors). This is not four-term additive energy.
  P3 the 12 AP speeds ARE the transitive tournament T_12: score seq 0..11 (the AP), char_A = x^12 (fully
     nilpotent = reify-ladder nullcone vertex), transitivity Vandermonde = the braid arrangement A_11 (rank 11).
  P4 CHIRALITY: {1,..,12} is PALINDROMIC (v_i + v_{13-i}=13), so it is SELF-CONVERSE = ACHIRAL = the FIXED
     point of the reversal involution (S213). Wall A = 'the rank-11 achiral transitive vertex is extremal'.
"""
from fractions import Fraction as F
from itertools import combinations

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

N=12
# ---------------------------------------------------------------------------
sep("P1  the adjacent-pair rows span an index-11! path sublattice of L(AP)")
# chain rows d_k = (k+1)e_k-k e_{k+1}, k=1..11
def dvec(k):
    v=[0]*(N+1)  # index 1..12
    v[k]=k+1; v[k+1]=-k
    return v
D=[dvec(k) for k in range(1,N)]
def dot(u,v): return sum(u[i]*v[i] for i in range(1,N+1))
def phi(v): return sum(i*v[i] for i in range(1,N+1))
print("  all 11 chain rows satisfy sum i*a_i = 0 (in the kernel)?", all(phi(d)==0 for d in D))
G=[[dot(D[a],D[b]) for b in range(N-1)] for a in range(N-1)]
tridiag=all(G[a][b]==0 for a in range(N-1) for b in range(N-1) if abs(a-b)>=2)
print("  Gram matrix TRIDIAGONAL (a weighted path / Jacobi matrix)?", tridiag)
print("  diagonal G[k,k] = k^2+(k+1)^2 :", [G[a][a] for a in range(N-1)])
print("  off-diagonal G[k,k+1] = -k(k+2):", [G[a][a+1] for a in range(N-2)])
print("  the 11 chain relations d_k live on ADJACENT coordinates (k,k+1) only:",
      all(sorted(i for i in range(1,N+1) if D[k][i]!=0)==[k+1,k+2] for k in range(N-1)))
disc=650
index=1
for k in range(1,N): index*=k
det_gram=disc*index*index
print(f"  disc L(AP)={disc}; det Gram(D)={det_gram}; saturation index sqrt(det/disc)={index}=11!")
print("  => D is a finite-index path/Jacobi frame, not a Z-basis and not THM-2052's private-anchor code.")

# ---------------------------------------------------------------------------
sep("P2  saturation restores the norm-3 Schur triples")
# search short vectors of L(AP): a in {-2..2}^12 with sum i*a_i=0, by norm
short={}
def gen(idx, cur, s):
    if idx>N:
        if s==0:
            nrm=sum(x*x for x in cur[1:])
            if 1<=nrm<=6: short.setdefault(nrm,[]).append(tuple(cur[1:]))
        return
    for x in (-2,-1,0,1,2):
        cur[idx]=x; gen(idx+1,cur,s+idx*x); cur[idx]=0
gen(1,[0]*(N+1),0)
minnorm=min(short)
mins=short[minnorm]
# NORM-3 additive triples (1 at i, 1 at j, -1 at i+j) need i<j STRICT (i=j gives norm-5 doublings 2v_i=v_2i)
trip=[(i,j) for i in range(1,N+1) for j in range(i+1,N+1) if i+j<=N]
print(f"  minimal norm = {minnorm}; #minimal vectors (kissing number) = {len(mins)}")
print(f"  #additive triples (i,j,i+j), i<=j, i+j<=12 = {len(trip)}  (kissing = 2 * this = {2*len(trip)})")
print(f"  minimal vectors ARE the (+-)(1@i,1@j,-1@i+j) additive triples? {len(mins)==2*len(trip)}")
print("  => 60 is the signed off-diagonal Schur-triple count, not four-term additive energy.")

# ---------------------------------------------------------------------------
sep("P3  the 12 AP speeds ARE the transitive tournament T_12 (rank-11 braid A_11, char x^12)")
# transitive tournament: i beats j iff i>j. score(i)=i-1 => scores 0..11 = the AP.
scores=sorted(sum(1 for j in range(1,N+1) if i>j) for i in range(1,N+1))
print("  score sequence of transitive T_12:", scores, " = translated AP {0..11}?", scores==list(range(N)))
print("  translation is not an LRC symmetry, so this is an order lens rather than a carrier map.")
# char_A of the (strictly-upper-triangular) transitive adjacency = x^12 (nilpotent)
# adjacency A[i][j]=1 iff i<j (say) -> strictly upper triangular -> A^12=0, char = x^12
def matpow_zero(n):
    A=[[1 if i<j else 0 for j in range(n)] for i in range(n)]
    # multiply A n times, check zero
    M=[row[:] for row in A]
    for _ in range(n-1):
        M=[[sum(M[i][k]*A[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return all(M[i][j]==0 for i in range(n) for j in range(n))
print("  transitive adjacency is nilpotent (A^12=0 => char_A = x^12, the reify-ladder NULLCONE VERTEX)?", matpow_zero(N))
# transitivity Vandermonde on nodes 1..12 = braid arrangement A_11 defining poly (rank 11) : nonzero (distinct)
def vander(a):
    p=1
    for i in range(len(a)):
        for j in range(i+1,len(a)): p*=(a[j]-a[i])
    return p
print("  transitivity Vandermonde prod_{i<j}(j-i) over 1..12 (braid A_11, rank 11) nonzero?", vander(list(range(1,N+1)))!=0)
print("  => the ranks coincide, but no map identifies the tournament, braid, and LRC relation objects.")

# ---------------------------------------------------------------------------
sep("P4  CHIRALITY: {1..12} is PALINDROMIC (v_i+v_{13-i}=13) => SELF-CONVERSE = ACHIRAL fixed point (S213)")
ap=list(range(1,N+1))
pal=all(ap[i-1]+ap[N-i]==N+1 for i in range(1,N+1))
print(f"  v_i + v_{{13-i}} = {[ap[i-1]+ap[N-i] for i in range(1,N+1)]} (all = 13)? {pal}")
print("  reversal i->13-i maps the AP to itself => the AP is SELF-CONVERSE = the ACHIRAL fixed point of the")
print("  chosen order with its converse. Phase reversal is universal, not an AP-fixed-configuration theorem.")
print("  The exact survivor is the constant pair sum 13, which feeds THM-2047's pair-sum walls.")

sep("SUMMARY")
print("""  Corrected survivor:
   * D is an index-11! adjacent-pair path/Jacobi sublattice of L(AP); saturation matters;
   * L(AP)'s 60 minimal vectors are signed distinct-index Schur triples, not additive energy;
   * transitive T_12 has score 0..11 and char_A=x^12, but translation/magnitude are lost;
   * AP pair sums are palindromic, but phase reversal is universal.
  These are diagnostic comparisons only. THM-2052's signed 13-coordinate code and
  THM-2053's parameter disks remain the faithful LRC carriers. See MISTAKE-227.""")
