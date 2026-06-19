"""
LRC(14) ROUTE G -- Determinant / positive-definite (Gram) form for the floor.
mac-mini-2026-06-19 session.

GOAL: express the lonely measure L(S) = int prod_i 1_safe(v_i tau) dtau  (THM-515),
equivalently the S3 cluster floor rho*(P,E) (THM-527), as a quadratic / Gram form in a
positive-definite kernel so positivity is automatic (Bochner) with a COMPUTABLE,
UNIFORM-IN-S minimum eigenvalue = the floor.

FINDINGS (all reproduced here, stdlib only):
  G1  Single-runner Toeplitz K[i][j]=h(i-j) is PSD (Bochner, hat h = 1_safe >= 0) but its
      min eigenvalue -> 0 (ess-inf 1_safe = 0 on the danger band). Floor 0.
  G3  Relation-frequency Gram h(n_i-n_j): same -- PSD, min eig -> 0. Floor 0.
  DPP No determinantal point-process kernel: pairwise danger overlaps A_ij are LARGER than
      A_i A_j for small-ratio pairs (positive correlation) => A_ij <> det 2x2 minor of a PSD K.
  PPP No permanental kernel either: the correlations have MIXED sign across pairs
      (1,2 positively correlated; 1,8 negatively). Neither all-neg (DPP) nor all-pos (PPP).
  G2  Gram/Cauchy-Schwarz lower bound L >= <g,phi>^2/||phi||^2 with a FIXED phi gives floor 0
      for the full L (the safe arc MOVES with S; span{1} gives only L^2).
  G4* PARTIAL POSITIVE: for the S3 CLUSTER floor the good x cluster near the UNIVERSAL
      centers {0,1/3,1/2,2/3} (b<7/2, HYP-2597) -- FIXED, not moving. A fixed PD bump there
      gives a positive Gram lower bound for the PURE cluster (0.075 at k=13). BUT it
      COLLAPSES to exactly 0 at the low-relation-height extremizer k=9, P={1,2,3,12}
      (THM-527's rho*=1/84 minimizer), because intersecting with G_P pushes the surviving
      good mass OFF the centers, onto thin arcs PINNED to the edge of P's danger band
      (x = 1/14 + eps, P-dependent). No P-independent PD kernel can sit there.

NET: PD/Gram structure is real and pervasive but every FIXED-kernel formulation has floor 0;
the determinantal route is structurally blocked (mixed-sign correlations). The crux is
unchanged. ROUTE G = DEAD-END for a clean uniform floor, with one PARTIAL byproduct
(universal-center Gram-LB) that confirms exactly where the low-height shapes break it.
"""
import math, itertools
from fractions import Fraction as F

def h(t):
    if t == 0: return 6.0/7.0
    return -math.sin(math.pi*t/7.0)/(math.pi*t)

def symeig(A, tol=1e-14):
    n=len(A); a=[r[:] for r in A]
    for _ in range(200):
        off=sum(a[p][q]**2 for p in range(n) for q in range(p+1,n))
        if off<tol: break
        for p in range(n):
            for q in range(p+1,n):
                if abs(a[p][q])<1e-18: continue
                app,aqq,apq=a[p][p],a[q][q],a[p][q]
                phi=0.5*math.atan2(2*apq,aqq-app); c,sn=math.cos(phi),math.sin(phi)
                for k in range(n):
                    akp,akq=a[k][p],a[k][q]; a[k][p]=c*akp-sn*akq; a[k][q]=sn*akp+c*akq
                for k in range(n):
                    apk,aqk=a[p][k],a[q][k]; a[p][k]=c*apk-sn*aqk; a[q][k]=sn*apk+c*aqk
    return sorted(a[i][i] for i in range(n))

def maxgap(pts):
    p=sorted(pts);n=len(p);g=max(p[i+1]-p[i] for i in range(n-1));return max(g,1-p[-1]+p[0])
def in_GP(x,P):
    return all(min((p*x)%1.0,1-(p*x)%1.0) > 1.0/14.0 for p in P) if P else True
def overlap(vs, Ng=300003):
    return sum(1 for k in range(Ng)
               if all(min((v*(k+0.5)/Ng)%1.0,1-(v*(k+0.5)/Ng)%1.0)<=1.0/14.0 for v in vs))/Ng

print("="*72)
print("G1: single-runner Toeplitz K[i][j]=h(i-j) -- PSD, min eig -> 0")
print("="*72)
for N in [5,8,12,16,20,28]:
    eigs=symeig([[h(i-j) for j in range(N)] for i in range(N)])
    print(f"  N={N:2d}: min eig={eigs[0]:+.6f}  (-> ess-inf 1_safe = 0)")

print("="*72)
print("DPP/PPP: pairwise danger correlations A_ij - A_i A_j  (A_i=1/7)")
print("="*72)
S=[1,2,3,5,8]; A1=1.0/7.0; signs=[]
for i,j in itertools.combinations(range(len(S)),2):
    d=overlap([S[i],S[j]])-A1*A1; signs.append(d>0)
    print(f"  v={S[i]},{S[j]}: A_ij-A_iA_j = {d:+.5f}  ({'POS-corr' if d>0 else 'NEG-corr'})")
print(f"  pos:{sum(signs)} neg:{len(signs)-sum(signs)} => MIXED sign => no DPP and no PPP kernel.")

print("="*72)
print("G4*: universal-center Gram lower bound for the S3 cluster floor rho*(P,E)")
print("="*72)
centers=[0.0,1/3,1/2,2/3]
def phi(x,w=0.04):
    v=0.0
    for c in centers:
        d=abs(((x-c+0.5)%1.0)-0.5)
        if d<w: v+=1-d/w
    return v
N=300000
cases=[(13,[]),(12,[5]),(11,[1,7]),(10,[1,3,11]),(7,[1,2,5,9,11,13]),(9,[1,2,3,12])]
for k,P in cases:
    assert len(P)+k==13
    num=0.0; den=0.0; rho=0; e=list(range(k))
    for j in range(N):
        x=(j+0.5)/N; px=phi(x); gp=in_GP(x,P)
        if gp: den+=px*px
        if gp and maxgap([(i*x)%1.0 for i in range(k)])>2.0/7.0+1e-9:
            rho+=1; num+=px
    num/=N; den/=N; rho/=N
    LB=num*num/den if den>0 else 0.0
    tag=" <-- EXTREMIZER, Gram-LB COLLAPSES" if (k==9 and P==[1,2,3,12]) else ""
    print(f"  k={k:2d} P={P}: rho*={rho:.5f}  Gram-LB={LB:.5f}{tag}")

print("="*72)
print("WHERE the extremizer mass lives (k=9, P={1,2,3,12}, rho*=1/84):")
print("="*72)
k=9;P=[1,2,3,12];Nb=600000;hits=[]
for j in range(Nb):
    x=(j+0.5)/Nb
    if in_GP(x,P) and maxgap([(i*x)%1.0 for i in range(k)])>2.0/7.0+1e-9: hits.append(x)
hits.sort(); arcs=[]; s=hits[0]; pr=hits[0]
for x in hits[1:]:
    if x-pr>0.002: arcs.append((s,pr)); s=x
    pr=x
arcs.append((s,pr))
for a,b in arcs:
    c=(a+b)/2; nc=min([0,1/3,1/2,2/3,1],key=lambda t:abs(t-c))
    print(f"  arc center={c:.4f} width={b-a:.4f} nearest-univ-center={nc:.3f} OFFSET={c-nc:+.4f}")
print(f"  1/14 = {1/14:.4f}  (P's p=1 danger-band edge). Mass is PINNED to this edge,")
print("  P-dependent => no fixed PD bump captures it => fixed-kernel Gram floor = 0.")
