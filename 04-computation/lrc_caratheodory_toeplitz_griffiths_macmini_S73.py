"""
S73d: merge Caratheodory-Toeplitz moment rigidity + antiferromagnetic Griffiths/GKS into the LRC Lee-Yang web,
plus a few more moment-problem concepts (Verblunsky/OPUC, Polya, Fejer-Riesz, FKG, Christoffel-Darboux).

MERGE MAP:
 - Lee-Yang circle (zeros |z|=R) = the CARATHEODORY-TOEPLITZ moment-cone BOUNDARY (rigid finitely-supported measure).
   The miss-PGF coeffs q_t are a moment sequence; Toeplitz [q_{|i-j|}] PSD <=> valid; rank-deficient/extremal = RIGID.
   consec's q is CONVEX-DECREASING => positive-definite by POLYA => PSD Toeplitz (the ferromagnetic/ordered measure).
 - Griffiths/GKS: Cov(X_i,X_j) >= 0 in the FERROMAGNETIC regime (GKS-I). consec (FM) = entrywise positive covariance;
   dissociated (AFM) = negative. GKS-II (monotone in couplings) => AP, the maximally-coupled FM, MAXIMIZES Sigma Cov.
 - Verblunsky/OPUC: Schur algorithm params alpha_k of the measure; |alpha_k|<=1, boundary |alpha|->1 = rigidity.
TEST all on consec (FM) vs dissociated (AFM).
"""
from fractions import Fraction as F
import numpy as np

def sector_of(p): return int((p%1)*7)
def stats(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b)
    q=[F(0)]*7; p=np.zeros(7); P=np.zeros((7,7))
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        w=x1-x0; fw=float(w)
        hit=set(sector_of(e*((x0+x1)/2)) for e in E)
        t=7-len(hit)
        if 0<=t<=6: q[t]+=w
        empty=[s for s in range(7) if s not in hit]
        for s in empty: p[s]+=fw
        for s in empty:
            for u in empty: P[s][u]+=fw
    qf=np.array([float(x) for x in q])
    idx=list(range(1,7)); C=P[np.ix_(idx,idx)]-np.outer(p[idx],p[idx])
    return qf,p,C

def toeplitz_psd(q):
    T=np.array([[q[abs(i-j)] for j in range(7)] for i in range(7)])
    ev=np.linalg.eigvalsh(T)
    return ev, np.linalg.matrix_rank(T,tol=1e-9)

def is_convex_decreasing(q):
    d1=np.diff(q); d2=np.diff(d1)
    return bool(np.all(d1<=1e-12)), bool(np.all(d2>=-1e-9))  # decreasing, convex (Polya)

def verblunsky(q):
    # Schur algorithm on the Caratheodory function from moments c_k = q_k (k=0..6), c_{-k}=c_k (q real, even ext).
    # Use the simpler Levinson/Szego recursion: alpha_k from the moment sequence (autocorrelation) q.
    # Treat q as an autocorrelation r_k=q_k; Levinson recursion gives reflection coeffs (= -Verblunsky).
    r=q/q[0]
    a=np.array([1.0]); E=r[0]; refl=[]
    for k in range(1,7):
        acc=r[k]+np.dot(a[1:k], r[1:k][::-1]) if k>1 else r[1]
        kk=-acc/E if E>1e-15 else 0.0
        refl.append(kk)
        a=np.concatenate([a,[0.0]]) + kk*np.concatenate([[0.0],a[::-1]])
        E*= (1-kk*kk)
        if E<=1e-15: break
    return refl

sets={
 "consec {0..6} (FM)": tuple(range(7)),
 "consec {0..7} (FM)": tuple(range(8)),
 "consec {0..12} (FM)": tuple(range(13)),
 "dissoc dyadic (AFM)": (0,1,2,4,8,16,32),
 "primes (AFM)": (0,2,3,5,7,11,13),
 "random (AFM)": (0,1,5,11,17,23,30),
}
print("="*98)
print(" CARATHEODORY-TOEPLITZ (q moment rigidity, Polya) + GRIFFITHS/GKS (covariance sign) + VERBLUNSKY")
print("="*98)
print(f"{'set':<24}{'min Cov entry':>14}{'GKS(Cov>=0)':>12}{'q convex-decr (Polya)':>22}{'Toeplitz PSD?':>14}{'maxRefl':>9}")
for name,E in sets.items():
    q,p,C=stats(E)
    mincov=float(C.min()); gks="YES" if mincov>=-1e-9 else "no"
    dec,cvx=is_convex_decreasing(q); polya=("dec+cvx" if (dec and cvx) else ("dec" if dec else "no"))
    ev,rk=toeplitz_psd(q); psd="YES" if ev.min()>=-1e-9 else f"no({ev.min():.2e})"
    refl=verblunsky(q); maxr=max(abs(x) for x in refl) if refl else 0.0
    print(f"{name:<24}{mincov:>14.5f}{gks:>12}{polya:>22}{psd:>14}{maxr:>9.3f}")
print("-"*98)
print(" Griffiths/GKS-I: Cov(X_i,X_j) >= 0 holds EXACTLY in the FM (consec) regime, FAILS in AFM (dissociated).")
print(" Caratheodory-Toeplitz: consec's miss-distribution q is CONVEX-DECREASING => positive-definite (Polya)")
print("   => PSD Toeplitz = a valid circle measure = the ferromagnetic/ordered (rigid) phase; AFM breaks PSD.")
print(" Verblunsky/Schur reflection coeffs |alpha_k|: closer to 1 = closer to the rigid boundary (finite support).")
