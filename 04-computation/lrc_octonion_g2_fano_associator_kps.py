"""
lrc_octonion_g2_fano_associator_kps.py  (kind-pasteur-2026-06-27-S31an)

THE OCTONION / G_2 / FANO SYNTHESIS. Abstract claim:
  14 runners = dim(G_2) = dim Aut(octonions);  7 sectors = the G_2 7-dim irrep
  = the 7 imaginary octonion units = the Fano plane (PG(2,2), |PSL(2,7)|=168).
  The project's doubling map z->2z (QR(7)={1,2,4}, order 3) IS the octonion
  multiplication's base Fano line. The associativity-compression failure (kappa_3,
  the 3-way joint-emptiness associator, my S31ai) should be the OCTONION
  NON-ASSOCIATIVITY, organized by the FANO PLANE.

TEST: compute kappa_3 over all C(6,3)=20 inner-sector triples. Group by the mod-7
multiplicative/Fano structure (the doubling orbits {1,2,4}=QR and {3,5,6}=NQR,
the Fano lines). Do the FANO-LINE triples have extremal/distinct associator?
Is kappa_3 invariant under the doubling map (the octonion automorphism)?
"""
import sys, itertools
from fractions import Fraction as F
INNER=list(range(1,7))
def sector_of(p): return int((p%1)*7)
def cells(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(F(mm,7*e))
    b=sorted(b); out=[]
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        cov=set(sector_of(e*((x0+x1)/2)) for e in E)
        out.append((cov,x1-x0))
    return out
def Pempty(cl,S):
    S=set(S); t=F(0)
    for cov,w in cl:
        if S.isdisjoint(cov): t+=w
    return t
def kappa3(cl,i,j,k):
    return (Pempty(cl,(i,j,k)) - Pempty(cl,(i,j))*Pempty(cl,(k,))
            - Pempty(cl,(i,k))*Pempty(cl,(j,)) - Pempty(cl,(j,k))*Pempty(cl,(i,))
            + 2*Pempty(cl,(i,))*Pempty(cl,(j,))*Pempty(cl,(k,)))

# Fano structure on Z/7: the doubling map d(x)=2x mod 7. Orbits: {1,2,4} (QR), {3,6,5} (NQR).
# Octonion base Fano lines (QR translates): {1,2,4} and its additive translates mod 7.
FANO_LINES=[]
base={1,2,4}
for s in range(7):
    line=frozenset((x+s)%7 for x in base)
    FANO_LINES.append(line)
FANO_LINES=set(FANO_LINES)  # 7 lines

def double(x): return (2*x)%7

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    E=tuple(range(8))  # consec_8
    cl=cells(E)
    print("=== kappa_3 (associator) over the 20 inner-sector triples of consec_8 ===")
    rows=[]
    for tr in itertools.combinations(INNER,3):
        k3=kappa3(cl,*tr)
        # is this triple a Fano line (mod 7)? (using sectors as Z/7 elements 1..6, none=0)
        is_fano = frozenset(tr) in FANO_LINES
        # is it a doubling orbit? {1,2,4} or {3,5,6}
        is_qr = set(tr)=={1,2,4}; is_nqr=set(tr)=={3,5,6}
        rows.append((float(k3),tr,is_fano,is_qr or is_nqr))
    rows.sort(reverse=True)
    for k3,tr,fano,orbit in rows:
        tag=[]
        if fano: tag.append("FANO-LINE")
        if orbit: tag.append("DOUBLING-ORBIT(QR/NQR)")
        print(f"  triple {tr}: kappa_3={k3:+.6f}  {' '.join(tag)}")
    # the two doubling orbits
    k_qr=float(kappa3(cl,1,2,4)); k_nqr=float(kappa3(cl,3,5,6))
    print(f"\n  QR orbit {{1,2,4}}: kappa_3={k_qr:+.6f}")
    print(f"  NQR orbit {{3,5,6}}: kappa_3={k_nqr:+.6f}")
    print(f"  (octonion: QR=the base Fano line e1 e2=e4; NQR=its 'conjugate')")
    # doubling-invariance test: kappa_3(i,j,k) vs kappa_3(2i,2j,2k mod 7)
    print("\n=== doubling-map (octonion automorphism z->2z) invariance of kappa_3? ===")
    ninv=0; ntot=0
    for tr in itertools.combinations(INNER,3):
        dtr=tuple(sorted(double(x) for x in tr))
        if 0 in dtr: continue  # maps out of inner
        if len(set(dtr))<3: continue
        ntot+=1
        a=kappa3(cl,*tr); b=kappa3(cl,*dtr)
        if abs(float(a-b))<1e-9: ninv+=1
    print(f"  kappa_3 doubling-invariant on {ninv}/{ntot} inner triples (1=full octonion-aut symmetry)")
    # the QR/NQR partition of single-sector emptiness
    print("\n=== single-sector emptiness p_j by QR/NQR class ===")
    for j in INNER:
        cls="QR" if j in {1,2,4} else "NQR"
        print(f"  sector {j} ({cls}): p={float(Pempty(cl,(j,))):.5f}")
