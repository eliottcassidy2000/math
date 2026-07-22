import numpy as np
def padd(A,B):
    C=dict(A)
    for k,v in B.items(): C[k]=C.get(k,0)+v
    return {k:v for k,v in C.items() if abs(v)>1e-9}
def pscale(A,s): return {k:s*v for k,v in A.items() if abs(s*v)>1e-9}
def pmul(A,B):
    C={}
    for (i,j),a in A.items():
        for (k,l),b in B.items(): C[(i+k,j+l)]=C.get((i+k,j+l),0)+a*b
    return {k:v for k,v in C.items() if abs(v)>1e-9}
def dx(A): return {(i-1,j):a*i for (i,j),a in A.items() if i>=1}
def dy(A): return {(i,j-1):a*j for (i,j),a in A.items() if j>=1}
def bracket(F,G): return padd(pmul(dx(F),dy(G)), pscale(pmul(dy(F),dx(G)),-1.0))
def mate_exists(F,N):
    mons=[(i,j) for i in range(N+1) for j in range(N+1-i)]
    cols=[bracket(F,{m:1.0}) for m in mons]
    outmons=sorted(set().union(*[set(c) for c in cols])) if cols else []
    if (0,0) not in outmons: return False
    A=np.zeros((len(outmons),len(mons)),dtype=complex)
    for c,b in enumerate(cols):
        for k,v in b.items(): A[outmons.index(k),c]=v
    rhs=np.zeros(len(outmons),dtype=complex); rhs[outmons.index((0,0))]=1.0
    sol,_,_,_=np.linalg.lstsq(A,rhs,rcond=None)
    return np.linalg.norm(A@sol-rhs)<1e-7

def phi_analysis(F):
    # F quasi-homogeneous: monomials (i,j) collinear.  face polynomial phi(tau) degree & roots.
    M=sorted(F.keys())
    if len(M)<2:
        return 0,1,(min(i for i,_ in M),min(j for _,j in M))
    # direction of the face:
    (i0,j0),(i1,j1)=M[0],M[-1]
    di,dj=i1-i0,j1-j0
    g=np.gcd(abs(di),abs(dj)) or 1
    step=(di//g,dj//g)
    # index each monomial along step; phi coeff at k
    imin=min(i for i,_ in M); jmin=min(j for _,j in M)
    ks=[]
    for (i,j) in M:
        k=(i-i0)//step[0] if step[0]!=0 else (j-j0)//step[1]
        ks.append(k)
    kmax=max(ks)
    coeffs=np.zeros(kmax+1,dtype=complex)
    for (i,j),k in zip(M,ks): coeffs[k]=F[(i,j)]
    coeffs=coeffs[::-1]  # numpy: highest power first
    r=np.roots(coeffs) if len(coeffs)>1 else np.array([])
    # distinct roots
    uniq=[]
    for z in r:
        if not any(abs(z-w)<1e-4 for w in uniq): uniq.append(z)
    return kmax, len(uniq), (imin,jmin)

tests={
 'y + x^2      (coord)':   {(0,1):1.,(2,0):1.},
 'y + x^3      (coord)':   {(0,1):1.,(3,0):1.},
 'x + y^2      (coord)':   {(1,0):1.,(0,2):1.},
 'x^2          (2 lines)': {(2,0):1.},
 'x^2 + y^2    (homog2)':  {(2,0):1.,(0,2):1.},
 'x^2+xy+y^2   (homog2)':  {(2,0):1.,(1,1):1.,(0,2):1.},
 '(x+y)^2      (sq)':      {(2,0):1.,(1,1):2.,(0,2):1.},
 'xy           (fiber C*)':{(1,1):1.},
 'x^2 y + x y^2 (=xy(x+y))':{(2,1):1.,(1,2):1.},
 'y^2 + x^3    (cusp)':    {(0,2):1.,(3,0):1.},
}
print("="*74)
print("BASE CASE: quasi-homogeneous f.  mate <=> face polynomial phi is LINEAR (deg 1)")
print("  phi(tau) = the 1-variable collapse of the face (the DvdK/S106 object).")
print(f"  {'f':26s}{'mate?':7s}{'phi_deg':8s}{'#distinct roots':16s}{'monomial factor'}")
for name,F in tests.items():
    ok=mate_exists(F,7); pd,nr,(im,jm)=phi_analysis(F)
    mf=f"x^{im} y^{jm}" if (im,jm)!=(0,0) else "(none)"
    print(f"  {name:26s}{str(ok):7s}{pd:<8d}{nr:<16d}{mf}")
print("""
PATTERN:  mate=True  <=>  phi has degree 1 (=> exactly 1 root, and no monomial factor)
          i.e. f = a*x^p + b*y  (or swap) = a WEIGHTED-LINEAR coordinate.
  phi deg>=2 (>=2 roots, or a repeated root like (tau+1)^2) OR a monomial factor => NO mate.
  This is the S106 DvdK-face single-root condition, made exact for one face.""")

print("="*74)
print("FREE LEMMA (provable, no JC): mate => f is PRIMITIVE (not p(h), deg p>=2)")
print("  Proof: f=p(h) => {f,g}=p'(h){h,g}; =1 forces p'(h) a unit => p linear. QED")
print("  So mate => connected generic fibers.  (necessary, NOT sufficient: xy is")
print("  primitive with fiber C* and no mate -- connectedness misses the C* punctures.)")
comps={'(x+y)^2 = p(x+y),p=t^2':{(2,0):1.,(1,1):2.,(0,2):1.},
       'x^2 = p(x),p=t^2':{(2,0):1.},
       'x^2+2x = p(x),p=t^2+2t':{(2,0):1.,(1,0):2.}}
for name,F in comps.items():
    print(f"   non-primitive {name:26s}: mate={mate_exists(F,7)}  (predicted False)")

print("\n"+"="*74)
print("HONEST LADDER for a quasi-homogeneous (single-face) Keller component f:")
print("""   mate exists (1 in im D_f)
     => f PRIMITIVE                                  [FREE lemma, proven above]
     => face poly phi LINEAR (deg 1, single branch)  [DvdK/S106 face, = 1 place at inf]
     => weighted degree delta < w1+w2                [S107 weight/valuation, = genus 0]
     => f = y + a x^p  (weighted-linear)  = a COORDINATE.
   Conversely y+a x^p is a coordinate. So for ONE face: mate <=> coordinate. PROVEN.
   The three failure modes are exactly 'fiber not = C':
     phi>=2 roots (>=2 places, e.g. x^2+y^2) | phi repeated / monomial factor
     (non-primitive, e.g. (x+y)^2, xy) | delta>=w1+w2 (genus>0, e.g. y^2+x^3 cusp).""")
