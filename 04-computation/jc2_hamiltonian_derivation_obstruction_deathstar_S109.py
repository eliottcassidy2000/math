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
def deg(A): return max((i+j for (i,j) in A), default=0)
def mate_exists(F,N):     # solve {F,g}=1 EXACTLY for g of deg<=N (reliable: full-output residual)
    mons=[(i,j) for i in range(N+1) for j in range(N+1-i)]
    cols=[bracket(F,{m:1.0}) for m in mons]
    outmons=sorted(set().union(*[set(c) for c in cols]))
    A=np.zeros((len(outmons),len(mons)),dtype=complex)
    for c,b in enumerate(cols):
        for k,v in b.items(): A[outmons.index(k),c]=v
    if (0,0) not in outmons: return False
    rhs=np.zeros(len(outmons),dtype=complex); rhs[outmons.index((0,0))]=1.0
    sol,_,_,_=np.linalg.lstsq(A,rhs,rcond=None)
    return np.linalg.norm(A@sol-rhs)<1e-7
def weight_obstructed(F):   # exists positive w with (i-1)w1+(j-1)w2>0 for all monomials
    M=list(F)
    for a in range(1,50):
        for b in range(1,50):
            w1,w2=a/17.0,b/17.0
            if all((i-1)*w1+(j-1)*w2>1e-9 for (i,j) in M): return (round(w1,2),round(w2,2))
    return None
tests={'x [coord, fiber C]':{(1,0):1.},'x^2 [fiber 2 lines]':{(2,0):1.},
 'x^2+y [coord]':{(2,0):1.,(0,1):1.},'xy [fiber C*]':{(1,1):1.},
 'x^2+y^3 [genus>0]':{(2,0):1.,(0,3):1.},'x+x^2y [fiber C*]':{(1,0):1.,(2,1):1.},
 'x^3 [3 lines]':{(3,0):1.},'x^2+xy+y^2 [homog2]':{(2,0):1.,(1,1):1.,(0,2):1.}}
print("="*72)
print("JC(2) as an operator condition:  f is a Keller component  <=>  1 in im(D_f),  D_f={f,.}")
print("  and (the conjecture)  mate exists  =>  f is a coordinate  <=>  generic fiber = C.")
print(f"  {'f':24s}{'mate?':7s}{'weight-obstruction (S107 valuation)':36s}")
for name,F in tests.items():
    ok=mate_exists(F,7); w=weight_obstructed(F)
    tag = "NO mate (valuation)" if w else "-- (dodges every weight)"
    print(f"  {name:24s}{str(ok):7s} w={str(w):14s} {tag}")
print("""
READING:
 * weight obstruction v_w(f)>w1+w2 (S107 manufactured valuation) catches x^2, x^2+y^3, x^3.
 * xy, homog2 sit ON a resonant w-face: there 1 in im(D_f) is a w-quasi-homogeneous
   CONSTANT-TERM condition (1-variable Laurent) = exactly the DvdK/S106 object.
 * x+x^2y has fiber C* (no mate) yet DODGES every weight => the true obstruction is the
   GLOBAL fiber topology (coker D_f = fiber cohomology), not one valuation = the JC(2) core.""")
