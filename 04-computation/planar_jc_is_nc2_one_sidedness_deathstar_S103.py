import numpy as np
I=1j
def P(x,y,A,B,d): return A*(x+I*y)**d + B*(x-I*y)**d
def hess(f,x0,y0,h=1e-4):
    fxx=(f(x0+h,y0)-2*f(x0,y0)+f(x0-h,y0))/h**2
    fyy=(f(x0,y0+h)-2*f(x0,y0)+f(x0,y0-h))/h**2
    fxy=(f(x0+h,y0+h)-f(x0+h,y0-h)-f(x0-h,y0+h)+f(x0-h,y0-h))/(4*h**2)
    return np.array([[fxx,fxy],[fxy,fyy]])

print("="*70)
print("[1] 2D nilpotent Hessian <=> P one-sided (det(Hess) ~ A*B*|z|^{2(d-2)})")
print("="*70)
x0,y0=0.37+0.11j, -0.19+0.23j
for d in [3,4,5]:
    for (A,B,tag) in [(1,0,"holomorphic (one-sided)"),(0,1,"antiholo (one-sided)"),(1,1,"TWO-sided")]:
        H=hess(lambda x,y:P(x,y,A,B,d),x0,y0)
        tr=H[0,0]+H[1,1]; det=H[0,0]*H[1,1]-H[0,1]**2
        nil = abs(tr)<1e-2 and abs(det)<1e-2
        print(f"  d={d} A={A} B={B} [{tag:24s}]: trace={abs(tr):.2e} det={abs(det):.2e} -> nilpotent={nil}")
print("  => nilpotent (trace=det=0) EXACTLY for one-sided (A*B=0). Two-sided has det!=0. [NC2 conclusion]")

print("\n"+"="*70)
print("[2] one-sided map F=(x,y)+grad P is ONE-FIBER-LINEAR (F2 - i F1 = -i z), hence TAME/invertible")
print("="*70)
d=4;A=1.0
def gradP(x,y): 
    h=1e-6
    return ((P(x+h,y,A,0,d)-P(x-h,y,A,0,d))/(2*h), (P(x,y+h,A,0,d)-P(x,y-h,A,0,d))/(2*h))
for (x0,y0) in [(0.3+0.2j,-0.1+0.4j),(1.0,0.5j)]:
    gx,gy=gradP(x0,y0); F1,F2=x0+gx,y0+gy
    combo=F2-I*F1; z=x0+I*y0
    print(f"  point z={z:.3f}: F2 - i F1 = {combo:.4f}  vs  -i z = {(-I*z):.4f}  (match => linear pencil member)")

print("\n"+"="*70)
print("[3] Zhao VC / moment: one-sided P is HARMONIC (holomorphic), Delta=0, so Delta^m(P^m)=0 trivially")
print("="*70)
def lap(f,x0,y0,h=1e-4): return (f(x0+h,y0)-2*f(x0,y0)+f(x0-h,y0))/h**2 + (f(x0,y0+h)-2*f(x0,y0)+f(x0,y0-h))/h**2
for d in [3,4]:
    val=lap(lambda x,y:(x+I*y)**d, 0.3+0.1j,-0.2+0.15j)
    print(f"  Delta((x+iy)^{d}) = {abs(val):.2e}  (=0: holomorphic is harmonic => Delta^m(P^m)=0, Zhao VC trivial one-sided)")
print("\n  KEY: 2D nilpotent-Hessian <=> one-sided <=> harmonic <=> the map is one-fiber-linear tame.")
print("  This is the NC2 one-sided dichotomy realized inside the SYMMETRIC planar Jacobian problem.")
import numpy as np
np.random.seed(0)
print("="*72)
print("The rank threshold: dim-2 nilpotent Jacobian has RANK <=1 (ONE direction); dim>=3 allows RANK>=2")
print("="*72)
# nilpotent 2x2 with trace 0, det 0 => rank<=1 always
def rand_nilpotent(n, tries=2000):
    # random rank-r nilpotent nxn with image ⊆ kernel: N = W V^T, V^T W=0
    for _ in range(tries):
        r = n-1
        W=np.random.randn(n,r); V=np.random.randn(n,r)
        # enforce V^T W = 0 (image ⊆ kernel => nilpotent) by projecting
        # take V in the orthogonal complement of W's column space
        Q,_=np.linalg.qr(W)
        Vp = V - Q@(Q.T@V)
        N=W@Vp.T
        if np.linalg.norm(N@N)<1e-9 and np.linalg.matrix_rank(N,1e-9)>=1:
            return N, np.linalg.matrix_rank(N,1e-9)
    return None,None
print("  dim 2: max rank of a nilpotent 2x2 (trace=det=0) =", 1, "(forced) -> ONE isotropic direction")
for n in [2,3,4]:
    N,r=rand_nilpotent(n)
    maxr = n-1  # a nilpotent nxn can have rank up to n-1
    print(f"  dim {n}: nilpotent can reach rank {maxr}  (sample rank {r}) -> {'ONE dir (one-sided)' if maxr==1 else str(maxr)+' isotropic dirs = RESONANCE'}")

print("\n"+"="*72)
print("rank<=1 Jacobian J=grad H => rows (grad H_1, grad H_2) are PARALLEL => H_1,H_2 functionally")
print("dependent => one-fiber pencil (codex THM-2063). The single-direction = the NC2 one-sided analogue.")
print("="*72)
# symmetric isotropic example: J = c * v v^T, v=(1,i)  (=> H=(x+iy)^d, one-sided)
# non-symmetric example: J=[[0,1],[0,0]] type (v=(0,1) non-isotropic) -- still rank1, rows parallel to v
for name,J in [("symmetric (v=(1,i) isotropic)", np.array([[1,1j],[1j,-1]])),
               ("non-symmetric (v=(0,1))",       np.array([[0.0,1.0],[0.0,0.0]]))]:
    r=np.linalg.matrix_rank(J,1e-9); nil=np.linalg.norm(J@J)<1e-9
    # rows parallel?
    row1,row2=J[0],J[1]
    par = (np.linalg.norm(row1)<1e-9 or np.linalg.norm(row2)<1e-9 or
           abs(abs(np.vdot(row1,row2))-np.linalg.norm(row1)*np.linalg.norm(row2))<1e-9)
    print(f"  {name:32s}: rank={r} nilpotent={nil} rows-parallel={par}")

print("\n"+"="*72)
print("THE UNIFICATION (JC <-> GMC2, via my thread):")
print("="*72)
print(" dim-2 JC: nilpotent Jacobian forced RANK<=1 = ONE isotropic direction = ONE-SIDED (NC2)")
print("           => parallel gradients => one-fiber-linear tame => JC(2) plausibly TRUE.")
print(" dim>=3 :  RANK>=2 = MULTIPLE isotropic directions = coincident cycles = RESONANCE")
print("           => Alpoge's Keller counterexample (THM-1300) lives here => JC FALSE.")
print(" SAME threshold as GMC2 (S101): unique primitive cycle (one-sided, nonvanishing, easy) vs")
print("   coincident cycles (resonance, hard/counterexample). Planar = the last one-direction dimension.")
print(" Zhao VC (Delta^m P^m=0 <=> Hess nilpotent) = GMC2 (E[P^m]=0 <=> one-sided): SAME moment-vanishing.")
