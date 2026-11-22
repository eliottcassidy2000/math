from subprocess import check_output
from pathlib import Path
C={4:(2465,2293760,1344270,712944,1740,986,423385),
   8:(37265,144834560,84839910,57602544,24360,14906,33838015)}
def f(y):
    z=0
    while y:
        y,j=divmod(y,2)
        if not j:z+=y+(y&1)
    return z
def H(n,x,y):
    h,q,c,t,e,r,k=C[n]
    return k+q*x*y+t*x+c*y+r*(x&1)+e*(y&1)+h*(f(x)+f(y))
for n,scale in ((8,10**12),(4,10**13)):
    for parityx in (0,1,4,7):
        for parityy in (0,1):
            for rev in (False,True):
                X,Y=13*scale+17,15*scale+123
                if rev:X,Y=Y,X
                x,y=10*X+parityx,2*Y+parityy
                L=min(X,Y)
                t=H(n,x,y)
                # Large controls in both normalized orientations, multiple blocks.
                out=check_output([str(Path(__file__).resolve().with_name('farey-native-slab')),str(n),str(L-99),str(L+200),'101',str(t)],text=True)
                assert f"HIT {x} {y}\n" in out,(n,x,y,out)
print("32 planted native tangent controls PASS")
from fractions import Fraction
def phase_read(n, width):
    rev=int(f"{n:0{width}b}"[::-1],2) if width else 0
    return Fraction((1<<width)-1-rev,1<<width)
for n in range(1,4096):
    k=n.bit_length()
    rho=phase_read(n,k); eps=f(n)-n*rho
    for b in (0,1):
        z=2*n+b
        rho1,eps1=(rho/2,eps-rho/2) if b else ((rho+1)/2,eps+(n&1))
        assert rho1==phase_read(z,k+1) and f(z)-z*rho1==eps1
    j=k//2;w=k-j;v=n>>j;suffix=n&((1<<j)-1)
    epsv=f(v)-v*phase_read(v,w)
    low=phase_read(suffix,j)
    assert low<=rho<low+Fraction(1,1<<j)
    assert n*low+epsv-Fraction(j,2)<=f(n)<=n*(low+Fraction(1,1<<j))+epsv+j
print("dyadic phase and prefix/suffix envelope controls PASS")
