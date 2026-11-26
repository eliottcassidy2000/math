"""Exact checks for determinant norms, balanced-face residues and shear descent."""
from math import factorial, prod
from collections import Counter
from itertools import permutations
import sympy as s

def norm(n,d):
    return prod(factorial(n+j)//factorial(j) for j in range(d))

def baryon(d,nmax):
    z=s.symbols(f'z0:{d*d}')
    perms=list(permutations(range(d)))
    signs=[(-1)**sum(p[i]>p[j] for i in range(d) for j in range(i+1,d)) for p in perms]
    u=s.Poly(s.Matrix(d,d,z).det(),*z)
    for n in range(nmax+1):
        F=u**n
        actual=sum(int(c)**2*prod(factorial(e) for e in ex) for ex,c in F.terms())
        assert actual==norm(n,d)
        if n:
            B=s.Poly(0,*z)
            for p,sgn in zip(perms,signs):
                f=F
                for i,j in enumerate(p):f=f.diff(z[d*i+j])
                B+=sgn*f
            assert B==prod(n+j for j in range(d))*(u**(n-1))
    print("dimension",d,"factorial norms through level",nmax, [norm(j,d) for j in range(nmax+1)])

def comps(k,t):
    if k==1:
        yield (t,)
    else:
        for j in range(t+1):
            for c in comps(k-1,t-j):yield (j,)+c

def multi(n,r):
    return factorial(n)//prod(factorial(x) for x in r)

def residues():
    # Lowest face has four charged supports, two different CT pairings and a separate bulk term.
    faces=[(3,0,1),(2,1,2),(1,2,3),(0,3,5),(2,2,7)]
    A0,m0=3,2
    Q=22
    for p in (5,7,13):
        m=p*m0
        for d in (1,2,3,5):
            total=0;counts=Counter()
            for r in comps(len(faces),m):
                A=sum(a*u for (a,b,c),u in zip(faces,r))
                B=sum(b*u for (a,b,c),u in zip(faces,r))
                if A!=B:continue
                assert A>=p*A0
                v=multi(m,r)*prod(c**u for (a,b,c),u in zip(faces,r))
                relative=prod(factorial(A+j)//factorial(p*A0+j) for j in range(d))
                w=v*relative%p
                if any(u%p for u in r):group="carry";assert w==0
                elif A>p*A0:group="off-face";assert w==0
                else:group="face"
                total=(total+w)%p
                counts[group]+=1
            assert total==pow(Q,p,p) and total!=0
        print("prime",p,"d=1,2,3,5 residue",total,"last-size allocations",dict(counts))

def shear_check(d):
    x,y,u,v,N=s.symbols('x y u v N')
    f=u+2*u**2
    bare=prod(N+j for j in range(d))
    for shift in range(1-d,0):
        assert s.rem(bare.subs(N,N+shift),bare,N)!=0
    # Classical lift of the marked shear to the plane.
    Y=y+x**(d-1)*(1+2*x**d)
    assert s.Matrix((x,Y)).jacobian((x,y)).det()==1
    assert s.expand(x*Y-(x*y+f.subs(u,x**d)))==0
    # The noncanonical radical partner can be corrected by an inverse number polynomial.
    for n in range(12):
        raw=(n+1)*prod(n+1+j for j in range(1,d))-n*prod(n+j for j in range(1,d))
        assert raw==d*prod(n+j for j in range(1,d))
    print("dimension",d,"checked fixed-u classical lift and quantum constant-shift rejection")
def borel():
    t,x,z=s.symbols("t x z",positive=True)
    w=(1-2*z*x)/(1+z*x)**4
    assert s.factor(s.diff(x**2/(2*(1+z*x)**2),x,2)-w)==0
    coeff=s.series((1-2*z)/(1+z)**4,z,0,8).removeO()
    for n in range(8):
        assert coeff.coeff(z,n)==(-1)**n*norm(n,3)//(factorial(n)**3)
    print("d=3 Borel transform verified to order seven:",(1-2*z)/(1+z)**4)

if __name__=="__main__":
    baryon(2,6);baryon(3,5)
    residues()
    for d in (2,3,4,5):shear_check(d)
    borel()
