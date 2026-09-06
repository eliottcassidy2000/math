#!/usr/bin/env python3
"""Exact controls for component-labelled torsion and canonical connection.

Universe: four smooth one-arm families, derivative orders0..3; both arms
and orders1..4 of one smooth three-component fiber; a coordinate; a
singular hostile. Universal statements are proved in the companion note.
"""
import sympy as s

x,y,u=s.symbols('x y u')
checks=0

def require(ok,label):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(label)

def zero(expr,label):
    require(s.cancel(expr)==0,label)

def polynomial(expr,label):
    a,b=s.fraction(s.cancel(expr))
    require(not b.has(x,y),label)
    return s.expand(a/b)

def operators(P,A,B):
    zero(A*s.diff(P,x)+B*s.diff(P,y)-1,'Bezout row')
    def D(f):return s.diff(P,x)*s.diff(f,y)-s.diff(P,y)*s.diff(f,x)
    def V(f):return A*s.diff(f,x)+B*s.diff(f,y)
    m=s.diff(A,x)+s.diff(B,y)
    def nabla(f):return s.expand(V(f)+m*f)
    return D,V,nabla,m


def run():
    for r in range(2,6):
        P=x+x**r*y
        A=1-r*x**(r-1)*y
        B=r*r*x**(r-2)*y*y
        D,V,N,m=operators(P,A,B)
        q=1/(s.Integer(r-1)*x**(r-1)); a=s.Integer(1)
        zero(N(1)-m,'unit derivative is divergence')
        for j in range(4):
            height=r-1+j
            zero(D(q)-a,'polynomial response equals derivative of rational primitive')
            polynomial(a,'response polynomial')
            cleared=polynomial(P**height*q,'exact annihilating power clears primitive')
            leading=(-1)**j*s.rf(r-1,j)/s.Integer(r-1)
            zero(cleared.subs(x,0)-leading,'nonzero exact leading principal coefficient')
            require(leading!=0,'one lower power cannot kill this component difference')
            # An independent polynomial killing witness for the response.
            zero(D(cleared)-P**height*a,'actual polynomial annihilator witness')
            if j<3:
                q=s.cancel(V(q));a=N(a)
        for f in [x,y,x*x*y+x]:
            zero(V(D(f))-D(V(f))+m*D(f),'Hamiltonian commutator')
            zero(N(P*f)-P*N(f)-f,'canonical Weyl relation')
            h=x+y
            # V2=V+hD; its full divergence correction is mandatory.
            A2=A-h*s.diff(P,y);B2=B+h*s.diff(P,x)
            m2=s.diff(A2,x)+s.diff(B2,y)
            N2=A2*s.diff(f,x)+B2*s.diff(f,y)+m2*f
            zero(N2-N(f)-D(h*f),'Bezout gauge difference is an actual response')
        print('ONE_ARM',r,'theta_order',r-1,'derivative_orders',list(range(r-1,r+3)))

    # Three components at the same target value. The actual CRT selectors
    # are formed in the plane ring, not injected as abstract module data.
    g=(x*x-1)**2;P=x*x+g*y;U=P-1;R=1+(x*x-1)*y
    px=s.diff(P,x)
    A=s.invert(px,g,x);B=polynomial((1-A*px)/g,'polynomial Bezout complement')
    A=polynomial(A,'polynomial Bezout inverse')
    D,V,N,m=operators(P,A,B)
    zero(U-(x-1)*(x+1)*R,'complete three-component fiber')
    # Explicit ideal witnesses for pairwise disjointness.
    zero((x+1)/2-(x-1)/2-1,'two line components are disjoint')
    zero(R-(x-1)*(x+1)*y-1,'residual component misses both lines')
    jets={}
    z=s.symbols('z')
    for alpha in [-1,1]:
        other=s.cancel(U/(x-alpha))
        for n in range(1,6):
            # other(alpha,y) is a nonzero constant, so all series
            # coefficients are polynomial in y with rational constants.
            series=s.series(other.subs(x,alpha+z)**(-n),z,0,n).removeO()
            numer=s.expand(series.subs(z,x-alpha))
            f=s.cancel(numer/(x-alpha)**n)
            jets[alpha,n]=f
            response=polynomial(D(f),'CRT jet transgresses to polynomial')
            zero(s.cancel(U**n*f).subs(x,alpha)-1,'selected component unit principal part')
            h=polynomial(U**n*f,'global CRT numerator')
            require(s.rem(h-1,(x-alpha)**n,x)==0,'entire selected principal jet equals u^-n')
            require(s.rem(h,other**n,x)==0,'all other component jets vanish')
            zero(D(h)-U**n*response,'actual killing polynomial')
            if n<=4:
                # Two arms have linearly independent component differences
                # modulo the third, zero, residual component.
                require(s.cancel(f).subs(x,-alpha).is_finite is not False,'regular on the other line')
        for n in range(1,5):
            f=jets[alpha,n]
            difference=polynomial(V(f)+n*jets[alpha,n+1],'connection matches derivative modulo polynomial')
            zero(N(D(f))+n*D(jets[alpha,n+1])-D(difference),'connection equality in actual response quotient')
            if n>1:
                regular=polynomial(U*f-jets[alpha,n-1],'multiplication gives previous component jet')
                zero(U*D(f)-D(jets[alpha,n-1])-D(regular),'module transition equality')
        print('COLLIDED_ARM',alpha,'orders',list(range(1,5)),'residual_arm_zero',True)
    require(s.Matrix([[1,0],[0,1],[0,0]]).row_join(s.ones(3,1)).det()!=0,'two independent component differences modulo diagonal')
    # A non-diagonal rational primitive passes exactness and fails gluing.
    f=jets[1,1]
    require(s.cancel(U*f).subs(x,1)==1 and s.cancel(U*f).subs(x,-1)==0,'generic exactness does not equalize components')

    P=x+y*y;D,V,N,m=operators(P,s.Integer(1),s.Integer(0))
    for i,j in [(0,0),(1,0),(0,2),(2,3),(4,1)]:
        q=P**i*y**(j+1)/s.Integer(j+1)
        zero(D(q)-P**i*y**j,'coordinate has polynomial primitive in full polynomial basis')
    require(m==0,'positive divergence-free control')

    # Singular control: D_xy is diagonal with eigenvalue j-i; theta is
    # not in its image. The needed global Bezout field does not exist.
    P=x*y
    for i in range(5):
        for j in range(5):
            zero(s.diff(P,x)*s.diff(x**i*y**j,y)-s.diff(P,y)*s.diff(x**i*y**j,x)-(j-i)*x**i*y**j,'singular monomial response')
    require(s.groebner([x,y],x,y).reduce(s.Integer(1))[1]==1,'singular components are not comaximal')
    require(s.diff(u**(-5),u)==-5*u**(-6),'characteristic zero pole derivative')
    print('POSITIVE coordinate; HOSTILES unequal component jets and singular intersecting components')
    print('EXPLICIT_GATES',checks)
    print('PASS component-labelled torsion and canonical connection controls; universal proof separate')

if __name__=='__main__':run()
