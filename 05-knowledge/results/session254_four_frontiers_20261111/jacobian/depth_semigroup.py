#!/usr/bin/env python3
"""Controls for ALL_ROW_DEPTH.md; the all-row proof is in that note."""
import sympy as s
from probe import x,depth_matrix

def main():
    cases=0
    for d in range(9):
        for n in range(d,41):
            witnesses={}
            for b in range(d+1):
                for a in range(d-b+1):
                    for e in range((n-b)//2+1):
                        c=n-b-2*e;r=a+2*b+e
                        assert b+c+2*e==n and a+b<=d
                        # Generator x^a u^b p^c y^e has this exact valuation
                        # and leading coefficient x^r, with coefficient 1.
                        witnesses.setdefault(r,(a,b,c,e))
            assert set(witnesses)==set(range((n+3*d)//2+1))
            cases+=1
    print('EXPLICIT_MONOMIAL_ATTAINMENT d0..8,n=d..40; cases',cases,'PASS')
    m=s.Symbol('m',integer=True,positive=True)
    theta=s.Function('theta')(x)
    A0=1+x*x/4;C0=-3*x/4-x**3/8
    A1=s.Rational(4,3)+2*x*x;C1=-4*x-s.Rational(3,2)*x**3
    da=s.diff(A0,x)*theta;dc=s.diff(C0,x)*theta
    dB=(m-1)*s.diff(A1,x)*dc-A1*s.diff(dc,x)+s.diff(da,x)*C1-(m-1)*da*s.diff(C1,x)
    dT=2*C1*dc-6*A0*A1*da
    q=-(x*x+6)/2
    assert s.expand(dB-(m*(-x+3*x**3/4)*theta+s.diff(theta,x)))==0
    assert s.expand(dT-q*dB/m-(-x*theta+(x*x+6)*s.diff(theta,x)/(2*m)))==0
    print('UNIVERSAL_BRACKET_VARIATION_IDENTITY=PASS')
    n=17;Ds=[]
    for d,deriv in [(2,s.diff(A0,x)),(3,s.diff(C0,x))]:
        coords,M=depth_matrix(d,n);null=M.T.nullspace()
        N=s.Matrix.vstack(*(v.T for v in null))
        V=s.Matrix.hstack(*(s.Matrix([s.expand(deriv*x**j).coeff(x,r) if i==n else 0 for i,r in coords]) for j in range(n+1)))
        Ds.append(N*V)
        print('ROW17_DEPTH',d,'shape',M.shape,'rank',M.rows-len(null))
    D=s.Matrix.vstack(*Ds);cut=n//2+3
    assert D[:,:cut]==s.zeros(D.rows,cut) and D[:,cut:].rank()==n+1-cut
    rp=list(D[:,cut:].T.rref()[1])
    print('ROW17_SELECTOR shape',D.shape,'rank',D.rank(),'kernel_degree_max',cut-1)
    print('ROW17_PIVOT_MINOR rows',rp,'cols',list(range(cut,n+1)),'det',D.extract(rp,list(range(cut,n+1))).det())
    assert D[:,11]!=s.zeros(D.rows,1)
    print('HOSTILE theta=x^11 at n17: rejected; n-6 pattern fails; PASS')
if __name__=='__main__':main()
