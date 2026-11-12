"""Exact algebra checks; the accompanying proof supplies universal statements."""
from pathlib import Path
import json
import sympy as s

HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def zero(x):return s.simplify(s.expand_complex(x))==0

def run():
    x,y,A,b,d,C,p0=s.symbols('x y A b d C p0',real=True)
    z=s.symbols('z'); I=s.I; zz=x+I*y; B=b+I*d
    lam=A*(x*x+y*y)+2*b*x-2*d*y+C
    checks=[]
    for name,h in [('cartesian',1),('square_net',z*z),('cubic_net',(1+z*z)**2)]:
        J=s.integrate(h*(A*z+s.conjugate(B)),z)
        p=s.expand_complex(p0-J.subs(z,zz)-s.conjugate(J.subs(z,zz)))
        q=s.expand_complex(s.sympify(h).subs(z,zz)*lam)
        xx=p+s.re(q); yy=p-s.re(q); xy=-s.im(q)
        residuals=[s.diff(xx,x)+s.diff(xy,y),s.diff(xy,x)+s.diff(yy,y),s.diff(p,x,2)+s.diff(p,y,2)]
        need(all(zero(v) for v in residuals),name)
        checks.append({'net':name,'equilibrium_and_local_compatibility':True})
    # A rational conformal coordinate, on a domain excluding poles and zeros
    # of Z'. Stress single-valuedness is weaker than displacement closure.
    examples=[('rank0',z,[]),('one_real_pole',z+1/z,[0]),('one_imaginary_pole',z+I/z,[0]),('two_poles',z+1/z+I/(z-1),[0,1]),('three_poles',z+1/z+I/(z-1)+(1+I)/(z-I),[0,1,I])]
    periods=[]
    for name,Z,poles in examples:
        h=s.diff(Z,z)**2
        # Exact residues by differentiation of the finite principal part.
        def res(expr,a):return s.simplify(s.residue(expr,z,a))
        R=[res(h,a) for a in poles]
        S=[res(z*h,a) for a in poles]
        V=[res(z*z*h,a) for a in poles]
        stress=[];elastic=[]
        for r,t,v in zip(R,S,V):
            f=s.expand(A*t+(b-I*d)*r)
            stress.append(s.im(s.expand_complex(f)))
            elastic.extend([s.re(s.expand_complex(f)),s.im(s.expand_complex(f))])
            # kappa=3 (Poisson ratio 0), after division by 2*pi*i:
            # kappa res(zF) - conjugate(res G) = 0.
            tr=s.expand(3*(A*v+(b-I*d)*t)-s.conjugate((b+I*d)*t+C*r))
            elastic.extend([s.re(s.expand_complex(tr)),s.im(s.expand_complex(tr))])
        def matrix(es):return s.Matrix([[s.expand(e).coeff(v) for v in (A,b,d,C)] for e in es])
        sr=matrix(stress).rank() if stress else 0
        er=matrix(elastic).rank() if elastic else 0
        periods.append({'name':name,'Z':str(Z),'res_h':list(map(str,R)),'res_zh':list(map(str,S)),'stress_period_matrix':matrix(stress).tolist() if stress else [],'stress_dimension':5-sr,'elastic_dimension_nu0':5-er,'elastic_period_matrix':matrix(elastic).tolist() if elastic else []})
    need([r['stress_dimension'] for r in periods]==[5,5,4,3,2],'period ranks')
    need([r['elastic_dimension_nu0'] for r in periods]==[5,2,2,1,1],'elastic ranks')
    # Cartesian sensor matrix in coordinates A,b,d,C,p0.
    xx=2*A*y*y-4*d*y+C+p0
    yy=-2*A*x*x-4*b*x-C+p0
    vars=(A,b,d,C,p0)
    def sensors(points):
        return s.Matrix([[f.subs({x:u,y:v}).expand().coeff(w) for w in vars] for u,v in points for f in (xx,yy)])
    sensing=[]
    for name,pts in [('two_generic',[(0,0),(1,2)]),('three_generic',[(0,0),(1,1),(2,3)]),('three_rectangle_corners',[(0,0),(1,0),(0,1)]),('four_rectangle_corners',[(0,0),(1,0),(0,1),(1,1)])]:
        mat=sensors(pts);sensing.append({'name':name,'points':pts,'rank':mat.rank(),'kernel':[list(v) for v in mat.nullspace()]})
    need([t['rank'] for t in sensing]==[4,5,4,4],'sensor ranks')
    # Circle evaluations, independent of h, and exact cocircular failure.
    circle=[]
    for pts in [[(0,0),(1,0),(0,1),(1,1)],[(0,0),(1,0),(0,1),(2,2)]]:
        M=s.Matrix([[u*u+v*v,2*u,-2*v,1] for u,v in pts])
        circle.append({'points':pts,'rank':M.rank(),'det':M.det()})
    need([v['rank'] for v in circle]==[3,4],'circle sampling')
    # Translation operators on the five-dimensional Cartesian stress space.
    basis=[s.Matrix([xx.diff(v),yy.diff(v)]) for v in vars]
    def coords(field):
        coeff=s.symbols('c:5')
        diff=field-sum((c*f for c,f in zip(coeff,basis)),s.zeros(2,1))
        eq=[a for f in diff for a in s.Poly(f,x,y).coeffs()]
        sol=s.solve(eq,coeff);return s.Matrix([sol[c] for c in coeff])
    Dx=s.Matrix.hstack(*(coords(v.diff(x)) for v in basis))
    Dy=s.Matrix.hstack(*(coords(v.diff(y)) for v in basis))
    need(Dx*Dy==s.zeros(5) and Dx**3==s.zeros(5) and Dy**3==s.zeros(5),'jet relations')
    alg=[s.eye(5),Dx,Dy,Dx**2,Dy**2]
    need(s.Matrix.hstack(*(m.reshape(25,1) for m in alg)).rank()==5,'jet algebra dimension')
    # Four sign changes are impossible on the positive parabola, but possible
    # already on the positive cubic, with exactly rational controls.
    t=s.symbols('t',positive=True)
    f=t**6-10*t**3+t**2-t/100+s.Rational(1,100000)
    pts=[s.Rational(1,10000),s.Rational(1,200),s.Rational(3,100),s.Integer(1),s.Integer(3)]
    signs=[s.sign(f.subs(t,u)) for u in pts]
    need(signs==[1,-1,1,-1,1],'cubic stress C5 witness')
    out={'pde_checks':checks,'periods':periods,'stress_sensing':sensing,'circle_sensing':circle,'translation_Dx':Dx.tolist(),'translation_Dy':Dy.tolist(),'translation_algebra_dimension':5,'cubic_curve_control':{'lambda':str(f),'parameters_A_b_d_C':[1,s.Rational(-1,200),5,s.Rational(1,100000)],'t_values':pts,'signs':signs}}
    dest=HERE/('stress_receipt.optimized.json' if not __debug__ else 'stress_receipt.json')
    dest.write_text(json.dumps(out,indent=2,default=str)+'\n',encoding='utf8')
    print(json.dumps({'pde_checks':checks,'period_dimensions':[(q['name'],q['stress_dimension'],q['elastic_dimension_nu0']) for q in periods],'sensor_ranks':[q['rank'] for q in sensing],'output':dest.name}))

if __name__=='__main__':run()
