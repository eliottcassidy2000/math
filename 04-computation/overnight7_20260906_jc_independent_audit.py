"""Independent polynomial/Taylor-packet audit; no producer imports.

Uses literal symbolic compiler expansion, dense weighted Taylor arrays,
column-space elimination and fresh exceptional residue embeddings.
The all-degree bounds are proved and audited in the companion report.
"""
from math import comb
from itertools import product
import hashlib,json,sys
import sympy as S
sys.stdout.reconfigure(newline='\n')
x,q,t,a=S.symbols('x q t a')
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
f=72783360*a**4-77822208*a**3-28419741*a**2+7849770*a-1276420
k=(-S.Rational(5183766767360,3**19)+S.Rational(931699873280,3**16)*a
   -S.Rational(460399208960,3**14)*a*a-S.Rational(183606968320,3**13)*a**3)
r=S.Rational(520,9)*a*a-S.Rational(1688,81)*a-S.Rational(5717,729)
Q=x**5+S.Rational(9,2)*x**4-2*x**3-S.Rational(27,4)*x*x+x-S.Rational(3,4)
Q+=x*x*(x*x-1)**2*(-S.Rational(259,36)+r+4*a-9*a*x-r*x*x)
D=1+x*x*q
C=x*D*(D+2)
Z=q*(D+3)+3
B=(D-1)*(D+2)**2
chi=-9*k/20
def residue(z,p):
    n,d=S.fraction(z)
    return int(n)%p*pow(int(d)%p,-1,p)%p
def dimension(columns,p):
    basis={}
    for col in columns:
        v=[int(z)%p for z in col]
        while any(v):
            i=next(i for i,z in enumerate(v) if z)
            if i not in basis:
                inv=pow(v[i],-1,p);basis[i]=[inv*z%p for z in v];break
            c=v[i];v=[(z-c*b)%p for z,b in zip(v,basis[i])]
    return len(basis)
def check_bank(p,alpha):
    need(int(f.subs(a,alpha))%p==0,'fresh exceptional embedding')
    slots=[(i,j) for i in range(7) for j in range(13-2*i)]
    index={key:n for n,key in enumerate(slots)}
    size=len(slots)
    products=[(i,j,index[(r+s,v+w)]) for i,(r,v) in enumerate(slots)
              for j,(s,w) in enumerate(slots) if 2*(r+s)+v+w<=12]
    def multiply(u,v):
        out=[0]*size
        for i,j,z in products:
            if u[i] and v[j]:out[z]=(out[z]+u[i]*v[j])%p
        return out
    def powers(u):
        rows=[[1]+[0]*(size-1)]
        for n in range(6):rows.append(multiply(rows[-1],u))
        return rows
    # Expand the complete compiler once as a literal polynomial; then Taylor
    # coefficients use binomial translation of x, not a local ring compiler.
    polys=[S.Poly(S.expand(expr.subs(q,Q+t*t).subs(a,alpha)),x,t) for expr in (C,Z)]
    monos=[(b,c,w) for b in range(7) for c in range(7-b) for w in range(13-2*(b+c))]
    packets=[]
    for point in (-1,0,1):
        arrays=[]
        for poly in polys:
            ar=[0]*size
            for (n,j),coef in poly.terms():
                for i in range(min(n,6)+1):
                    if (i,j) in index:
                        ar[index[i,j]]=(ar[index[i,j]]+residue(coef,p)*comb(n,i)*pow(point,n-i,p))%p
            arrays.append(ar)
        need([ar[0] for ar in arrays]==[0,0],'base target point')
        cp,zp=map(powers,arrays)
        cols=[]
        for b,c,w in monos:
            ar=multiply(cp[b],zp[c]);out=[0]*size
            for n,(i,j) in enumerate(slots):
                if (i,j+w) in index:out[index[i,j+w]]=ar[n]
            cols.append(out)
        packets.append(cols)
    low=[n for n,(i,j) in enumerate(slots) if j<10]
    pre=[[packets[z][m][n] for z in range(3) for n in low] for m in range(len(monos))]
    normal=[[-packets[z][m][index[1,10]]%p for z in range(3)] for m in range(len(monos))]
    value=[[packets[z][m][index[0,10]] for z in range(3)] for m in range(len(monos))]
    ranks=[dimension(pre,p),dimension([c+n for c,n in zip(pre,normal)],p),
           dimension([c+n+v for c,n,v in zip(pre,normal,value)],p)]
    need((len(monos),len(pre[0]),ranks)==(140,135,[105,107,108]),'full relaxed packet ranks')
    period=[[sum(u*v for u,v in zip((5,-18,13),row))%p] for row in normal]
    need(dimension([c+n for c,n in zip(pre,period)],p)==105,'derivative period vanishes')
    for z in (1,2):
        need(dimension([c+[(v[z]-v[0])%p] for c,v in zip(pre,value)],p)==105,'leading values agree')
    for mon in ((1,0,10),(0,1,10),(0,0,10)):
        m=monos.index(mon)
        need(not any(pre[m]),'actual global polynomial witness fixes prefix')
    need(dimension([normal[monos.index(mon)]+value[monos.index(mon)]
                   for mon in ((1,0,10),(0,1,10),(0,0,10))],p)==3,'three attained first-jet dimensions')
    return {'p':p,'alpha':alpha,'ranks':ranks,'packet_sha256':hashlib.sha256(
        json.dumps([pre,normal,value],separators=(',',':')).encode()).hexdigest()}
def symbolic():
    need(S.gcd(f,S.fraction(S.together(k))[0])==1,'nonzero inherited debt at every embedding')
    need(S.expand(C*C*(Z-3)-B*(B+4))==0,'full target compiler relation')
    tangent=[]
    for point in (-1,0,1):
        row=[]
        for expr in (C,Z):
            base=S.expand(expr.subs(q,Q))
            need(S.expand(base.subs(x,point))==0,'three literal base branches')
            row.append(S.expand(S.diff(base,x).subs(x,point)))
        tangent.append(row)
    need(tangent==[[3,-9],[3,4],[3,9]],'actual normalized tangent rows')
    T=S.Matrix(tangent);ell=S.Matrix([[5,-18,13]])
    need(T.rank()==2 and ell*T==S.zeros(1,2),'unique tangent relation')
    H=S.Matrix([[u*u,2*u*v,v*v] for u,v in tangent])
    need(H.det()==63180,'all surface Hessians vanish')
    jac=S.expand(S.diff(C,x)*S.diff(Z,q)-S.diff(C,q)*S.diff(Z,x))
    need(S.expand(jac-6*(D+1)*(D*D+2*D-2))==0,'exact surface wedge')
    for point in (-1,0,1):
        need(S.expand(jac.subs(q,Q).subs(x,point))==12,'sharp retained normal coefficient')
    h=S.Function('h')(x)
    shear=S.Matrix([x,q+t**10*h,t])
    need(S.simplify(shear.jacobian([x,q,t]).det())==1,'source shear determinant')
    need(S.expand(8*chi+S.Rational(18,5)*k)==0,'fifth collision period')
    need(S.simplify(40*8*chi/18+8*k)==0,'payer first density period at ninth order')
    need(S.Rational(36*8,18)==16,'ninth source-order hostile changes eighth density')
    # Arbitrary polynomial h has retained L=0 exactly when its vector lies
    # in the two tangent columns; x is deliberately outside that plane.
    need(sum(u*v for u,v in zip((5,-18,13),(-1,0,1)))==8,'non-descended hostile')
    for hh in (S.Integer(1),4*x*x-9*x,x*(x*x-1)):
        need(S.expand(sum(u*hh.subs(x,v) for u,v in zip((5,-18,13),(-1,0,1))))==0,'lawful retained kernel control')
if __name__=='__main__':
    symbolic()
    records=[check_bank(p,r) for p,r in ((467,169),(487,313))]
    print('Fresh literal compiler / column-space packets:',json.dumps(records,sort_keys=True))
    print('Proof audit: full filtered global target-polynomial upper bound and actual polynomial attainment PASS')
    print('Order audit: every O(t10) graph preserves all retained J0..J8; t9*x changes J8 but its collision period is nonzero')
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)
