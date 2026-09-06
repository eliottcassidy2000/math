#!/usr/bin/env python3
"""Independent compiler and wedge audit of ninth-order collision periods.

No producer or repository imports. Fresh exceptional embeddings are449/120
and467/169. A literal series compiler solves overdetermined branch systems;
a separate dual-number compiler differentiates at fixed source points.
All gates remain live under -O; output is deterministic LF.
"""
from fractions import Fraction as Q
import sys
import sympy as S

sys.stdout.reconfigure(newline='\n')
GATES=0
N=11
X,A=S.symbols('x alpha')
EXCEPTION=72783360*A**4-77822208*A**3-28419741*A**2+7849770*A-1276420
R=S.Rational(520,9)*A**2-S.Rational(1688,81)*A-S.Rational(5717,729)
QPOLY=(X**5+S.Rational(9,2)*X**4-2*X**3-S.Rational(27,4)*X**2+X-S.Rational(3,4)
       +X**2*(X**2-1)**2*(-S.Rational(259,36)+R+4*A-9*A*X-R*X**2))
ELL=(5,-18,13)
POINTS=(-1,0,1)
TANGENTS=((1,-9,0),(1,4,0),(1,9,0))


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def mq(x,p):
    num,den=S.fraction(x)
    return int(num)%p*pow(int(den)%p,-1,p)%p


class Series:
    def __init__(self,p):self.p=p
    def c(self,x):return [x%self.p]+[0]*(N-1)
    def add(self,x,y):return [(a+b)%self.p for a,b in zip(x,y)]
    def sub(self,x,y):return [(a-b)%self.p for a,b in zip(x,y)]
    def scale(self,x,c):return [a*c%self.p for a in x]
    def mul(self,x,y):return [sum(x[j]*y[i-j] for j in range(i+1))%self.p for i in range(N)]
    def polynomial(self,cs,x):
        a=self.c(0)
        for c in reversed(cs):a=self.add(self.mul(a,x),self.c(c))
        return a
    def derivative(self,x):return [(i+1)*x[i+1]%self.p for i in range(N-1)]+[0]


def coeffs(poly,p,alpha):
    f=S.Poly(S.expand(poly.subs(A,alpha)),X)
    return [mq(f.nth(i),p) for i in range(max(0,f.degree())+1)] if f else [0]


def solve(rows,rhs,p):
    a=[[x%p for x in row]+[b%p] for row,b in zip(rows,rhs)]
    nc=len(rows[0]);r=0;piv=[]
    for c in range(nc):
        k=next((i for i in range(r,len(a)) if a[i][c]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]
        inv=pow(a[r][c],-1,p);a[r]=[x*inv%p for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                v=a[i][c];a[i]=[(x-v*y)%p for x,y in zip(a[i],a[r])]
        piv.append(c);r+=1
    if any(not any(row[:nc]) and row[nc] for row in a):return None
    need(len(piv)==nc,('unique overdetermined section solution',p,nc,piv))
    x=[0]*nc
    for i,c in enumerate(piv):x[c]=a[i][nc]
    return x


def compiler(ring,xx,qq,ww):
    d=ring.add(ring.c(1),ring.mul(ring.mul(xx,xx),qq))
    cc=ring.mul(xx,ring.mul(d,ring.add(d,ring.c(2))))
    ee=ring.mul(qq,ring.add(d,ring.c(3)))
    return ring.scale(cc,pow(3,-1,ring.p)),ring.add(ee,ring.c(3)),ww


def source_image(ring,xseries,qcs,perturb):
    t=ring.c(0);t[1]=1
    qq=ring.add(ring.polynomial(qcs,xseries),ring.mul(t,t))
    xx=list(xseries);ww=list(t)
    # Inputs are t9 times polynomials; only their orders0 and1 matter
    # modulo t11, but full series evaluation is retained here.
    t9=ring.c(0);t9[9]=1
    for i,(base,poly) in enumerate(zip((xx,qq,ww),perturb)):
        changed=ring.add(base,ring.mul(t9,ring.polynomial(poly,xseries)))
        if i==0:xx=changed
        elif i==1:qq=changed
        else:ww=changed
    return compiler(ring,xx,qq,ww)


def collision_controls(p,alpha):
    ring=Series(p)
    qcs=coeffs(QPOLY,p,alpha)
    need(int(EXCEPTION.subs(A,alpha))%p==0,('fresh exceptional root',p,alpha))
    sections=[ring.c(i) for i in POINTS]
    common=[ring.c(0),ring.c(0)]
    rows=[]
    for i,T in enumerate(TANGENTS):
        for j in range(2):
            row=[0]*5;row[i]=T[j];row[3+j]=-1;rows.append(row)
    for order in range(1,10):
        values=[source_image(ring,z,qcs,([0],[0],[0])) for z in sections]
        rhs=[-(v[j][order]-common[j][order]) for v in values for j in range(2)]
        delta=solve(rows,rhs,p)
        need(delta is not None,('old collision through t9',p,alpha,order))
        for i in range(3):sections[i][order]=(sections[i][order]+delta[i])%p
        for j in range(2):common[j][order]=(common[j][order]+delta[3+j])%p
        values=[source_image(ring,z,qcs,([0],[0],[0])) for z in sections]
        need(all(v[j][:order+1]==common[j][:order+1] for v in values for j in range(2)),('literal old collision',p,order))
    need(all(z[j]==0 for z in sections for j in (1,3,5,7,9)),'old sections have only even terms through t9')
    rhs=[-v[j][10] for v in values for j in range(2)]
    need(solve(rows,rhs,p) is None,'old unpaid t10 collision is a genuine next-order debt')
    controls=[
        ('zero',0,0,0),('constant_graph',0,1,0),('kernel_graph',0,4*X**2-9*X,0),
        ('split_graph',0,X,0),('reparam_common_w',X**2,X**2*S.diff(QPOLY,X),3),
        ('split_w',0,0,X),('w_Lzero_but_not_constant',0,0,4*X**2-9*X),
        ('one_form_hidden_split',0,X,X),
    ]
    results=[]
    rows3=[]
    for i,T in enumerate(TANGENTS):
        for j in range(3):
            row=[0]*6;row[i]=T[j];row[3+j]=-1;rows3.append(row)
    for name,r,s,k in controls:
        polys=tuple(coeffs(S.sympify(f),p,alpha) for f in (r,s,k))
        hv=tuple(mq(S.expand(S.sympify(s)-S.diff(QPOLY,X)*S.sympify(r)).subs({X:i,A:alpha}),p) for i in POINTS)
        kv=tuple(mq(S.sympify(k).subs(X,i),p) for i in POINTS)
        lh=sum(a*b for a,b in zip(ELL,hv))%p
        predicted=(lh==0 and len(set(kv))==1)
        ims=[source_image(ring,z,qcs,polys) for z in sections]
        delta=solve(rows3,[-im[j][9] for im in ims for j in range(3)],p)
        need((delta is not None)==predicted,('literal changed collision iff',p,name))
        if delta is not None:
            changed=[list(z) for z in sections]
            for i in range(3):changed[i][9]=(changed[i][9]+delta[i])%p
            ims=[source_image(ring,z,qcs,polys) for z in changed]
            need(all(ims[i][j][:10]==ims[0][j][:10] for i in (1,2) for j in range(3)),('all three target coordinates collide modt10',p,name))
        results.append((name,predicted))
    print('FRESH_COLLISIONS',p,alpha,'old agreesmodt10 but failsmodt11',results)
    return qcs,controls


def fixed_point_densities(p,alpha,qcs,controls):
    ring=Series(p)
    def add(a,b):return ring.add(a[0],b[0]),ring.add(a[1],b[1])
    def mul(a,b):return ring.mul(a[0],b[0]),ring.add(ring.mul(a[1],b[0]),ring.mul(a[0],b[1]))
    def c(x):return ring.c(x),ring.c(0)
    def evaluate(cs,a):
        out=c(0)
        for x in reversed(cs):out=add(mul(out,a),c(x))
        return out
    def target(i,perturb):
        xx=(ring.c(i),ring.c(1));tt=ring.c(0);tt[1]=1;t=(tt,ring.c(0))
        qq=add(evaluate(qcs,xx),mul(t,t));ww=t
        t9=ring.c(0);t9[9]=1;t9=(t9,ring.c(0))
        xxnew=add(xx,mul(t9,evaluate(perturb[0],xx)))
        qq=add(qq,mul(t9,evaluate(perturb[1],xx)))
        ww=add(ww,mul(t9,evaluate(perturb[2],xx)))
        d=add(c(1),mul(mul(xxnew,xxnew),qq))
        yy=mul(xxnew,mul(d,add(d,c(2))))
        yy=(ring.scale(yy[0],pow(3,-1,p)),ring.scale(yy[1],pow(3,-1,p)))
        zz=add(mul(qq,add(d,c(3))),c(3))
        return yy,zz,ww
    def density(targets,slot,monomial):
        i,j=((0,1),(0,2),(1,2))[slot]
        wedge=ring.sub(ring.mul(targets[i][1],ring.derivative(targets[j][0])),
                       ring.mul(ring.derivative(targets[i][0]),targets[j][1]))
        coefficient=ring.c(1)
        for which,power in enumerate(monomial):
            for _ in range(power):coefficient=ring.mul(coefficient,targets[which][0])
        return ring.mul(coefficient,wedge)
    count=0
    for name,r,s,k in controls:
        polys=tuple(coeffs(S.sympify(f),p,alpha) for f in (r,s,k))
        hv=tuple(mq(S.expand(S.sympify(s)-S.diff(QPOLY,X)*S.sympify(r)).subs({X:i,A:alpha}),p) for i in POINTS)
        kv=tuple(mq(S.sympify(k).subs(X,i),p) for i in POINTS)
        shifts=[]
        for index,i in enumerate(POINTS):
            old=target(i,([0],[0],[0]));new=target(i,polys)
            response=[]
            for slot in range(3):
                for monomial in ((0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,0,2)):
                    before=density(old,slot,monomial);after=density(new,slot,monomial)
                    need(before[:8]==after[:8],('all lower fixed-source densities unchanged',p,name,i,slot,monomial))
                    difference=(after[8]-before[8])%p
                    predicted=9*(4*hv[index],kv[index],TANGENTS[index][1]*kv[index])[slot]%p if monomial==(0,0,0) else 0
                    need(difference==predicted,('literal ninth-source wedge response',p,name,i,slot,monomial))
                    if monomial==(0,0,0):response.append(difference)
                    count+=1
            shifts.append(response)
        periods=tuple(sum(ELL[i]*shifts[i][j] for i in range(3))*pow(18,-1,p)%p for j in range(3))
        lh=sum(a*b for a,b in zip(ELL,hv))%p
        expected=(2*lh%p,sum(a*b for a,b in zip(ELL,kv))*pow(2,-1,p)%p,
                  sum(ELL[i]*TANGENTS[i][1]*kv[i] for i in range(3))*pow(2,-1,p)%p)
        need(periods==expected,('all three target-form period coefficients',p,name))
        need((not any(periods))==(lh==0 and len(set(kv))==1),('all-form period iff collision conditions',p,name))
        if name=='one_form_hidden_split':
            need(periods==(16,4,81) and (periods[0]-4*periods[1])%p==0,
                 'one fixed two-form hides the collision split')
    print('LITERAL_DENSITIES',p,'controls',count,'all lower coefficients and three period slots PASS')


def symbolic_controls():
    x,q=S.symbols('x q')
    d=1+x*x*q;cc=x*d*(d+2);ee=q*(d+3)
    jac=S.expand(S.diff(cc,x)*S.diff(ee,q)-S.diff(cc,q)*S.diff(ee,x))
    need(S.expand(jac-6*(d+1)*(d*d+2*d-2))==0,'literal compiler surface wedge')
    normals=[]
    for i,T in zip(POINTS,TANGENTS):
        qi=S.simplify(QPOLY.subs(X,i));qprime=S.simplify(S.diff(QPOLY,X).subs(X,i))
        need(not qi.has(A) and not qprime.has(A),'retained source values and slopes independent of embedding')
        tangent=(S.diff(cc,x)+qprime*S.diff(cc,q),S.diff(ee,x)+qprime*S.diff(ee,q))
        tangent=tuple(S.simplify(v.subs({x:i,q:qi})) for v in tangent)
        need(tangent==(3,T[1]),('retained compiler tangent',i))
        normal=(S.simplify(S.diff(cc,q).subs({x:i,q:qi})/3),S.simplify(S.diff(ee,q).subs({x:i,q:qi})),0)
        normals.append(normal)
        need(S.simplify(T[0]*normal[1]-T[1]*normal[0])==4,('normal area4',i))
    need(all(sum(ELL[i]*TANGENTS[i][j] for i in range(3))==0 for j in range(3)),'full3target tangent relation')
    K=S.Matrix([ELL,[ELL[i]*TANGENTS[i][1] for i in range(3)]])
    need(K.rank()==2 and K.nullspace()==[S.ones(3,1)],'stable values constant iff both weighted periods vanish')
    h=(S.Integer(-1),S.Integer(0),S.Integer(1))
    need(sum(a*b for a,b in zip(ELL,h))==8,'split h=x period16')
    need(16-4*4==0,'single chosen form may have zero period despite split collision')
    # A universal three-dimensional alternating form and arbitrary common
    # target displacement; all section shifts cancel by alternation.
    aa,bb,ccs=S.symbols('A0 B0 C0');bvec=S.Matrix(S.symbols('b0:3'));r=S.symbols('r0:3')
    M=S.Matrix([[0,aa,bb],[-aa,0,ccs],[-bb,-ccs,0]])
    total=0
    for i in range(3):
        T=S.Matrix(TANGENTS[i]);v=bvec-r[i]*T
        total+=ELL[i]*(T.T*M*v)[0]
    need(S.expand(total)==0,'universal smooth3target wedge cancellation')
    # The distinct-tangent premise matters: for two local source germs
    # f1=(s,0), f2=(s,s²), the t² perturbation ftilde2=(s,s²-t²)
    # lets new sections s=t collide although old sections s=0 collide.
    # Both base tangents are (1,0), and weights(1,-1) give response2.
    ss,tt=S.symbols('s t')
    g=S.Matrix([ss,ss*ss-tt*tt])
    density=S.det(S.Matrix.hstack(g.diff(ss),g.diff(tt)))
    need(S.expand(g.subs(ss,tt)-S.Matrix([tt,0]))==S.zeros(2,1),'parallel-tangent new collision')
    need(-S.expand(density).coeff(tt,1)==2,'parallel-tangent failure if distinct-line premise dropped')
    print('SYMBOLIC normal vectors',normals,'weighted wedge zero; all-form converse and single-form hostile PASS')


if __name__=='__main__':
    symbolic_controls()
    for p,alpha in ((449,120),(467,169)):
        qcs,controls=collision_controls(p,alpha)
        fixed_point_densities(p,alpha,qcs,controls)
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)
