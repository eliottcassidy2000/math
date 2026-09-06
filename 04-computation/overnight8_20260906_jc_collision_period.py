"""Exact universal tangent/collision/wedge identities, no repo imports."""
import sympy as S
import sys,hashlib,json
sys.stdout.reconfigure(newline='\n')
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
x,q,t=S.symbols('x q t');aa,bb,cc=S.symbols('A B C')
D=1+x*x*q;Y=x*D*(D+2)/3;Z=q*(D+3)+3
baseq=[-3,S.Rational(-3,4),-3]
baseqp=[-S.Rational(9,2),S.Integer(1),S.Rational(9,2)]
points=[-1,0,1];ell=[5,-18,13]
T=[];U=[]
for i,point in enumerate(points):
    at={x:point,q:baseq[i]}
    row=[S.expand((S.diff(F,x)+baseqp[i]*S.diff(F,q)).subs(at)) for F in (Y,Z)]+[0]
    normal=[S.expand(S.diff(F,q).subs(at)) for F in (Y,Z)]+[0]
    T.append(S.Matrix(row));U.append(S.Matrix(normal))
    need([S.expand(F.subs(at)) for F in (Y,Z)]==[0,0],'same literal target point')
    need(row==[1,[-9,4,9][i],0],'literal source tangent')
    need(row[0]*normal[1]-row[1]*normal[0]==4,'actual normal wedge coefficient')
need(sum((ell[i]*T[i] for i in range(3)),S.zeros(3,1))==S.zeros(3,1),'tangent relation')
def wedge(u,v):
    return aa*(u[0]*v[1]-u[1]*v[0])+bb*(u[0]*v[2]-u[2]*v[0])+cc*(u[1]*v[2]-u[2]*v[1])
rv=S.symbols('r:3');hv=S.symbols('h:3');kv=S.symbols('k:3')
V=[rv[i]*T[i]+hv[i]*U[i]+S.Matrix([0,0,kv[i]]) for i in range(3)]
L=lambda v:sum(u*w for u,w in zip(ell,v))
period=S.expand(sum(S.Rational(ell[i],2)*wedge(T[i],V[i]) for i in range(3)))
formula=2*aa*L(hv)+(bb*L(kv)+cc*L([T[i][1]*kv[i] for i in range(3)]))/2
need(S.expand(period-formula)==0,'all source components and arbitrary form coefficients')
b=S.Matrix(S.symbols('b:3'));N=S.symbols('N')
need(S.expand(sum(ell[i]*N*wedge(T[i],b-rv[i]*T[i]) for i in range(3)))==0,'universal common-motion wedge identity')
# Explicit complete nine-equation collision matrix: three section variables
# and one common target vector. Its full left nullspace has dimension three.
M=S.zeros(9,6)
for i in range(3):
    for j in range(3):M[3*i+j,i]=T[i][j];M[3*i+j,3+j]=-1
need(M.rank()==6,'all ninth collision unknowns, rank six')
left=M.T.nullspace();need(len(left)==3,'exactly three collision conditions')
forcing=S.Matrix([v for row in V for v in row])
vars=list(hv)+list(kv)+list(rv)
conditions=[S.expand((v.T*forcing)[0]) for v in left]
actual=S.Matrix([[s.coeff(v) for v in vars] for s in conditions])
desired=[L(hv),kv[1]-kv[0],kv[2]-kv[0]]
expected=S.Matrix([[s.coeff(v) for v in vars] for s in desired])
need(actual.rank()==expected.rank()==actual.col_join(expected).rank()==3,'necessary and sufficient h-period / stable equality')
normalM=M.extract([0,1,3,4,6,7],[0,1,2,3,4])
need(normalM.rank()==5,'graph rank-five repair')
for h in (S.Integer(1),4*x*x-9*x,x*(x*x-1),S.Integer(0)):
    v=[h.subs(x,z) for z in points]
    need(L(v)==0,'lawful graph kernel control')
    vals=dict(zip(hv,v))|dict(zip(kv,[0,0,0]))
    need(S.expand(period.subs(vals))==0,'actual leading period of lawful graph')
need(L(points)==8,'collision-splitting x hostile')
need(S.expand(period.subs(dict(zip(hv,points))|dict(zip(kv,[0,0,0])))-16*aa)==0,'sharp h=x period')
km=S.Matrix([ell,[ell[i]*T[i][1] for i in range(3)]])
need(km.rank()==2 and km*S.ones(3,1)==S.zeros(2,1),'all-form converse detects stable equality')
allform=S.Matrix([[S.expand(period).coeff(c).coeff(v) for v in vars] for c in (aa,bb,cc)])
need(allform.rank()==allform.col_join(expected).rank()==3,'complete period and collision cokernels coincide')
single=S.expand(period.subs(dict(zip(hv,points))|dict(zip(kv,points))))
need(single==16*aa+4*bb+81*cc and single.subs({aa:1,bb:-4,cc:0})==0,'one chosen form hides a collision split')
# Independent literal local compiler expansion of the first density response.
# A complete base 1-jet in x is enough for the retained values. Use a local
# coordinate u and an arbitrary ninth source displacement r,s,k.
u=S.symbols('u')
records=[]
for i,point in enumerate(points):
    for r,s,k in ((0,1,0),(0,point,0),(1,0,0),(0,0,1),(2,3,-1)):
        xx=point+u+t**9*r
        qq=baseq[i]+baseqp[i]*u+t*t+t**9*s
        ww=t+t**9*k
        yy=Y.subs({x:xx,q:qq});zz=Z.subs({x:xx,q:qq})
        d=[S.diff(F,u).subs(u,0) for F in (yy,zz,ww)]
        e=[S.diff(F,t).subs(u,0) for F in (yy,zz,ww)]
        full=S.expand(wedge(d,e));j8=full.coeff(t,8)
        # The unchanged constant-pencil density has no t8 term in this
        # base1-jet test; subtract it explicitly rather than assume parity.
        xx0=point+u;qq0=baseq[i]+baseqp[i]*u+t*t
        ff=[Y.subs({x:xx0,q:qq0}),Z.subs({x:xx0,q:qq0}),t]
        old=S.expand(wedge([S.diff(F,u).subs(u,0) for F in ff],[S.diff(F,t).subs(u,0) for F in ff]))
        delta=S.expand(j8-old.coeff(t,8))
        hh=s-baseqp[i]*r
        prediction=36*aa*hh+9*k*(bb+cc*T[i][1])
        need(S.expand(delta-prediction)==0,'literal compiler first visible density, all slots')
        need(all(S.expand(full.coeff(t,j)-old.coeff(t,j))==0 for j in range(8)),'every lower density coefficient unchanged')
        records.append([i,r,s,k,str(delta)])
print('Complete collision matrix rank6, cokernel3; graph rank5, cokernel1')
print('Universal wedge-period identity: common ninth collision motion annihilates delta Lambda(J8)')
print('Source normal L(h)=0 and stable retained values equal are necessary and sufficient at order9')
print('Literal all-slot density controls',len(records),'sha256',hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
print('PASS_OPTIMIZATION_LIVE_GATES',GATES)
