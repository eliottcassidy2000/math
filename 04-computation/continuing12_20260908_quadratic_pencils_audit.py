"""Independent exact referee: full W2 quadratic boxes and three entire pencils.

No mathematical producer is imported.  Always-active gates survive -O.
The analytic proof, including unbounded PBW equality, is in the paired report.
"""
from pathlib import Path
import hashlib,json,sys
import sympy as S
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
DEST=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
u,t,x,r,b,a,d,h,k,s,z,Y,g,T,v,p,q=S.symbols('u t x r b a d h k s z Y g T v p q')
gates=0
def need(ok,label):
    global gates
    gates+=1
    if not ok: raise ArithmeticError(label)
def eq(A,B,label): need(S.cancel(A-B)==0,label)
def jac(A,B):return S.diff(A,u)*S.diff(B,t)-S.diff(A,t)*S.diff(B,u)
def global_Q(N,P):
    nx=S.Poly(N.subs(u,x-h),x);px=S.Poly(P.subs(u,x-h),x)
    ans=S.expand((nx.nth(8)*x**4+nx.nth(7)*x**3+(px.nth(4)-nx.nth(6))*x*x+(px.nth(3)-nx.nth(5))*x).subs(x,u+h))
    return S.expand(ans-ans.subs(u,0)+s)
def chart(A):return S.cancel(A.subs({u:1/r-h,t:-r*r-r**4*b},simultaneous=True))
def poly_coeff(E,X,n):return S.expand(E).coeff(X,n)

# Recover the complete global coefficient kernel directly from negative Laurent terms.
nc=S.symbols('n0:9');pc=S.symbols('p0:7');qc=S.symbols('q0:5')
AA=sum(nc[i]*x**i for i in range(9));BB=sum(pc[i]*x**i for i in range(7));CC=sum(qc[i]*x**i for i in range(5))
raw=S.expand((AA*t*t+BB*t+CC).subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
negative=[]
for e in range(-4,0):
    for j in range(3):
        E=raw.coeff(r,e).coeff(b,j)
        if E!=0:negative.append(E)
M,_=S.linear_eq_to_matrix(negative,nc+pc+qc)
need(M.rank()==6 and len(nc+pc+qc)-M.rank()==15,'entire 15-dimensional Laurent kernel')
forced={pc[6]:2*nc[8],pc[5]:2*nc[7],qc[4]:nc[8],qc[3]:nc[7],qc[2]:pc[4]-nc[6],qc[1]:pc[3]-nc[5]}
for E in negative:eq(E.subs(forced),0,'every negative Laurent coefficient paid')
need(len(S.solve(negative,list(forced),dict=True))==1,'complete linear kernel has a unique solved chart')
for V,E in S.solve(negative,list(forced),dict=True)[0].items():eq(E,forced[V],'raw kernel matches recovered global box')

# Third pencil: divisibility bounds all possible P before coefficient elimination.
N3=u**3*(a*u+d);P3=u*u*(p*u*u+q*u+v);Q3=global_Q(N3,P3)
D3=S.expand(P3*P3-4*N3*Q3)
eq(poly_coeff(D3,u,8),p*p,'u3 first forced square')
eq(poly_coeff(D3.subs(p,0),u,6),q*q,'u3 second forced square')
F3=(N3*t*t+P3*t+Q3).subs({p:0,q:0})
eq(S.diff(F3,u).subs(u,0),0,'u3 actual entire affine critical line')
eq(S.diff(F3,t).subs(u,0),0,'u3 actual second gradient entry')
disc3=S.expand(D3.subs({p:0,q:0}))
for j in [0,1,2,5,6,7,8]:eq(disc3.coeff(u,j),0,'u3 complete discriminant support')
eq(S.det(S.Matrix([[d,a],[disc3.coeff(u,3),disc3.coeff(u,4)]])),d*v*v,'u3 exact independence determinant')

# Fourth pencil: global boxes followed by the whole discriminant support.
N4=u**4*(a*u*u+d);P4=u*u*(p*u*u+q*u+v);Q4=global_Q(N4,P4)
D4=S.expand(P4*P4-4*N4*Q4)
eq(D4.coeff(u,8),(p-2*a)**2,'u4 forced square')
disc4=S.expand(D4.subs(p,2*a))
eq(disc4.coeff(u,5),2*q*(v-2*d),'u4 exhaustive branch equation')
for j in [0,1,2,3,7,8]:eq(disc4.coeff(u,j),0,'u4 all other off-pencil coefficients')
eq(a*disc4.coeff(u,4)-d*disc4.coeff(u,6),a*(v-2*d)**2-d*q*q,'u4 exact rank condition')
F4=S.expand((N4*t*t+P4*t+Q4).subs(p,2*a))
F4a=F4.subs({q:0,v:2*d+k});F4b=F4.subs(v,2*d)
for diffvar in [u,t]:eq(S.diff(F4a,diffvar).subs(u,0),0,'u4 branch A omitted source line paid')
eq(S.diff(F4b,u).subs(u,0),q,'u4 branch B source line derivative')
Z=S.symbols('Z')
open4=s-d+(a*u*u+d)*Z*Z+q*u*Z
eq(F4b,open4.subs(Z,1+u*u*t),'u4 genuine localized coordinate change')
eq(S.diff(open4,u),Z*(2*a*u*Z+q),'u4 open first derivative')
eq(S.diff(open4,Z).subs(Z,-q/(2*a*u)),-d*q/(a*u),'u4 open second derivative on nonzero first-zero branch')
inf4=chart(F4b)
eq(S.diff(inf4,b).subs(r,0),0,'u4 whole boundary tangential derivative')
eq(S.diff(inf4,r).subs(r,0),-(4*a*h+q)*(3*h*h+b),'u4 boundary normal derivative')
for V in [r,b]:eq(S.diff(inf4,V).subs({r:0,b:-3*h*h}),0,'u4 actual boundary critical point for every parameter')
eq(S.diff(inf4,r).subs({r:0,q:-4*a*h}),0,'u4 matching wall is whole critical divisor')

# Seventh pencil: derive P's two coefficients from the recovered linear box.
c6,c5,c4=S.symbols('c6 c5 c4')
N7=u**7*(a*u+d);P7trial=u**4*(c6*u*u+c5*u+c4)
nx=S.Poly(N7.subs(u,x-h),x);px=S.Poly(P7trial.subs(u,x-h),x)
sol=S.solve([px.nth(6)-2*nx.nth(8),px.nth(5)-2*nx.nth(7)],[c6,c5],dict=True)
need(len(sol)==1,'u7 unique mandatory top matching')
eq(sol[0][c6],2*a,'u7 degree-six coefficient')
eq(sol[0][c5],2*d-4*a*h,'u7 degree-five coefficient')
P7=P7trial.subs(sol[0]).subs(c4,k-4*d*h);Q7=global_Q(N7,P7)
F7=S.expand(N7*t*t+P7*t+Q7)
zs=u-2*h+u**3*t;L=(a*u+d)*z+k;f=u*z*L
eq(F7,s+f.subs(z,zs),'u7 exhausted family factors in original coordinates')
eq(P7*P7-4*N7*(Q7-s-g),u**7*((k*k+4*a*g)*u+4*d*g),'u7 whole discriminant identity')
eq(S.det(S.Matrix([[d,a],[0,k*k]])),d*k*k,'u7 exact whole-pencil independence')
need(S.Poly(chart(F7),r,b).is_multivariate,'u7 complete second chart is polynomial')
R=3*h*h+b
inf7=chart(F7)
eq(inf7.subs(r,0),s+a*R*R-k*R,'u7 boundary restriction from literal substitution')
eq(S.diff(inf7,b).subs(r,0),2*a*R-k,'u7 all boundary tangencies')
eq(S.diff(inf7,r).subs(r,0).subs(k,2*a*R),d*R*R,'u7 boundary normal on tangency locus')
eq(S.diff(inf7,b).subs({r:0,a:0}),-k,'u7 a=0 degree drop remains smooth')
eq(S.diff(F7,u).subs(u,0),-2*h*(k-2*d*h),'u7 necessary and sufficient source-line condition')
eq(S.diff(F7,t).subs(u,0),0,'u7 source-line tangency paid')
fu=s+(a+d/u)*Y*Y+k*Y
eq(F7,fu.subs(Y,u*zs),'u7 actual open coordinates')
eq(jac(u,u*zs),u**4,'u7 coordinate Jacobian, no missing line')
eq(S.diff(fu,u),-d*Y*Y/u**2,'u7 first open derivative forces Y0')
eq(S.diff(fu,Y).subs(Y,0),k,'u7 second open derivative excludes open critical points')
for sub in [{h:0},{k:2*d*h}]:
    for V in [u,t]:eq(S.diff(F7,V).subs(sub).subs(u,0),0,'u7 excluded parameter has whole source critical line')

# Derive a rational mate by elementary integration in the recovered field.
inv=d*Y*Y/(g-a*Y*Y-k*Y)
eq(fu.subs(u,inv),s+g,'u7 field inverse over target')
dy=jac(F7,u*zs)
eq(dy,-d*u*u*(u*zs)**2,'u7 actual field derivation')
integrand=-(g-a*Y*Y-k*Y)**2/(d**3*Y**6)
primitive=S.integrate(S.expand(integrand),Y)
eq(S.diff(primitive,Y),integrand,'u7 independent primitive by integration')
mate=S.cancel(primitive.subs({g:f,Y:u*z},simultaneous=True))
eq(jac(F7,mate.subs(z,zs)),1,'u7 literal rational mate')
H=S.cancel(f**3*mate)
need(S.denom(H)==30*d**3,'u7 cubic witness has no source denominators')
eq(jac(F7,H.subs(z,zs)),(F7-s)**3,'u7 literal cubic polynomial response')
eq(zs.subs(u,0),-2*h,'u7 Eu/Ez separation')
eq(L.subs({u:0,z:-2*h}),k-2*d*h,'u7 Eu/EL separation')
eq(L.subs(z,0),k,'u7 Ez/EL separation')
eq(L.subs(z,zs).subs(u,-d/a),k,'u7 all possible linear-factor content paid')
eq(P7.subs(u,-d/a),k*d**4/a**4,'u7 only possible N/P gcd support is u0')
eq((F7-s-g).subs(u,0),-g,'u7 no vertical factor on any other source fibre')
for sub,depth in [({a:0,d:1,h:1,k:4},5),({a:1,d:2,h:1,k:8},6)]:
    eq(S.Poly(S.gcd(N7.subs(sub),P7.subs(sub)),u).monic().as_expr(),u**depth,'admitted gcd exponent jump retains only u0 support')
    need(S.cancel((d*k*h*(k-2*d*h)).subs(sub))!=0,'gcd-jump control remains globally smooth')

# Independently extract full component coefficients by triangular Taylor jets.
def extract(local,ff,HH):
    f1=S.diff(ff,local).subs(local,0);f2=S.diff(ff,local,2).subs(local,0)/2
    H0=HH.subs(local,0);H1=S.diff(HH,local).subs(local,0);H2=S.diff(HH,local,2).subs(local,0)/2
    A3=S.cancel(H0);A2=S.cancel(H1/f1);A1=S.cancel((H2-A2*f2)/f1**2)
    return [A1,A2,A3]
Cu=extract(u,f.subs(z,zs),H.subs(z,zs));Cz=extract(z,f,H)
expected_u=[a*a*k/d**3,(k-2*d*h)*(4*a*d*d*h*h+2*a*d*h*k+a*k*k-6*d**3*h)/(3*d**3),(k-2*d*h)**3*(k*k+6*d*h*k+24*d*d*h*h)/(30*d**3)]
expected_z=[a*a*k/d**3,a*k**3/(3*d**3),k**5/(30*d**3)]
for i in range(3):
    eq(Cu[i],expected_u[i],'u7 actual Eu jet recovered independently')
    eq(Cz[i],expected_z[i],'u7 complete Ez jet recovered independently')
    need(not Cu[i].has(t),'u7 Eu residue is scalar on its full component')
for local,ff,HH,cs in [(u,f.subs(z,zs),H.subs(z,zs),Cu),(z,f,H,Cz)]:
    rem=HH-cs[2]-cs[1]*ff-cs[0]*ff*ff
    for j in range(3):eq(S.diff(rem,local,j).subs(local,0),0,'u7 entire local cubic remainder')
need(Cz[2]!=0,'u7 nonzero highest pole coefficient remains on Ez')
eq((f*f*mate*z).subs(z,0),k**4/(30*d**3*u),'u7 strict order-three lower witness is a genuine Ez pole')
wrong=extract(u,f.subs(z,-2*h),H.subs(z,-2*h))
need(S.cancel((wrong[1]-Cu[1]).subs({a:1,d:1,h:1,k:1}))!=0,'hostile: dropping actual +u jet changes scalar second coefficient')

# Complete rank wall and smoothness, including both complex branches and a=0.
C=S.Matrix([Cu,Cz])
eq(Cu[1]-Cz[1],-2*h*(4*a*h*h-6*d*h+3*k)/3,'u7 second-coefficient alignment equation')
eq(Cu[2]-Cz[2],-8*h**3*(12*d*d*h*h-15*d*h*k+5*k*k)/15,'u7 third-coefficient alignment equation')
eq(S.det(C[:,1:3]).subs(a,0),-h*k**5*(k-2*d*h)/(15*d**3),'u7 entire a0 locus remains rank two')
wall=5*q*q-15*q+12
wallsub={d:1,h:1,k:q,a:3*(2-q)/4}
def quotient_zero(E,label):
    num,den=S.fraction(S.cancel(E.subs(wallsub)))
    eq(S.rem(num,wall,q),0,label)
    need(S.gcd(den,wall)==1,'wall quotient denominator is unit on both branches')
for i in range(3):quotient_zero(Cu[i]-Cz[i],'both exact complex branches align')
for E in [a,d,h,k,k-2*d*h,Cz[0],Cz[2]]:
    nn,dd=S.fraction(S.cancel(E.subs(wallsub)))
    need(S.gcd(nn*dd,wall)==1,'both complex wall branches satisfy every smooth/nonzero condition')
eq(S.discriminant(12*d*d*h*h-15*d*h*k+5*k*k,k),-15*d*d*h*h,'no real rank-one wall for nonzero real dh')
drop=q*q+6*q+24
drop_sub={a:1,d:1,h:1,k:q}
nn,dd=S.fraction(S.cancel(Cu[2].subs(drop_sub)))
eq(S.rem(nn,drop,q),0,'complex Eu leading pole-drop is retained')
for E in [Cu[1],Cz[2],d*k*h*(k-2*d*h)]:
    nn,dd=S.fraction(S.cancel(E.subs(drop_sub)))
    need(S.gcd(nn*dd,drop)==1,'Eu pole-drop keeps order2 there, order3 on Ez, and global smoothness')
need(S.gcd(wall,drop)==1,'Eu pole-drop is distinct from coefficient-alignment wall')

# Formal Weyl actions use actual finite Laurent arrays, independent of PBW formulas.
def norm(A):return {j:S.simplify(S.expand(vv)) for j,vv in A.items() if S.simplify(S.expand(vv))!=0}
def plus(*AA):
    ans={}
    for A0 in AA:
        for j,vv in A0.items():ans[j]=ans.get(j,0)+vv
    return norm(ans)
def scale(A0,c):return norm({j:c*vv for j,vv in A0.items()})
def mult(A0,n):return norm({j-n:vv for j,vv in A0.items() if j>n})
def derivative(A0,n):
    return norm({j+n:vv*(-1)**n*S.rf(j,n) for j,vv in A0.items()})
def polyapply(poly,A0):
    return plus(*(scale(derivative(A0,mon[0]),coeff) for mon,coeff in S.Poly(poly,T).terms()))
def op(pp,A0):return plus(*(polyapply(poly,mult(A0,j)) for j,poly in enumerate(pp)))
cc1,cc2,cc3=S.symbols('cc1 cc2 cc3');theta={1:cc1,2:cc2,3:cc3}
for j in range(3):
    for degree in range(7):
        pp=[S.Integer(0)]*3;pp[j]=T**degree
        q0,q1,q2=pp[0],pp[1]-T*pp[0],pp[2]-T*pp[1]+T*T*pp[0]/2
        predicted=plus(polyapply(q0,{1:cc1}),polyapply(q1,{1:cc2}),polyapply(q2,{1:cc3}))
        need(op(pp,theta)==predicted,'independent formal Laurent action equals full PBW remainder map')
        back=[q0,q1+T*q0,q2+T*q1+T*T*q0/2]
        for i in range(3):eq(back[i],pp[i],'PBW triangular inverse has no missing coefficient')
ll=S.symbols('l1:4')
compiled=[ll[0],ll[1]+ll[0]*T,ll[2]+ll[1]*T+ll[0]*T*T/2]
need(op(compiled,theta)=={1:sum(ll[i]*[cc1,cc2,cc3][i] for i in range(3))},'universal compiled relation sends theta to its exact coefficient relation')
need(mult(theta,3)=={},'g3 always annihilates cubic Laurent packet')
need(op([0,0,1],theta)!={},'g2 does not annihilate cubic packet')
CW=C.subs(wallsub)
wallrels=[[1,0,-CW[0,0]/CW[0,2]],[0,1,-CW[0,1]/CW[0,2]]]
for l in wallrels:
    rr=[l[0],l[1]+l[0]*T,l[2]+l[1]*T+l[0]*T*T/2]
    for row in range(2):
        packet={j+1:CW[row,j] for j in range(3)}
        for n in range(4):
            actual=op([T**n*entry for entry in rr],packet)
            for E in actual.values():quotient_zero(E,'two compiled wall generators and their derivatives annihilate both components')

controls=[]
for aa,dd,hh,kk in [(0,1,1,1),(1,1,1,1),(1,2,1,-1),(2,-1,2,3),(3,1,-1,2)]:
    sub={a:aa,d:dd,h:hh,k:kk};CM=C.subs(sub)
    need(dd*kk*hh*(kk-2*dd*hh)!=0,'independent real control lies in full smooth locus')
    need(CM.rank()==2,'independent real control has two coefficient directions')
    rel=CM.nullspace();need(len(rel)==1,'rank-two control has exactly one constant relation')
    l=list(rel[0]);RR=[l[0],l[1]+l[0]*T,l[2]+l[1]*T+l[0]*T*T/2]
    for row in range(2):
        packet={j+1:CM[row,j] for j in range(3)}
        for n in range(5):need(op([T**n*rr for rr in RR],packet)=={},'all chosen differentiated compiled generators annihilate actual component')
    for n in range(5):
        leading=derivative({j+1:CM[1,j] for j in range(3)},n)
        need(max(leading)==3+n,'actual canonical derivative has exact shifted scalar order')
    controls.append([aa,dd,hh,kk])

cert={'status':'INDEPENDENT EXACT CONTROLS; analytic proof in paired report',
      'complete_global_dimension':15,'pencils':['u3 span(1,u)','u4 span(1,u2)','u7 span(1,u)'],
      'seventh_smooth_locus':'d*k*h*(k-2*d*h)!=0; a,s arbitrary',
      'unit_order':3,'full_source_arms':2,'rank_one_locus':['4*a*h*h-6*d*h+3*k=0','12*d*d*h*h-15*d*h*k+5*k*k=0'],
      'controls':controls,'gates':gates,'imports_producer':False,
      'scope':'fixed W2; original source response ring C[x,t]; no other pencil or global mate claim'}
raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('INDEPENDENT PASS: complete 15-dimensional global coefficient kernel and all three declared pencils.')
print('INDEPENDENT PASS: all source and boundary points, field primitive, all affine fibre components and full cubic jets.')
print('INDEPENDENT PASS: two complex rank-one branches, exact coefficient relation compiler, and canonical derivative controls.')
print('CERTIFICATE_SHA256',hashlib.sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
