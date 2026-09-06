"""Exact all-phase rectangle certificate; standard library, raw LF output.

The original Laurent convolutions are reconstructed before eliminating z.
No producer imports, root approximations, or disabled asserts are used.
"""
from fractions import Fraction as F
from math import comb, factorial
from itertools import product
from pathlib import Path
import hashlib
import json
import sys
sys.stdout.reconfigure(newline='\n')
N=0


def need(ok,label):
    global N
    if not ok:raise RuntimeError(label)
    N+=1


# Sparse polynomials in x,y,z,t, with Laurent exponents allowed in t.
def clean(p):return {k:F(c) for k,c in p.items() if c}
def add(p,q):
    out=p.copy()
    for k,c in q.items():out[k]=out.get(k,F(0))+c
    return clean(out)
def scale(p,c):return clean({k:c*v for k,v in p.items()})
def mul(p,q):
    out={}
    for a,c in p.items():
        for b,d in q.items():
            k=tuple(x+y for x,y in zip(a,b));out[k]=out.get(k,F(0))+c*d
    return clean(out)
def power(p,n):
    out={(0,0,0,0):F(1)}
    for _ in range(n):out=mul(out,p)
    return out
def mono(c=1,x=0,y=0,z=0,t=0):return {(x,y,z,t):F(c)}
def evaluate(p,x,y,z,t):return sum(c*x**a*y**b*z**d*t**e for (a,b,d,e),c in p.items())
def specialize(p,x,y,z):
    out={}
    for (a,b,d,e),c in p.items():out[e]=out.get(e,F(0))+c*x**a*y**b*z**d
    return {k:v for k,v in out.items() if v}


beta={};cc={};dd={}
for j,c in enumerate([1,13,55]):beta=add(beta,mono(c,t=j-1))
for key in [mono(x=1,t=2),mono(y=1,t=3),mono(z=1,t=4)]:beta=add(beta,key)
for j,c in enumerate([1,12,45]):cc=add(cc,mono(c,t=j-1))
cc=add(add(cc,mono(F(2,3),x=1,t=2)),mono(F(3,7),y=1,t=3))
for j,c in enumerate([1,11,36]):dd=add(dd,mono(c,t=j-1))
dd=add(add(dd,mono(F(5,12),x=1,t=2)),mono(F(1,7),y=1,t=3))
O={j:F(comb(14,2*j+1)) for j in range(7)}
E={j:F(comb(14,2*j)) for j in range(8)}
carrier={}
for a in O:
    for b in O:carrier[a+b]=carrier.get(a+b,F(0))+O[a]*O[b]
for a in E:
    for b in E:carrier[a+b-1]=carrier.get(a+b-1,F(0))+E[a]*E[b]
raw=add(mul(beta,beta),scale(mul(mul(mono(t=1),cc),dd),2))
Q=clean({k:c*carrier.get(k[3],0) for k,c in raw.items()})
P=clean({k:c*O.get(k[3],0) for k,c in beta.items()})
need(set(k[3] for k in Q)==set(range(-1,9)),'complete Laurent support')
need(Q[(0,0,0,-1)]==28 and all(c>0 for c in Q.values()),'all raw coefficients and lower carry')
expectedP=add(add(add(mono(182),mono(20020,t=1)),mono(2002,x=1,t=2)),add(mono(3432,y=1,t=3),mono(2002,z=1,t=4)))
need(P==expectedP,'full original first polynomial')

# In the next polynomials t denotes positive s, with the same original root.
T=clean({(a,b,d,j+1):c*(-1 if j%2 else 1) for (a,b,d,j),c in Q.items()})
z_at_root=add(add(mono(F(12,7),y=1,t=-1),mono(-1,x=1,t=-2)),add(mono(10,t=-3),mono(F(-1,11),t=-4)))
eliminated={}
for (a,b,d,e),c in T.items():eliminated=add(eliminated,mul(mono(c,x=a,y=b,t=e),power(z_at_root,d)))
h=scale(eliminated,F(-11,14))
need(all(d==0 and 0<=e<=8 for a,b,d,e in h),'same-zero elimination is a degree-eight polynomial')
fixed=[-443993,73031400,-3278853435,46232902140,-234760993560,343030019640,-83518139925,-26087589000,3585421125]
need(specialize(h,F(84),F(35),F(0))==dict(enumerate(map(F,fixed))),'incoming fibre identity independently recovered')


def normalized_phase(s,x,y,z):return z*s**4-F(12,7)*y*s**3+x*s*s-10*s+F(1,11)


endpoints=[F(1,102),F(1,100),F(1,9),F(13,100),F(1),F(8,5),F(19,2)]
signs=[1,-1,-1,1,1,-1,-1]
endpoint_record=[]
for s,sgn in zip(endpoints,signs):
    vals=[normalized_phase(s,F(x),F(y),F(z)) for x,y,z in product((83,85),(34,36),(0,5))]
    for v in vals:need(sgn*v>0,'uniform original-phase endpoint sign')
    endpoint_record.append({'s':str(s),'min':str(min(vals)),'max':str(max(vals)),'sign':sgn})

# A direct signed termwise bound handles the first phase, without eliminating z.
first_upper=F(0)
for (a,b,d,e),c in T.items():
    x,y,z,s=(F(85),F(36),F(5),F(1,100)) if c>0 else (F(83),F(34),F(0),F(1,102))
    first_upper+=c*x**a*y**b*z**d*s**e
need(first_upper==F(-728550005046322718853208807,1704830652000000000000000),'exact raw first-branch bound')
need(first_upper < -427,'strict first response margin')


def phase_transform(a,b):
    # Return nine bivariate coefficient polynomials of the transformed h.
    out=[{} for _ in range(9)]
    for (i,j,d,e),c in h.items():
        need(d==0,'no lost coefficient survives elimination')
        if b is None:
            terms={k:F(comb(e,k))*a**(e-k) for k in range(e+1)}
        else:
            terms={}
            for k in range(e+1):
                for l in range(8-e+1):
                    terms[k+l]=terms.get(k+l,F(0))+comb(e,k)*a**(e-k)*b**k*comb(8-e,l)
        for k,c0 in terms.items():out[k][(i,j)]=out[k].get((i,j),F(0))+c*c0
    return out


def rectangle_bernstein(p):
    # x=83+2X,y=34+2Y; tensor Bernstein degree(2,2).
    translated={}
    for (a,b),c in p.items():
        for i in range(a+1):
            for j in range(b+1):
                translated[i,j]=translated.get((i,j),F(0))+c*comb(a,i)*83**(a-i)*2**i*comb(b,j)*34**(b-j)*2**j
    return [[sum(translated.get((i,j),F(0))*F(comb(k,i),comb(2,i))*F(comb(l,j),comb(2,j))
                 for i in range(k+1) for j in range(l+1)) for l in range(3)] for k in range(3)]


windows=[(F(1,9),F(13,100)),(F(1),F(8,5)),(F(19,2),None)]
cert=[]
for a,b in windows:
    transformed=phase_transform(a,b)
    blocks=[rectangle_bernstein(p) for p in transformed]
    allvals=[v for block in blocks for row in block for v in row]
    for v in allvals:need(v>0,'strict tensor Bernstein coefficient')
    need(len(allvals)==81,'complete coefficient block')
    # Independent evaluation of the identity at three parameter values.
    for X,Y,u in [(F(0),F(1),F(2)),(F(1,3),F(2,5),F(3,7)),(F(1),F(0),F(0))]:
        xx,yy=83+2*X,34+2*Y
        coeff=[]
        for block in blocks:
            coeff.append(sum(block[i][j]*comb(2,i)*X**i*(1-X)**(2-i)*comb(2,j)*Y**j*(1-Y)**(2-j) for i in range(3) for j in range(3)))
        result=sum(v*u**k for k,v in enumerate(coeff))
        s=a+u if b is None else (a+b*u)/(1+u)
        direct=evaluate(h,xx,yy,F(0),s)*(1 if b is None else (1+u)**8)
        need(result==direct,'independent transform evaluation')
    cert.append({'left':str(a),'right':None if b is None else str(b),
                 'minimum':str(min(allvals)),'coefficients':[[list(map(str,row)) for row in block] for block in blocks]})

# Nonnegative-root domain cap, including weak/discriminant boundaries.
v=F(8,25)
derivative_upper=5*v**4-52*v**3+165*v*v-166*v+36
need(derivative_upper==F(-146524,78125)<0,'uniform first-critical-point cutoff')
need(83-55*v==F(327,5),'quadratic cap coefficient')
need(F(36*36*5,4*327)==F(540,109)<5,'uniform z cap')

# Every x,y in the rectangle has a positive-root B at z=1.
b_nodes=[F(0),F(1,10),F(1),F(3),F(5),F(7)]
for w,sign in zip(b_nodes,[-1,1,-1,1,-1,1]):
    for x,y in product((83,85),(34,36)):
        need(sign*(w**5-13*w**4+55*w**3-x*w*w+y*w-1)>0,'nonvacuity over the entire coefficient rectangle')

# Literal factorial fibres verify the actual source at its named base point.
actual_counts=[]
for mass,model in [(14,P),(28,Q)]:
    literal={}
    for a in range(mass+1):
        for c in range(mass-a+1):
            b=mass-a-c
            if -27*a+b+15*c==0:
                j=a-mass//14
                literal[j]=literal.get(j,0)+factorial(mass)//(factorial(a)*factorial(b)*factorial(c))
    need(literal==specialize(model,F(84),F(35),F(1)),'literal -27,1,15 factorial fibre')
    actual_counts.append(len(literal))
need(actual_counts==[5,10],'all fifteen actual channels retained')


def moments(x,y,z,isD=True):
    den=[F(1),F(-13),F(55),-x,y,-z]
    num=[F(1),F(-11),F(36),-5*x/12,y/7] if isD else [F(1),F(-12),F(45),-2*x/3,3*y/7]
    out=[]
    for j in range(9):out.append((num[j] if j<5 else F(0))-sum(den[k]*out[j-k] for k in range(1,min(j,5)+1)))
    return out


def determinant(a):
    a=[row[:] for row in a];out=F(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None:return F(0)
        if j!=i:a[i],a[j]=a[j],a[i];out=-out
        q=a[i][i];out*=q
        for j in range(i+1,len(a)):
            c=a[j][i]/q
            for k in range(i+1,len(a)):a[j][k]-=c*a[i][k]
    return out


# The prior positive degree-seven surrogate is rejected by the exact domain.
hx,hy,hs=F(78071,1000),F(601,50),F(57,2)
hz=12*hy/(7*hs)-hx/hs**2+10/hs**3-1/(11*hs**4)
need(evaluate(P,hx,hy,hz,-hs)==0 and evaluate(Q,hx,hy,hz,-hs)>0,'inherited original positive surrogate')
hostile_dets=[]
for isD in [False,True]:
    mm=moments(hx,hy,hz,isD);d=determinant([[mm[i+j] for j in range(5)] for i in range(5)])
    need(d<0,'exact domain rejects the degree-seven hostile')
    hostile_dets.append(str(d))

# Dropping the fifth-coefficient cap is false inside this x,y rectangle.
lo,hi=F(16693,2000),F(41733,5000)
need(normalized_phase(lo,F(84),F(35),F(6))<0<normalized_phase(hi,F(84),F(35),F(6)),'z6 original hostile root bracket')
fixed_transform=phase_transform(lo,hi)
for coefficient in fixed_transform:
    need(sum(c*84**i*35**j for (i,j),c in coefficient.items())<0,'negative h on the entire z6 hostile bracket')
mm=moments(F(84),F(35),F(6))
need(determinant([[mm[i+j] for j in range(5)] for i in range(5)])==-25668,'z6 is outside exact D domain')

record={'rectangle':{'x':['83','85'],'y':['34','36'],'z':['0','5']},
        'Q_raw':[{'powers':list(k),'coefficient':str(c)} for k,c in sorted(Q.items())],
        'eliminated_h':[{'powers':list(k),'coefficient':str(c)} for k,c in sorted(h.items())],
        'phase_endpoints':endpoint_record,'raw_first_upper':str(first_upper),
        'positive_transforms':cert,'critical_derivative_upper':str(derivative_upper),
        'beta_geometry_z_cap':'540/109','degree7_hostile_H5_determinants':hostile_dets}
here=Path(__file__).resolve().parent
if here.name=='04-computation':here=here.parent/'05-knowledge'/'results'
path=here/(Path(__file__).stem+'_certificate.json')
path.write_bytes((json.dumps(record,indent=2,sort_keys=True)+'\n').encode())
print('Uniform coefficient prism: 83<=x<=85,34<=y<=36,0<=z<=5; every original positive phase has full Q<0.')
print('Nonnegative beta geometry gives z<540/109<5, with no interlacer hypothesis.')
print('Original phase endpoint extrema:',json.dumps(endpoint_record,sort_keys=True))
print('Raw first-branch upper:',first_upper)
for c in cert:print('Positive transform',c['left'],c['right'],'81 coefficients; minimum',c['minimum'])
print('All15 literal actual channels PASS; exact degree8 domain rejects both positive-response hostiles.')
print('Certificate SHA256',hashlib.sha256(path.read_bytes()).hexdigest())
print('PASS',N,'always-active exact gates; actual LF stdout')
