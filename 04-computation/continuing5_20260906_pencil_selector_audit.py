"""Independent rational pencil referee; no SymPy or producer imports.

Five degree-bounded evaluations reconstruct each determinant quartic twice:
from the literal series Hankel and from a 9x9 Sylvester matrix. Sturm chains
and separate leading minors replay the root-order and nonnesting controls.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import json,sys,hashlib
sys.stdout.reconfigure(newline='\n')
G=0
def need(ok,label):
    global G
    if not ok:raise RuntimeError(label)
    G+=1
def trim(a):
    a=list(map(F,a))
    while len(a)>1 and not a[-1]:a.pop()
    return a
def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b):out[i+j]+=c*d
    return trim(out)
def val(a,t):
    out=F(0)
    for c in reversed(a):out=out*t+c
    return out
def rem(a,b):
    a,b=trim(a),trim(b)
    while a!=[0] and len(a)>=len(b):
        a=add(a,[F(0)]*(len(a)-len(b))+scale(b,-a[-1]/b[-1]))
    return a
def sturm(a):
    a=trim(a);out=[a,[i*a[i] for i in range(1,len(a))]]
    while True:
        r=scale(rem(out[-2],out[-1]),-1)
        if r==[0]:break
        out.append(scale(r,1/abs(r[-1])))
    return out
def changes(seq,t):
    signs=[]
    for p in seq:
        a=p[-1] if t=='inf' else (p[-1]*(-1)**(len(p)-1) if t=='-inf' else val(p,t))
        if a:signs.append(a>0)
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def roots(p,a,b):
    seq=sturm(p)
    return changes(seq,a)-changes(seq,b)
def det(a):
    a=[list(map(F,row)) for row in a];d=F(1)
    for i in range(len(a)):
        k=next((k for k in range(i,len(a)) if a[k][i]),None)
        if k is None:return F(0)
        if i!=k:a[i],a[k]=a[k],a[i];d=-d
        q=a[i][i];d*=q
        for k in range(i+1,len(a)):
            r=a[k][i]/q
            for j in range(i+1,len(a)):a[k][j]-=r*a[i][j]
    return d
def leading(a):return [det([row[:n] for row in a[:n]]) for n in range(1,len(a)+1)]
def multiply(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def transpose(a):return [list(x) for x in zip(*a)]

def polynomials(x,y,z,A):
    B=[-z,y,-x,F(55),F(-13),F(1)]
    a=([3*y/7,-2*x/3,F(45),F(-12),F(1)] if A=='C' else
       [y/7,-5*x/12,F(36),F(-11),F(1)])
    return B,a
def moments(x,y,z,A):
    B,a=polynomials(x,y,z,A);den=list(reversed(B));num=list(reversed(a));m=[]
    for j in range(9):m.append((num[j] if j<5 else F(0))-sum(den[k]*m[j-k] for k in range(1,min(j,5)+1)))
    return m
def hankel(x,y,z,A):
    m=moments(x,y,z,A)
    return [[m[i+j] for j in range(5)] for i in range(5)]
def sylvester(B,A):
    b,a=list(reversed(B)),list(reversed(A))
    return [[F(0)]*i+b+[F(0)]*(3-i) for i in range(4)]+[[F(0)]*i+a+[F(0)]*(4-i) for i in range(5)]
def interpolate(vals):
    out=[F(0)]
    for i,c in enumerate(vals):
        p=[F(1)];den=F(1)
        for j in range(len(vals)):
            if i!=j:p=mul(p,[-j,1]);den*=i-j
        out=add(out,scale(p,c/den))
    return out

here=Path(__file__).resolve().parent
certificate=here/'continuing5_20260906_pencil_selector_certificate.json'
if not certificate.exists():certificate=here.parent/'05-knowledge'/'results'/certificate.name
raw=certificate.read_bytes();data=json.loads(raw)
need(hashlib.sha256(raw).hexdigest()=='4f4ea7bdcb823f05b6d74b818597a8e4a4a007cc87f8bcf0fc5bfa9033761b8d','frozen complete producer certificate')
need(len(data['fixtures'])==3 and len(data['channel_hostiles'])==2,'declared fixture universe')

# An explicit exact hyperbolic congruence supplements the all-parameter proof.
# S=[[0,A],[A^T,D]], A=[[0,1],[1,a]]. With u=u'-A^{-T}Dv/2,
# the transformed form is [[0,A],[A^T,0]], of inertia(2,2).
for a,b,c in [(F(0),F(0),F(0)),(F(14),F(130),F(904)),(F(15),F(147),F(1067)),(F(3,2),F(-7,3),F(11,5))]:
    S=[[0,0,0,1],[0,0,1,a],[0,1,a,b],[1,a,b,c]]
    need(det(S)==1,'slope determinant and rank')
    ainv=[[-a,F(1)],[F(1),F(0)]];D=[[a,b],[b,c]]
    L=multiply(ainv,D)
    U=[[F(int(i==j)) for j in range(4)] for i in range(4)]
    for i in range(2):
        for j in range(2):U[i][j+2]=-L[i][j]/2
    transformed=multiply(transpose(U),multiply(S,U))
    target=[[0,0,0,1],[0,0,1,a],[0,1,0,0],[1,a,0,0]]
    need(transformed==target,'exact hyperbolic congruence')

# Degree-bounded affine slope controls, using direct rational series division.
for x,y in [(F(0),F(0)),(F(1),F(1)),(F(2),F(3)),(F(155,2),F(9)),(F(311,4),F(21,2))]:
    for A in ('C','D'):
        m0,m1,m2=(moments(x,y,z,A) for z in [F(0),F(1),F(2)])
        slopes=([F(0)]*5+[F(1),F(14),F(130),904+4*x/3] if A=='C' else
                [F(0)]*5+[F(1),F(15),F(147),1067+19*x/12])
        need([b-a for a,b in zip(m0,m1)]==slopes,'native full slope')
        need([c-2*b+a for a,b,c in zip(m0,m1,m2)]==[0]*9,'affine slope control')

summary=[]
for fixture in data['fixtures']:
    need(set(fixture['channels'])=={'C','D'},'both channels retained')
    x,y,z0=map(F,[fixture['x'],fixture['y'],fixture['strict_z']]);boxes={}
    for A,rec in fixture['channels'].items():
        hvals=[];rvals=[]
        for z in map(F,range(5)):
            hvals.append(det(hankel(x,y,z,A)))
            B,a=polynomials(x,y,z,A);rvals.append(det(sylvester(B,a)))
            need(hvals[-1]==rvals[-1],'independent literal resultant evaluation')
        hpoly=interpolate(hvals);rpoly=interpolate(rvals)
        expected=list(reversed(list(map(F,rec['quartic_descending']))))
        need(hpoly==rpoly==expected,'complete degree-bounded polynomial reconstruction')
        need(len(hpoly)==5 and hpoly[-1]==1,'monic quartic')
        ab=[list(map(F,pair)) for pair in rec['isolating_intervals']]
        need(roots(hpoly,'-inf','inf')==4,'all four real roots')
        for lo,hi in ab:need(lo<hi and roots(hpoly,lo,hi)==1,'independent rational Sturm box')
        need(all(ab[i][1]<ab[i+1][0] for i in range(3)) and ab[1][1]<z0<ab[2][0],'middle-root ordering at strict reference')
        ls=leading(hankel(x,y,z0,A))
        need(ls==list(map(F,rec['reference_leading_minors'])),'all reference minors')
        for d in ls:need(d>0,'strict reference PD')
        B,a=polynomials(x,y,z0,A)
        need(roots(B,0,'inf')==5 and roots(a,0,'inf')==4,'native positive beta and interlacer root counts')
        for at in [ab[0][0]-1,ab[3][1]+1]:
            need(val(hpoly,at)>0 and not all(d>0 for d in leading(hankel(x,y,at,A))),'exterior positive determinant is not PSD')
        if 'z0_leading_minors' in rec:need(leading(hankel(x,y,F(0),A))==list(map(F,rec['z0_leading_minors'])),'zero-endpoint matrix data')
        boxes[A]=ab
    name=fixture['name']
    if name=='C_upper':need(boxes['C'][2][1]<boxes['D'][2][0],'C upper active')
    elif name=='D_upper':need(boxes['D'][2][1]<boxes['C'][2][0],'D upper active')
    else:
        need(boxes['D'][1]==[F(163,2101),F(9,116)] and boxes['D'][1][0]>0,'strict positive lower endpoint')
        need(boxes['D'][2]==[F(177,530),F(176,527)] and boxes['D'][2][1]<boxes['C'][2][0],'exact full upper interval')
        need(boxes['C'][1][1]<0 and det(hankel(x,y,F(0),'D'))<0,'zero genuinely outside full fibre')
    summary.append((name,{A:ab[1:3] for A,ab in boxes.items()}))

for control in data['channel_hostiles']:
    A=control['survives'];other='D' if A=='C' else 'C'
    x,y,z=map(F,[control['x'],control['y'],control['z']])
    good,bad=leading(hankel(x,y,z,A)),leading(hankel(x,y,z,other))
    need(good==list(map(F,control['surviving_leading_minors'])) and bad==list(map(F,control['failed_leading_minors'])),'nonnesting exact minors')
    need(all(d>0 for d in good) and bad[-1]<0,'nonnesting strict one-channel model')
    B,a=polynomials(x,y,z,A)
    need(roots(B,0,'inf')==5 and roots(a,0,'inf')==4,'nonnesting actual root geometry')

for A in ('C','D'):
    need(all(d>0 for d in leading(hankel(F(84),F(35),F(0),A))),'zero is a permitted strict matrix reference')
    B,a=polynomials(F(84),F(35),F(0),A)
    need(B[0]==0 and B[1]>0 and roots(B[1:],0,'inf')==4,'zero reference has one simple beta zero and four positive roots')

print('PASS exact selector proof audit: strict reference, inertia(2,2,1zero), middle quartic roots, and weak endpoints.')
for name,bs in summary:print(name,{A:[[str(v) for v in pair] for pair in ab] for A,ab in bs.items()})
print('All six quartics independently reconstructed from five Hankel and five 9x9 Sylvester determinants each.')
print('Full two-channel nonnesting and a positive lower endpoint retained; no global reference supplier inferred.')
print('Producer certificate SHA256',hashlib.sha256(raw).hexdigest())
print('PASS',G,'always-active exact gates; no producer imports or SymPy')
