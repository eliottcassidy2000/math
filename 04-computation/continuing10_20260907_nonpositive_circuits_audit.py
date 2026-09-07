"""Independent native determinant and inverse-Bernstein referee."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations
from math import comb
from pathlib import Path
import json
import sys
import sympy as s
sys.stdout.reconfigure(newline="\n")
HERE=Path(__file__).resolve().parent
G=0
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def trim(f):
    f=list(map(Q,f))
    while len(f)>1 and not f[-1]: f.pop()
    return f
def ev(f,x):
    out=Q(0)
    for a in reversed(f): out=out*x+a
    return out
def deriv(f): return trim([i*f[i] for i in range(1,len(f))] or [0])
def rem(f,g):
    f=trim(f); g=trim(g)
    while f!=[0] and len(f)>=len(g):
        n=len(f)-len(g); factor=f[-1]/g[-1]
        for j,a in enumerate(g): f[n+j]-=factor*a
        f=trim(f)
    return f
def sturm(f):
    a=trim(f); b=deriv(a); out=[a,b]
    while b!=[0]:
        c=[-q for q in rem(a,b)]
        if c==[0]: break
        c=[q/abs(c[-1]) for q in c]
        out.append(c); a,b=b,c
    return out
def variation(values):
    signs=[1 if v>0 else -1 for v in values if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def var_at(chain,x):
    return variation([ev(f,x) for f in chain])
def count(chain,a,b):
    va=var_at(chain,a)
    vb=variation([f[-1] for f in chain]) if b is None else var_at(chain,b)
    return va-vb
def reflected(f): return [(-1)**k*a for k,a in enumerate(f)]
def center(d,q): return [Q(comb(d,k))*q**(k*(k-1)//2) for k in range(d+1)]
def ratios(e):
    d=len(e)-1; h=[e[k]/comb(d,k) for k in range(d+1)]
    R=[h[k]*h[k]/(h[k-1]*h[k+1]) for k in range(1,d)]
    C=[R[k]/R[k-1] for k in range(1,len(R))]
    return R,C
def reconstruct(d,q,cs):
    h=[Q(1),Q(1)]; rk=1/q
    for k in range(1,d):
        if k>=2: rk*=cs[k-2]
        h.append(h[-1]*h[-1]/(rk*h[-2]))
    return [comb(d,k)*h[k] for k in range(d+1)]
def refine(f,chain,interval):
    lo,hi=interval
    if ev(f,hi)==0: return hi,hi
    mid=(lo+hi)/2
    return (lo,mid) if count(chain,lo,mid) else (mid,hi)
def isolate_positive(f):
    chain=sturm(f); degree=len(f)-1
    gate(count(chain,Q(0),None)==degree,'fresh exact positive root count')
    upper=Q(1)
    while count(chain,Q(0),upper)<degree: upper*=2
    pending=[(Q(0),upper)]; intervals=[]
    while pending:
        lo,hi=pending.pop(); n=count(chain,lo,hi)
        if n==0: continue
        if n==1:
            for _ in range(34): lo,hi=refine(f,chain,(lo,hi))
            intervals.append((lo,hi)); continue
        mid=(lo+hi)/2; pending.extend(((mid,hi),(lo,mid)))
    return sorted(intervals)

def model(x,y,z):
    return ([-z,y,-x,Q(55),Q(-13),Q(1)],
            [Q(3,7)*y,-Q(2,3)*x,Q(45),Q(-12),Q(1)],
            [Q(1,7)*y,-Q(5,12)*x,Q(36),Q(-11),Q(1)])

def imul(a,b):
    values=[x*y for x in a for y in b]
    return min(values),max(values)

def iev(f,interval):
    result=(Q(0),Q(0))
    for c in reversed(f):
        lo,hi=imul(result,interval)
        result=(lo+c,hi+c)
    return result

def sign_interval(f,interval):
    lo,hi=iev(f,interval)
    return 1 if lo>0 else -1 if hi<0 else 0

def verify_C_geometry(x,y,z):
    B,C,D=model(x,y,z)
    br=isolate_positive(B); cr=isolate_positive(C)
    gate(len(br)==5 and len(cr)==4,'five B/four C positive roots')
    gate(all(br[i][1]<cr[i][0] and cr[i][1]<br[i+1][0] for i in range(4)),
         'strict C interlacing by independently isolated polynomial roots')
    gate(all(sign_interval(C,br[i])==(-1)**i for i in range(5)),
         'positive C/B residues at every native root')
    return B,C,D,br,cr

def locate(name):
    for p in (HERE/name,HERE.parent/'05-knowledge/results'/name,
              Path('C:/w/continuing10_20260907_moments')/name):
        if p.is_file(): return p
    raise FileNotFoundError(name)
STEM='continuing10_20260907_nonpositive_circuits'
for suffix,pin in (('.py','d392980046171db2ea1df6cb759c04882fd57305ae811fe3d7f5746df6b6a973'),
                   ('.out','36e5c21f2113d4a3dfbc825af03ddb0d3570a3e359c58ac5a3a78f164db9908b'),
                   ('_certificate.json','be4762cf527f1add2887ee2b29a98440a83f0c68454571b54d182b7395b58092')):
    raw=locate(STEM+suffix).read_bytes()
    gate(sha256(raw).hexdigest()==pin,'frozen primary '+suffix)
    if suffix=='.out': gate(b'\r' not in raw,'actual producer LF')
data=json.loads(locate(STEM+'_certificate.json').read_text(encoding='utf-8'))
x,y,z,u,v,w=s.symbols('x y z u v w')
def expr(value): return s.sympify(value,locals={str(t):t for t in (x,y,z,u,v,w)})
def zero(value,label): gate(s.cancel(value)==0,label)
def moments(x,y,z,which):
    den=[1,-13,55,-x,y,-z]
    num=[1,-12,45,-Q(2,3)*x,Q(3,7)*y] if which=='C' else [1,-11,36,-Q(5,12)*x,Q(1,7)*y]
    out=[]
    for k in range(9):
        out.append(s.expand((num[k] if k<5 else 0)-sum(den[j]*out[k-j] for j in range(1,min(k,5)+1))))
    return out
def determinant(M):
    # Literal permutation expansion, independent of a symbolic matrix engine.
    n=len(M); answer=0
    for p in permutations(range(n)):
        term=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        for i in range(n): term*=M[i][p[i]]
        answer+=term
    return s.expand(answer)
dm=moments(x,y,z,'D')
for a,b in zip(dm[:7],data['D_moments']): zero(a-expr(b),'native D series moment')
gate(s.Poly(dm[6],y).coeff_monomial(y)==-s.Rational(915,7),'corrected prose coefficient minus915over7')
H3=determinant([[dm[i+j] for j in range(3)] for i in range(3)])
H4=determinant([[dm[i+j] for j in range(4)] for i in range(4)])
zero(H3-expr(data['D_H3']),'native third determinant')
zero(H4-expr(data['D_H4']),'native fourth determinant')
zero(s.diff(H4,z,2)+6,'full unbounded-z concavity')
x0=s.Rational(831875,8788); yl=13*x**3/s.Integer(166375)
yh=(-343*x*x+67788*x-3157056)/s.Integer(2592)
zl=44*y**3/x**3
for a,key in ((x0,'x0'),(yl,'y_lower'),(yh,'y_upper'),(zl,'z_lower')):
    zero(a-expr(data[key]),'native inequality chart '+key)
gate(x0>0 and x0<99,'positive x interval and denominator')
zero(H3+s.Rational(18,7)*(y-yh),'moving third determinant boundary')
gg=s.expand(343*x*x-67788*x+3157056+2592*yl)
zero(gg-expr(data['g']),'feasibility univariate bound')
zero(gg.subs(x,99)-expr(data['g_at99']),'strict99 endpoint obstruction')
shift=s.Poly(s.expand(s.diff(gg,x).subs(x,99+u)),u)
zero(shift.as_expr()-expr(data['g_derivative_shift']),'entire shifted derivative')
gate(gg.subs(x,99)>0 and all(c>0 for c in shift.all_coeffs()),'all-height exclusion x>=99')
X=x0+(99-x0)*u; Y=yl.subs(x,X)+(yh-yl).subs(x,X)*v
polys=[s.cancel(-7112448*x**6*H4.subs(z,zl)),s.cancel(-1008*x**3*s.diff(H4,z).subs(z,zl))]
gate([Q(data['scale0']),Q(data['scale1'])]==[7112448,1008],'positive normalization factors')
for i,P in enumerate(polys):
    zero(P-expr(data['P'+str(i)]),'complete primitive sign polynomial')
    gate(s.Poly(P,x,y).domain==s.ZZ and s.Poly(P,x,y).content()==1,'actual primitive integer normalization')
zero(H4.subs(z,zl+w)-H4.subs(z,zl)-w*s.diff(H4,z).subs(z,zl)+3*w*w,
     'exact all-z tail identity')

# Invert the complete Bernstein expansion to power coefficients. This does
# not compute Bernstein coefficients by the producer's triangular formula.
counts=[]
for index,degree in enumerate(((18,6),(9,3))):
    rec=data['chart'+str(index)]; d,e=degree
    bdict={tuple(k):Q(value) for k,value in rec['bernstein']}
    gate(len(bdict)==(d+1)*(e+1) and set(bdict)=={(r,k) for r in range(d+1) for k in range(e+1)},
         'entire tensor index universe')
    gate(all(value>0 for value in bdict.values()),'every tensor coefficient strictly positive')
    expanded={}
    for (r,k),value in bdict.items():
        for i in range(r,d+1):
            first=comb(d,r)*comb(d-r,i-r)*(-1)**(i-r)
            for j in range(k,e+1):
                term=value*first*comb(e,k)*comb(e-k,j-k)*(-1)**(j-k)
                expanded[i,j]=expanded.get((i,j),Q(0))+term
    expanded={k:value for k,value in expanded.items() if value}
    direct=s.Poly(s.expand(polys[index].subs({x:X,y:Y},simultaneous=True)),u,v)
    gate(direct.degree(u)==d and direct.degree(v)==e,'actual chart bidegree')
    expected={tuple(k):Q(value) for k,value in direct.terms()}
    gate(expanded==expected,'full inverse Bernstein identity in native determinant chart')
    gate(expected=={tuple(k):Q(value) for k,value in rec['power_terms']},'complete saved monomial packet')
    gate(tuple(rec['degrees'])==degree,'saved tensor degree bounds')
    minimum=min(bdict.values())
    claimed=(Q(38259412467502741144816725,262144),Q(4049307478755,256))[index]
    gate(minimum==claimed,'exact positive global tensor minimum')
    gate(minimum==Q(rec['minimum']),'saved exact minimum')
    # Also compare the stored monomial packet when its shape is decoded below.
    gate(all(sum(Q(comb(d,r))*U**r*(1-U)**(d-r) for r in range(d+1))==1
             for U in (Q(0),Q(1,3),Q(1))), 'closed-interval partition-of-unity controls')
    counts.append(len(bdict))

def actual_control(rec):
    xyz=list(map(Q,rec['xyz'])); es=[Q(1),Q(13),Q(55)]+xyz
    rr,cc=ratios(es)
    gate(rr==list(map(Q,rec['R'])) and cc==list(map(Q,rec['circuits'])),'literal native Newton/circuit ratios')
    minors={}
    for name in ('C','D'):
        mm=moments(*xyz,name)
        ds=[determinant([[mm[i+j] for j in range(k)] for i in range(k)]) for k in range(1,6)]
        gate(ds==list(map(Q,rec[name+'_minors'])),'all actual '+name+' leading minors')
        minors[name]=ds
    return xyz,rr,cc,minors
bad3,rr,cc,M=actual_control(data['H4_only_hostile'])
gate(all(r>1 for r in rr) and all(c<1 for c in cc) and M['D'][2]<0<M['D'][3],
     'fourth-alone false implication hostile')
good3,rr,cc,M=actual_control(data['H3_only_hostile'])
gate(all(c<1 for c in cc) and all(d>0 for d in M['C']) and M['D'][2]>0>M['D'][3],
     'third-alone false implication hostile')
verify_C_geometry(*good3)
center,rr,cc,M=actual_control(data['tie_center'])
gate(cc==[1,1,1] and M['D'][2]>0>M['D'][3],'all-tie boundary control')
gate({tuple(map(Q,rec['xyz'])) for rec in data['two_negative_survivors']}=={(Q(95),Q(68),Q(1)),(Q(86),Q(50),Q(9))},
     'complete two-negative survivor universe')
for rec in data['two_negative_survivors']:
    xyz,rr,cc,M=actual_control(rec)
    gate(all(r>1 for r in rr) and sum(c<1 for c in cc)==2,'exact two-negative strict Newton survivor')
    gate(all(d>0 for d in M['D'][:4]) and M['D'][4]<0,'through-six D positivity misses native degree-eight rejection')
    if xyz==[86,50,9]:
        gate(all(d>0 for d in M['C']),'actual C geometry for plus-minus-minus')
        B,C,D,br,cr=verify_C_geometry(*xyz)
        gate(sign_interval(D,br[0])!=0,'native D root evaluation retained')
full_words=set()
for rec in data['positive_controls']:
    xyz,rr,cc,M=actual_control(rec)
    gate(all(d>0 for d in M['C']+M['D']),'actual full-model positive control')
    B,C,D,br,cr=verify_C_geometry(*xyz)
    dr=isolate_positive(D)
    gate(all(br[i][1]<dr[i][0] and dr[i][1]<br[i+1][0] for i in range(4)),
         'independent full-model D root ordering')
    full_words.add(tuple(1 if c>1 else -1 if c<1 else 0 for c in cc))
gate(full_words=={(1,1,1),(-1,1,1),(1,-1,1),(1,1,-1)},'complete four positive-control word universe')
print('PASS independent closed-nonpositive-octant audit')
print('Native D moments, literal determinants, positive denominators, full curved chart and all-z tail independently reconstructed.')
print('Complete inverse Bernstein arrays:',counts,'all coefficients strictly positive including the boundary.')
print('Exact controls retain both one-determinant failures, all ties, two-negative through-six survivors and four full-model points.')
print('Only the closed three-nonpositive octant is excluded; general two-negative model exclusion remains OPEN.')
print('Always-active exact gates',G)
