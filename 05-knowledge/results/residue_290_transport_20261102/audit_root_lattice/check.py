#!/usr/bin/env python3
"""Independent mathematical controls: no producer or repo imports."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import product,combinations
from math import comb,gcd
import hashlib,json

G=0
def req(ok,why):
    global G
    G+=1
    if not ok:raise RuntimeError(why)
def v2(n):
    n=abs(n)
    if not n:raise ValueError('v2 zero')
    count=0
    while n%2==0:n//=2;count+=1
    return count

def add(a,b):
    r=[0]*max(len(a),len(b))
    for i,v in enumerate(a):r[i]+=v
    for i,v in enumerate(b):r[i]+=v
    return r
def mul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
      for j,y in enumerate(b):r[i+j]+=x*y
    return r
def neg(a):return [-x for x in a]
def clean(a):
    while a and a[-1]==0:a.pop()
    return a

# Independently clear the generation denominator symbolically.
D=[1,1,1];e=[0,1,-1];An=[0,1,0,1];Bn=[0,1,1];kn=[1,1]
B2n=clean(add(mul(mul(e,Bn),D),mul(add(An,neg(mul(e,D))),kn)))
req(B2n==[0,0,1,2,2,0,-1],'cleared B2 numerator')
slack=clean(add(mul([0,0,1],mul(D,D)),neg(B2n)))
req(slack==[0,0,0,0,1,2,2] and all(x>=0 for x in slack),'B2 upper polynomial')
req(F(4,3)**3>2 and F(25,4)**3>3**5 and 1+F(4,3)+6*F(25,4)==F(239,6),'M13 constants')

# U mode is built by backward edge propagation, not the closed formula.
cache={1:((),(),F(1))}
def urecord(n):
    cur=n; path=[];vals=[];seen=set()
    while cur not in cache:
        if cur in seen or len(path)>10000:raise RuntimeError('U sample not authenticated')
        seen.add(cur);path.append(cur);a=v2(3*cur+1);vals.append(a);cur=(3*cur+1)//2**a
    st,aa,hh=cache[cur]
    for x,a in reversed(list(zip(path,vals))):
        hh=hh/F((-1)**a*3*x,3*x+1)
        st=(x,)+st;aa=(a,)+aa;cache[x]=(st,aa,hh)
    return cache[n]
def w_from_counters(L,K):
    return F(2*(L+1),(L+K+1)*(L+K+2)*comb(L+K,L))

mass=F(0);norm=F(1);passive_cases=0
for s in range(3,2048,2):
    st,aa,h=urecord(s);tau=len(st);S=sum(aa);K=sum((a-1)//2 for a in aa)
    W=w_from_counters(tau-1,K)
    req(h==F((-1)**S*2**S,3**tau*s),'closed h versus propagated h')
    req(abs(h)**6<=tau+1,'distinct odd vertex bound')
    mass+=W;norm+=W*h*h
    # Verify every prefix's energy with the original squared phase mode.
    loss=F(0)
    for i,x in enumerate(st):
        hx=urecord(x)[2];b=F(3*x,3*x+1);d=1-b*b
        loss+=d*hx*hx
        y=1 if i==tau-1 else st[i+1]
        req(urecord(y)[2]**2+loss==h*h,'native passive prefix')
        passive_cases+=1
req(mass<F(13,3) and norm<55,'partial W controls')
for L in range(8):
 for K in range(1,12):
    req(w_from_counters(L,K)==w_from_counters(L+1,K)+w_from_counters(L,K+1),'original W Pascal')
    req(sum(w_from_counters(L+1,K+j) for j in range(9))+w_from_counters(L,K+9)==w_from_counters(L,K),'W tail payment')
ray=[]
for k in range(13):
    z=(64**(k+1)-1)//3;st,aa,h=urecord(z)
    W=w_from_counters(0,3*k+2);d=1-F(3*z,3*z+1)**2
    req(st==(z,) and aa==(6*k+6,) and z%3==0,'actual omitted source')
    req(W==F(2,(3*k+3)*(3*k+4)),'root-ray weight')
    ray.append(W/d)
req(all(a<b for a,b in zip(ray,ray[1:])) and ray[-1]>10**16,'frontier growth control')

# P guard/carry obtained directly by repeated inverse B stripping; compare
# every pair in each short word class, not only consecutive source pairs.
def p(s):
    if s<=0 or s%2==0 or s%8==5:raise ValueError('P source guard')
    a=2 if s%8==1 else 1; z=F(3*s+1,2**a)
    req(z.denominator==1,'initial integer P quotient')
    z=int(z);j=0
    while z%8==5:z=(z-1)//4;j+=1
    req(z>0 and z%2 and z%8!=5,'P target guard')
    return z,(a,j),a+2*j
pairs=0;word_counts=[]
for depth in (1,2,3):
    groups=defaultdict(list)
    for s in range(1,1024,2):
        if s%8==5:continue
        y=s;w=[];T=0
        for _ in range(depth):
            previous=y;y,letter,t=p(y);w.append(letter);T+=t
            req(2**t*(3*y+1)-3*(3*previous+1)==2**letter[0],'P carry')
        groups[tuple(w)].append((s,y,T))
    word_counts.append(len(groups))
    for group in groups.values():
      for (s,y,T),(z,w,T2) in combinations(group,2):
        req(T==T2 and (y-w)*2**T==3**depth*(s-z),'same-word affine difference')
        req(v2((y-w)//2)==v2((s-z)//2)-T>=0,'pair owner exact clock')
        pairs+=1
req(p(17)[0]==p(35)[0]==3 and (p(17)[2],p(35)[2])==(4,5),'cross-word collision hostile')

# Eisenstein multiplication via integral 2x2 multiplication matrices.
def emat(x):
    m,n=x;return ((m,-n),(n,m-n))
def mv(A,x):return tuple(sum(A[i][j]*x[j] for j in range(2)) for i in range(2))
def mm(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def det(A):return A[0][0]*A[1][1]-A[0][1]*A[1][0]
def Q(x):return det(emat(x))
def epow(x,p,mod):
    A=emat(x);ans=((1,0),(0,1))
    while p:
        if p&1:ans=tuple(tuple(v%mod for v in row) for row in mm(ans,A))
        A=tuple(tuple(v%mod for v in row) for row in mm(A,A));p//=2
    return mv(ans,(1,0))

split=[]
for k in (1,2,3):
    q=7**k;roots=[r for r in range(q) if (r*r+r+1)%q==0]
    al=next(r for r in roots if r%7==2);be=next(r for r in roots if r%7==4)
    req(len(roots)==2 and (al+be+1)%q==0 and (al*be-1)%q==0,'split roots')
    # Literal span of pi^k modulo q, compared to the evaluation kernel.
    M=((1,0),(0,1))
    for _ in range(k):M=mm(M,emat((2,-1)))
    req(det(M)==q,'ideal index')
    image={tuple(v%q for v in mv(M,(a,b))) for a in range(q) for b in range(q)}
    kernel={(m,n) for m in range(q) for n in range(q) if (m+al*n)%q==0}
    req(image==kernel and len(image)==q,'ideal kernel via image enumeration')
    observed={};primitive_colors=Counter()
    for u,v in product(range(q),repeat=2):
        # Inverse pair observation is the primary loop, independently of source enumeration.
        n=((v-u)*pow(be-al,-1,q))%q;m=(u-al*n)%q
        req((m+be*n)%q==v and Q((m,n))%q==u*v%q,'pair inverse and norm')
        du=(m+(al+2)*n)%q;dv=(m+(be+2)*n)%q
        req((du==v)==(k==1 or n%7**(k-1)==0),'exact exchange locus')
        req((du-u)%q==2*n%q,'successive observation inverse')
        observed[(u,du)]=(m,n)
        req(Q((m+2*n,n))-Q((m,n))==2*n*(2*m+n),'real norm shear defect')
        if m%7 or n%7:primitive_colors[u]+=1
    req(len(observed)==q*q,'two-sample observer bijective')
    req(sum(primitive_colors.values())==48*7**(2*k-2),'primitive state cardinality')
    req(all(primitive_colors[u]==(q if u%7 else 6*q//7) for u in range(q)),'actual residue fibre weights')
    split.append({'q':q,'ideal_image_size':len(image),'observer_size':len(observed),'primitive_states':sum(primitive_colors.values())})
# Exhaustively realize each primitive residue at q49 by genuine positive Euclid columns.
lifts=0
q=49
for m0,n0 in product(range(q),repeat=2):
    if m0%7==0 and n0%7==0:continue
    n=n0 if n0 else q
    found=next((m0+q*t for t in range(1,5000) if m0+q*t>n and gcd(m0+q*t,n)==1 and (m0+q*t+n)%2),None)
    req(found is not None,'primitive lift missing')
    req(((found+2*n)%q,n%q)==((m0+2*n0)%q,n0),'source D guard lift')
    lifts+=1

# Full branch matrices retain the determinant that endpoint residues forget.
Um=((2,-1),(1,0)); Am=((2,1),(1,0)); Dm=((1,2),(0,1))
for q0 in (7,49,343):
    transport=Am
    for _ in range(q0-1):transport=mm(Dm,transport)
    transport=mm(Dm,mm(Um,transport))
    endpoint=mv(transport,(2,1))
    req(endpoint==(16*q0+2,4*q0+1) and det(transport)==-1,'full word and determinant')
    req(tuple(z%q0 for z in endpoint)==(2,1),'empty word residue hostile')
req(det(Dm)==1 and Dm[0][0]+Dm[1][1]==2,'same spectral polynomial as identity')
for q0 in (7,49):
    rr=[z for z in range(q0) if (z*z+z+1)%q0==0]
    al=next(z for z in rr if z%7==2);be=next(z for z in rr if z%7==4)
    for m,n in product(range(q0),repeat=2):
        u=(m+al*n)%q0;v=(m+be*n)%q0
        rot=(-n,m-n);con=(m-n,-n)
        req(((rot[0]+al*rot[1])%q0,(rot[0]+be*rot[1])%q0)==(al*u%q0,be*v%q0),'rotation ports')
        req(((con[0]+al*con[1])%q0,(con[0]+be*con[1])%q0)==(v,u),'conjugate ports')
        req(Q(rot)==Q((m,n))==Q(con),'norm symmetries')

# Inert geometry is checked by full local multiplication-table inverses and
# direct norm distributions, then a separate full CRT target count.
inert=[]
for p0 in (2,5,29):
    histogram=Counter();nonzero=0
    for m,n in product(range(p0),repeat=2):
        x=(m,n);N=Q(x)%p0;histogram[N]+=1
        conjugate=((m-n)%p0,-n%p0)
        req(epow(x,p0,p0)==conjugate,'Frobenius conjugation')
        req(tuple(t%p0 for t in mv(emat(x),conjugate))==(N,0),'field norm')
        if x!=(0,0):
            req(N!=0,'no zero divisor norm')
            inv=tuple(t*pow(N,-1,p0)%p0 for t in conjugate)
            req(tuple(t%p0 for t in mv(emat(x),inv))==(1,0),'local explicit inverse')
            nonzero+=1
    req(histogram[0]==1 and all(histogram[r]==p0+1 for r in range(1,p0)),'local norm circle')
    inert.append({'p':p0,'units':nonzero,'norm_fibre':p0+1})
full=Counter(Q((m,n))%290 for m,n in product(range(290),repeat=2))
unit_bins={r:c for r,c in full.items() if gcd(r,290)==1}
req(len(unit_bins)==112 and set(unit_bins.values())=={540} and sum(unit_bins.values())==60480,'full mod290 product norm')
req(not any((x*x+x+1)%290==0 for x in range(290)),'no rational evaluation at290')
req(Q((6,2))==28,'Eisenstein28')
# Positive definite bound: Q(a,b)>=(a*a+b*b)/2 makes this window complete.
for target in (145,290):
    bound=25
    req(2*target<bound*bound,'complete norm window radius')
    req(not any(Q((a,b))==target for a,b in product(range(-bound,bound+1),repeat=2)),'no integral inert norm')
def gaussian(x,y):a,b=x;c,d=y;return a*c-b*d,a*d+b*c
req(gaussian((2,1),(5,-2))==(12,1) and 12**2+1==145,'Gaussian145')
req(gaussian((1,1),(12,1))==(11,13) and 11**2+13**2==290,'Gaussian290')
req(gaussian((1,1),gaussian((2,1),(5,2)))==(-1,17) and 1+17**2==290,'second Gaussian290')
req(F(112,580)==F(28,145) and 540!=112,'density and norm-fibre distinction')
# Gaussian 145 port inversion and the exact nonzero nilradical at 290.
req([r for r in range(145) if (r*r+1)%145==0]==[12,17,128,133],'four Gaussian roots')
gaussian_ports=set()
for u,v in product(range(145),repeat=2):
    n=(u-v)*pow(24,-1,145)%145;m=(u-12*n)%145
    req((m-12*n)%145==v and (m*m+n*n)%145==u*v%145,'Gaussian145 inverse norm')
    gaussian_ports.add((m,n))
req(len(gaussian_ports)==145**2,'Gaussian145 ring observer')
gaussian_bins=Counter();nilpotents=[]
for m,n in product(range(290),repeat=2):
    gnorm=(m*m+n*n)%290
    if gcd(gnorm,290)==1:gaussian_bins[gnorm]+=1
    square=gaussian((m,n),(m,n))
    if all(z%290==0 for z in square):nilpotents.append((m,n))
req(nilpotents==[(0,0),(145,145)],'exact Gaussian290 nilpotents')
req([z for z in range(290) if z*z%290==0]==[0],'scalar290 reduced')
# A nonfaithful observer still supports the shear: this is the hostile to
# accidentally widening a failed ring-isomorphism claim to every observer.
observed_D={};port_kernel=[]
for m,n in product(range(290),repeat=2):
    uv=((m+17*n)%290,(m-17*n)%290)
    nxt=((m+19*n)%290,(m-15*n)%290)
    req(uv not in observed_D or observed_D[uv]==nxt,'quotient D well-defined')
    observed_D[uv]=nxt
    if uv==(0,0):port_kernel.append((m,n))
req(port_kernel==[(0,0),(145,145)],'paired17 kernel is nilradical')
req(len(observed_D)==290**2//2 and all((u-v)%2==0 for u,v in observed_D),'paired17 parity image')
req(len(gaussian_bins)==112 and set(gaussian_bins.values())=={224} and sum(gaussian_bins.values())==25088,'Gaussian290 norm unit fibres')

print(json.dumps({'status':'PASS','gates':G,'U_initial_sources':1023,'U_passive_prefixes':passive_cases,
 'U_partial_mass':str(mass),'U_partial_norm':str(norm),'root_ray_checks':len(ray),
 'P_source_limit_exclusive':1024,'P_depths':[1,2,3],'P_word_classes':word_counts,'P_all_common_word_pairs':pairs,
 'split_quotients':split,'primitive_lifts_mod49':lifts,'inert_fields':inert,
 'norm290_unit_fibres':len(unit_bins),'norm290_fibre_size':540,'norm290_units':60480,
 'Gaussian290_unit_fibres':len(gaussian_bins),'Gaussian290_fibre_size':224,'Gaussian290_units':25088,
 'Gaussian290_nilpotents':nilpotents,
 'scope':'independent finite controls and source-sensitive hostiles; no producer imports'},indent=2,sort_keys=True))
