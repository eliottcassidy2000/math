#!/usr/bin/env python3
"""Exact, dependency-free checks. No assertion is used as a validation gate."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

counts = {}
def check(ok, label):
    if not ok:
        raise RuntimeError(label)
    counts[label] = counts.get(label, 0) + 1

def add(x,y): return (x[0]+y[0],x[1]+y[1])
def mul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
def scale(k,x): return (k*x[0],k*x[1])
def Q(x): return x[0]**2-x[0]*x[1]+x[1]**2
ONE=(1,0); OMEGA=(0,1); POW=[ONE,OMEGA,(-1,-1)]
def rows(x):
    a,b=x; return (2*a-b,-a-b,-a+2*b)
def R(v): return (v[1],v[2],v[0])
def K(v):
    r=R(v);r2=R(r);return tuple(a-b for a,b in zip(r,r2))
def primitive(s):return s>0 and s%2==1 and s%8!=5
def edge(s):
    if not primitive(s):raise ValueError(s)
    a=2 if s%8==1 else 1; v=(3*s+1)//2**a;j=0
    while v%8==5:
        v=(v-1)//4;j+=1
    return v,a,a+2*j

def source(N,s):
    z=(0,0);v=s
    while v<4*N:
        z=add(z,mul(POW[v%3], (1,2) if v<2*N else OMEGA))
        v=4*v+1
    return z

def source_rows(N,s):
    v=s;r=[0,0,0]
    while v<4*N:
        k=(0,-3,3) if v<2*N else (-1,-1,2)
        for phase in range(3):r[phase]+=k[(v+phase)%3]
        v=4*v+1
    return tuple(r)

# Exact arithmetic lattice and all three marked observations.
for a in range(-12,13):
 for b in range(-12,13):
    z=(a,b);v=rows(z)
    check(sum(v)==0 and len({x%3 for x in v})==1,'image_lattice')
    check(sum(x*x for x in v)==6*Q(z),'three_row_norm')
    check(v[0]**2+3*b*b==4*Q(z),'marked_row_identity')
    check(R(v)==rows(mul(OMEGA,z)),'cyclic_intertwiner')
    check(K(v)==rows(mul((1,2),z)),'skew_intertwiner')
    check(K(K(v))==tuple(-3*x for x in v),'skew_square')
    check(Q(mul((1,2),z))==3*Q(z),'skew_norm_multiplier')
for u in range(-12,13):
 for v in range(-12,13):
    w=-u-v
    if (u-v)%3==0:
        z=((u-v)//3,(-u-2*v)//3)
        check(rows(z)==(u,v,w),'lattice_inverse')
check((1-(-1))%3!=0,'full_A2_is_not_image')
check(rows((1,0))[0]!=rows(mul(OMEGA,(1,0)))[0],'marked_observer_rotation_hostile')

for N in range(1,129):
 for s in range(1,4*N,2):
    if not primitive(s):continue
    z=source(N,s)
    check(rows(z)==source_rows(N,s),'independent_original_source_rows')
    check(Q(z) in (0,1,3),'native_source_norm_shells')

# All-scale quotient obstruction sampled independently of its algebraic proof.
cutoff_examples=[]
for N in range(2,8193):
    s=4*N-(5 if N%4==1 else 1);y,a,t=edge(s)
    check(2*N<=s<4*N and primitive(s),'cutoff_source_domain')
    check(Q(source(N,s))==1,'cutoff_source_positive_mass')
    check(y>4*N and a==t==1,'cutoff_target_escapes')
    check(source(N,y)==(0,0),'cutoff_target_zero_mass')
    if N in (2,3,5,4096,4097,8192):cutoff_examples.append({'N':N,'s':s,'P(s)':y})

inverse_examples=[]
for k in range(16):
    s=(32*64**k-5)//9;y,a,t=edge(s)
    b=F(3*(3*s+1),3*(3*s+1)+2**a)
    check((y,a,t)==(1,1,3+6*k),'infinite_ROOT_inverse_family')
    check(b==1-F(1,16*64**k),'inverse_gain_exact')
    h=-1/b
    check(F(-1)*b*h==1,'inverse_boundary_mode')
    check(h*h>=1,'inverse_counting_norm_diverges')
    if k<4:inverse_examples.append({'k':k,'s':s,'clock':t,'b':str(b),'h':str(h)})

# Authenticate finite first-hit words; no global hitting claim is used.
path_examples=[]
max_depth=0
q_sorted=sorted([F(18*s+6+(4 if s%8==1 else 2),4 if s%8==1 else 2) for s in range(3,10001,2) if primitive(s)])
for k,q in enumerate(q_sorted[:1000],1):
    check(q==18*k+10+3*(-1)**(k+1),'reciprocal_loss_rank')
for s0 in range(1,20001,2):
    if not primitive(s0):continue
    s=s0;visited=set();B=F(1);eps=1;clock=0;h=0;carry=0
    while s!=1:
        if s in visited or h>=10000:raise RuntimeError(('unauthenticated finite universe',s0,s,h))
        visited.add(s);y,a,t=edge(s);X=3*s+1;Y=3*y+1
        check(2**t*Y==3*X+2**a,'native_affine_carry')
        b=F(3*X,3*X+2**a);d=1-b*b
        check(F(0)<b<F(1) and d>0,'positive_branch')
        carry=3*carry+2**(clock+a)
        B*=b;eps*=(-1)**t;clock+=t;h+=1;s=y
    max_depth=max(max_depth,h)
    X0=3*s0+1
    check(2**clock*4==3**h*X0+carry,'word_carry')
    check(B==F(3**h*X0,4*2**clock),'product_telescope')
    mode=eps/B
    check(mode==eps*F(4*2**clock,3**h*X0),'root_conditioned_mode')
    check((1/B)**9<=1+3*h,'polynomial_depth_inflation')
    # Source charge=ONE: exact homogeneous terminal, inflated initial account.
    e=1/B;spent=F(0);s=s0
    while s!=1:
        y,a,t=edge(s);b=F(3*(3*s+1),3*(3*s+1)+2**a);d=1-b*b;c=(-1)**t*b
        r=d*e;spent+=r*r/d;e=c*e;s=y
    check(e==eps and spent==1/(B*B)-1,'zero_forcing_boundary_invoice')
    if s0 in (1,3,7,11,17,113,227,4095,19999):
        path_examples.append({'s':s0,'depth':h,'clock':clock,'B':str(B),'mode':str(mode),'inflation':str(1/(B*B))})

# ROOT self-edge is excluded, not silently made homogeneous.
y,a,t=edge(1);b=F(3*(3*1+1),3*(3*1+1)+2**a)
check((y,a,t,b)==(1,2,2,F(3,4)),'ROOT_loop_obstruction')

# Polarized C3 transfer: every pair interaction scales by 3, no debt gain.
Z=[(a,b) for a in range(-4,5) for b in range(-4,5)]
for z in Z:
 for w in Z:
    beta=Q(add(z,w))-Q(z)-Q(w)
    check(sum(a*b for a,b in zip(rows(z),rows(w)))==3*beta,'pair_row_polarization')
    kz=mul((1,2),z);kw=mul((1,2),w)
    check(Q(add(kz,kw))-Q(kz)-Q(kw)==3*beta,'pair_debt_scaling')


# Authenticate the inherited frontier and its exact internal boundary factor.
frontier_path=Path(__file__).parent/'sources'/'frontier.txt'
frontier={int(l.split()[0]):int(l.split()[1]) for l in frontier_path.read_text().splitlines() if not l.startswith('#')}
check(len(frontier)==627 and max(frontier)==63665 and max(frontier.values())==64,'frontier_size_height_rank')
check(min(s for s in range(1,1010,2) if primitive(s) and s not in frontier)==1009,'frontier_minimum_excluded')
internal={1:F(1)}
for s in sorted(frontier,key=frontier.get):
    y,a,t=edge(s)
    check(y in frontier and (s==1 or frontier[y]<frontier[s]),'frontier_closure_descent')
    if s!=1:internal[s]=F(3*(3*s+1),3*(3*s+1)+2**a)*internal[y]
worst=min(internal,key=internal.get);C=1/internal[worst]**2
check(worst==993 and C==F(332306998946228968225951765070086144,211753373595316860199691410086465225),'frontier_exact_internal_factor')
for s0 in range(1,20001,2):
    if not primitive(s0):continue
    s=s0;B=F(1);h=0;lam=F(0);seen=set()
    while s not in frontier:
        if s in seen:raise RuntimeError(('repeated exterior source',s0,s))
        seen.add(s);y,a,t=edge(s);b=F(3*(3*s+1),3*(3*s+1)+2**a)
        B*=b;lam+=(1-b)/(1+b);h+=1;s=y
    check((1/B)**9<=1+F(3*h,754),'exterior_polynomial_inflation')
    qbound=sum((F(1,4543+18*k) for k in range(h)),F(0))
    check(lam<=qbound,'exterior_ranked_forcing')
    check((1/(B*internal[s])**2/C)**9<=(1+F(3*h,754))**2,'combined_boundary_invoice_bound')
# Authenticated actual-source marked-observer hostile from the dyadic worker.
z=source(512,799);w=source(512,1065)
check(z==(-2,-1) and w==(0,1),'actual_marked_hostile_charges')
check(799 in frontier and 1065 not in frontier and edge(1065)==(799,2,2),'actual_marked_hostile_shared_terminal')
check(Q(add(z,w))-Q(z)-Q(w)==0 and 2*rows(z)[0]*rows(w)[0]==6,'actual_marked_hostile_cross_terms')
frontier_result={'file_sha256':sha256(frontier_path.read_bytes()).hexdigest(),'maximum_internal_inflation':str(C),'maximizer':worst,'maximizer_rank':frontier[worst]}

receipt={'frontier':frontier_result,'status':'PASS','counts':counts,'finite_first_hit_universe':'all positive primitive seeds s <= 20000; each authenticated to first ROOT with cap 10000','maximum_native_depth':max_depth,'cutoff_hostiles':cutoff_examples,'ROOT_inverse_family':inverse_examples,'path_controls':path_examples,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'claims_not_made':['all-source ROOT arrival','a fixed all-source native Hilbert measure','uniform spatial ROOT energy','dyadic source-debt payment']}
print(json.dumps(receipt,indent=2,sort_keys=True))
