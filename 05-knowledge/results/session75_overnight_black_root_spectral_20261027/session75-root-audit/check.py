#!/usr/bin/env python3
"""Independent Session 75 native-ROOT audit. Stdlib only; no producer imports.

The source compiler traverses literal odd inputs and assigns their primitive
ancestor, independently of the producer's seed/B-chain implementation.
Native edges use a 2-adic valuation formula and are checked by literal peeling.
All validation gates are explicit and survive python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

F = Fraction
checks = Counter()

def require(label, condition, context=None):
    checks[label] += 1
    if not condition:
        raise RuntimeError(f'{label}: {context!r}')

def plus(z, w):
    return tuple(a+b for a, b in zip(z, w))

def times(z, w):
    a, b = z
    c, d = w
    return (a*c-b*d, a*d+b*c-b*d)

def scale(k, z):
    return tuple(k*a for a in z)

def norm(z):
    a, b = z
    return a*a-a*b+b*b

def pair(z, w):
    return norm(plus(z,w))-norm(z)-norm(w)

def row(z):
    a, b = z
    return (2*a-b, -a-b, -a+2*b)

def rotate(v):
    return (v[1], v[2], v[0])

def skew(v):
    x,y,z = v
    return (y-z,z-x,x-y)

ZERO = (0,0)
ONE = (1,0)
OM = (0,1)
OMEGA = (ONE,OM,(-1,-1))
ALPHA = (1,2)

def is_seed(s):
    return s > 0 and s % 8 in (1,3,7)

def peel(m):
    j = 0
    while (m-1) % 8 == 4:
        m = (m-1)//4
        j += 1
    return m, j

def native(s):
    if not is_seed(s):
        raise ValueError(('nonprimitive source',s))
    a = 2 if s % 8 == 1 else 1
    odd = (3*s+1)//(2**a)
    n = 3*odd+1
    v = (n & -n).bit_length()-1
    j = (v-1)//2
    target = (n//4**j-1)//3
    require('native_valuation_vs_peeling', (target,j)==peel(odd), s)
    require('native_affine_identity', 2**(a+2*j)*(3*target+1)==3*(3*s+1)+2**a,s)
    require('native_target_primitive',is_seed(target),s)
    return target,a,a+2*j

def gain(s,a):
    return F(9*s+3,9*s+3+2**a)

def literal_source(N):
    out = defaultdict(lambda: ZERO)
    original_rows = defaultdict(lambda: (0,0,0))
    for m in range(1,4*N,2):
        s,_ = peel(m)
        increment = times(OMEGA[m%3], ALPHA if m<2*N else OM)
        out[s] = plus(out[s],increment)
        k = (0,-3,3) if m<2*N else (-1,-1,2)
        original_rows[s] = plus(original_rows[s],tuple(k[(m+r)%3] for r in range(3)))
    return dict(out),dict(original_rows)

def chain_source(N,s):
    total = ZERO
    for _ in range(10000):
        if s >= 4*N:
            return total
        total = plus(total,times(OMEGA[s%3],ALPHA if s<2*N else OM))
        s = 4*s+1
    raise RuntimeError('source-chain cap')

# Source formula and the marked first row, including literal input ownership.
source_hash = sha256()
for N in range(1,257):
    charges,rows = literal_source(N)
    require('literal_primitive_support',set(charges)=={s for s in range(1,4*N,2) if is_seed(s)},N)
    for s,A in charges.items():
        require('literal_vs_chain_source',A==chain_source(N,s),(N,s))
        require('all_original_source_rows',row(A)==rows[s],(N,s))
        require('source_shells',norm(A) in (0,1,3),(N,s))
        source_hash.update(f'{N},{s},{A[0]},{A[1]};'.encode())

# Cutoff obstruction: all residue classes, exact edge, positive vs zero mass.
examples=[]
Ns=list(range(2,20001))
Ns += [4*k+r for k in (2**60,3**100,10**150) for r in range(4)]
for N in Ns:
    s = 4*N-(5 if N%4==1 else 1)
    y,a,t = native(s)
    A = chain_source(N,s)
    require('cutoff_source_in_band',2*N<=s<4*N and is_seed(s),N)
    require('cutoff_single_chain_term',4*s+1>=4*N,N)
    require('cutoff_exact_charge',A==OMEGA[(s+1)%3] and norm(A)==1,N)
    require('cutoff_j_zero',a==t==1,N)
    require('cutoff_exact_target',y==(6*N-7 if N%4==1 else 6*N-1),N)
    require('cutoff_target_outside',y>4*N and chain_source(N,y)==ZERO,N)
    # delta_y and zero are the same input quotient class. At s their
    # composition and (composition-c*identity) values differ by one.
    require('cutoff_distinct_null_representatives',s!=y and norm(A)*1==1,N)
    if N in (2,3,4,5,6,9,20000):
        examples.append({'N':N,'s':s,'P(s)':y,'A':A,'mu_s':norm(A),'mu_y':0})
charges,_=literal_source(1)
support={s for s,A in charges.items() if norm(A)>0}
require('N1_exact_support',support=={1,3},support)
require('N1_forward_closed',all(native(s)[0] in support for s in support),support)
# The finite quotient criterion can be checked directly, not just inferred.
null_examples=[]
for N in range(2,65):
    charges,_=literal_source(N)
    support={s for s,A in charges.items() if norm(A)>0}
    s=4*N-(5 if N%4==1 else 1)
    y=native(s)[0]
    norm_C_delta=sum(norm(charges[x]) for x in support if native(x)[0]==y)
    require('explicit_quotient_failure',y not in support and norm_C_delta>=1,N)
    if N<=5:
        null_examples.append({'N':N,'null_coordinate':y,'input_norm_sq':0,'pullback_norm_sq':norm_C_delta})

# ROOT's actual loop and forcing in both conventions.
y,a,t=native(1)
b=gain(1,a)
d=1-b*b
require('ROOT_self_edge',(y,a,t,b)==(1,2,2,F(3,4)),(y,a,t,b))
require('ROOT_boundary_forcing_price',(1-b)**2/d==F(1,7))

# Fixed path versus root-conditioned mode, using backward recursion for modes.
max_depth=max_clock=0
path_count=0
path_examples=[]
path_profiles={}
for initial in range(1,24001,2):
    if not is_seed(initial):
        continue
    path=[]
    seen=set()
    s=initial
    while s!=1:
        if s in seen or len(path)>=10000:
            raise RuntimeError(('finite first-hit universe failed',initial,s))
        seen.add(s)
        y,a,t=native(s)
        b=gain(s,a)
        require('strict_edge_contraction',0<b<1,s)
        path.append((s,y,a,t,b))
        s=y
    depth=len(path)
    clock=sum(edge[3] for edge in path)
    eps=(-1)**clock
    product=F(1)
    for _,_,_,_,b in path:
        product*=b
    require('path_product_telescope',product==F(3**depth*(3*initial+1),4*2**clock),initial)
    require('global_first_hit_polynomial_bound',1/product**9<=1+3*depth,initial)
    path_profiles[initial]=(product,depth)
    mode=F(1)
    for s,y,a,t,b in reversed(path):
        old=mode
        mode/=(-1)**t*b
        require('backward_homogeneous_mode',(-1)**t*b*mode==old,s)
    require('first_hit_boundary_mode',mode==eps/product,initial)
    # Test actual source charges for an independent varying native cutoff.
    N=initial//4+2
    A=chain_source(N,initial)
    mu=norm(A)
    e_original=A
    e_inflated=scale(1/product,A)
    dissipation=F(0)
    for s,y,a,t,b in path:
        c=(-1)**t*b
        d=1-b*b
        r=scale(d,e_inflated)
        dissipation+=norm(r)/d
        e_inflated=scale(c,e_inflated)
        e_original=scale(c,e_original)
    require('unchanged_initial_attenuation',e_original==scale(eps*product,A),initial)
    require('inflated_initial_terminal',e_inflated==scale(eps,A),initial)
    require('initial_energy_inflation',norm(scale(1/product,A))==mu/product**2,initial)
    require('zero_signed_forcing_invoice',dissipation==mu*(1/product**2-1),initial)
    max_depth=max(max_depth,depth)
    max_clock=max(max_clock,clock)
    path_count+=1
    if initial in (1,3,7,11,113,227,4095,23999):
        path_examples.append({'s':initial,'depth':depth,'clock':clock,'B':str(product),'h_ROOT_1':str(mode),'mu':mu})

# Direct inverse family excludes counting-l2, with exact branch authentication.
inverse_examples=[]
for k in range(32):
    numerator=32*64**k-5
    require('direct_inverse_integrality',numerator%9==0,k)
    s=numerator//9
    y,a,t=native(s)
    b=gain(s,a)
    require('direct_ROOT_inverse',(y,a,t)==(1,1,3+6*k),k)
    require('direct_ROOT_gain',b==1-F(1,16*64**k),k)
    require('counting_l2_nondecay',1/b**2>1,k)
    if k<4:
        inverse_examples.append({'k':k,'s':s,'clock':t,'gain':str(b)})

# Full two-port polarization and literal positive-b versus signed-c controls.
forcing_examples=[]
for s in (1,3,7,9,11,15,17,113,227):
    _,a,t=native(s)
    b=gain(s,a)
    d=1-b*b
    eps=(-1)**t
    c=eps*b
    lam=(1-b)/(1+b)
    for e in ((1,0),(0,1),(1,2),(F(1,2),F(-2,3))):
        for forcing in ((0,0),(1,-2),(F(-1,3),F(2,5))):
            e_next=plus(scale(c,e),forcing)
            residual=plus(scale(d,e),scale(-c,forcing))
            require('signed_two_port_identity',norm(e_next)+norm(residual)/d==norm(e)+norm(forcing)/d,(s,e,forcing))
        physical=scale(eps,e)
        f_positive=plus(physical,scale(-b,e))
        f_signed=plus(physical,scale(-c,e))
        require('original_positive_b_invoice',norm(f_positive)/d==(lam if eps==1 else 1/lam)*norm(e),(s,e))
        require('signed_c_invoice',norm(f_signed)/d==lam*norm(e),(s,e))
    forcing_examples.append({'s':s,'epsilon':eps,'b':str(b),'positive_b_unit_price':str(lam if eps==1 else 1/lam),'signed_c_unit_price':str(lam)})
require('s3_conventions_differ',(-gain(3,1))!=gain(3,1))

# Lattice, exact K image, orientation, and retained marked observer.
for a in range(-20,21):
    for b in range(-20,21):
        z=(a,b)
        v=row(z)
        kz=times(ALPHA,z)
        require('E_centered_congruence',sum(v)==0 and len({c%3 for c in v})==1,z)
        require('E_inverse',((v[0]-v[1])//3,(-v[0]-2*v[1])//3)==z,z)
        require('E_norm',sum(c*c for c in v)==6*norm(z),z)
        require('marked_observer',v[0]**2+3*b*b==4*norm(z),z)
        require('cyclic_orientation',row(times(OM,z))==rotate(v),z)
        require('K_intertwines_alpha',row(kz)==skew(v),z)
        require('K_square',skew(skew(v))==scale(-3,v),z)
        require('alpha_norm_multiplier',norm(kz)==3*norm(z),z)
        kv=skew(v)
        require('K_image_divisible_by_3',all(x%3==0 for x in kv),z)
        require('conjugated_marked_observer',kv[2]-kv[1]==3*v[0],z)
# Every point 3(u,-u+v,-v) comes from a unique input in L.
for u in range(-20,21):
    for v in range(-20,21):
        x=(3*u,3*(-u+v),-3*v)
        z=(x[2]//3,-x[0]//3)
        require('K_image_surjects_3A2',skew(row(z))==x,(u,v))
        a2=(u,-u+v,-v)
        image=(a2[0]-a2[1])%3==0
        if image:
            require('L_surjects_congruence_plane',row(((a2[0]-a2[1])//3,(-a2[0]-2*a2[1])//3))==a2,a2)
require('lattice_index_3',abs(2*(-2)-(-1)*1)==3)
require('centered_not_in_L',(1-(-1))%3!=0)
require('K_not_L_automorphism',row(ONE)[0]%3!=0)
require('marked_rotation_hostile',row(ONE)[0]==2 and row(OM)[0]==-1)
require('alpha_square',times(ALPHA,ALPHA)==(-3,0))
Z=[(a,b) for a in range(-4,5) for b in range(-4,5)]
for z in Z:
    for w in Z:
        require('bilinear_lift',sum(a*b for a,b in zip(row(z),row(w)))==3*pair(z,w),(z,w))
        require('bilinear_skew_scale',pair(times(ALPHA,z),times(ALPHA,w))==3*pair(z,w),(z,w))
require('marked_cross_hostile',pair(OM,(-2,-1))==0 and 2*row(OM)[0]*row((-2,-1))[0]==6)

# Actual authenticated terminal records, dyadic pair ownership and debt scaling.
# This does NOT claim a uniform debt bound. All sources at these finite N arrive.
ledger_examples=[]
for N in (1,2,3,4,5,17,64,129,256):
    charges,_=literal_source(N)
    records=[]
    for s,A in sorted(charges.items()):
        y=s
        T=0
        seen=set()
        while y!=1:
            if y in seen or len(seen)>10000:
                raise RuntimeError(('ledger finite arrival failed',N,s,y))
            seen.add(y)
            y,_,t=native(y)
            T+=t
        records.append((s,(-1)**T,A))
    def levels(transform):
        es=[]
        q=1
        while q<=4*N:
            bins=defaultdict(lambda:ZERO)
            for s,sign,A in records:
                key=((s-1)//2)%q
                bins[key]=plus(bins[key],scale(sign,transform(A)))
            es.append(sum(norm(A) for A in bins.values()))
            q*=2
        return es
    E=levels(lambda x:x)
    KE=levels(lambda x:times(ALPHA,x))
    Delta=[x-y for x,y in zip(E,E[1:])]
    KDelta=[x-y for x,y in zip(KE,KE[1:])]
    S=sum(norm(A) for _,_,A in records)
    D=sum(max(x,0) for x in Delta)
    KD=sum(max(x,0) for x in KDelta)
    require('source_is_final_refinement',E[-1]==S,N)
    require('actual_energy_K_scaling',KE==[3*x for x in E],N)
    require('signed_level_debt_K_scaling',KDelta==[3*x for x in Delta],N)
    require('positive_after_level_debt_K_scaling',KD==3*D,N)
    require('normalized_debt_unchanged',F(KD,3*S)==F(D,S),N)
    require('energy_debt_telescope',E[0]==S+sum(Delta),N)
    # Pair owner is label-only; transforming charges cannot change it.
    for i,(s,sign,A) in enumerate(records):
        for other,other_sign,B in records[:i]:
            gap=abs((s-other)//2)
            q=gap & -gap
            require('unique_pair_owner',((s-1)//2)%q==((other-1)//2)%q and ((s-1)//2)%(2*q)!=((other-1)//2)%(2*q),(N,s,other))
    ledger_examples.append({'N':N,'S':S,'D_positive':D,'D_over_S':str(F(D,S)),'levels':len(E)})

# Phase 1: exact ranked-loss bounds and independently traversed frontier.
def reciprocal_loss(s):
    return F(9*s+5,2) if s%8==1 else F(9*s+4)

def ranked_q(k):
    return 18*k+10+3*(-1)**(k+1)

ordered=sorted(reciprocal_loss(s) for s in range(3,100001,2) if is_seed(s))
for k in range(1,10001):
    require('ranked_q_formula',ordered[k-1]==ranked_q(k),k)
    require('ranked_q_simple_lower_bound',ranked_q(k)-1>=18*k+6,k)
    require('exterior_tail_lower_bound',ranked_q(251+k)>=4543+18*(k-1),k)
rank_product=F(1)
for h in range(1,513):
    q=ranked_q(h)
    rank_product*=F(q+1,q-1)
    require('ranked_global_polynomial',rank_product**9<=1+3*h,h)
require('root_excluded_from_ranked_q',reciprocal_loss(1)==7 and ordered[0]==31)

base=Path(__file__).resolve().parent
frontier_lines=[l for l in (base/'sources'/'frontier.txt').read_text().splitlines() if l.strip() and not l.startswith('#')]
frontier_pairs=[tuple(map(int,l.split())) for l in frontier_lines]
V=dict(frontier_pairs)
require('frontier_no_duplicate_labels',len(V)==len(frontier_pairs))
require('frontier_size_maximum',len(V)==627 and max(V)==63665 and max(V.values())==64)
require('frontier_root_rank',V[1]==0)
require('frontier_minimum_exterior',min(s for s in range(1,1010,2) if is_seed(s) and s not in V)==1009)
require('frontier_minimum_q_exterior',reciprocal_loss(1009)==4543)
internal={}
for s0,rank in sorted(V.items()):
    require('frontier_primitive_labels',is_seed(s0),s0)
    s=s0
    depth=0
    product=F(1)
    visited=set()
    while s!=1:
        if s in visited:
            raise RuntimeError(('frontier repeat',s0,s))
        visited.add(s)
        y,a,t=native(s)
        require('frontier_closed',y in V,s)
        require('frontier_rank_exact_step',V[y]+1==V[s],s)
        product*=gain(s,a)
        s=y
        depth+=1
    require('frontier_rank_matches_literal_word',rank==depth,s0)
    internal[s0]=product
C=max(1/b**2 for b in internal.values())
maximizers=[s for s,b in internal.items() if 1/b**2==C]
expected_C=F(332306998946228968225951765070086144,211753373595316860199691410086465225)
require('frontier_exact_C',C==expected_C)
require('frontier_unique_C_maximizer',maximizers==[993],maximizers)
require('frontier_C_maximizer_rank_gain',V[993]==31 and internal[993]==F(460166680231540515,576460752303423488))
require('frontier_C_lt_1_57',C<F(157,100))

# Certified rational lower bound to log(1+x), from the positive atanh series.
# Checking our sum <= this LOWER bound is a valid exact numerical confirmation.
def log1p_lower(x):
    z=x/(2+x)
    total=F(0)
    power=z
    for j in range(32):
        total+=2*power/(2*j+1)
        power*=z*z
    return total

log_bounds={h:log1p_lower(F(18*h,4525))/18 for h in range(257)}
ranked_harmonic=F(0)
for h in range(257):
    if h:
        ranked_harmonic+=F(1,4543+18*(h-1))
    require('log_integral_direction_exact_control',ranked_harmonic<=log_bounds[h],h)

max_exterior=0
for initial,(full_product,depth) in path_profiles.items():
    s=initial
    seen=set()
    h=0
    exterior_product=F(1)
    paid=F(0)
    while s not in V:
        if s in seen:
            raise RuntimeError(('finite exterior repetition',initial,s))
        seen.add(s)
        require('exterior_q_at_least_4543',reciprocal_loss(s)>=4543,s)
        y,a,t=native(s)
        b=gain(s,a)
        exterior_product*=b
        paid+=(1-b)/(1+b)
        s=y
        h+=1
    require('exterior_product_factorization',full_product==exterior_product*internal[s],initial)
    require('exterior_polynomial_bound',1/exterior_product**9<=1+F(3*h,754),initial)
    require('combined_frontier_bound',(1/full_product**2/C)**9<=(1+F(3*h,754))**2,initial)
    require('exterior_logarithmic_forcing',paid<=log_bounds[h],initial)
    max_exterior=max(max_exterior,h)

# Synthetic visit multisets are not asserted to be legal native path words;
# they test the algebraic repeated-site repair independently of arrival.
exterior_sites=[s for s in range(1009,1300,2) if is_seed(s) and s not in V][:24]
for m in range(1,len(exterior_sites)+1):
    for phase in range(4):
        counts=[1+(i+phase)%5 for i in range(m)]
        M=max(counts)
        product=F(1)
        paid=F(0)
        for s,n in zip(exterior_sites[:m],counts):
            q=reciprocal_loss(s)
            product*=((q+1)/(q-1))**n
            paid+=F(n)/q
        require('repetition_inflation_repair',product**9<=(1+F(3*m,754))**M,(m,M))
        require('repetition_forcing_repair',paid<=M*log_bounds[m],(m,M))

# Actual marked-observer hostile, not merely an arbitrary two-vector example.
charges512,_=literal_source(512)
z,w=charges512[799],charges512[1065]
require('actual_marked_hostile_charges',(z,w)==((-2,-1),(0,1)))
require('actual_marked_hostile_arrival',799 in V and 1065 not in V and native(1065)==(799,2,2))
require('actual_marked_hostile_pair',pair(z,w)==0 and 2*row(z)[0]*row(w)[0]==6)
require('actual_marked_hostile_owner',abs((1065-799)//2)==133 and 133%2==1)
frontier_receipt={'size':len(V),'maximum_internal_inflation':str(C),'unique_maximizer':993,'maximizer_rank':V[993],'maximizer_gain':str(internal[993]),'maximum_exterior_depth_in_finite_universe':max_exterior}

expected={
 'frontier.txt':'df88f03b576b13c98c0c5378fc1cf61685eeb6817073681266d2a781dc329b26',
 'paid_root_README.md':'21b5845d0fea0aa2a5b1f347d208686d05f30b8f410fbf6f894f8d5a0fb16b3c',
 'root_eisenstein_first_arrival.md':'accf35a22c7fc6557e55a9098134eb6e6a74092abc372cc3f24d4f5a41556233'
}
base=Path(__file__).resolve().parent
for filename,digest in expected.items():
    require('pinned_source_bytes',sha256((base/'sources'/filename).read_bytes()).hexdigest()==digest,filename)
result={
 'status':'PASS',
 'source_commit':'96e35642110104cc547516a63f37eebcf930f707',
 'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
 'source_digest':source_hash.hexdigest(),
 'checks':dict(sorted(checks.items())),
 'source_universe':'all literal odd inputs m<4N for N=1..256',
 'cutoff_universe':'N=2..20000 plus all four residues for bases 2^60, 3^100, 10^150',
 'finite_first_hit_universe':'every primitive s<=24000; no cap hit, no repeats before ROOT',
 'path_count':path_count,'maximum_depth':max_depth,'maximum_clock':max_clock,
 'cutoff_examples':examples,'null_representative_examples':null_examples,
 'ROOT_inverse_examples':inverse_examples,'path_examples':path_examples,
 'forcing_convention_examples':forcing_examples,'actual_debt_ledgers':ledger_examples,
 'frontier':frontier_receipt,
 'claims_not_made':['universal ROOT arrival','a fixed native all-source Hilbert weight','uniform source debt bound','a positive-b and signed-c forcing identification','unmarked C3 energy as the original marked readout']
}
print(json.dumps(result,sort_keys=True,indent=2))
