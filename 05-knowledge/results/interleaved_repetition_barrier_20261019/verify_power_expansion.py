#!/usr/bin/env python3
"""Exact controls for power-expansion zero-head and endpoint boundaries.

No network, file writes, floating-point tests, or removable assertions.
This checks finite controls and proof constants, not Matveev's theorem.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
import hashlib
import json

checks=0

def require(value,label):
    global checks
    checks+=1
    if not value:
        raise RuntimeError(label)

def branch(word):
    a,b,m=1,0,1
    for bit in word:
        if bit=='1': a,b=3*a,3*b+m
        elif bit!='0': raise ValueError(bit)
        m*=2
    return a,b,m

def compose(first,second):
    a,b,m=first
    c,d,q=second
    return c*a,c*b+d*m,q*m

def repeat(triple,k):
    a,b,m=triple
    return a**k,b*(a**k-m**k)//(a-m),m**k

def residue(triple):
    a,b,m=triple
    return (-b*pow(a,-1,m))%m

def replay(x,word):
    initial=x
    for bit in word:
        require(x%2==int(bit),'literal parity guard')
        x=(3*x+1)//2 if x%2 else x//2
    a,b,m=branch(word)
    require(m*x==a*initial+b,'literal affine endpoint')
    return x

def log_interval(z,terms=160):
    require(0<z<1 and terms>=1,'atanh parameters')
    lower=2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    tail=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return lower,lower+tail

K='110110110110110101011011101010'
Fword='10'
H='110101000'
V=K*4+Fword
bank={'K':K,'V':V,'F':Fword,'H':H}
data={name:branch(word) for name,word in bank.items()}
P,C0,Q=data['K']
D=P-Q
av,bv,mv=data['V']
C=F(bv,av-mv)
E=C+1
B=10**10
Nstar=918644030413102796640970006656455
L0=2**43
require((P,Q,C0,D)==(3**19,2**30,4912087879,88519643),'K branch')
require(C==F(83829391903560839836419813590128166431,157489106280555890461927692448373459),'V center')
require(0<F(C0,D)<C<533 and E<534,'center and carry envelopes')
for name in ('K','V'):
    a,b,m=data[name]
    require(a>m and F(b,a-m)<=C,'expansion center')
    require(3*a<4*m,'single expansion slope below four thirds')
    for j in range(1,len(bank[name])+1):
        aj,bj,mj=branch(bank[name][:j])
        require(aj>mj,'every expansion prefix grows')
for name in ('F','H'):
    a,b,m=data[name]
    require(0<a<m and F(a,m)<=F(3,4) and b<=m-a,'reset envelope')
    require(len(bank[name])<=9,'reset clock')
require(F(3,2)**2<F(4,3)**3,'reset count ratio')
require(2*30**5*32*2<B,'Matveev coefficient majorant')
require(F(2**40)-1>10+48*B,'long-range base strictness')
require(L0>8*B,'long-range positive derivative')
require(15*L0<2**47 and 15*L0<Nstar,'short-long overlap')
require(F(L0,8)+1>1068 and F(8560,8)+1>1068,'half-exponential threshold')
require(534<2**10,'log carry bound')
require(Q**2>2*C0 and 2*D<2**28,'K prefix boundary constants')

# Reproduce the inherited short-length certificate independently.
r1=residue(data['K'])
r4=residue(repeat(data['K'],4))
require(r1==815464443,'K source')
require(r4==588943162011782273059357097228697595,'K4 source')
small=compose(repeat(data['K'],3),data['F'])
a,b,m=small
require(a<m and (m-a)*r1>b,'K3 F descent')
require(F(81,512)<F(3,4) and F(101,512)<F(1,4),'H below F')
a0=253438635377865553122983825848491
b0=401690733307859070263635245069941
c0=326161216350563845551964374972875
d0=516953297105243726377334761586514
l2,u2=log_interval(F(1,3))
l3,u3=log_interval(F(1,2))
require(F(1,2)<l2<u2<1 and 1<l3<u3<2,'log height bounds')
lo,hi=l2/u3,u2/l3
eta=(E/r4)/(1-E/r4)
require(c0*b0-a0*d0==1,'neighbor determinant')
require(F(a0,b0)<lo<hi<F(c0,d0),'neighbor brackets')
require(b0*lo-a0>eta,'small-form exclusion')
require(b0+d0==Nstar,'certified inherited length')
ll,ul=log_interval(F(1,7),32)
require(ll>F(1,4),'log four thirds lower bound')

# Complete all F/H reset trees for each declared finite control.
# Roots: all 14 nonempty words U over K,V with 1..3 letters, k=2..5.
# Also include single-generator powers K^k,V^k for k=1..10.
cases=set()
for nletters in range(1,4):
    for letters in product('KV',repeat=nletters):
        for k in range(2,6): cases.add((''.join(letters),k))
for name in ('K','V'):
    for k in range(1,11):cases.add((name,k))

hash_rows=hashlib.sha256()
leaf_count=0
branch_count=0
literal_steps=0
endpoint_filters=0
hist=Counter()
max_length=0
max_reset_letters=0
samples=[]
for name,k in sorted(cases):
    U=''.join(bank[z] for z in name)
    ud=(1,0,1)
    for z in name:ud=compose(ud,data[z])
    require(ud==branch(U),'macro and literal expansion data')
    au,bu,mu=ud
    require(F(bu,au-mu)<=C,'composite expansion center')
    expansion=U*k
    ell=len(expansion)
    base=repeat(ud,k)
    require(base==branch(expansion),'geometric repeated carry')
    stack=[('',base)]
    while stack:
        reset,cur=stack.pop()
        ac,bc,mc=cur
        require(ac>=mc,'nonterminal reset-tree multiplier')
        branch_count+=1
        for z in ('F','H'):
            new=compose(cur,data[z])
            aa,bb,mm=new
            tail=reset+z
            if aa>=mm:
                stack.append((tail,new))
                continue
            full=expansion+''.join(bank[t] for t in tail)
            N=len(full)
            require(new==branch(full),'full macro and literal coefficients')
            require(F(bb,mm)<E,'single-episode carry')
            require(N<=15*ell,'linear reset completion bound')
            r=residue(new)
            require(r>0 and (mm-aa)*r>bb,'all-lift endpoint descent')
            y=replay(r,full)
            literal_steps+=N
            require(y<r,'direct least-source descent')
            # Lift source once and check the exact monotone inequality.
            lifted=r+mm
            yl=(aa*lifted+bb)//mm
            require((aa*lifted+bb)%mm==0 and yl<lifted,'positive next-lift control')
            if k>=2:
                fuel=(au-mu)*r+bu
                require(fuel>0 and fuel%mu**k==0,'repetition fuel divisibility')
                require(fuel>=mu**k,'positive fuel integer')
                require(F(r)+C>F(mu**k,au),'source-height lower bound')
                require(au<=3**len(U) and 2*len(U)<=ell,'power-density input')
            for d in range(9):
                num=bb-mm*d
                den=mm-aa
                require(num<=0 or num%den!=0,'endpoint displacement filter')
                endpoint_filters+=1
            leaf_count+=1
            hist[len(tail)]+=1
            max_length=max(max_length,N)
            max_reset_letters=max(max_reset_letters,len(tail))
            row=[name,k,tail,N,str(r),str(y)]
            hash_rows.update(json.dumps(row,separators=(',',':')).encode()+b'\n')
            if len(samples)<6:samples.append(row)

# Hostile: a contracting general word can overshoot; it is not a bank
# first crossing and has an earlier descent. Odd valuations from old compiler.
vals=[4,1,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3]
hostile=''.join('1'+'0'*(v-1) for v in vals)
aha,bha,mha=branch(hostile)
require(aha<mha,'hostile contracting multiplier')
require(residue((aha,bha,mha))==165,'hostile least source')
require(replay(165,hostile)==167,'hostile overshoot endpoint')
require(replay(165,'1000')==31,'hostile earlier descent')
# Formal slope and carry without guard cannot replace the source.
short=compose(data['K'],data['F'])
sa,sb,sm=short
require(sa<sm and sa+sb>sm,'formal K F source one overshoots')
require(residue(short)!=1,'formal source one lacks the guard')
# Generic contracting squares retain exact cycles, outside expansion bank.
require(replay(1,'1010')==1,'contracting-square equality hostile')

report={
 'status':'PROVED with cited Matveev input; finite exact controls only in declared universe',
 'theorem':'Every lawful U^k R with U in {K,V}+, k>=2, R in {F,H}+, and contracting total multiplier has endpoint below source.',
 'external_input':{'name':'Matveev bound, positive-rational specialization',
   'checked_display':'Languasco--Luca--Moree--Togbe, Theorem 2.1',
   'url':'https://link.springer.com/article/10.1007/s12188-025-00293-9#Sec2',
   'safe_coefficient':B},
 'proof_constants':{'L0':L0,'N_star':Nstar,'small_long_overlap_max':15*L0,
    'long_base_lhs':2**40-1,'long_base_rhs_majorant':10+48*B,
    'C':str(C),'E':str(E)},
 'finite_universe':{'expansion_roots':len(cases),'reset_tree_nonterminal_nodes':branch_count,
   'first_crossing_leaves':leaf_count,'reset_letter_histogram':dict(sorted(hist.items())),
   'endpoint_offsets_checked_in_each_leaf':'0 through 8',
   'endpoint_filter_checks':endpoint_filters,'literal_T_steps':literal_steps,
   'max_total_T_steps':max_length,'max_reset_letters':max_reset_letters,
   'leaf_data_sha256':hash_rows.hexdigest(),'sample_rows':samples},
 'hostiles':['165 -> 167 contracting general word, earlier descent to 31',
   'formal K F at 1 overshoots but violates the source guard',
   'contracting square (10)^2 returns 1 and is not an expanding-bank power'],
 'checks':checks,
 'scope_limits':['Full Lambda_infinity(S_d) membership remains open.',
   'No arbitrary-source bank entry or global Collatz theorem.',
   'No numerical tests are presented as a proof of Matveev.']}
print(json.dumps(report,indent=2,sort_keys=True))
