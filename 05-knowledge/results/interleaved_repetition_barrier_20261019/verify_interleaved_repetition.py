#!/usr/bin/env python3
"""Standalone exact controls; theorems and Matveev input are in the note.

Complete macro-prefix universe: roots K,V, all four children through depth 8,
stop at the first coefficient crossing. Separate family controls include
unbounded-episode construction indices 2,4,10,100,4096. No network/file writes.
"""
from fractions import Fraction as F
from collections import Counter
from math import gcd
import hashlib
import json

checks=0

def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def branch(w):
    a,b,m=1,0,1
    for bit in w:
        if bit=='1':a,b=3*a,3*b+m
        elif bit!='0':raise ValueError(bit)
        m*=2
    return a,b,m

def comp(u,v):
    a,b,m=u;c,d,q=v
    return c*a,c*b+d*m,q*m

def power(u,k):
    a,b,m=u
    return a**k,b*(a**k-m**k)//(a-m),m**k

def source(u):
    a,b,m=u
    return -b*pow(a,-1,m)%m

def replay(x,w):
    old=x
    for bit in w:
        need(x%2==int(bit),'literal chronological guard')
        x=(3*x+1)//2 if x%2 else x//2
    a,b,m=branch(w)
    need(m*x==a*old+b,'literal endpoint identity')
    return x

def logs(z,m=160):
    lo=2*sum((z**(2*i+1)/F(2*i+1) for i in range(m)),F(0))
    hi=lo+2*z**(2*m+1)/(F(2*m+1)*(1-z*z))
    return lo,hi

def int_hash(x):
    return hashlib.sha256(x.to_bytes((x.bit_length()+7)//8,'big')).hexdigest()

bits={'K':'110110110110110101011011101010','F':'10','H':'110101000'}
bits['V']=bits['K']*4+bits['F']
data={z:branch(w) for z,w in bits.items()}
P,C0,Q=data['K'];D=P-Q
av,bv,mv=data['V']
C=F(bv,av-mv)
B=10**10;L0=2**48;H0=2**52
p=5750934602875680;q=9115015689657667
r4=source(power(data['K'],4))
need(r4==588943162011782273059357097228697595,'inherited K4 cylinder')
need(C==F(83829391903560839836419813590128166431,157489106280555890461927692448373459),'center constant')
need(C<533 and r4>C,'source exceeds center')
need(Q<P and 12*P<13*Q,'slow expansion rate')
need(F(av,mv)**30<F(P,Q)**122,'V rate per bit')
for z in 'KVFH':
    a,b,m=data[z]
    need(F(b,a)<=F(len(bits[z]),6),'step-aware carry tariff')
    if z!='F':need(F(b,a)<F(len(bits[z]),6),'strict non-F tariff')
    if z in 'KV':
        need(a>m and F(b,a-m)<=C,'expansion center envelope')
        for j in range(1,len(bits[z])+1):
            aa,bb,mm=branch(bits[z][:j]);need(aa>mm,'prefix expansion')
    else:need(F(a,m)<=F(3,4) and b<=m-a,'reset multiplier and positive monotonicity')
need(F(3,2)**2<F(4,3)**3,'linear reset clock')
# Forced K4 or immediate descent.
r1=source(data['K'])
for k in (1,2,3):
    for z in 'FH':
        a,b,m=comp(power(data['K'],k),data[z])
        need(a<m and (m-a)*r1>b,'short initial run descends uniformly')
# Rational short-range certificate.
l2,u2=logs(F(1,3));l3,u3=logs(F(1,2))
lo,hi=l2/u3,u2/l3
need(F(1,2)<l2<u2<1 and 1<l3<u3<2,'log enclosures')
need(gcd(p,q)==1 and q>2*H0,'coprime bracket and denominator horizon')
need(F(p,q)-F(1,q*q)<lo<hi<F(p,q)+F(1,q*q),'strict rational approximation')
need(F(H0,6*r4)<F(1,2),'small relative slope defect')
need(2*H0*q<3*r4,'short-range contradiction')
need(2136*H0*q<r4,'coarse-envelope independent margin')
# Long-range proof constants, no giant exponential evaluation.
need(2*30**5*32*2<B,'Matveev safe coefficient')
need(15*L0<H0,'short-long overlap')
need(F(L0,450)>53*B+52,'long-range base strictness')
need(L0>450*(B+1),'positive long-range derivative')
need(F(1,200)-F(1,360)==F(1,450),'repetition height rate')

leaf_count=0;node_count=0;frontier_count=0;factor_count=0;qualifying=0
literal_steps=0;max_episodes=0;hist=Counter();digest=hashlib.sha256()

def parsed_record(word):
    u=(1,0,1);prefix=[u];lengths=[0]
    for z in word:
        u=comp(u,data[z]);prefix.append(u);lengths.append(lengths[-1]+len(bits[z]))
    return u,prefix,lengths

def inspect_leaf(word,do_replay=True):
    global factor_count,qualifying,literal_steps,max_episodes
    u,prefix,lengths=parsed_record(word)
    a,b,m=u;N=lengths[-1]
    L=sum(len(bits[z]) for z in word if z in 'KV')
    need(a<m,'first-crossing leaf contracts')
    for aa,bb,mm in prefix[:-1]:need(aa>=mm,'no earlier crossing')
    need(F(b,a)<F(N,6),'normalized first-crossing carry')
    need(F(b,m)<F(N,6),'step-aware carry bound')
    need(N<=15*L,'interleaved time bound')
    x=source(u)
    need(x>0 and (m-a)*x>b,'all positive lifts descend in finite control')
    if do_replay:
        actual=''.join(bits[z] for z in word)
        need(branch(actual)==u,'literal and macro carrier agree')
        need(replay(x,actual)<x,'direct least-source descent')
        literal_steps+=N
    # Verify every prefix's carry debt and every expansion-only repeated factor.
    for i,(aa,bb,mm) in enumerate(prefix):
        need(F(bb,aa)<=F(lengths[i],6),'normalized prefix debt')
    found=False
    for i in range(len(word)):
        for nmac in range(1,(len(word)-i)//2+1):
            name=word[i:i+nmac]
            if any(z not in 'KV' for z in name):continue
            k=2
            while i+k*nmac<=len(word) and word[i:i+k*nmac]==name*k:
                factor_count+=1
                U=''.join(bits[z] for z in name);au,bu,mu=branch(U)
                ap,bp,mp=prefix[i];sp=F(ap,mp)
                y=(ap*x+bp)//mp
                need((ap*x+bp)%mp==0,'context endpoint integer')
                fuel=(au-mu)*y+bu
                need(fuel>0 and fuel%mu**k==0,'interior power fuel')
                r=(k-1)*len(U)
                Z=F(x)+F(N,6)+C
                need(F(y)+C<=sp*(F(x)+F(lengths[i],6)+C),'source-preserving context transport')
                need(F(y)+C>F(mu**k,au),'local fuel height')
                need(Z>F(2**r)/(sp*F(au,mu)),'original augmented-source height')
                if 100*r>=L:
                    qualifying+=1;found=True
                k+=1
    episodes=sum(z in 'KV' and (i==0 or word[i-1] in 'FH') for i,z in enumerate(word))
    max_episodes=max(max_episodes,episodes)
    hist[episodes]+=1
    return {'word':word,'N':N,'L_exp':L,'episodes':episodes,'qualifies':found,
            'source_bits':x.bit_length(),'source_sha256':int_hash(x)}

# Complete bounded forest; all four children are retained at every open node.
stack=[z for z in 'KV']
while stack:
    word=stack.pop();u,_,_=parsed_record(word);a,b,m=u
    if a<m:
        need(word[-1] in 'FH','crossing is reset endpoint')
        rec=inspect_leaf(word)
        digest.update(json.dumps(rec,sort_keys=True,separators=(',',':')).encode()+b'\n')
        leaf_count+=1
    else:
        node_count+=1
        if len(word)==8:frontier_count+=1
        else:stack.extend(word+z for z in 'KVFH')

# Closed-form family with genuinely unbounded episode count.
family=[]
for ell in (2,4,10,100,4096):
    ac,bc,mc=power(data['K'],ell)
    j=0
    while True:
        need(ac>=mc,'family pre-reset not crossed')
        ac,bc,mc=comp((ac,bc,mc),data['F']);j+=1
        if ac<mc:break
        ac,bc,mc=comp((ac,bc,mc),data['K'])
    word='K'*ell+'FK'*(j-1)+'F'
    N=30*(ell+j-1)+2*j
    L=30*(ell+j-1);r=30*(ell-1)
    need(j<=ell and 100*r>=L,'family macroscopic repeated factor')
    need(N<=15*L,'family clock')
    # Exact coefficients compiled from the full literal bit string provide a
    # second route even for the long control. Full direct replay is limited
    # to ell<=100 and is labelled below.
    literal=''.join(bits[z] for z in word)
    need(branch(literal)==(ac,bc,mc),'family independent literal coefficients')
    x=source((ac,bc,mc));y=(ac*x+bc)//mc
    need((ac*x+bc)%mc==0 and y<x,'family endpoint and all-lift descent')
    need(((P-Q)*x+C0)%Q**ell==0,'family initial repetition guard')
    need(F(bc,ac)<F(N,6),'family step-aware carry')
    if ell<=100:
        need(replay(x,literal)==y,'family direct parity replay');literal_steps+=N
    family.append({'initial_K_power':ell,'episodes':j,'N':N,'L_exp':L,
       'duplicated_bits':r,'source_bits':x.bit_length(),
       'source_sha256':int_hash(x),'endpoint_sha256':int_hash(y),
       'direct_parity_replay':ell<=100,'independent_literal_coefficients':True})
need(family[-1]['episodes']>1000,'nonvacuity beyond inherited episode cap')
need(F(P,Q)*F(3,4)<1 and F(P,Q)**2*F(3,4)<1,'unbounded family symbolic crossings')

# Hostile controls: the first-crossing and guard assumptions matter.
vals=[4,1,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3]
hostile=''.join('1'+'0'*(v-1) for v in vals)
aa,bb,mm=branch(hostile)
need(aa<mm and source((aa,bb,mm))==165,'generic contracting overshoot source')
need(replay(165,hostile)==167 and replay(165,'1000')==31,'generic word has earlier crossing/descent')
need(replay(1,'1010')==1,'contracting-square equality outside expansion-bank theorem')
aa,bb,mm=comp(data['K'],data['F'])
need(aa<mm and aa+bb>mm and source((aa,bb,mm))!=1,'formal wrong-guard overshoot')

report={'status':'PROVED scoped theorem with cited logarithmic-form input; finite controls reported separately',
 'theorem':'Every lawful first-crossing bank word with an expansion-bank factor U^k and 100(k-1)|U|>=L_exp descends at its endpoint; episode count unrestricted.',
 'new_carry_bound':'e/S<N/6 and e<N/6 at every first crossing',
 'short_range':{'H0':H0,'p':p,'q':q,'r4':r4,'series_terms_each_log':160,
   'gcd':gcd(p,q),'denominator_margin':q-2*H0,'spacing_margin':3*r4-2*H0*q,
   'coarse_envelope_margin':r4-2136*H0*q},
 'long_range':{'B':B,'L0':L0,'N_upper_at_join':15*L0,
   'base_lhs':str(F(L0,450)),'base_rhs_upper':53*B+52,
   'derivative_margin':L0-450*(B+1)},
 'finite_forest':{'max_macro_depth':8,'open_nodes':node_count,
   'uncrossed_depth_8_frontier':frontier_count,'first_crossing_leaves':leaf_count,
   'episode_histogram':dict(sorted(hist.items())),'max_episodes':max_episodes,
   'repeated_factors_tested':factor_count,'qualifying_factors':qualifying,
   'leaf_digest_sha256':digest.hexdigest()},
 'unbounded_episode_family_controls':family,
 'literal_replay_T_steps':literal_steps,'checks':checks,
 'external_source':{'display':'Languasco--Luca--Moree--Togbe, Theorem 2.1 (Matveev)',
   'url':'https://link.springer.com/article/10.1007/s12188-025-00293-9#Sec2'},
 'scope_limits':['Only exact first-crossing bank witnesses are covered.',
  'Uncrossed finite frontier is retained, not claimed empty.',
  'Full realized 2-adic boundary and arbitrary bank entry remain open.',
  'The large family control uses independent literal coefficients, not full parity replay.']}
print(json.dumps(report,indent=2,sort_keys=True))
