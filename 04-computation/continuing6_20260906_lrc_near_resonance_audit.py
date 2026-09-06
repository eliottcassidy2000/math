"""Independent referee: native spatial walls, subset masks, complete words.

No producer code is imported. Every gate remains active under -O.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd, comb, isqrt
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline='\n')
BASE=Path(__file__).resolve().parent
G=0
def need(ok,label):
    global G
    G+=1
    if not ok:raise ArithmeticError(label)
def locate(name):
    for path in (BASE/name,BASE.parent/'05-knowledge/results'/name,
                 Path('C:/w/s0905/04-computation')/name,
                 Path('C:/w/s0905/05-knowledge/results')/name):
        if path.is_file():return path
    raise FileNotFoundError(name)
def digest(data):
    return hashlib.sha256(json.dumps(data,separators=(',',':')).encode()).hexdigest()
def ceildiv(a,b):return -(-a//b)

raw=locate('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
need(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','raw full-profile supplier pin')
data=json.loads(raw)
profiles={int(k):{(a,tuple(w)) for a,w in value['profiles']} for k,value in data['levels'].items()}
certificate=json.loads(locate('continuing6_20260906_lrc_near_resonance_certificate.json').read_bytes())
need([c['clock'] for c in certificate['clocks']]==[14904,14898,14892],'exact new clock universe')

primes=[p for p in range(2,357) if all(p%d for d in range(2,isqrt(p)+1))]
def allowed_sum(total):
    for p in primes:
        e=0
        while total%p==0:
            total//=p;e+=1
        if e and (p%3!=2 or e>2):return False
    return total==1
def native(n,alpha,p,q):
    den=n*alpha.denominator
    return sum(all(14*min(v*(alpha.numerator+j*alpha.denominator)%den,
                         (-v*(alpha.numerator+j*alpha.denominator))%den)<den
                   for v in (p,q)) for j in range(n))
def spatial(p,q):
    L=14*p*q
    walls={0,L}
    for v in (p,q):
        d=L//(14*v)
        for j in range(v):
            walls.add((14*j-1)*d%L);walls.add((14*j+1)*d%L)
    cuts=sorted(walls)
    good=[]
    for a,b in zip(cuts,cuts[1:]):
        if all(14*min(v*(a+b)%(2*L),(-v*(a+b))%(2*L))<2*L for v in (p,q)):
            good.append((a,b))
    need(good[0][0]==0 and good[-1][1]==L,'literal origin component')
    return L,[(good[-1][0]-L,good[0][1])]+good[1:-1]
atlas=[]
geometry={}
for total in range(3,357):
    if not allowed_sum(total):continue
    for p in range(1,(total+1)//2):
        q=total-p
        if p<q and gcd(p,q)==1:
            atlas.append((p,q));geometry[p,q]=spatial(p,q)
need(len(atlas)==5855 and sorted(atlas)==[tuple(a) for a in certificate['atlas']],'independent whole strict pair atlas')

def valid(word):
    mask_gcd=[0]*128
    for mask in range(1,128):
        bit=mask&-mask
        mask_gcd[mask]=gcd(mask_gcd[mask^bit],word[bit.bit_length()-1])
    if mask_gcd[127]!=1:return False
    for mask in range(1,127):
        a=mask_gcd[mask]
        complement=tuple(sorted(gcd(a,word[j]) for j in range(7) if not(mask>>j&1)))
        if (a,complement) not in profiles[7-mask.bit_count()]:return False
    return True

total_words=total_candidates=0
for record in certificate['clocks']:
    T=record['clock']
    values=[d for d in data['levels']['6']['gcds'] if T%d==0]
    need(values==record['allowed_divisors'],'entire clock-specific state domain')
    costs={d:d*ceildiv(T,7*d) for d in values}
    words=[]
    maximum={pair:None for pair in combinations_with_replacement(values,2)}
    count=0
    for word in combinations_with_replacement(values,7):
        count+=1
        if not valid(word):continue
        E=sum(costs[d] for d in word)-T
        words.append([list(word),E])
        for i in range(7):
            for j in range(i+1,7):
                key=word[i],word[j]
                old=maximum[key]
                if old is None or E>old[0]:maximum[key]=(E,list(word))
    need(count==comb(len(values)+6,7)==record['word_universe'],'complete word universe including repeated occurrences')
    need(words==record['valid_words_and_costs'] and digest(words)==record['valid_word_sha256'],'every valid full word and native ceiling cost')
    M=max(E for w,E in words)
    need(M==record['global_maximum'] and len(words)==record['valid_word_count'],'attained full-word global maximum')
    need([[a,b,*maximum[a,b]] if maximum[a,b] is not None else [a,b,None,None]
          for a,b in maximum]==record['conditional_maxima'],'entire independently optimized conditional table and first owners')
    lookup={(e,p,q):(dp,dq,Csep,bound,mode) for e,p,q,dp,dq,Csep,bound,mode in record['all_candidates']}
    es=[e for e in range(1,7) if T%e==0]
    need(es==record['edge_gcds'] and len(lookup)==len(es)*5855==record['base_count'],'all unique sheet/edge candidates')
    modes=Counter();survivors=[]
    for e in es:
        n=T//e
        for p,q in atlas:
            L,intervals=geometry[p,q]
            credit=e*sum(ceildiv(n*(b-a),L)-1 for a,b in intervals)
            margins=e*gcd(n,p),e*gcd(n,q)
            bound=maximum.get(tuple(sorted(margins)))
            cost=None if bound is None else bound[0]
            mode='infeasible' if cost is None else 'separate' if credit>cost else 'located'
            need(lookup[e,p,q]==(*margins,credit,cost,mode),'every native candidate agrees, including infeasible margins')
            modes[mode]+=1
            if credit<=M:survivors.append([e,p,q,credit])
    need(modes==record['conditional_first_counts'] and survivors==record['global_survivors']==[[6,1,355,0]],'complete global and conditional exclusion accounting')
    need(record['located_credit']==300>M and record['global_surplus']==300-M,'actual located credit closes global cost')
    print('CLOCK',T,'WORDS',count,len(words),'GLOBAL',M,'NATIVE_CANDIDATES',len(lookup),'UNIQUE_LOCATION_SURVIVOR',survivors)
    total_words+=count;total_candidates+=len(lookup)

profiles_to_check=certificate['strict_lemma_controls']+certificate['equality_hostiles']+[c['located_profile'] for c in certificate['clocks']]
native_profiles=0
for profile in profiles_to_check:
    n,p,q=profile['n'],profile['p'],profile['q']
    L,I=geometry.get((p,q)) or spatial(p,q)
    walls=sorted({F(n*x%L,L) for interval in I for x in interval})
    phases=walls+[(walls[j]+walls[(j+1)%len(walls)]+(1 if j==len(walls)-1 else 0))/2%1 for j in range(len(walls))]
    actual=[[str(alpha),native(n,alpha,p,q)] for alpha in phases]
    need(actual==profile['walls_and_cells'],'literal entire translated-grid wall/cell profile')
    counts=[c for a,c in actual]
    need(min(counts)==profile['minimum'] and max(counts)==profile['maximum'],'exact native extrema')
    R,s=q//14,7*q-n
    centers=[F(1,2)-F(s*k,q) for k in range(-R,R+1)]
    half=F(s,14*q)
    for alpha,count in actual:
        x=F(alpha)
        absent=sum(min((x-c)%1,(c-x)%1)<=half for c in centers)
        need(count==2*R+1-absent,'native count equals number of closed absence gaps avoided')
    if profile['endpoint']:
        need(s==7 and q==14*R+1 and native(n,F(0),1,q)==2*R-1,'open danger endpoints cause two missing teeth at contact')
    else:
        need(s*(14*R+1)<7*q and min(counts)==2*R and max(counts)==2*R+1,'strict lemma controls')
    native_profiles+=1

for T,h,beta,value in certificate['multiplicity_controls']:
    alpha=F(beta)
    # Full physical grid uses the actual scaled pair, not the quotient pair.
    need(native(T,alpha,6*h,6*h*355)==value==6*native(T//6,(h*alpha)%1,1,355),'all physical sheet-multiplicity controls')

old=json.loads(locate('continuing5_20260906_lrc_clock16704_complete_certificate.json').read_bytes())['new_scales']
need(digest(old)==certificate['old_scales_sha256']=='ddbe0e091d36d54c8f6a7c8ea631bbf363d799b9adaa2ec4fe4cf56250d11a76','prior exact necessary-set semantic pin')
new=[t for t in old if t not in (14904,14898,14892)]
need(new==certificate['new_scales'] and len(new)==8198 and max(new)==14886,'exact three-clock deletion')
need(digest(new)==certificate['new_scales_sha256']=='c146d73c149c448e744e138aa3bb6f8c286748be81353b8a6605adc47a129117','new complete-set semantic pin')
print('COMPLETE_WORDS',total_words,'ALL_CANDIDATES',total_candidates,'NATIVE_PROFILES',native_profiles)
print('NEW_NECESSARY_SET 8198 clocks; maximum14886')
print('PASS',G,'always-active exact gates; actual LF')
