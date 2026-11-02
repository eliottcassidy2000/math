"""Exact wheels, marked-port Burnside sums and residue-coset jet controls."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import product,combinations,permutations
from math import gcd,prod,factorial
from pathlib import Path
import importlib.util,json,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def units(Q):return [a for a in range(Q) if gcd(a,Q)==1]
def words(Q,g):
    u=units(Q); out=Counter()
    for i,a in enumerate(u):
        offsets=[0];j=1
        while True:
            b=u[(i+j)%len(u)]+((i+j)//len(u))*Q
            if b-a>g:break
            offsets.append(b-a);j+=1
        if offsets[-1]==g:out[tuple(b-a for a,b in zip(offsets,offsets[1:]))]+=1
    return out
def contract(w,B):
    s=[0]
    for g in w:s.append(s[-1]+g)
    kept=[v for i,v in enumerate(s) if i not in B]
    return tuple(b-a for a,b in zip(kept,kept[1:]))
def transfer(w,p):
    s=[0]
    for g in w:s.append(s[-1]+g)
    blocks=defaultdict(list)
    for i,a in enumerate(s):blocks[a%p].append(i)
    out=Counter({w:p-len(blocks)})
    for B in blocks.values():
        if 0 not in B and len(w) not in B:out[contract(w,B)]+=1
    return +out
def oddparts(n,lo=1):
    if not n:yield ()
    for j in range(lo,n+1,2):
        for rest in oddparts(n-j,j):yield (j,)+rest
def weight(parts):
    counts=Counter(parts)
    e=sum((j-1)//2 for j in parts)+sum(gcd(i,j) for i,j in combinations(parts,2))
    z=prod(j**m*factorial(m) for j,m in counts.items())
    return e,z
def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def literal_port_burnside(m,r):
    edges={e for e in combinations(range(m+r),2) if e[1]>=r};total=0
    for body in permutations(range(r,r+m)):
        perm=tuple(range(r))+body;todo=set(edges);orbits=0;valid=True
        while todo:
            first=next(iter(todo));cur=first;sign=0
            while True:
                todo.remove(cur);a,b=perm[cur[0]],perm[cur[1]];sign^=a>b;cur=tuple(sorted((a,b)))
                if cur==first:break
            if sign:valid=False;break
            orbits+=1
        if valid:total+=1<<orbits
    need(total%factorial(m)==0,'literal port Burnside integer')
    return total//factorial(m)
def main():
    wheel_checks=0; short_checks=0; collision=None
    for Q,p in [(2,3),(6,5),(30,7),(210,11),(2310,13)]:
        for g in range(2,min(Q,62)+1,2):
            old=words(Q,g); actual=words(Q*p,g); predicted=Counter()
            for w,c in old.items():
                for v,a in transfer(w,p).items():predicted[v]+=c*a
            need(+predicted==actual,'full collision-block transfer')
            wheel_checks+=1
            before=Counter();after=Counter()
            for w,c in old.items():before[len(w)]+=c
            for w,c in actual.items():after[len(w)]+=c
            naive={j:(p-j-1)*before[j]+j*before[j+1] for j in range(1,g//2+1)}
            if g<2*p:
                need(all(naive[j]==after[j] for j in naive),'collision-free length recurrence');short_checks+=1
            elif collision is None and any(naive[j]!=after[j] for j in naive):
                collision={'Q':Q,'p':p,'g':g,'old_lengths':dict(before),'actual_lengths':dict(after),'unguarded_prediction':naive}
    need(collision is not None,'collision hostile exists')
    pair=[(2,4,2,4,6),(4,2,4,6,2)];profiles=[]
    for w in pair:
        need(w in words(30,18),'actual word quotient hostile')
        profile=Counter()
        for v,c in transfer(w,7).items():profile[len(v)]+=c
        profiles.append(dict(profile))
    need(profiles==[{5:1,4:4},{5:2,4:3}],'same-length distinct future')
    prime=[n for n in range(2,1001) if all(n%d for d in range(2,int(n**.5)+1))]
    records={};ranges=[];active=None
    for x in range(3,1001):
        P=[p for p in prime if p<=x];C=Counter(b-a for a,b in zip(P,P[1:]));champ=sorted(d for d,c in C.items() if c==max(C.values()))
        if champ==[4] and active is None:active=x
        if champ!=[4] and active is not None:ranges.append([active,x]);active=None
        if x in (131,138,139,491,546,547):records[x]={'counts':dict(C),'champions':champ,'prime_count':len(P)}
    need(ranges==[[131,139],[491,547]],'unique4 intervals through1000')
    pair_counts={}
    for d in [2,4,10,30,14,42,54,162]:
        Q=210;a=sum(gcd(x,Q)==gcd(x+d,Q)==1 for x in range(Q))
        b=prod(p-(1 if d%p==0 else 2) for p in [2,3,5,7]);need(a==b,'CRT pair count');pair_counts[d]=a
    need(pair_counts[2]==pair_counts[4] and pair_counts[30]==2*pair_counts[10] and pair_counts[42]==2*pair_counts[14] and pair_counts[54]==pair_counts[162],'named prime-support controls')
    patterns={}
    for S,expected in [((0,19,54),0),((0,38,108),120),((0,3,5,15,17),0),((0,6,10,30,34),36)]:
        literal=sum(all(gcd(x+s,2310)==1 for s in S) for x in range(2310))
        predicted=prod(p-len({s%p for s in S}) for p in [2,3,5,7,11])
        need(literal==predicted==expected,'supplied arithmetic patterns')
        patterns[str(S)]=literal
    counts={r:[] for r in range(3)}; partitions_seen=0
    for n in range(56):
        totals=[0,0,0]; fac=factorial(n)
        for parts in oddparts(n):
            e,z=weight(parts);need(fac%z==0,'permutation-class integrality')
            term=(fac//z)<<e
            for r in range(3):totals[r]+=term<<(r*len(parts))
            partitions_seen+=1
        for r in range(3):
            need(totals[r]%fac==0,'marked-port Burnside integrality');counts[r].append(totals[r]//fac)
    known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168]
    need(counts[0][:13]==known,'A000568 prefix')
    engine=ROOT/'05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py'
    inherited=load_module(engine,'local_scale_port_control')
    a,stats=inherited.local_dp(55);need(a==counts[0],'independent packed-response DP through55')
    for r in range(3):
        for m in range(7):need(literal_port_burnside(m,r)==counts[r][m],'literal marked edge-action control')
    # Exact functional differential identity, including a nonconstant marking.
    differential=0
    for n in range(1,18):
        y=lambda j:j+1
        lhs=sum(F(1<<weight(t)[0],weight(t)[1])*prod(y(j) for j in t) for t in oddparts(n))*n
        rhs=F(0)
        for i in range(1,n+1,2):
            rhs+=(1<<((i-1)//2))*y(i)*sum(F(1<<(weight(t)[0]+sum(gcd(i,j) for j in t)),weight(t)[1])*prod(y(j) for j in t) for t in oddparts(n-i))
        need(lhs==rhs,'marked functional differential equation');differential+=1
    # Coset norm controls: arbitrary modulus, fixed leading residue, all proper coefficients divisible by m.
    jet_checks=0
    for m in range(2,10):
        for a in range(1,m):
            delta=min(a,m-a)
            for b0 in [a-m,a,a+m]:
                for tail in product([-m,0,m],repeat=3):
                    b=[b0,*tail]
                    for x,y in [(F(0),F(0)),(F(1,3),F(-2,3)),(F(1),F(-1))]:
                        energy=b0*b0+(b[1]+x*b0)**2+(b[2]+y*b0)**2+(b[3]+x*b[2]+y*b[1]+x*y*b0)**2
                        need(energy>=delta*delta*(1+x*x)*(1+y*y),'all-modulus coset jet');jet_checks+=1
    hostile=load_module(ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py','symmetry_hostile_ports').hostile_checks()
    norm_checks=0
    for m in range(2,13):
        for a in range(1,m):
            delta=min(a,m-a)
            for b,c in product(range(-3,4),repeat=2):
                x=a+m*b;y=m*c;need(x*x+x*y+y*y>=delta*delta,'Eisenstein coset norm');norm_checks+=1
    out={'status':'PASS','wheel_transitions':wheel_checks,'collision_free_transitions':short_checks,'collision_hostile':collision,'length_quotient_hostile':{'words':pair,'profiles':profiles},'prime4_half_open_intervals':ranges,'prime_endpoint_records':records,'CRT_pair_counts_mod210':pair_counts,'odd_partitions_seen':partitions_seen,'A000568_through55':counts[0], 'fixed_one_port_body_counts_through55':counts[1],'fixed_two_port_body_counts_through55':counts[2], 'functional_identity_checks':differential,'literal_port_Burnside_checks':21,'jet_checks':jet_checks,'Eisenstein_norm_checks':norm_checks,'inherited_DP_stats':stats,'hostile_81_27':hostile}
    out['pattern_counts_mod2310']=patterns
    print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
