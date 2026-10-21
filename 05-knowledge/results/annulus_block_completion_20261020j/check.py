#!/usr/bin/env python3
"""Exact all-subset envelopes and source-owned temporal completion. No global bound."""
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, replace
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from itertools import combinations, product
from math import gcd
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('completion_owned',HERE.parent/'owned_pair_blocks_20261020i/check.py')
owned=module_from_spec(spec);sys.modules[spec.name]=owned;spec.loader.exec_module(owned)
native=owned.native
V=((1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),(-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1))
CHECKS=0

def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def integer(n,lo=0):
    if type(n) is not int or n<lo:raise ValueError('exact integer domain')
    return n

def typed_counts(c):
    if type(c) is not tuple or len(c)!=12:raise ValueError('immutable twelve-type multiplicities')
    for n in c:integer(n)
    return c

def add(v,w):return v[0]+w[0],v[1]+w[1]
def scale(n,v):return n*v[0],n*v[1]
def total(vs):
    a=b=0
    for x,y in vs:a+=x;b+=y
    return a,b
def Q(v):return v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
def beta(v,w):return Q(add(v,w))-Q(v)-Q(w)
def encode(obj):return (json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()

@dataclass(frozen=True)
class Envelope:
    counts:tuple
    source_norm:int
    maximum:int
    cut:int
    vector:tuple

def envelope(c):
    typed_counts(c)
    choices=[]
    for i in range(12):
        w=total(scale(c[j%12],V[j%12]) for j in range(i,i+6))
        choices.append((Q(w),w))
    k=max(q for q,w in choices);i=next(j for j,(q,w) in enumerate(choices) if q==k)
    return Envelope(c,sum(n*Q(v) for n,v in zip(c,V)),k,i,choices[i][1])

def audit_envelope(e):
    if type(e) is not Envelope:raise ValueError('exact envelope receipt')
    typed_counts(e.counts)
    for n in (e.source_norm,e.maximum,e.cut):integer(n)
    if type(e.vector) is not tuple or len(e.vector)!=2 or any(type(x) is not int for x in e.vector):
        raise ValueError('exact integral vector')
    ck(e==envelope(e.counts),'recomputed immutable envelope')
    return e

def counts(vs):
    out=[0]*12
    for v in vs:
        if type(v) is not tuple or len(v)!=2 or any(type(x) is not int for x in v) or v not in V:
            raise ValueError('literal signed norm-one/norm-three alphabet')
        out[V.index(v)]+=1
    return tuple(out)

def circuits():
    ans=[]
    for i,j in combinations(range(12),2):
        if add(V[i],V[j])==(0,0):
            c=[0]*12;c[i]=c[j]=1;ans.append(tuple(c))
    def det(v,w):return v[0]*w[1]-v[1]*w[0]
    for ids in combinations(range(12),3):
        a,b,c=(V[i] for i in ids);cs=(det(b,c),det(c,a),det(a,b))
        if all(n<0 for n in cs):cs=tuple(-n for n in cs)
        if not all(n>0 for n in cs):continue
        g=gcd(gcd(cs[0],cs[1]),cs[2]);out=[0]*12
        for i,n in zip(ids,cs):out[i]=n//g
        ans.append(tuple(out))
    return tuple(sorted(ans))

def brute(c):
    ids=[i for i,n in enumerate(c) if n]
    return max(Q(total(scale(n,V[i]) for i,n in zip(ids,ns)))
               for ns in product(*(range(c[i]+1) for i in ids)))

def reject(fn):
    try:fn()
    except (ValueError,TypeError):return True
    return False

def controls():
    # Every binary multiset and every subset: 3^12 independent subset choices.
    vec=[(0,0)]*4096
    for mask in range(1,4096):
        bit=mask&-mask;vec[mask]=add(vec[mask^bit],V[bit.bit_length()-1])
    for mask in range(4096):
        c=tuple((mask>>i)&1 for i in range(12));e=envelope(c);best=0;sub=mask
        while True:
            q=Q(vec[sub]);best=max(best,q);ck(q<=e.maximum,'all binary subset energies')
            if sub==0:break
            sub=(sub-1)&mask
        ck(best==e.maximum,'exact binary envelope')
        U=sum(c[i] for i in range(0,12,2));W=sum(c[i] for i in range(1,12,2))
        ck(27*e.maximum>=3*U*U+9*U*W+9*W*W>=e.source_norm**2,'sharp rotation bound')
    cs=circuits();ck(len(cs)==46,'six pairs and forty primitive triples')
    hist=Counter()
    for c in cs:
        e=envelope(c);active=[n for n in c if n]
        ck(total(scale(n,v) for n,v in zip(c,V))==(0,0),'circuit zero')
        ck(e.maximum==brute(c)<=e.source_norm,'primitive all-order safety')
        ck(e.maximum==max(n*n*Q(v) for n,v in zip(c,V)),'circuit corner formula')
        hist[(tuple(sorted(active)),e.source_norm,e.maximum)]+=1
    sharp=tuple(0 if i%2==0 else (2 if (i//2)%2==0 else 1) for i in range(12))
    e=envelope(sharp);ck((sum(sharp),e.source_norm,e.maximum)==(9,27,27),'sharp integer zero packet')
    ck(total(scale(n,v) for n,v in zip(sharp,V))==(0,0) and brute(sharp)==27,'sharp independent subset test')
    hostile=[0]*12;hostile[1]=hostile[7]=3;h=envelope(tuple(hostile))
    ck((h.source_norm,h.maximum)==(18,27),'three overlapping opposite-pair hostile')
    for u,w in product(range(8),repeat=2):
        c=tuple(u if i%2==0 else w for i in range(12));p=envelope(c)
        ck(p.maximum==4*u*u+12*u*w+12*w*w,'rotation-symmetric form')
    ck(reject(lambda:envelope((True,)+(0,)*11)),'reject bool multiplicity')
    ck(reject(lambda:envelope((1.0,)+(0,)*11)),'reject float multiplicity')
    ck(reject(lambda:audit_envelope(replace(e,source_norm=27.0))),'reject equal float receipt')
    ck(reject(lambda:counts(((True,0),))),'reject bool vector')
    return {'binary_count_vectors':4096,'subset_choices':3**12,
            'circuit_histogram':[[list(k[0]),k[1],k[2],v] for k,v in sorted(hist.items())],
            'sharp_counts':sharp,'sharp_envelope':e.__dict__,'hostile_envelope':h.__dict__}

def packet(N,rows,routes):
    integer(N,1)
    if type(rows) is not tuple or not rows:raise ValueError('nonempty immutable physical packet')
    seen=set();term=None;vs=[]
    for r in rows:
        owned.audit_row(N,r,routes)
        ck(4*N<=r.source<256*N and r.source not in seen,'distinct new physical source')
        seen.add(r.source)
        if term is None:term=r.terminal
        ck(term==r.terminal,'same actual terminal')
        v=owned.signed(r);ck(v in V,'nonzero signed native alphabet');vs.append(v)
        ck(r.unsigned==native.independent_charge(64*N,r.source),'independent selected charge')
    e=audit_envelope(envelope(counts(tuple(vs))))
    profile=[]
    for d in sorted({0,*[r.depth for r in rows]}):
        w=total(v for r,v in zip(rows,vs) if r.depth<=d)
        profile.append([d,*w,Q(w),sum(Q(v) for r,v in zip(rows,vs) if r.depth<=d)])
    return {'N':N,'terminal':term,'rows':[owned.payload(r) for r in rows],
            'envelope':e.__dict__,'actual_profile':profile,
            'actual_maximum':max(r[3] for r in profile),'final_vector':profile[-1][1:3]}

def source_controls(N,rows,routes,sharp):
    groups=defaultdict(list)
    for r in rows:
        if r.source>=4*N and r.terminal==1 and owned.signed(r) in V:
            groups[owned.signed(r)].append(r)
    selected=[]
    for v,c in zip(V,sharp):
        ck(len(groups[v])>=c,'actual sharp-template availability')
        selected.extend(sorted(groups[v],key=lambda r:r.source)[:c])
    s=packet(N,tuple(selected),routes)
    ck(s['envelope']['source_norm']==s['envelope']['maximum']==27 and s['final_vector']==[0,0],
       'native sharp nine-source completion')
    pairs=owned.make_pairs(N,rows);table={r.source:r for r in rows}
    candidates=[]
    for p in pairs:
        x,y=table[p.positive],table[p.negative]
        if x.terminal!=711 or owned.signed(x)!=(2,1) or x.depth==y.depth:continue
        lo,hi=sorted((x.depth,y.depth));candidates.append((hi,lo,p.positive,p.negative))
    selected=[];end=-1;intervals=[]
    for hi,lo,p,m in sorted(candidates):
        if lo<end:continue
        selected.extend((table[p],table[m]));intervals.append([lo,hi,p,m]);end=hi
    ck(len(intervals)>=5,'native temporal block exceeds timing-blind capacity')
    t=packet(N,tuple(selected),routes);t['exposure_intervals']=intervals
    ck(t['envelope']['source_norm']>27 and t['actual_maximum']==3,'large clock-aware safe packet')
    ck(t['envelope']['maximum']>t['envelope']['source_norm'] and t['final_vector']==[0,0],
       'clock sidecar genuinely necessary')
    ck(reject(lambda:packet(N,(selected[0],selected[0]),routes)),'duplicate source rejection')
    ck(reject(lambda:packet(N,(replace(selected[0],clock=selected[0].clock+2),),routes)),
       'same parity wrong full clock rejected')
    return s,t

def sibling_blocks(N,routes):
    """Reconstruct the inherited deterministic consecutive near-sibling bank."""
    frontier=routes.V;blocks=[];used=set()
    for b,(lo,hi) in enumerate(zip((4,8,16,32,64),(8,16,32,64,128))):
        run=[];runs=[]
        for x in range(lo*N+1,hi*N,8):
            if x in frontier or 2*x+1 in frontier:
                if run:runs.append(run);run=[]
            else:run.append(x)
        if run:runs.append(run)
        for run in runs:
            for j in range(0,len(run)-2,3):
                terms=[];ids=[]
                for x in run[j:j+3]:
                    y=2*x+1;u,a,d=routes.get(x);z,c,e=routes.get(y)
                    ck(u==z and d==e and c==a+1,'actual sibling clock and arrival')
                    w=-beta(native.charge(64*N,x),native.charge(64*N,y))
                    terms.append((w,d,x,y));ids.extend((x,y))
                ck(not used.intersection(ids),'physical disjointness of inherited triple bank');used.update(ids)
                positive=[d for w,d,x,y in terms if w==3];negative=[d for w,d,x,y in terms if w==-3]
                ck((not positive and not negative) or (len(positive)==len(negative)==1),'native signed triple pattern')
                if not positive:sign=0;start=end=0
                else:
                    start,end=sorted((positive[0],negative[0]))
                    sign=0 if start==end else (1 if positive[0]<negative[0] else -1)
                blocks.append({'id':len(blocks),'band':b,'terms':tuple(terms),'sign':sign,'start':start,'end':end})
    ck(len(blocks)==21152,'inherited fixed 4096 bank')
    return blocks

def audit_block(N,b,routes):
    integer(N,1)
    if type(b) is not dict or set(b)!={'id','band','terms','sign','start','end'}:
        raise ValueError('exact native block fields')
    for key in ('id','band','start','end'):integer(b[key])
    if type(b['sign']) is not int or b['sign'] not in (-1,0,1):raise ValueError('signed interval type')
    if type(b['terms']) is not tuple or len(b['terms'])!=3:raise ValueError('three immutable native terms')
    ck(b['band']<5,'native sibling band')
    xs=[]
    for term in b['terms']:
        if type(term) is not tuple or len(term)!=4 or any(type(x) is not int for x in term):
            raise ValueError('exact native pair term')
        w,d,x,y=term;integer(d);integer(x,1);integer(y,1)
        ck(x%8==1 and y==2*x+1 and x not in routes.V and y not in routes.V,'native sibling source guards')
        ck((4<<b['band'])*N<=x<(8<<b['band'])*N and y<256*N,'source annulus band')
        u,a,dx=routes.get(x);v,c,dy=routes.get(y)
        ck(u==v and dx==dy==d and c==a+1,'reauthenticated actual first arrivals')
        ck(w==-beta(native.charge(64*N,x),native.charge(64*N,y)),'source-bound pair interaction')
        xs.append(x)
    ck(xs[1]==xs[0]+8 and xs[2]==xs[1]+8,'consecutive native triple')
    pos=[d for w,d,x,y in b['terms'] if w==3];neg=[d for w,d,x,y in b['terms'] if w==-3]
    if not pos and not neg:sign=start=end=0
    else:
        ck(len(pos)==len(neg)==1,'native signed term alphabet');start,end=sorted((pos[0],neg[0]))
        sign=0 if start==end else (1 if pos[0]<neg[0] else -1)
    ck((b['sign'],b['start'],b['end'])==(sign,start,end),'exact interval source readout')

def interval_guard(container,jobs):
    if type(container) is not tuple or len(container)!=2 or type(jobs) is not tuple:
        raise ValueError('immutable interval data')
    lo,hi=container;integer(lo);integer(hi);ck(lo<hi,'nonempty negative interval')
    for p in jobs:
        if type(p) is not tuple or len(p)!=2:raise ValueError('immutable positive interval')
        integer(p[0]);integer(p[1]);ck(lo<=p[0]<p[1]<=hi,'positive interval containment')
    ck(bool(jobs),'nonempty completion')
    ordered=sorted(jobs)
    ck(all(a[1]<=b[0] for a,b in zip(ordered,ordered[1:])),'no simultaneous reuse of negative amplitude')

def audit_interval_certificate(N,blocks,negative,positive,routes):
    """Reauthenticate every selected native term; all omitted interactions remain."""
    if type(blocks) is not tuple or type(positive) is not tuple:raise ValueError('immutable block bank/certificate')
    integer(negative);ids=(negative,)+positive
    for i in ids:integer(i)
    ck(len(set(ids))==len(ids) and all(i<len(blocks) for i in ids),'distinct retained block IDs')
    for i in ids:
        audit_block(N,blocks[i],routes);ck(blocks[i]['id']==i,'canonical block address')
    n=blocks[negative];ck(n['sign']==-1,'negative container')
    jobs=[blocks[i] for i in positive]
    ck(all(p['sign']==1 for p in jobs),'positive jobs')
    interval_guard((n['start'],n['end']),tuple((p['start'],p['end']) for p in jobs))
    addresses=[tuple(sorted((x,y))) for i in ids for w,d,x,y in blocks[i]['terms']]
    ck(len(set(addresses))==len(addresses),'pair-address ownership')
    for d in sorted({0,*[t for i in ids for w,t,x,y in blocks[i]['terms']]}):
        ck(sum(w for i in ids for w,t,x,y in blocks[i]['terms'] if t<=d)<=0,'all-horizon combined selected debit')
    return ids

def interval_completion(N,routes):
    blocks=tuple(sibling_blocks(N,routes));jobs=defaultdict(deque)
    for b in blocks:
        if b['sign']==1:jobs[b['start'],b['end']].append(b['id'])
    keys=sorted(jobs,key=lambda k:(k[1],k[0]));matches=[];used=set()
    negatives=sorted((b for b in blocks if b['sign']==-1),key=lambda b:(b['start']-b['end'],b['start'],b['id']))
    for n in negatives:
        end=n['start'];chosen=[]
        for lo,hi in keys:
            if end<=lo<hi<=n['end'] and jobs[lo,hi]:
                chosen.append(jobs[lo,hi].popleft());end=hi
        if chosen:
            ids=audit_interval_certificate(N,blocks,n['id'],tuple(chosen),routes)
            ck(not used.intersection(ids),'no block term reused across completion packets');used.update(ids)
            matches.append([n['id'],chosen])
    ck(matches and max(len(p) for n,p in matches)>=2,'actual larger multi-petal completion')
    before=[0]*257;after=[0]*257;selected_neg=[0]*257
    for n,ps in matches:
        for d in range(257):
            nv=sum(w for w,t,x,y in blocks[n]['terms'] if t<=d)
            pv=sum(w for i in ps for w,t,x,y in blocks[i]['terms'] if t<=d)
            ck(nv+pv<=0,'combined interaction certificate at every finite cut')
            before[d]+=pv;after[d]+=nv+pv;selected_neg[d]+=nv
    representative=max(matches,key=lambda p:(len(p[1]),-p[0]));repids=[representative[0],*representative[1]]
    repsources=sorted({s for i in repids for w,d,x,y in blocks[i]['terms'] for s in (x,y)})
    literal=[]
    for s in repsources:
        u,t,d=routes.get(s);r=owned.Row(s,s%15,s%35,u,t,d,native.charge(64*N,s));owned.audit_row(N,r,routes)
        literal.append(owned.payload(r))
    # A single negative interval cannot be spent twice at the same horizon.
    ck(reject(lambda:interval_guard((1,5),((2,4),(3,4)))),'reject overlapping positive claims')
    ck(reject(lambda:audit_interval_certificate(N,blocks,True,(0,),routes)),'reject bool block ID')
    ck(reject(lambda:audit_interval_certificate(N,blocks,matches[0][0],(matches[0][1][0],)*2,routes)),'reject reused positive block')
    bank_events=[0]*257;positive_events=[0]*257
    for b in blocks:
        for w,d,x,y in b['terms']:bank_events[d]+=w
        if b['sign']==1:
            positive_events[b['start']]+=3;positive_events[b['end']]-=3
    bank=[];positive=[];p=q=0
    for d in range(257):
        p+=bank_events[d];q+=positive_events[d];bank.append(p);positive.append(q)
    unmatched=[p-c for p,c in zip(positive,before)]
    ck(all(x>=0 and x%3==0 for x in unmatched),'unmatched positive interval count')
    return {'scale':N,'inherited_blocks':len(blocks),'block_sha256':sha256(encode(blocks)).hexdigest(),
            'sign_counts':dict(Counter(b['sign'] for b in blocks)),
            'matching_rule':'longest negative container; earliest-ending available disjoint positive intervals',
            'matches':matches,'containers':len(matches),'completed_positive_blocks':sum(len(ps) for n,ps in matches),
            'max_positive_blocks_per_container':max(len(ps) for n,ps in matches),
            'selected_positive_debit_by_depth':before,'selected_negative_debit_by_depth':selected_neg,
            'combined_debit_by_depth':after,'representative_blocks':[blocks[i] for i in repids],
            'complete_sibling_bank_debit':bank,'all_positive_debit':positive,'unmatched_positive_debit':unmatched,
            'maximum_upper_bound_gain':max(before),'gain_maximizing_depths':[d for d,x in enumerate(before) if x==max(before)],
            'maximum_unmatched_debit':max(unmatched),'unmatched_maximizing_depths':[d for d,x in enumerate(unmatched) if x==max(unmatched)],
            'representative_source_rows':literal,
            'scope':'only listed pair-address interactions; no source-norm credit and all other cross terms remain'}

def full_source_ledger(N,rows,temporal,intervals):
    """One final insertion, retaining every other original source and clock."""
    ids={r[0] for r in temporal['rows']};terminal=temporal['terminal'];v=(2,1)
    batches=defaultdict(list);digest=sha256()
    for r in rows:batches[r.depth].append(r);digest.update(encode(owned.payload(r)))
    B=E=(0,0);W=defaultdict(lambda:(0,0));H=defaultdict(lambda:(0,0));A=0
    insertion=[];ledger=[];T=temporal['envelope']['source_norm']
    for d in range(257):
        for r in batches[d]:
            z=owned.signed(r)
            if r.source<4*N:W[r.terminal]=add(W[r.terminal],z)
            else:H[r.terminal]=add(H[r.terminal],z);A+=Q(z)
            if r.terminal==terminal:
                if r.source in ids:B=add(B,z)
                else:E=add(E,z)
        ck(B in ((0,0),v,scale(-1,v)),'literal time-disjoint block vector')
        delta=Q(B)+beta(E,B);sigma=0 if B==(0,0) else (1 if B==v else -1)
        ck(beta(E,v)==3*E[0],'directional native projection')
        ck(delta==(0 if sigma==0 else 3+3*sigma*E[0]),'exact directional insertion invoice')
        ck(delta<=T,'all-horizon full last-insertion payment')
        insertion.append([d,*B,*E,Q(E),None if sigma==0 else sigma*E[0],delta])
        energy=sum(Q(z) for z in H.values());mixed=sum(beta(W[u],z) for u,z in H.items())
        I=energy-A+mixed;M=153*N-A;R=I-intervals['complete_sibling_bank_debit'][d]
        old=R+intervals['all_positive_debit'][d];new=R+intervals['unmatched_positive_debit'][d]
        ck(I<=new<=old and old-new==intervals['selected_positive_debit_by_depth'][d],
           'source-owned monotone residual upper bound')
        ledger.append([d,A,M,energy,mixed,I,R,old,new])
    ck(A==153*N,'full finite annulus source norm')
    maximum=max(r[-1] for r in insertion)
    ck(maximum==156 and [r[0] for r in insertion if r[-1]==maximum]==[33],
       'fixed native full-interaction positive control')
    ck(max(r[6] for r in insertion if r[6] is not None)==51,'active directional maximum')
    return {'source_count':len(rows),'source_rows_sha256':digest.hexdigest(),
            'last_insertion_scope':'one block inserted after ALL other new sources; old baseline retained',
            'last_insertion_rows':insertion,'insertion_budget':T,'maximum_insertion':maximum,
            'maximizing_depths':[r[0] for r in insertion if r[-1]==maximum],
            'minimum_insertion':min(r[-1] for r in insertion),'maximum_active_signed_external_a':51,
            'ledger_columns':['d','arriving_new_norm','moving','new_energy','old_new_cross','I','R','old_upper','new_upper'],
            'full_ledger':ledger,'old_upper_max_excess':max(r[7]-r[2] for r in ledger),
            'new_upper_max_excess':max(r[8]-r[2] for r in ledger),
            'I_max_excess':max(r[5]-r[2] for r in ledger)}

def branch(s):
    integer(s,1)
    if s%2==0 or s%8==5:raise ValueError('original primitive source')
    return 2 if s%8==1 else 1

def lift_mod(s,k,modulus):
    """Exact source L_a^k(s) modulo any positive modulus; retain division by63."""
    a=branch(s);integer(k);integer(modulus,1);c=35 if a==1 else 49
    numerator=((63*s+c)*pow(64,k,63*modulus)-c)%(63*modulus)
    ck(numerator%63==0,'full-denominator lift residue')
    return numerator//63

def rotate(v,k):
    for _ in range(k%3):v=(-v[1],v[0]-v[1])
    return v

def phase_lift_controls(N,rows,temporal,routes):
    table={r.source:r for r in rows};result=[]
    for a in (1,2):
        chosen=[p for p in temporal['exposure_intervals'] if branch(p[2])==branch(p[3])==a]
        base=tuple(table[s] for lo,hi,x,y in chosen for s in (x,y));c=35 if a==1 else 49
        ck(len(base)==(48 if a==1 else 20),'retained same-branch subpacket')
        initial=packet(N,base,routes)
        ck(initial['actual_maximum']==3 and initial['envelope']['source_norm']>27,'infinite-lift seed capacity')
        samples=[]
        for k in range(7):
            factor=64**k;lifted=[]
            for r in base:
                s=factor*r.source+c*(factor-1)//63
                ck(s not in routes.V and branch(s)==a,'lift outside stopped frontier on same branch')
                u,t,d=routes.get(s)
                ck((u,t,d)==(r.terminal,r.clock+6*k,r.depth),'actual lifted terminal clock and depth')
                z=native.charge(64*factor*N,s)
                ck(z==rotate(r.unsigned,c*k),'exact annulus band rotation')
                rr=owned.Row(s,s%15,s%35,u,t,d,z);lifted.append(rr)
                for modulus in (1,15,35,63,105,315,1024):
                    ck(lift_mod(r.source,k,modulus)==s%modulus,'modular reader versus literal lift')
            p=packet(factor*N,tuple(lifted),routes)
            ck(p['actual_maximum']==3 and p['envelope']['source_norm']==initial['envelope']['source_norm'],
               'literal all-cut profile preservation')
            samples.append({'k':k,'N':factor*N,'row_sha256':sha256(encode(p['rows'])).hexdigest(),
                            'source_mod105':[r.source%105 for r in lifted],'actual_maximum':p['actual_maximum']})
        registers=[[lift_mod(r.source,k,105) for r in base] for k in range(6)]
        for r in base:
            for k in range(6):ck(lift_mod(r.source,k+6,105)==lift_mod(r.source,k,105),'joint-register six-period')
        result.append({'branch':a,'translation':c,'base_packet':initial,'exposure_intervals':chosen,
                       'six_joint_mod105_rows':registers,'literal_controls':samples,
                       'external_field_scope':'internal profile lifts for all k; later full insertion still needs its directional guard'})
    for args in ((True,0,15),(3,False,15),(3,1,35.0)):
        ck(reject(lambda args=args:lift_mod(*args)),'typed symbolic reader')
    return result

def run():
    finite=controls()
    f=json.loads((HERE.parent/'session18_native_annulus/frontier.json').read_text())
    routes=native.Routes(f['frontier'],dict(f['frontier_ranks']))
    rows=owned.enumerate_rows(4096,routes)
    sharp,temporal=source_controls(4096,rows,routes,tuple(finite['sharp_counts']))
    intervals=interval_completion(4096,routes)
    full=full_source_ledger(4096,rows,temporal,intervals)
    lifts=phase_lift_controls(4096,rows,temporal,routes)
    # Small complete hostile ledger; this bound fails, so no unconditional upgrade.
    ledger=native.decomposition(7,routes)[0]
    row=ledger[8]
    ck(row[4:]==[1054,17,32,30,-13],'literal N7 finite-horizon obstruction retained')
    return {'status':'PROVED scoped envelopes and temporal guards; FINITE-EXACT native controls; uniform 153N OPEN',
            'finite_controls':finite,'native_sharp_packet':sharp,'native_temporal_packet':temporal,
            'interval_completion':intervals,'full_native_control':full,'same_branch_lifts':lifts,
            'N7_depth8_ledger':row,'route_audit':routes.audit(),
            'checks':CHECKS}

if __name__=='__main__':
    obj=run();data=(json.dumps(obj,sort_keys=True,indent=2)+'\n').encode();path=HERE/'receipt.json'
    if '--write' in sys.argv:path.write_bytes(data)
    else:ck(path.read_bytes().replace(b'\r\n',b'\n')==data,'saved receipt agrees')
    print(json.dumps({'checks':CHECKS,'sha256':sha256(data).hexdigest(),
                     'native_sharp_source_ids':[r[0] for r in obj['native_sharp_packet']['rows']],
                     'native_temporal_mass':obj['native_temporal_packet']['envelope']['source_norm'],
                     'native_temporal_max':obj['native_temporal_packet']['actual_maximum'],
                     'interval_completed':obj['interval_completion']['completed_positive_blocks'],
                     'interval_containers':obj['interval_completion']['containers'],
                     'maximum_petals':obj['interval_completion']['max_positive_blocks_per_container']}))
