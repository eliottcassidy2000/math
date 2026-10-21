#!/usr/bin/env python3
"""Source-authenticated temporal opposite-pair normal form; no universal claim."""
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from itertools import combinations_with_replacement, permutations
import json
import heapq
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'session18_native_annulus'
spec=spec_from_file_location('owned_pair_native',BASE/'verify_annulus.py')
native=module_from_spec(spec);spec.loader.exec_module(native)
CHECKS=0

def check(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(message)

def integer(x,minimum=0):
    if type(x) is not int or x<minimum:raise ValueError('exact integer domain')
    return x

def vector(v):
    if type(v) is not tuple or len(v)!=2 or any(type(x) is not int for x in v):
        raise ValueError('exact vector required')
    return v

def add(v,w):return v[0]+w[0],v[1]+w[1]
def neg(v):return -v[0],-v[1]
def norm(v):return v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
def beta(v,w):return 2*v[0]*w[0]-v[0]*w[1]-v[1]*w[0]+2*v[1]*w[1]
def encode(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()

@dataclass(frozen=True)
class Row:
    source:int
    mod15:int
    mod35:int
    terminal:int
    clock:int
    depth:int
    unsigned:tuple

def typed(r):
    if type(r) is not Row:raise ValueError('source row required')
    for x in (r.source,r.terminal):
        integer(x,1)
        if x%2!=1 or x%8==5:raise ValueError('primitive source/terminal')
    for x in (r.mod15,r.mod35,r.clock,r.depth):integer(x)
    vector(r.unsigned)
    if r.mod15!=r.source%15 or r.mod35!=r.source%35:raise ValueError('joint source registers')
    return r

def signed(r):
    typed(r)
    return r.unsigned if r.clock%2==0 else neg(r.unsigned)

def payload(r):
    typed(r)
    return [r.source,r.mod15,r.mod35,r.terminal,r.clock,r.depth,*r.unsigned,norm(r.unsigned)]

def audit_row(N,r,routes):
    integer(N,1);typed(r)
    check(r.source<256*N,'original source domain')
    check(r.unsigned==native.charge(64*N,r.source),'full enlarged source charge')
    check((r.terminal,r.clock,r.depth)==routes.get(r.source),'actual first-arrival source binding')
    check(r.mod15%5==r.mod35%5,'compatible fifteen/thirty-five registers')
    return r

def axis(v):
    vector(v)
    if v==(0,0):raise ValueError('zero has no opposite axis')
    return max(v,neg(v))

@dataclass(frozen=True)
class Pair:
    positive:int
    negative:int

def make_pairs(N,rows,chronological=True):
    integer(N,1)
    if type(rows) is not tuple:raise ValueError('immutable source rows')
    if type(chronological) is not bool:raise ValueError('explicit matching policy')
    seen=set();groups=defaultdict(lambda:[[],[]])
    for r in rows:
        typed(r)
        if r.source in seen:raise ValueError('duplicate physical source')
        seen.add(r.source)
        v=signed(r)
        if r.source<4*N or v==(0,0):continue
        a=axis(v);groups[r.terminal,a][v!=a].append(r)
    pairs=[]
    key=(lambda r:(r.depth,r.source)) if chronological else (lambda r:r.source)
    for _,(p,m) in sorted(groups.items()):
        p.sort(key=key);m.sort(key=key)
        pairs.extend(Pair(x.source,y.source) for x,y in zip(p,m))
    return tuple(pairs)

def audit_pairs(N,rows,pairs,canonical=True):
    integer(N,1)
    if type(rows) is not tuple or type(pairs) is not tuple or type(canonical) is not bool:
        raise ValueError('immutable source/pair bank')
    table={}
    for r in rows:
        typed(r)
        if r.source in table:raise ValueError('duplicate source')
        table[r.source]=r
    used=set();result=[]
    for p in pairs:
        if type(p) is not Pair:raise ValueError('exact pair type')
        integer(p.positive,1);integer(p.negative,1)
        check(p.positive in table and p.negative in table,'pair endpoints retained')
        x,y=table[p.positive],table[p.negative]
        check(x.source!=y.source and min(x.source,y.source)>=4*N,'distinct new sources')
        check(x.source not in used and y.source not in used,'no reused source')
        used.update((x.source,y.source))
        v=signed(x);check(v!=(0,0) and v==axis(v) and signed(y)==neg(v),'clocked opposite pair')
        check(x.terminal==y.terminal,'actual terminal retained')
        result.append((max(x.depth,y.depth),2*norm(v)))
    if canonical:check(pairs==make_pairs(N,rows),'canonical arrival-rank matching')
    return tuple(result)

def exposed(p,table,d):
    integer(d)
    x,y=table[p.positive],table[p.negative]
    return add(signed(x) if x.depth<=d else (0,0),signed(y) if y.depth<=d else (0,0))

def aligned_block_hostile(pairs,table):
    """Find three disjoint real zero pairs with overlapping equal exposures."""
    groups=defaultdict(list)
    for p in pairs:
        x,y=table[p.positive],table[p.negative]
        lo,hi=sorted((x.depth,y.depth))
        if lo<hi:
            v=exposed(p,table,lo)
            groups[x.terminal,v].append((lo,hi,p.positive,p.negative))
    for (u,v),intervals in sorted(groups.items()):
        active=[]
        for lo,hi,p,m in sorted(intervals):
            while active and active[0][0]<=lo:heapq.heappop(active)
            heapq.heappush(active,(hi,lo,p,m))
            if len(active)<3:continue
            chosen=heapq.nsmallest(3,active);selected=[Pair(p,m) for _,_,p,m in chosen]
            check(all(exposed(p,table,lo)==v for p in selected),'three actual aligned exposures')
            members=[table[s] for p in selected for s in (p.positive,p.negative)]
            check(len({r.source for r in members})==6,'six disjoint source owners')
            check(all(r.terminal==u for r in members),'one literal terminal')
            total=(0,0)
            for r in members:total=add(total,signed(r))
            check(total==(0,0),'complete six-source packet zero')
            arrived=(0,0)
            for r in members:
                if r.depth<=lo:arrived=add(arrived,signed(r))
            budget=sum(norm(r.unsigned) for r in members)
            check(norm(arrived)==9*norm(v) and budget==6*norm(v),'sharp temporal aligned-petal failure')
            return {'terminal':u,'horizon':lo,'source_rows':[payload(r) for r in members],
                    'pair_source_ids':[[p.positive,p.negative] for p in selected],
                    'exposed_each':list(v),'full_packet_norm_budget':budget,
                    'arrived_packet_energy':norm(arrived),'deficit':norm(arrived)-budget}
    return None

def rejects(fn):
    try:fn()
    except (ValueError,TypeError):return True
    return False

def enumerate_rows(N,routes):
    rows=[]
    for s in range(1,256*N,2):
        if s%8==5:continue
        u,t,d=routes.get(s)
        r=Row(s,s%15,s%35,u,t,d,native.charge(64*N,s))
        audit_row(N,r,routes)
        check(r.unsigned==native.independent_charge(64*N,s),'independent charge coordinates')
        if s<4*N:check(r.unsigned==native.charge(N,s),'fixed-source scale transport')
        rows.append(r)
    return tuple(rows)

def build(N,routes):
    rows=enumerate_rows(N,routes);table={r.source:r for r in rows}
    pairs=make_pairs(N,rows);events=audit_pairs(N,rows,pairs)
    historical=make_pairs(N,rows,False);hevents=audit_pairs(N,rows,historical,False)
    maturity=Counter();old_maturity=Counter()
    for d,q in events:maturity[d]+=q
    for d,q in hevents:old_maturity[d]+=q
    batches=defaultdict(list)
    for r in rows:batches[r.depth].append(r)
    W=defaultdict(lambda:(0,0));H=defaultdict(lambda:(0,0));counts=Counter()
    arrived=credit=old_credit=0;total=sum(norm(r.unsigned) for r in rows if r.source>=4*N)
    check(total==153*N,'complete outer-annulus source budget')
    ledger=[];strict=[]
    for d in range(max(batches)+1):
        for r in batches[d]:
            v=signed(r)
            if r.source<4*N:W[r.terminal]=add(W[r.terminal],v)
            else:
                H[r.terminal]=add(H[r.terminal],v);arrived+=norm(v)
                if v!=(0,0):counts[r.terminal,v]+=1
        credit+=maturity[d];old_credit+=old_maturity[d]
        direct=0;remaining_norm=0;rebuilt=defaultdict(lambda:(0,0))
        for (u,v),c in counts.items():
            a=axis(v)
            if v==a:
                other=counts[u,neg(a)];direct+=2*norm(a)*min(c,other)
            diff=max(c-counts[u,neg(v)],0)
            remaining_norm+=norm(v)*diff
            rebuilt[u]=add(rebuilt[u],(diff*v[0],diff*v[1]))
        check(credit==direct and credit>=old_credit,'all-cut maximum mature cancellation')
        check(remaining_norm+credit==arrived,'residual source-norm decomposition')
        check(all(rebuilt[u]==v for u,v in H.items()),'signed residual vector unchanged')
        energy=sum(norm(v) for v in H.values());mixed=sum(beta(W[u],v) for u,v in H.items())
        residual_coherence=energy-remaining_norm+mixed
        interaction=energy-arrived+mixed;moving=total-arrived
        check(residual_coherence-credit==interaction,'owned pair credit identity')
        check((interaction<=moving)==(residual_coherence<=moving+credit),'exact residual payment criterion')
        ledger.append([d,arrived,moving,credit,remaining_norm,energy,mixed,residual_coherence,interaction,total-energy-mixed,old_credit])
        if credit>old_credit:strict.append([d,credit-old_credit])
    first_examples=[]
    for p in pairs:
        x,y=table[p.positive],table[p.negative]
        if x.depth!=y.depth and len(first_examples)<6:
            lo,hi=sorted((x.depth,y.depth))
            check(exposed(p,table,lo)!=(0,0) and exposed(p,table,hi)==(0,0),'actual temporal exposed interval')
            first_examples.append({'rows':[payload(x),payload(y)],'birth':lo,'maturity':hi,'exposed':list(exposed(p,table,lo))})
    return {'N':N,'sources':len(rows),'nonzero_new_sources':sum(r.source>=4*N and norm(r.unsigned)>0 for r in rows),
            'source_rows_sha256':sha256(b''.join(encode(payload(r)) for r in rows)).hexdigest(),
            'pairs':len(pairs),'pair_rows_sha256':sha256(b''.join(encode([p.positive,p.negative]) for p in pairs)).hexdigest(),
            'total_paired_source_norm':credit,'simultaneous_pairs':sum(table[p.positive].depth==table[p.negative].depth for p in pairs),
            'ledger':ledger,'strictly_better_horizons':strict,'interval_examples':first_examples,
            'aligned_three_pair_hostile':aligned_block_hostile(pairs,table),
            'all_horizons_paid':all(r[9]>=0 for r in ledger)},rows,pairs

def small_controls():
    # Independent exhaustive arbitrary arrival schedules and all pairings.
    cases=0
    for m in range(1,5):
        schedules=tuple(combinations_with_replacement(range(4),m))
        for p in schedules:
            for q in schedules:
                canonical=[sum(max(p[i],q[i])<=d for i in range(m)) for d in range(4)]
                check(canonical==[min(sum(x<=d for x in p),sum(x<=d for x in q)) for d in range(4)],'rank matching exact min formula')
                for perm in permutations(range(m)):
                    check(all(sum(max(p[i],q[perm[i]])<=d for i in range(m))<=canonical[d] for d in range(4)),'all pairing optimum at every cut')
                    cases+=1
    # Source-order permutation can hide already available cancellation.
    p=(0,2);q=(2,0)
    check(sum(max(p[i],q[i])<=0 for i in range(2))==0,'arbitrary pairing hostile')
    check(min(sum(x<=0 for x in p),sum(x<=0 for x in q))==1,'chronological repair')
    # A shared old core is not a license to discard petal-petal interactions.
    v=(1,0);core=(-1,0)
    check(norm(v)+beta(core,v)==-1,'every single petal separately profitable')
    check(norm((4,0))+beta(core,(4,0))==8 and 8>4*norm(v),'four-petal shared-core failure')
    for a in range(6):
        for r in range(13):
            for v in ((1,0),(2,1),(0,1)):
                x=(-a*v[0],-a*v[1]);y=(r*v[0],r*v[1])
                cost=norm(y)+beta(x,y)
                check(cost-r*norm(v)==norm(v)*r*(r-2*a-1),'exact collinear petal capacity')
                check((cost<=r*norm(v))==(r==0 or r<=2*a+1),'sharp finite petal guard')
    triple=((1,0),(0,1),(-1,-1))
    check(add(add(*triple[:2]),triple[2])==(0,0) and all(neg(v) not in triple for v in triple),'opposite packing not all zero-sum packing')
    return cases

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=json.loads((BASE/'frontier.json').read_text());routes=native.Routes(data['frontier'],dict(data['frontier_ranks']))
    cases=small_controls();outputs=[]
    for N in (7,4096):
        output,rows,pairs=build(N,routes);outputs.append(output)
        r=next(r for r in rows if r.clock>0)
        for bad in (replace(r,clock=r.clock+1),replace(r,depth=r.depth+1),replace(r,mod15=(r.mod15+1)%15),replace(r,source=float(r.source))):
            check(rejects(lambda bad=bad:audit_row(N,bad,routes)),'tampered immutable source row rejected')
        check(rejects(lambda:audit_pairs(N,rows,pairs+(pairs[0],))), 'duplicate source spending rejected')
        check(rejects(lambda:audit_pairs(N,rows,(Pair(True,pairs[0].negative),))), 'boolean pair source rejected')
        check(rejects(lambda:make_pairs(N,rows+(rows[0],))), 'duplicate physical source rejected')
    check(outputs[0]['ledger'][8][8:10]==[30,-13],'N7 exact failure retained')
    check(outputs[1]['all_horizons_paid'],'N4096 every horizon positive control')
    check(outputs[1]['aligned_three_pair_hostile'] is not None,'real native zero-block transient obstruction')
    dag=routes.audit()
    out={'status':'PROVED temporal pair normal form and sharp shared-core guard; FINITE-EXACT two inherited scales; unbounded renewal OPEN',
         'source_schema':['s','mod15','mod35','terminal','clock','depth','unsigned_a','unsigned_b','norm'],
         'ledger_schema':['depth','arrived_new_norm','moving_new_norm','mature_pair_norm_credit','residual_source_norm','annulus_energy','old_new_cross','residual_coherence','owned_interaction','remaining_full_bank','source_order_pair_credit'],
         'universe':{'N':[7,4096],'sources':'all primitive positive s<256N; old iff s<4N','route_cap':native.ROUTE_CAP,'arrival_depth_cap':native.DEPTH_CAP,'all_later_horizons':'constant after authenticated complete arrival'},
         'exhaustive_pairing_cases':cases,'records':outputs,'route_dag':dag,
         'frontier_lf_sha256':sha256((BASE/'frontier.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),'checks':CHECKS}
    raw=encode(out);path=HERE/'receipt.json'
    if args.verify:
        if path.read_bytes().replace(b'\r\n',b'\n')!=raw:raise ValueError('saved receipt mismatch')
    else:path.write_bytes(raw)
    print(json.dumps({'checks':CHECKS,'receipt_lf_sha256':sha256(raw).hexdigest(),'pairs':[(r['N'],r['pairs'],r['total_paired_source_norm']) for r in outputs]},sort_keys=True))

if __name__=='__main__':main()
