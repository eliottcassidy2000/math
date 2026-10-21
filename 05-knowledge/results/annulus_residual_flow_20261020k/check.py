#!/usr/bin/env python3
"""Source-owned interval capacity, exact primal/dual receipts, no uniform claim."""
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from itertools import product
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('residual_previous',HERE.parent/'annulus_block_completion_20261020j/check.py')
old=module_from_spec(spec);sys.modules[spec.name]=old;spec.loader.exec_module(old)
native=old.native
CHECKS=0

def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def integer(n,lo=0):
    if type(n) is not int or n<lo:raise ValueError('exact integer domain')
    return n

def encode(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()

@dataclass(frozen=True)
class Interval:
    identity:int
    start:int
    end:int

def typed_interval(x,D=None):
    if D is not None:integer(D,1)
    if type(x) is not Interval:raise ValueError('exact interval object')
    for n in (x.identity,x.start,x.end):integer(n)
    if x.start>=x.end or (D is not None and x.end>D):raise ValueError('nonempty half-open interval')
    return x

def interval_bank(xs,D=None):
    if D is not None:integer(D,1)
    if type(xs) is not tuple:raise ValueError('immutable interval bank')
    seen=set()
    for x in xs:
        typed_interval(x,D)
        if x.identity in seen:raise ValueError('duplicate interval source address')
        seen.add(x.identity)
    return {x.identity:x for x in xs}

def capacity(xs,D):
    integer(D,1);interval_bank(xs,D);delta=[0]*(D+1)
    for x in xs:delta[x.start]+=1;delta[x.end]-=1
    ans=[];n=0
    for d in range(D):n+=delta[d];ans.append(n)
    return tuple(ans)

def typed_problem(jobs,b):
    if type(b) is not tuple or not b:raise ValueError('nonempty immutable horizon capacity')
    for n in b:integer(n)
    return interval_bank(jobs,len(b))

def propose(jobs,b):
    """Deterministic proposal only; optimality is accepted via a separate dual."""
    table=typed_problem(jobs,b);starts=defaultdict(set);ends=defaultdict(set)
    for x in jobs:starts[x.start].add(x.identity);ends[x.end].add(x.identity)
    chosen=set(table);active=set()
    for d,cap in enumerate(b):
        active.difference_update(ends[d]);active.update(starts[d])
        while len(active)>cap:
            i=max(active,key=lambda i:(table[i].end,i));active.remove(i);chosen.remove(i)
    return tuple(sorted(chosen))

def dual(jobs,b):
    """Exact maximizer of hit(D)-sum b(d), retaining its horizon set."""
    typed_problem(jobs,b);score=[];paths=[]
    key=lambda x:(x[0],-len(x[1]),tuple(-d for d in x[1]))
    for d,cap in enumerate(b):
        hist=[0]*(d+1)
        for x in jobs:
            if x.start<=d<x.end:hist[x.start]+=1
        active=sum(hist);candidates=[(active-cap,(d,))];prefix=0
        for h in range(d):
            prefix+=hist[h]
            candidates.append((score[h]+active-prefix-cap,paths[h]+(d,)))
        best=max(candidates,key=key);score.append(best[0]);paths.append(best[1])
    return max([(0,())]+list(zip(score,paths)),key=key)

def audit_solution(jobs,b,selected,cuts):
    table=typed_problem(jobs,b)
    if type(selected) is not tuple or type(cuts) is not tuple:raise ValueError('immutable source/dual lists')
    for i in selected:integer(i)
    for d in cuts:integer(d)
    ck(tuple(sorted(set(selected)))==selected and all(i in table for i in selected),'distinct original selected jobs')
    ck(tuple(sorted(set(cuts)))==cuts and all(d<len(b) for d in cuts),'distinct actual dual horizons')
    load=capacity(tuple(table[i] for i in selected),len(b))
    ck(all(x<=c for x,c in zip(load,b)),'complete horizon feasibility')
    hit=tuple(x.identity for x in jobs if any(x.start<=d<x.end for d in cuts))
    lower=len(hit)-sum(b[d] for d in cuts);removed=len(jobs)-len(selected)
    ck(lower==removed,'weak-dual lower bound meets actual deletion count')
    # This equality proves optimality without trusting either discovery algorithm.
    retained=set(selected)
    return {'selected':selected,'removed':tuple(i for i in sorted(table) if i not in retained),
            'cuts':cuts,'cut_capacities':[b[d] for d in cuts],'hit_ids':hit,
            'minimum_deletions':removed,'load':load}

def allocate(selected,negative,D):
    """Persistent horizon matching; a positive interval may migrate at a seam."""
    ps=interval_bank(selected,D);ns=interval_bank(negative,D)
    sp=defaultdict(set);ep=defaultdict(set);sn=defaultdict(set);en=defaultdict(set)
    for x in selected:sp[x.start].add(x.identity);ep[x.end].add(x.identity)
    for x in negative:sn[x.start].add(x.identity);en[x.end].add(x.identity)
    activep=set();activen=set();assigned={};segments=[]
    for d in range(D+1):
        for p,(n,start) in list(assigned.items()):
            if ps[p].end==d or ns[n].end==d:
                segments.append((p,n,start,d));del assigned[p]
        activep.difference_update(ep[d]);activen.difference_update(en[d])
        activep.update(sp[d]);activen.update(sn[d])
        free=sorted(activen-{n for n,a in assigned.values()})
        missing=sorted(activep-set(assigned))
        ck(len(missing)<=len(free),'legal temporal capacity assignment')
        for p,n in zip(missing,free):assigned[p]=(n,d)
    ck(not assigned,'finite allocation closes')
    return tuple(sorted(segments))

def audit_allocation(selected,negative,D,segments):
    ps=interval_bank(selected,D);ns=interval_bank(negative,D)
    if type(segments) is not tuple:raise ValueError('immutable source-owned segments')
    byp=defaultdict(list);byn=defaultdict(list)
    for s in segments:
        if type(s) is not tuple or len(s)!=4:raise ValueError('exact ownership segment')
        for n in s:integer(n)
        p,n,a,b=s;ck(p in ps and n in ns,'actual positive and negative owners')
        ck(ps[p].start<=a<b<=ps[p].end and ns[n].start<=a<b<=ns[n].end,'contained temporal ownership')
        byp[p].append((a,b,n));byn[n].append((a,b,p))
    ck(set(byp)==set(ps),'all selected positive intervals covered')
    for p,parts in byp.items():
        parts.sort();ck(parts[0][0]==ps[p].start and parts[-1][1]==ps[p].end,'full positive endpoints')
        ck(all(a[1]==b[0] for a,b in zip(parts,parts[1:])),'positive coverage without gaps or duplication')
    for n,parts in byn.items():
        parts.sort();ck(all(a[1]<=b[0] for a,b in zip(parts,parts[1:])),'no simultaneous negative ownership reuse')
    return max((len(v) for v in byp.values()),default=0)

def reject(f):
    try:f()
    except (ValueError,TypeError):return True
    return False

def small_controls():
    universe=tuple((a,b) for a in range(3) for b in range(a+1,4));cases=0
    for mask in range(64):
        jobs=tuple(Interval(i,*universe[i]) for i in range(6) if mask>>i&1)
        for b in product(range(3),repeat=3):
            selected=propose(jobs,b);lower,cuts=dual(jobs,b);audit_solution(jobs,b,selected,cuts)
            best=0
            for bits in range(1<<len(jobs)):
                sub=tuple(x for j,x in enumerate(jobs) if bits>>j&1)
                if all(x<=c for x,c in zip(capacity(sub,3),b)):best=max(best,len(sub))
            ck(len(selected)==best and len(jobs)-best==lower,'independent complete small subset optimum')
            brute=max(sum(any(x.start<=d<x.end for d in range(3) if bits>>d&1) for x in jobs)
                      -sum(b[d] for d in range(3) if bits>>d&1) for bits in range(8))
            ck(brute==lower,'independent all dual horizon subsets');cases+=1
    p=(Interval(10,0,2),);n=(Interval(20,0,1),Interval(21,1,2))
    flow=allocate(p,n,2);ck(audit_allocation(p,n,2,flow)==2,'handoff defeats whole-container restriction')
    ck(all(not(x.start<=0 and 2<=x.end) for x in n),'no whole-container certificate exists')
    gap=(Interval(30,0,1),Interval(31,0,1),Interval(32,2,3),Interval(33,2,3))
    gcap=capacity(gap,3);gjobs=(Interval(40,0,3),)
    ck(sum(gcap)>3 and propose(gjobs,gcap)==(),'total capacity area does not fill a missing horizon')
    gp=propose(gjobs,gcap);_,gd=dual(gjobs,gcap);audit_solution(gjobs,gcap,gp,gd)
    duplicate=(Interval(0,0,2),Interval(1,0,2));audit_solution(duplicate,(1,1),propose(duplicate,(1,1)),(0,))
    ck(reject(lambda:typed_interval(Interval(True,0,1))),'bool source ID rejected')
    ck(reject(lambda:typed_problem(p,(1.0,1))),'float capacity rejected')
    ck(reject(lambda:typed_interval(Interval(0,0,1),True)),'bool horizon rejected')
    ck(reject(lambda:interval_bank((),2.0)),'float empty-bank horizon rejected')
    ck(reject(lambda:audit_allocation(p,n,2.0,flow)),'float allocation horizon rejected')
    ck(reject(lambda:audit_allocation((),(),True,())),'bool empty allocation horizon rejected')
    ck(reject(lambda:audit_allocation(p,n,2,flow+(flow[0],))),'duplicated temporal ownership rejected')
    ck(reject(lambda:audit_allocation(p,n,2,((10,20,0,2),))),'out-of-owner segment rejected')
    return {'complete_small_cases':cases,'source_intervals':universe,'handoff_segments':flow,
            'gap_capacity':gcap,'gap_dual':gd}

def source_bank(N,routes):
    """Same inherited canonical runs, now with the exact x=1 mod8 start for all N."""
    integer(N,1);blocks=[];rows={};used=set()
    for band in range(5):
        lo,hi=(4<<band)*N,(8<<band)*N;runs=[];run=[]
        for x in range(lo+(1-lo)%8,hi,8):
            if x in routes.V or 2*x+1 in routes.V:
                if run:runs.append(run);run=[]
            else:run.append(x)
        if run:runs.append(run)
        for run in runs:
            for j in range(0,len(run)-2,3):
                terms=[];ids=[]
                for x in run[j:j+3]:
                    y=2*x+1
                    for s in (x,y):
                        u,t,d=routes.get(s);a=native.charge(64*N,s)
                        row=old.owned.Row(s,s%15,s%35,u,t,d,a);old.owned.audit_row(N,row,routes)
                        ck(a==native.independent_charge(64*N,s),'independent source charge')
                        rows[s]=old.owned.payload(row)
                    p,q=rows[x],rows[y]
                    ck(p[3]==q[3] and p[5]==q[5] and q[4]==p[4]+1,'literal native pair terminal and clocks')
                    w=-old.beta(tuple(p[6:8]),tuple(q[6:8]));terms.append((w,p[5],x,y));ids.extend((x,y))
                ck(not used.intersection(ids),'disjoint source owners in canonical bank');used.update(ids)
                pos=[d for w,d,x,y in terms if w==3];neg=[d for w,d,x,y in terms if w==-3]
                if not pos and not neg:sign=start=end=0
                else:
                    ck(len(pos)==len(neg)==1,'native interval weight alphabet');start,end=sorted((pos[0],neg[0]))
                    sign=0 if start==end else (1 if pos[0]<neg[0] else -1)
                blocks.append({'id':len(blocks),'band':band,'terms':tuple(terms),'sign':sign,'start':start,'end':end})
    # Reuse the independent native compiler/auditor, not a guessed residue graph.
    for b in blocks:old.audit_block(N,b,routes)
    return tuple(blocks),rows

def solve_native(N,routes,prior=None):
    blocks,rows=source_bank(N,routes);D=native.DEPTH_CAP+1
    positive=tuple(Interval(b['id'],b['start'],b['end']) for b in blocks if b['sign']==1)
    negative=tuple(Interval(b['id'],b['start'],b['end']) for b in blocks if b['sign']==-1)
    b=capacity(negative,D);fixed=()
    if prior is not None:
        ck(sha256(encode(blocks)).hexdigest()==prior['interval_completion']['block_sha256'],'unchanged native bank source')
        fixed=tuple(sorted({i for n,ps in prior['interval_completion']['matches'] for i in ps}))
        ck(len(fixed)==6088,'previous completed positive owners preserved')
    table=interval_bank(positive,D);fixed_set=set(fixed)
    residual=tuple(x for x in positive if x.identity not in fixed_set)
    fixedload=capacity(tuple(table[i] for i in fixed),D)
    leftcap=tuple(c-f for c,f in zip(b,fixedload));ck(all(c>=0 for c in leftcap),'fixed prior proof remains feasible')
    proposed=propose(residual,leftcap);number,cuts=dual(residual,leftcap)
    nested=audit_solution(residual,leftcap,proposed,cuts)
    selected=tuple(sorted(fixed+proposed));global_receipt=audit_solution(positive,b,selected,cuts)
    ck(number==global_receipt['minimum_deletions'],'nested refinement is globally optimal')
    chosen=tuple(table[i] for i in selected);segments=allocate(chosen,negative,D)
    migrations=audit_allocation(chosen,negative,D,segments)
    # Check the selected term sum directly from the original signed source interactions.
    selected_load=capacity(chosen,D);total_load=capacity(positive,D)
    removed_set=set(global_receipt['removed']);remaining=tuple(x for x in positive if x.identity in removed_set)
    residual_load=capacity(remaining,D);bank_debit=tuple(3*(p-c) for p,c in zip(total_load,b))
    upper=tuple(3*x for x in residual_load)
    cutwise=tuple(max(j,0) for j in bank_debit)
    ck(all(c<=u for c,u in zip(cutwise,upper)),'marked-segment guard dominates whole-interval guard')
    # Canonical per-cut matching is source-indexed. Hash all actual owner pairs;
    # a matched piece alone does not certify its entire positive interval.
    cut_match_hash=sha256()
    for d in range(D):
        pp=sorted(x.identity for x in positive if x.start<=d<x.end)
        nn=sorted(x.identity for x in negative if x.start<=d<x.end)
        matched=tuple(zip(pp,nn));residual=pp[len(matched):]
        ck(3*len(residual)==cutwise[d],'actual cutwise pair-owner matching')
        cut_match_hash.update(encode([d,matched,residual]))
    ck(all(j<=u for j,u in zip(bank_debit,upper)),'all-horizon retained original interaction upper bound')
    grouped=tuple(3*(p-c) for p,c in zip(selected_load,b))
    ck(max(grouped)<=0,'joint completion retires only nonpositive selected debit')
    byp=defaultdict(list)
    for p,n,a,z in segments:byp[p].append((p,n,a,z))
    example_id=max(byp,key=lambda i:(len(byp[i]),-i));example=byp[example_id]
    involved={example_id,*[n for p,n,a,z in example]}
    source_ids=sorted({s for i in involved for w,d,x,y in blocks[i]['terms'] for s in (x,y)})
    removed_sources=sorted({s for i in global_receipt['removed'] for w,d,x,y in blocks[i]['terms'] for s in (x,y)})
    result={'N':N,'block_count':len(blocks),'source_row_count':len(rows),
            'block_sha256':sha256(encode(blocks)).hexdigest(),
            'source_rows_sha256':sha256(encode([rows[s] for s in sorted(rows)])).hexdigest(),
            'positive_count':len(positive),'negative_count':len(negative),'previous_fixed':fixed,
            'extra_selected_count':len(proposed),'nested_certificate':nested,'global_certificate':global_receipt,
            'negative_capacity':b,'joint_debit':grouped,'original_bank_debit':bank_debit,
            'residual_debit':upper,'flow_segments':segments,'max_owners_for_one_positive':migrations,
            'marked_segment_residual_debit':cutwise,'cutwise_matching_sha256':cut_match_hash.hexdigest(),
            'marked_segment_peak':max(cutwise),'marked_segment_peak_horizons':[d for d,z in enumerate(cutwise) if z==max(cutwise)],
            'removed_blocks':[blocks[i] for i in global_receipt['removed']],
            'removed_source_rows':[rows[s] for s in removed_sources],
            'handoff_example':{'positive_id':example_id,'segments':example,
                               'blocks':[blocks[i] for i in sorted(involved)],'source_rows':[rows[s] for s in source_ids]},
            'cut_sources':[[d,[x.identity for x in positive if x.start<=d<x.end],b[d]] for d in cuts]}
    if prior is not None:
        oldupper=prior['interval_completion']['unmatched_positive_debit'][:D]
        gain=[a-z for a,z in zip(oldupper,upper)];ck(all(g>=0 for g in gain),'pointwise monotone refinement of prior guard')
        prev=prior['full_native_control']['full_ledger'][:D]
        absolute=[r[6]+u for r,u in zip(prev,upper)]
        cut_absolute=[r[6]+u for r,u in zip(prev,cutwise)]
        ck(all(r[5]<=a for r,a in zip(prev,absolute)),'original physical I and R retained')
        ck(all(r[5]<=c<=a for r,c,a in zip(prev,cut_absolute,absolute)),'marked segment bound retains full actual ledger')
        result.update({'prior_residual_debit':oldupper,'pointwise_gain':gain,'maximum_gain':max(gain),
                       'gain_horizons':[d for d,g in enumerate(gain) if g==max(gain)],
                       'max_residual_debit':max(upper),'residual_peak_horizons':[d for d,g in enumerate(upper) if g==max(upper)],
                       'actual_R_plus_residual':absolute,'actual_R_plus_marked_segments':cut_absolute,
                       'maximum_excess_over_moving':max(a-r[2] for a,r in zip(absolute,prev))})
    return result

def run():
    controls=small_controls()
    f=json.loads((HERE.parent/'session18_native_annulus/frontier.json').read_text())
    routes=native.Routes(f['frontier'],dict(f['frontier_ranks']))
    prior_path=old.HERE/'receipt.json';prior_bytes=prior_path.read_bytes().replace(b'\r\n',b'\n');prior=json.loads(prior_bytes)
    ck(sha256(prior_bytes).hexdigest()=='b85cfd8a95962c827d9e21d1418a8fe7f4086c6cca6203f042f547d32ad3b609',
       'frozen original source ledger')
    cases=[solve_native(4096,routes,prior),solve_native(8192,routes)]
    ck(cases[0]['global_certificate']['minimum_deletions']==94 and cases[1]['global_certificate']['minimum_deletions']==159,
       'two exact finite optima; no uniform residual constant inferred')
    return {'scope':'source-owned signed pair terms only; all other interactions remain in R; uniform153N OPEN',
            'controls':controls,'native_cases':cases,'route_audit':routes.audit(),'checks':CHECKS}

if __name__=='__main__':
    result=run();data=(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n').encode();path=HERE/'receipt.json'
    if '--write' in sys.argv:path.write_bytes(data)
    else:ck(path.read_bytes().replace(b'\r\n',b'\n')==data,'saved exact receipt')
    print(json.dumps({'checks':CHECKS,'sha256':sha256(data).hexdigest(),
                      'cases':[[c['N'],c['positive_count'],len(c['global_certificate']['selected']),
                                c['global_certificate']['minimum_deletions'],c['global_certificate']['cuts']]
                               for c in result['native_cases']]}))
