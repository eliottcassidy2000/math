"""Bounded inverse-selector replay and independent full-count target controls."""
from pathlib import Path
from math import isqrt
import importlib.util
import hashlib
import json
import sys
import check as core

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
INVERSE=ROOT/'05-knowledge/results/changing_core_returns_20261113/returns/inverse.py'
spec=importlib.util.spec_from_file_location('inherited_tail_inverse',INVERSE)
inv=importlib.util.module_from_spec(spec);spec.loader.exec_module(inv)

def induced(rows,keep):
    return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(keep)) for v in keep)

def direct_start(rows,start):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)];dp[1<<start][start]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            if not value:continue
            for w in range(n):
                if not(mask>>w&1) and rows[v]>>w&1:dp[mask|(1<<w)][w]+=value
    return sum(dp[-1])

def one_apex(rows,apex,force_first_spine=False):
    rest=[v for v in range(len(rows)) if v!=apex]
    spine=sorted(rest,key=lambda v:-sum((rows[v]>>w)&1 for w in rest))
    core.need(all((rows[spine[i]]>>spine[j])&1 for i in range(len(spine)) for j in range(i+1,len(spine))),
              'literal transitive spine')
    bits=[(rows[apex]>>v)&1 for v in spine]
    state={(0,0):1}
    for index,bit in enumerate(bits):
        new={}
        for (left,right),value in state.items():
            key=(1 if bit else 2,right)
            new[key]=new.get(key,0)+value
            if not(force_first_spine and index==0):
                key=(left,right if right else (1 if bit else 2))
                new[key]=new.get(key,0)+value
        state=new
    H=sum(v for (left,right),v in state.items() if left in (0,2) and right in (0,1))
    if not force_first_spine:
        inversion=1+sum(2 if j-i==1 else 1<<(j-i-1)
                        for i in range(len(bits)) for j in range(i+1,len(bits)) if bits[i] and not bits[j])
        core.need(H==inversion,'independent two-slot path assignment versus weighted inversions')
    return H,spine,''.join(map(str,bits))

def main():
    candidates=[]
    for r in range(1,7):
        source=15*32**r-5**r
        for digit in (2,4):
            n=3*source+digit;f=core.factor(n)
            candidates.append({'r':r,'source':source,'digit':digit,'target':n,'factorization':f,
                               'prime':f==[n],'seven_prime':len(f)==2 and f[0]==7})
    core.need(core.factor(47184047)==[47184047],'critical r4 target is prime by complete trial division')
    bank=inv.cores(10)
    records=[]
    for target in (1427,1429,47184047,317):
        found=None;nodes=0
        for index,row in enumerate(bank):
            h,e,s,p,ee,t=row['state'];G=s+p
            addresses,meta=inv.inverse(target,G,h,e,G-e-ee-1);nodes+=meta['nodes']
            if row['order']==1:addresses=[x for x in addresses if x]
            if not addresses:continue
            address=min(addresses,key=lambda x:(x.bit_length(),x))
            T,z,V,x,word=inv.physical(row,address)
            state=inv.f.advance(row['state'],word)
            core.need(state[0]+state[4]-1==target,'independent six-state observer')
            core.need(core.valid(T) and core.strong(T),'literal strong target')
            predecessors=[v for v in range(len(T)) if (T[v]>>z)&1]
            core.need(predecessors==[x],'actual indegree-one root')
            if row['order']==1:
                actual,spine,apex_word=one_apex(T,x)
                core.need(spine[0]==z,'new root is first transitive-spine vertex')
                starts,_,_=one_apex(T,x,True)
                deletion,_,_=one_apex(V,x)
                verification='Independent actual one-apex two-slot path DP and inversion formula'
            else:
                core.need(len(T)<=18,'bounded independent subset control')
                actual=core.hp(T);starts=direct_start(T,z);deletion=core.hp(V)
                spine=apex_word=None
                verification='Independent full subset path DP, including rooted starts and deletion'
            core.need(actual==target and deletion-starts==1,'full target and unit-ready defect')
            found={'target':target,'core_order':row['order'],'core_rank':row['rank'],'core_words':row['words'],
                   'tail_address':address,'tail_word':''.join(map(str,word)),'order':len(T),'rows':T,
                   'root':z,'sole_predecessor':x,'H_starting_at_root':starts,'H_after_root_deletion':deletion,
                   'ready_defect':deletion-starts,'parent_six_state':state,'cores_tried':index+1,
                   'inverse_nodes':nodes,'verification':verification,'transitive_spine':spine,'apex_word':apex_word}
            break
        if found is None:
            found={'target':target,'status':'NO HIT IN DECLARED BANK','core_order_cap':10,
                   'distinct_core_observers':len(bank),'inverse_nodes':nodes,
                   'scope':'Complete inherited arbitrary-tail inverse for each listed core; not a global absence claim'}
        records.append(found)
    core.need(len(bank)==1596 and records[-1]['target']==317 and 'rows' not in records[-1],'bounded hostile scope')
    result={'status':'PASS','Q_power_unit_target_factors':candidates,'source_replacements':records,
            'prime47184047_trial_division_bound':isqrt(47184047),
            'inverse_source_pin':{'path':str(INVERSE.relative_to(ROOT)).replace('\\','/'),
                                  'sha256_lf':hashlib.sha256(INVERSE.read_bytes().replace(b'\r\n',b'\n')).hexdigest()},
            'scope':'Successful target witnesses with ready roots; no retained-Q-module repair or source-to-target transition is claimed'}
    suffix='.optimized' if sys.flags.optimize else ''
    (HERE/('target_receipt'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','prime_r4_target':47184047,'witnesses':[(r['target'],r.get('order'),r.get('ready_defect')) for r in records],
                      'bounded_bank_size':len(bank)}))


if __name__=='__main__':main()
