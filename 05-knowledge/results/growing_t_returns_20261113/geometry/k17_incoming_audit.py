"""Review controls for the incoming intrinsic complete k17 exclusion."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE=HERE.parent.parent/'target_directed_cores_20261112'
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def table(rows):
    n=len(rows); full=(1<<n)-1
    e=[[0]*n for _ in range(full+1)]; h=[1]+[0]*full
    for mask in range(1,full+1):
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                e[mask][v]=sum(e[rest][u] for u in range(n) if rows[u]>>v&1) if rest else 1
        h[mask]=sum(e[mask])
    return e,h
def p2table(h):
    out=[]
    for mask in range(len(h)):
        total=0; sub=mask
        while True:
            total+=h[sub]*h[mask^sub]
            if not sub: break
            sub=(sub-1)&mask
        out.append(total)
    return out
parents=roots=births=0
for n in range(1,5):
    full=(1<<n)-1
    for code in range(1<<(n*(n-1)//2)):
        rows=[0]*n
        for k,(i,j) in enumerate(combinations(range(n),2)):
            a,b=(i,j) if code>>k&1 else(j,i)
            rows[a]|=1<<b
        end,h=table(rows); p=p2table(h)
        start,_=table([full^(1<<v)^rows[v] for v in range(n)])
        need(p[full]>=1<<n,'P2 assignment lower bound'); parents+=1
        for r in range(n):
            f=sum(end[a][r]*h[full^a] for a in range(1,full+1))
            s=sum(start[a][r]*h[full^a] for a in range(1,full+1))
            t=sum(start[a][r]*p[full^a] for a in range(1,full+1))
            q0=h[full]+s+p[full]-end[full][r]-f-1
            need(2*f<=p[full] and 2*t>=p[full] and 2*q0>=p[full],'intrinsic cover bounds')
            roots+=1
        for out in range(1<<n):
            child=[row | (0 if out>>v&1 else 1<<n) for v,row in enumerate(rows)]+[out]
            ch=table(child)[1]; pp=p2table(ch)[-1]
            need(pp>=2*p[full],'every birth doubles P2'); births+=1
word_checks=0
for length in range(1,11):
    for word in range(1<<(length-1),1<<length):
        L=2**length+sum(3**i*2**(length-i-1) for i in range(length) if word>>i&1)
        D=1+sum(2**(length-i-1) for i in range(length) if not word>>i&1)
        v=[1,0,0,0,0,0]
        for i in range(length):
            bit=word>>i&1; w=[0]*6
            for old,a in enumerate(v):
                last,started=divmod(old,2)
                w[old]+=a
                w[2*(1 if bit else 2)+started]+=a
                if started or not bit: w[2*last+1]+=a
            v=w
        B=sum(v[:4])
        need(L>=2**length+3**(length-1) and B>=2**(length-1) and D<=2**length-1,'word floors')
        need(2*B-D>=1,'positive tail remainder'); word_checks+=1
producer=SOURCE/'complete_53084231.json'
auditor=SOURCE/'ports/complete_audit_53084231.json'
a=json.loads(producer.read_text()); b=json.loads(auditor.read_text())
pa={x['order']:x['retained_cores'] for x in a['layers']}
pb={x['order']:x['retained'] for x in b['layers']}
need(a['target']==b['target']==53084231 and a['core_p2_bound']==106168462 and a['core_order_bound']==26,'complete universe')
need(all(pa.get(i,0)==pb.get(i,0) for i in range(1,27)),'full core layers')
need(sum(pa.values())==b['total_cores']==5640026,'all core count')
need(not a['found'] and a['complete_core_universe'] and not b['found'] and not b['hits'],'no target')
need(b['direct_outer_tests']==74845805,'direct address universe')
sources=[SOURCE/'TARGET_INVERSE.md',SOURCE/'ports/COMPLETE_INVERSE_AUDIT.md',
         SOURCE/'complete_ready_inverse.cpp',SOURCE/'target_inverse.cpp',SOURCE/'ports/complete_audit.cpp',producer,auditor,
         HERE.parent.parent/'clone_clocks_negative_collatz_20261113/CLOCKS.md']
pins=[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
out={'status':'VERIFIED intrinsic proof/code review and FINITE-EXACT local controls; large census inherited and replayed by parent',
     'consequence':'REFUTED all-k indegree-one unit-ready405 interface; ordinary count and Chan-Pak unaffected',
     'small_controls':{'tournaments':parents,'roots':roots,'births':births,'legal_words':word_checks},
     'target':53084231,'core_p2_bound':106168462,'core_order_bound':26,
     'retained_cores':5640026,'direct_outer_tests':74845805,'found':False,'source_pins':pins}
suffix='.optimized' if sys.flags.optimize else ''
(HERE/('K17_INCOMING_AUDIT'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out['small_controls']))
