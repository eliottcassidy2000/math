"""Independent physical three-return audit; no imported path-count implementation."""
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
REPAIRS=HERE.parent/'repairs'
def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def table(rows):
    n=len(rows); full=(1<<n)-1
    for v,row in enumerate(rows):
        need(0<=row<=full and not(row>>v&1),'vertex row range/loop')
        for u in range(v):
            need(((row>>u)&1)+((rows[u]>>v)&1)==1,'tournament orientation')
    incoming=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
    e=[[0]*n for _ in range(full+1)]; h=[1]+[0]*full
    for mask in range(1,full+1):
        verts=mask
        while verts:
            bit=verts&-verts; verts-=bit; v=bit.bit_length()-1
            rest=mask^bit
            if not rest: e[mask][v]=1; continue
            prior=rest&incoming[v]
            while prior:
                b=prior&-prior; prior-=b
                e[mask][v]+=e[rest][b.bit_length()-1]
        h[mask]=sum(e[mask])
    return e,h

def strong(rows):
    full=(1<<len(rows))-1
    for r in range(len(rows)):
        reached=1<<r
        while True:
            newer=reached
            for v,row in enumerate(rows):
                if reached>>v&1: newer|=row
            if newer==reached: break
            reached=newer
        if reached!=full: return False
    return True

def state(rows,r):
    full=(1<<len(rows))-1
    e,h=table(rows)
    s,_=table([full^(1<<v)^rows[v] for v in range(len(rows))])
    f=sum(e[m][r]*h[full^m] for m in range(1,full+1))
    return [h[full],s[full][r],h[full^(1<<r)],f]

def birth(rows,r):
    n=len(rows)
    return [row|((1<<n) if v==r else 0) for v,row in enumerate(rows)]+[((1<<n)-1)^(1<<r)]

def patch(rows,edges):
    need(len({tuple(sorted(e)) for e in edges})==len(edges),'distinct patch arcs')
    result=rows[:]
    for a,b in edges:
        need(0<=a<len(rows) and 0<=b<len(rows) and a!=b,'old patch vertices')
        result[a]^=1<<b; result[b]^=1<<a
    return result

def induced(rows,names):
    need(len(set(names))==len(names),'injective name map')
    return [sum(1<<j for j,v in enumerate(names) if rows[u]>>v&1) for u in names]

def first_record(name):
    return json.loads((REPAIRS/name).read_text(encoding='utf-8-sig').splitlines()[0])

normal=REPAIRS/'chain_881_1691_3311_6551.json'
optimized=REPAIRS/'chain_881_1691_3311_6551.optimized.json'
d=json.loads(normal.read_text(encoding='utf-8'))
need(d==json.loads(optimized.read_text(encoding='utf-8')),'original replay pair')
a=d['source_rows']; helper=birth(a,d['helper_port'])
need(helper==d['helper_rows'],'exact helper birth')
stages=[]
def check(label,rows,root,expected,is_strong):
    actual=state(rows,root)
    connected=strong(rows)
    need(actual==expected,(label,'state',actual,expected))
    need(connected==is_strong,(label,'strongness',connected))
    stages.append({'stage':label,'order':len(rows),'root':root,'state':actual,'strong':connected})

need(state(a,3)[0]==881 and strong(a),'literal A2 source')
check('helper',helper,len(a),[1589,698,881,1772],True)
p=patch(helper,d['first_patch']); need(p==d['first_parent_rows'],'first parent arcs')
need(len(d['first_patch'])==5,'first patch width')
raw1=first_record('helper_alt303_two_star.raw.jsonl')
need(raw1['rows']==p and raw1['changed_arcs']==d['first_patch'],'raw first repair')
check('parent351',p,0,[351,1,3,1341],True)
t1=birth(p,0); need(t1==d['first_ready_rows'],'first target arcs')
check('ready1691',t1,d['first_ready_root'],[1691,350,351,1692],True)
names=d['canonical_to_actual_names']
need(sorted(names)==list(range(len(t1))),'canonical map permutation')
need(names[0]==0 and names[-1]==d['first_ready_root'],'path map endpoints')
need(all(p[u]>>v&1 for u,v in zip(names[:-2],names[1:-1])),'named unique parent Hamilton path')
canonical=json.loads((REPAIRS/'alternate_1691_inputs.json').read_text(encoding='utf-8'))[2]['rows']
need(induced(t1,names)==canonical,'full canonical target isomorphism')
raw2=first_record('alternate303_two_star.raw.jsonl')
need(d['second_patch']==[[names[u],names[v]] for u,v in raw2['changed_arcs']],'second patch transport')
need(d['second_parent_root']==names[raw2['root']],'second port transport')
q=patch(t1,d['second_patch']); need(q==d['second_parent_rows'],'second parent arcs')
need(len(d['second_patch'])==6 and induced(q,names)==raw2['rows'],'second canonical repair')
check('parent853',q,d['second_parent_root'],[853,1,5,2459],True)
t2=birth(q,d['second_parent_root']); need(t2==d['second_ready_rows'],'second target arcs')
check('ready3311',t2,d['second_ready_root'],[3311,852,853,3312],True)
names2=names+[d['second_ready_root']]
need(sorted(names2)==list(range(len(t2))),'second canonical map permutation')
raw3=first_record('alternate6551_two_star.raw.jsonl')
need(d['third_patch']==[[names2[u],names2[v]] for u,v in raw3['changed_arcs']],'third patch transport')
need(d['third_parent_root']==names2[raw3['root']],'third port transport')
r=patch(t2,d['third_patch']); need(r==d['third_parent_rows'],'third parent arcs')
need(len(d['third_patch'])==12 and induced(r,names2)==raw3['rows'],'third canonical repair')
check('parent495',r,d['third_parent_root'],[495,1,33,6057],False)
t3=birth(r,d['third_parent_root']); need(t3==d['third_ready_rows'],'third target arcs')
check('ready6551_before_deletion',t3,d['third_ready_root'],[6551,494,495,6552],False)
deleted=d['deleted_old_name']; need(deleted==8==names2[10] and t3[deleted]==0,'named global sink')
keep=[v for v in range(len(t3)) if v!=deleted]
need(keep==d['retained_names'] and keep.index(d['third_ready_root'])==d['final_root'],'retained root map')
final=induced(t3,keep); need(final==d['final_rows'],'final deletion arcs')
check('ready6551_final',final,d['final_root'],[6551,494,495,6552],True)
for pkey,tkey in [('first_parent_state','first_ready_state'),('second_parent_state','second_ready_state'),('third_parent_state','third_ready_state')]:
    pp,tt=d[pkey],d[tkey]
    need(tt[0]==pp[0]+pp[3]-pp[1] and tt[3]-tt[0]==pp[1]==1,'birth full-target and readiness identity')
sources=[normal,optimized,REPAIRS/'chain_881_1691_3311_6551.py',REPAIRS/'alternate_1691_inputs.json',
         REPAIRS/'helper_alt303_two_star.raw.jsonl',REPAIRS/'alternate303_two_star.raw.jsonl',REPAIRS/'alternate6551_two_star.raw.jsonl']
pins=[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
out={'status':'FINITE-EXACT independent full physical audit accepted',
     'scope':'one helper and three target returns with old-arc repairs; no all-k claim',
     'stages':stages,'patch_widths':[5,6,12],'canonical_to_actual_names':names,
     'deleted_old_name':deleted,'source_pins':pins}
suffix='.optimized' if sys.flags.optimize else ''
(HERE/('THREE_STEP_AUDIT'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(stages))
