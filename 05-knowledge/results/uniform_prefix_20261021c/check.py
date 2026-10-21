"""Independent literal-source replay of a native one-step payment obstruction."""
from collections import defaultdict
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse
import json

HERE=Path(__file__).resolve().parent
COUNT=0
ZERO=(0,0,0,0)


def need(ok,msg):
    global COUNT
    COUNT+=1
    if not ok: raise ValueError(msg)


def encode(obj): return (json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()
def Q(v): return v[0]**2-v[0]*v[1]+v[1]**2
def delta(v): return Q((v[0]+v[2],v[1]+v[3]))-Q(v[:2])
def add(v,w): return tuple(a+b for a,b in zip(v,w))


def edge(u):
    a=2 if u%8==1 else 1
    z=(3*u+1)//(1<<a)
    while z%8==5: z=(z-1)//4; a+=2
    return z,a


def charge(N,s):
    a=b=0
    lower=((1,2),(-2,-1),(1,-1)); upper=((0,1),(-1,-1),(1,0))
    while s<4*N:
        x,y=(lower if s<2*N else upper)[s%3]
        a+=x; b+=y; s=4*s+1
    return a,b


def summary(values):
    total=peak=0
    for x in values: total+=x; peak=max(peak,total)
    return total,peak


def compose(a,b): return a[0]+b[0],max(a[1],a[0]+b[1])


def run():
    for n in range(7):
        for word in product((-2,0,3),repeat=n):
            for k in range(n+1):
                need(compose(summary(word[:k]),summary(word[k:]))==summary(word),'excursion composition')
    need(summary((4,-5))==( -1,4),'negative net block can require positive reserve')
    need(summary((-5,4))==(-1,0),'same net debt can have different reserve')
    raw=(HERE.parent/'session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n',b'\n')
    need(sha256(raw).hexdigest()=='391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f','frontier pin')
    V=set(json.loads(raw)['frontier'])
    N=4096; D=17
    states=[defaultdict(lambda:ZERO),defaultdict(lambda:ZERO)]
    source_hash=sha256(); local_hash=sha256(); local_counts=defaultdict(int); clock_ranges={}
    source_count=norm=0; local_witness={}
    # Direct independent source routing: no early zero coalescence or native cache.
    for s in range(1,256*N,2):
        if s%8==5: continue
        source_count+=1; c=charge(64*N,s); old=s<4*N
        if old: need(c==charge(N,s),'old named charge stable')
        else: norm+=Q(c)
        need(s%15%5==s%35%5,'joint registers compatible')
        source_hash.update(encode([s,s%15,s%35,int(old),*c]))
        if not any(c): continue
        u=s; clock=0; word=[]
        for _ in range(D):
            if u in V: break
            z,a=edge(u); word.append(a); clock+=a; u=z
        sg=1-2*(clock%2)
        v=tuple(sg*x for x in (c+(0,0) if old else (0,0)+c))
        states[0][u]=add(states[0][u],v)
        if u in (811,1081,2163):
            row=[s,s%15,s%35,int(old),*c,u,len(word),clock,word]
            local_hash.update(encode(row)); local_counts[u]+=1
            if u not in clock_ranges: clock_ranges[u]=[clock,clock]
            else: clock_ranges[u]=[min(clock_ranges[u][0],clock),max(clock_ranges[u][1],clock)]
            local_witness.setdefault(u,row)
        if u not in V:
            u,a=edge(u); v=tuple((1-2*(a%2))*x for x in v)
        states[1][u]=add(states[1][u],v)
    need(source_count==393216 and norm==153*N,'complete source universe and unchanged budget')
    vals=[]
    for state in states:
        arrived=sum(delta(v) for u,v in state.items() if u in V)
        moving=sum(max(0,delta(v)) for u,v in state.items() if u not in V)
        vals.append([arrived,moving,arrived+moving])
    need(vals==[[166652,268827,435479],[184903,254791,439694]],'direct original-source two-cut replay')
    base=states[0][811]; a=states[0][1081]; b=states[0][2163]
    need(base==(-6,-10,7,-36) and a==(-2,-2,-9,-25) and b==(0,0,-1,0),'literal local fibre charges')
    need(edge(1081)==(811,2) and edge(2163)==(811,3),'actual local child clocks')
    incoming=add(a,tuple(-x for x in b))
    increase=delta(add(base,incoming))-delta(base)-max(0,delta(a))-max(0,delta(b))
    need(increase==2004,'positive local charge cannot be paid pointwise')
    need(vals[1][2]-vals[0][2]==4215,'positive global one-step change')
    buckets=defaultdict(list)
    for u,v in states[0].items():
        if u in V or not any(v): continue
        z,c=edge(u); buckets[z].append((v,tuple((1-2*(c%2))*x for x in v)))
    debit=0
    for z,rows in buckets.items():
        total=ZERO
        for _,v in rows: total=add(total,v)
        before=sum(max(0,delta(v)) for v,_ in rows)
        d=(delta(add(states[0].get(z,ZERO),total))-delta(states[0].get(z,ZERO)) if z in V else max(0,delta(total)))-before
        debit+=d
    need(debit==4215,'grouped exact debit identity matches independent source replay')
    probe=json.loads((HERE/'probe_4096.json').read_text())
    need(probe['rows'][17][1:4]==vals[0] and probe['rows'][18][1:4]==vals[1],'discovery replay agrees')
    need(summary([probe['rows'][i+1][3]-probe['rows'][i][3] for i in range(150)])[1]==0,'declared finite producer keeps initial reserve')
    return {'status':'PROVED debit and excursion algebra; FINITE-EXACT native local obstruction; uniform block payment OPEN',
            'N':N,'depths':[17,18],'budget':153*N,'source_count':source_count,'source_hash':source_hash.hexdigest(),
            'local_source_hash':local_hash.hexdigest(),'local_source_counts':dict(local_counts),
            'local_clock_ranges':clock_ranges,'local_source_witnesses':local_witness,'states':vals,
            'local_target':811,'local_before':base,'children':[[1081,2,a],[2163,3,b]],
            'local_increase':increase,'global_increase':debit,'checks':COUNT}


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--write',action='store_true'); a=p.parse_args()
    data=encode(run()); path=HERE/'receipt.json'
    if a.write: path.write_bytes(data)
    else: need(path.read_bytes().replace(b'\r\n',b'\n')==data,'saved receipt comparison')
    print(json.dumps({'checks':COUNT,'sha256_lf':sha256(data).hexdigest()}))
