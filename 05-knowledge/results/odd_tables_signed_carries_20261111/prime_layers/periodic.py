"""Finite synchronized carry graph and an actual multiplier-seven rebase."""
from collections import deque
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import types

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.setrecursionlimit(10000)
VS=(-3,0,3,10,13,16)


def g(n):
    answer=1
    while n>=2:
        q=n//2
        if n%2==0:answer+=q+(q&1)
        n=q
    return answer


def edge(t,state,bit):
    j,u,v=state
    J,r=divmod(t*bit+j,2)
    delta=J&1
    successor=(J,2*u+13*(1-bit)-3*t*(1-r),
               13*(1-bit)-3*(1-r)*(1-2*delta))
    drift=(u+v)*bit-3*(1-r)*(J+delta)
    return successor,drift


def graph(t):
    vertices=[(j,u,v) for j in range(t) for u in range(-13,3*t+1) for v in VS]
    edges={x:[] for x in vertices}
    for x in vertices:
        for bit in(0,1):
            y,drift=edge(t,x,bit)
            if -13<=y[1]<=3*t:edges[x].append((y,bit,drift))
    return edges


def sccs(edges):
    serial=0;stack=[];active=set();index={};low={};parts=[]
    def visit(x):
        nonlocal serial
        index[x]=low[x]=serial;serial+=1;stack.append(x);active.add(x)
        for y,_,_ in edges[x]:
            if y not in index:visit(y);low[x]=min(low[x],low[y])
            elif y in active:low[x]=min(low[x],index[y])
        if low[x]==index[x]:
            part=[]
            while True:
                y=stack.pop();active.remove(y);part.append(y)
                if y==x:break
            parts.append(part)
    for x in edges:
        if x not in index:visit(x)
    return parts


def zero_cycle_certificate(edges,part):
    allowed=set(part);potential={part[0]:0};queue=deque([part[0]])
    while queue:
        x=queue.popleft()
        for y,_,weight in edges[x]:
            if y not in allowed:continue
            value=potential[x]+weight
            if y in potential:
                if potential[y]!=value:return None
            else:potential[y]=value;queue.append(y)
    return [dict(state=x,potential=potential[x]) for x in sorted(part)]


def reachable(edges,start):
    seen={start};queue=deque([start])
    while queue:
        for y,_,_ in edges[queue.popleft()]:
            if y not in seen:seen.add(y);queue.append(y)
    return seen


def walk(t,state,word):
    drift=0
    for bit in word:
        state,increment=edge(t,state,int(bit));drift+=increment
    return state,drift


def L(B):return 13*g(B)+32*B+6*(B&1)+18
def N(C):return 3*g(C)+12*C+10
def E(n):return 13*g(n)-3*g(7*n+4)-n-3*(n&1)-80
def pad(q,m):return (16**m)*q+9*(16**m-1)//15


def main():
    five=graph(5);parts=sccs(five);certificates=[]
    for part in parts:
        cyclic=len(part)>1 or any(y==part[0] for y,_,_ in five[part[0]])
        if not cyclic:continue
        cert=zero_cycle_certificate(five,part)
        if cert is None:raise RuntimeError(('nonzero five-cycle',part))
        certificates.append(cert)
    ranks={}
    def rank(x):
        if x not in ranks:
            ranks[x]=1+max((rank(y) for y,_,_ in five[x] if y!=x),default=-1)
        return ranks[x]
    for x in five:rank(x)
    for x,arcs in five.items():
        for y,_,weight in arcs:
            if (x==y and weight!=0) or (x!=y and ranks[x]<=ranks[y]):
                raise RuntimeError(('rank certificate',x,y))
    rank_data=json.dumps([dict(state=x,rank=ranks[x]) for x in sorted(ranks)],indent=2)+'\n'
    (HERE/'five_rank_certificate.json').write_text(rank_data)
    actual_channels=[]
    for t,u,v in((5,4,0),(7,-7,13)):
        gr=graph(t);components=sccs(gr)
        cyclic={x for part in components if len(part)>1 or any(y==part[0] for y,_,_ in gr[part[0]]) for x in part}
        for j in range(t):
            seen=reachable(gr,(j,u,v))
            if seen&cyclic:raise RuntimeError(('direct channel reaches loop',t,j))
            actual_channels.append(dict(multiplier=t,initial=[j,u,v],reachable_states=len(seen),cycle_reachable=False))
    cycle=walk(7,(4,-1,-3),'1001')
    if cycle!=((4,-1,-3),4):raise RuntimeError(('seven cycle',cycle))
    for n in range(256):
        if L(64*n+36)-N(224*n+135)!=E(n):raise RuntimeError(('physical defect',n))
        if E(16*n+9)!=E(n)+4:raise RuntimeError(('actuator',n))
    preparations=[]
    for p in range(64):
        q1,q2=16*p+7,32*p+15
        if E(q1)>=0 or E(q2)>=0 or E(q2)-E(q1)!=-100*p-50:
            raise RuntimeError(('preparation',p))
        options=[q for q in(q1,q2) if E(q)%4==0]
        if len(options)!=1:raise RuntimeError(('phase choice',p))
        q=options[0];m=-E(q)//4;n=pad(q,m)
        B,C=64*n+36,224*n+135
        if E(n)!=0 or L(B)!=N(C):raise RuntimeError(('rebase',p))
        preparations.append(dict(prefix=p,q=q,loops=m,n_bits=n.bit_length(),H_bits=L(B).bit_length()))
    # Small literal graphs authenticate the four source charts and constants.
    path='05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
    source=subprocess.check_output(['git','show','HEAD:'+path],cwd=ROOT)
    module=types.ModuleType('inherited_compiler');exec(compile(source,path,'exec'),module.__dict__)
    literal=[]
    for n in(0,1):
        B,C=64*n+36,224*n+135
        charts=[([18,8,3,5,14],6,B,L(B)),([18,8,3,5,14],3,B,L(B)+2),
                ([2,4,1],3,2*C+1,N(C)),([2,4,1],1,4*C+3,N(C)+2)]
        rows=[]
        for core,cut,address,expected in charts:
            tournament=module.construct(core,cut,0,address)
            actual=module.hp_subsets(tournament)[-1]
            if actual!=expected:raise RuntimeError(('literal',n,cut,address))
            rows.append(dict(core=core,cut=cut,address=address,vertices=len(tournament),H=actual))
        literal.append(dict(n=n,charts=rows))
    q=13;m=-E(q)//4;n=pad(q,m);B,C=64*n+36,224*n+135
    seed=dict(q=q,defect=E(q),loops=m,n=n,B=B,C=C,H_pair=[L(B),L(B)+2],
              C3_addresses=[2*C+1,4*C+3])
    report=dict(status='PROVED constructions and finite complete carry-graph certificate; no arbitrary-target coverage',
        five_graph=dict(vertices=len(five),edges=sum(map(len,five.values())),cyclic_components=len(certificates),
                        cyclic_vertices=sum(map(len,certificates)),zero_drift_potentials=certificates,
                        rank_max=max(ranks.values()),rank_certificate_sha256=hashlib.sha256(rank_data.encode()).hexdigest()),
        direct_channels=actual_channels,seven_cycle=dict(state=[4,-1,-3],lsb_word='1001',suffix_value=9,drift=4),
        preparation_controls=preparations,literal_controls=literal,seed=seed,
        inherited_compiler_sha256=hashlib.sha256(source).hexdigest())
    (HERE/'periodic_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in('preparation_controls','literal_controls','inherited_compiler_sha256')},indent=2))


if __name__=='__main__':main()
