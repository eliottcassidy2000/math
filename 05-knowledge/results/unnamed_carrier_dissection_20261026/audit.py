"""Independent Python checks: rooted recovery, actual multiset/probability closure,
and corrected atom pushforwards. Canonicalization/DFS do not import the old
producer; its repaired aggregation helper alone is extracted for a regression check.
The frozen atlas supplies the complete strong-class universe through seven only.
"""
from collections import Counter,defaultdict
from functools import lru_cache,reduce
from itertools import permutations,product,chain
from pathlib import Path
from math import gcd,isqrt
import json,random,ast
ROOT=Path(__file__).resolve().parents[3]
def need(v,s):
    if not v:raise RuntimeError(s)
def arc(t,u,v):return t[u]>>v&1
def relabel(t,p):return tuple(sum(arc(t,p[i],p[j])<<j for j in range(len(p))) for i in range(len(p)))
@lru_cache(None)
def canon(t,roots=()):
    fixed=tuple(dict.fromkeys(roots));d=defaultdict(list)
    for i in range(len(t)):
        if i not in fixed:d[t[i].bit_count()].append(i)
    return min(relabel(t,fixed+tuple(chain.from_iterable(q))) for q in product(*(permutations(d[k]) for k in sorted(d))))
def strong(t):
    n=len(t)
    for reverse in (False,True):
        seen={0};old=set()
        while old!=seen:
            old=set(seen)
            seen|={j for i in old for j in range(n) if (arc(t,j,i) if reverse else arc(t,i,j))}
        if len(seen)!=n:return False
    return True
def induced(t,keep):return relabel(t,keep)
def paths(t):
    full=(1<<len(t))-1
    def dfs(p,m):
        if m==full:yield p;return
        q=(full^m) if not p else t[p[-1]]&(full^m)
        while q:
            b=q&-q;q-=b;yield from dfs(p+(b.bit_length()-1,),m|b)
    yield from dfs((),0)
def partner(t,p):
    r=[((1<<len(t))-1)^x^(1<<i) for i,x in enumerate(t)]
    for a,b in zip(p,p[1:]):r[a]^=1<<b;r[b]^=1<<a
    return tuple(r)
def pres(n,m):
    r=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            f=j==i+1 or m>>k&1
            if j>i+1:k+=1
            if f:r[i]|=1<<j
            else:r[j]|=1<<i
    return tuple(r)
def sandwich(t,s,e):
    n=len(t);r=list(t)+[1<<s,1<<n]
    for i in range(n):
        if i!=s:r[i]|=1<<n
        if i==e:r[i]|=1<<(n+1)
        else:r[n+1]|=1<<i
    return tuple(r)
def decoders(t):
    n=len(t);out=[]
    for a in range(n):
        if t[a].bit_count()!=1:continue
        for b in range(n):
            if a==b or t[b].bit_count()!=n-2:continue
            keep=tuple(v for v in range(n) if v not in (a,b))
            if not strong(induced(t,keep)):continue
            s=next(i for i in range(n) if arc(t,a,i));e=next(i for i in range(n) if i!=b and arc(t,i,b))
            if s in keep and e in keep:out.append((a,b,s,e))
    return out
def refine(states,succ,init,prob=False):
    cell={x:init[x] for x in states};rounds=0
    def norm(d):
        ids={};return {x:ids.setdefault(d[x],len(ids)) for x in states}
    cell=norm(cell)
    while True:
        sig={}
        for x in states:
            c=Counter(cell[y] for y in succ[x]);g=reduce(gcd,c.values(),0) if prob else 1
            row=tuple(sorted((k,v//g) for k,v in c.items())) if c else ()
            sig[x]=(cell[x],row)
        new=norm(sig)
        if len(set(new.values()))==len(set(cell.values())):return new,rounds
        rounds+=1;cell=new
def main():
    legacy=ROOT/'04-computation/experiments/joint_port_bisimulation_20261026/atom_responses.py'
    node=next(x for x in ast.parse(legacy.read_text()).body if isinstance(x,ast.FunctionDef) and x.name=='forget_ports')
    scope={'Counter':Counter};exec(compile(ast.Module(body=[node],type_ignores=[]),str(legacy),'exec'),scope)
    need(scope['forget_ports']([('a',2),('b',3)],{'a':'c','b':'c'})==Counter({'c':5}),'repaired aggregation hostile')
    atlas=json.loads((ROOT/'05-knowledge/results/strong_atom_continuation_vitali_20261026/strong-atlas.json').read_text())
    result={'sandwich':[],'atom_pushforward':[],'closure':[]}
    for n in range(3,8):
        atoms=[a for a in atlas if a['n']==n];states=set();fs=[];js=[];lost=[]
        for a in atoms:
            t=tuple(a['source']);need(strong(t),'atlas source strong')
            for s in range(n):
                for e in range(n):
                    if s==e:continue
                    states.add(canon(t,(s,e)))
                    u=sandwich(t,s,e);need(decoders(u)==[(n,n+1,s,e)],'unique strong-core sandwich decoder')
            c=Counter();j=Counter();wrong={};legacy_response=[];legacy_fc={}
            for h in a['response']:
                v=canon(tuple(h['rows']));c[v]+=h['total'];j[canon(tuple(h['rows']),tuple(h['roots']))]+=h['total'];wrong[v]=h['total']
                key=canon(tuple(h['rows']),tuple(h['roots']));legacy_response.append((key,h['total']));legacy_fc[key]=v
            need(sum(c.values())==a['presentations'],'pushforward mass')
            need(scope['forget_ports'](legacy_response,legacy_fc)==c,'repaired legacy pushforward matches independent sum')
            fs.append(tuple(sorted(c.items())));js.append(tuple(sorted(j.items())));lost.append(a['presentations']-sum(wrong.values()))
        need(len(set(fs))==len(atoms),'free signatures classify finite strong atlas')
        need(len(set(js))==len(atoms),'joint signatures classify finite strong atlas')
        result['sandwich'].append({'n':n,'atoms':len(atoms),'joint_states':len(states),'decoded_labelled_root_pairs':len(atoms)*n*(n-1)})
        result['atom_pushforward'].append({'n':n,'free_signatures':len(set(fs)),'joint_signatures':len(set(js)),'old_overwrite_mass_loss_max':max(lost)})
    rnd=random.Random(2026102602)
    for n in range(8,16):
        for _ in range(20):
            while True:
                t=pres(n,rnd.getrandbits((n-1)*(n-2)//2))
                if strong(t):break
            s,e=rnd.sample(range(n),2);u=sandwich(t,s,e);p=rnd.sample(range(n+2),n+2);v=relabel(u,p)
            need(decoders(v)==[(p.index(n),p.index(n+1),p.index(s),p.index(e))],'random relabelled decoder')
    result['random_sandwich']={'seed':2026102602,'orders':list(range(8,16)),'samples_each':20}
    for n in range(3,8):
        classes=sorted({canon(pres(n,m)) for m in range(1<<((n-1)*(n-2)//2))});succ={};free={};one={}
        for t in classes:
            ps=list(paths(t));by=defaultdict(list)
            for p in ps:by[p[0],p[-1]].append(canon(partner(t,p),(p[0],p[-1])))
            for s in range(n):
                for e in range(n):
                    if s==e:continue
                    x=canon(t,(s,e));new=sorted(by[s,e])
                    if x in succ:need(succ[x]==new,'root orbit successor equality')
                    succ[x]=new;free[x]=t;one[x]=(canon(t,(s,)),canon(t,(e,)))
        states=sorted(succ);need(all(y in succ for x in states for y in succ[x]),'dynamics closed')
        st={x for x in states if strong(x)};rr={'n':n,'states':len(states),'strong_states':len(st),'observations':{}}
        for name,obs in [('free',free),('one_root_pair',one)]:
            for prob in (False,True):
                cell,depth=refine(states,succ,obs,prob);rr['observations'][name+('/probability' if prob else '/multiset')]={'cells':len(set(cell.values())),'strong_cells':len({cell[x] for x in st}),'strict_refinements':depth}
        result['closure'].append(rr)
    # Exact Gaussian and cyclic-phase transport controls.
    a,b=1,0;ga=[]
    for k in range(1,9):
        a,b=2*a-b,a+2*b;x,y,c=abs(a*a-b*b),abs(2*a*b),5**k
        need(x*x+y*y==c*c and gcd(x,y)==1,'Gaussian primitive ray');ga.append([k,a,b,x,y,c])
    need(ga[2][3:]==[117,44,125],'44/5 norm-power identity')
    result['gaussian_ray']=ga
    result['two_square_44']=any(a*a+b*b==44 for a in range(7) for b in range(7));need(not result['two_square_44'],'44 inert obstruction')
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
