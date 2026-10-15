"""Exact heavy-child gold defects; exploratory selector tests, not a census."""
from functools import lru_cache
from itertools import combinations, permutations
from fractions import Fraction
from pathlib import Path
from random import Random
import json


def closure(n, edges):
    if type(n) is not int or not 1 <= n <= 18:
        raise ValueError('exact size in 1..18 required')
    if type(edges) is not tuple: raise ValueError('tuple edges required')
    pred=[0]*n
    for edge in edges:
        if type(edge) is not tuple or len(edge)!=2: raise ValueError('edge pair required')
        a,b=edge
        if any(type(z) is not int or not 0<=z<n for z in edge) or a==b:
            raise ValueError('distinct vertices required')
        pred[b]|=1<<a
    for k in range(n):
        for v in range(n):
            if pred[v]>>k&1:pred[v]|=pred[k]
    if any(pred[i]>>i&1 for i in range(n)):raise ValueError('cyclic constraints')
    return tuple(pred)


def relations(p):
    return tuple((a,b) for b in range(len(p)) for a in range(len(p)) if p[b]>>a&1)


def extend(p,a,b):
    return closure(len(p),relations(p)+((a,b),))


def pair_counts(p):
    """Classical forward/backward order-ideal DP, exact Python integers."""
    if type(p) is not tuple or not p or any(type(x) is not int or not 0<=x<(1<<len(p)) for x in p):
        raise ValueError('exact predecessor-mask tuple required')
    n=len(p)
    if closure(n,relations(p))!=p:raise ValueError('closed poset required')
    return _pair_counts(p)


@lru_cache(maxsize=4096)
def _pair_counts(p):
    n=len(p);full=(1<<n)-1
    f=[0]*(full+1);back=[0]*(full+1);f[0]=1;back[full]=1
    transitions={}
    for mask in range(full+1):
        if not f[mask]:continue
        eligible=[v for v in range(n) if not mask>>v&1 and p[v]&mask==p[v]]
        transitions[mask]=eligible
        for v in eligible:f[mask|1<<v]+=f[mask]
    for mask in sorted(transitions,reverse=True):
        if mask!=full:back[mask]=sum(back[mask|1<<v] for v in transitions[mask])
    c=[[0]*n for _ in range(n)]
    for mask,vs in transitions.items():
        xs=[x for x in range(n) if mask>>x&1]
        for y in vs:
            mass=f[mask]*back[mask|1<<y]
            for x in xs:c[x][y]+=mass
    return f[full],tuple(tuple(row) for row in c)


def balanced(p):
    e,c=pair_counts(p);best=-1;pairs=[]
    for x,y in combinations(range(len(p)),2):
        minor=min(c[x][y],c[y][x])
        if minor>best:best=minor;pairs=[(x,y)]
        elif minor==best:pairs.append((x,y))
    return max(best,0),tuple(pairs)


def first_receipt(p,x,y):
    e,c=pair_counts(p)
    if any(type(v) is not int or not 0 <= v < len(p) for v in (x,y)) or x==y:
        raise ValueError('distinct exact integer vertices required')
    if not 0<c[x][y]<e:raise ValueError('incomparable first pair required')
    if c[x][y]<c[y][x]:x,y=y,x
    q=extend(p,x,y);m=c[x][y];light=e-m
    s,pairs=balanced(q)
    second=pairs[0] if pairs else None
    spin=None
    if second:
        u,v=second;_,child=pair_counts(q);d=child[u][v]
        an=m-light;bn=c[u][v]-c[v][u];rn=4*d-2*m-bn
        spin={'mean_first_numerator':an,'mean_second_numerator':bn,
              'correlation_numerator':rn,'denominator':e,
              'cone_slack_numerator':e-3*an-abs(bn+rn),
              'cells':[[d,m-d],[c[u][v]-d,light-(c[u][v]-d)]]}
        if spin['cone_slack_numerator'] != -4*(m-light-s):raise ArithmeticError('spin identity')
    return {'first':[x,y],'total':e,'heavy':m,'light':light,
            'child_minor':s,'second':list(second) if second else None,
            'defect':m-light-s,'valid':m-light-s<=0,'spin':spin}


def literal_extensions(p):
    n=len(p)
    def walk(mask,prefix):
        if len(prefix)==n:
            yield prefix;return
        for v in range(n):
            if not mask>>v&1 and p[v]&mask==p[v]:
                yield from walk(mask|1<<v,prefix+(v,))
    return tuple(walk(0,()))


def delete_vertex(p,v):
    vertices=[x for x in range(len(p)) if x!=v];index={x:i for i,x in enumerate(vertices)}
    return closure(len(vertices),tuple((index[a],index[b]) for a,b in relations(p) if a!=v and b!=v))


def inspect(p):
    e,c=pair_counts(p);s,pairs=balanced(p)
    if e==1:return {'chain':True}
    greedy=[first_receipt(p,*pair) for pair in pairs]
    all_results=[first_receipt(p,x,y) for x,y in combinations(range(len(p)),2)
                 if 0<c[x][y]<e]
    return {'n':len(p),'relations':[list(e) for e in relations(p)],'extensions':e,
            'balance':str(Fraction(s,e)),'greedy':greedy,
            'best':min(all_results,key=lambda r:r['defect']),
            'all_best_work':all(r['valid'] for r in greedy),
            'some_best_works':any(r['valid'] for r in greedy)}


def main():
    rng=Random(20261014);universe=set()
    # Deliberately varied finite probes of a NEW greedy amplification claim.
    for n in range(3,12):
        for numerator in (1,2,3,5,7,9):
            for _ in range(18):
                universe.add(closure(n,tuple((a,b) for a,b in combinations(range(n),2)
                                             if rng.randrange(10)<numerator)))
    input_path=Path(__file__).resolve().parents[2]/'05-knowledge/results/gold_low_balance_inputs_20261014j.json'
    hard_inputs=json.loads(input_path.read_text(encoding='utf-8'))
    hard_posets=[]
    for h in hard_inputs['predecessor_mask_hex']:
        p=tuple(int(h[i:i+4],16) for i in range(0,len(h),4))
        if len(p)!=14 or closure(14,relations(p))!=p:raise ValueError('invalid archived input')
        hard_posets.append(p);universe.add(p)
        covers=tuple((a,b) for a,b in relations(p) if not any(
            p[z]>>a&1 and p[b]>>z&1 for z in range(14)))
        # New perturbation probes: remove each cover, and add eight selected relations.
        for e in covers:universe.add(closure(14,tuple(x for x in covers if x!=e)))
        inc=[(a,b) for a,b in combinations(range(14),2) if not(p[a]>>b&1 or p[b]>>a&1)]
        for a,b in inc[:8]:
            universe.add(extend(p,a,b));universe.add(extend(p,b,a))
    checks=0;fail_all=[];fail_some=[];minimum=None;tested=0;hard_rows=[];global_failures=[]
    for p in sorted(universe,key=lambda p:(len(p),p)):
        row=inspect(p)
        if row.get('chain'):continue
        tested+=1
        if not row['all_best_work'] and len(fail_all)<6:fail_all.append(row)
        if not row['some_best_works'] and len(fail_some)<6:fail_some.append(row)
        if p in hard_posets:hard_rows.append(row)
        if minimum is None or Fraction(row['balance'])<Fraction(minimum['balance']):minimum=row
        checks+=1
        if not row['best']['valid']:global_failures.append(row)
    # Independent permutations check both the DP and literal four-cell criterion.
    for p in sorted(universe,key=lambda p:(len(p),p)):
        if len(p)>6:break
        n=len(p);extensions=[]
        for order in permutations(range(n)):
            pos={x:i for i,x in enumerate(order)}
            if all(pos[a]<pos[b] for a,b in relations(p)):extensions.append(pos)
        e,c=pair_counts(p)
        checks+=1
        if e!=len(extensions):raise ArithmeticError('count disagreement')
        for x,y in combinations(range(n),2):
            checks+=1
            if c[x][y]!=sum(pos[x]<pos[y] for pos in extensions):raise ArithmeticError('pair disagreement')
            if not 0<c[x][y]<e:continue
            r=first_receipt(p,x,y);a,b=r['first'];u,v=r['second']
            valid=True
            for first in (False,True):
                subset=[pos for pos in extensions if (pos[a]<pos[b])==first]
                for second in (False,True):
                    cell=sum((pos[u]<pos[v])==second for pos in subset)
                    if e<len(subset)+cell:valid=False
            checks+=1
            if valid!=r['valid']:raise ArithmeticError('literal two-comparison mismatch')
    hostile_audits=[]
    for row in fail_some:
        p=closure(row['n'],tuple(map(tuple,row['relations'])))
        orders=literal_extensions(p);positions=[{x:i for i,x in enumerate(o)} for o in orders]
        e,c=pair_counts(p)
        for a,b in combinations(range(len(p)),2):
            checks+=1
            if c[a][b]!=sum(x[a]<x[b] for x in positions):raise ArithmeticError('hostile pair count')
        receipts=row['greedy']+[row['best']]
        for r in receipts:
            a,b=r['first'];selected=[x for x in positions if x[a]<x[b]]
            s=max(min(sum(x[u]<x[v] for x in selected),sum(x[v]<x[u] for x in selected))
                  for u,v in combinations(range(len(p)),2))
            checks+=1
            if s!=r['child_minor'] or len(selected)!=r['heavy']:raise ArithmeticError('hostile heavy child')
            u,v=r['second']
            cells=[[sum(x[u]<x[v] for x in selected),sum(x[v]<x[u] for x in selected)],
                   [sum(x[a]>x[b] and x[u]<x[v] for x in positions),
                    sum(x[a]>x[b] and x[v]<x[u] for x in positions)]]
            checks+=1
            if cells!=r['spin']['cells']:raise ArithmeticError('hostile four cells')
        deletion_results=[inspect(delete_vertex(p,v))['some_best_works'] for v in range(len(p))]
        checks+=len(deletion_results)
        # Retain every tied best reply, not merely the lexicographic one.
        component={(4,7),(6,7),(6,9)}
        reply_rows=[]
        for pair in sorted(component):
            receipt=first_receipt(p,*pair)
            replies=balanced(extend(p,*receipt['first']))[1]
            checks+=1
            if receipt['defect']<=0 or not replies or not set(replies)<=component:
                raise ArithmeticError('claimed closed hostile reply component')
            reply_rows.append({'pair':list(pair),'defect':receipt['defect'],
                               'all_best_replies':[list(x) for x in replies]})
        hostile_audits.append({'n':len(p),'literal_extensions':len(orders),
                              'every_single_vertex_deletion_repairs_greedy':all(deletion_results),
                              'deletion_results':deletion_results,
                              'closed_best_response_component':reply_rows,
                              'hasse_edges':[list((a,b)) for a,b in relations(p) if not any(
                                  p[z]>>a&1 and p[b]>>z&1 for z in range(len(p)))]})
    type_source=closure(3,())
    pair_counts((0,0))
    for bad_source in ((False,0),(0.0,0),[0,0]):
        checks+=1
        try:pair_counts(bad_source)
        except ValueError:pass
        else:raise ArithmeticError('cached malformed source accepted')
    for bad_pair in ((True,2),(-1,2),(0,3),(1,1),(1.0,2)):
        checks+=1
        try:first_receipt(type_source,*bad_pair)
        except ValueError:pass
        else:raise ArithmeticError('malformed comparison accepted')
    out={'status':'FINITE-EXACT exploratory selector probe; universal conjectures OPEN',
         'universe':{'seed':20261014,'random_sizes':[3,11],'raw_draws':972,
                     'hard_inputs':42,'perturbation':'delete each cover or add either orientation of first eight incomparable pairs',
                     'unique':len(universe),'nonchains':tested},
         'checks':checks,'all_best_pair_counterexamples':fail_all,
         'some_best_pair_counterexamples':fail_some,'least_balance_example':minimum,
         'hard_input_results':hard_rows,'global_gold_failures':global_failures,
         'hostile_independent_audits':hostile_audits}
    path=Path(__file__).resolve().parents[2]/'05-knowledge/results/gold_heavy_child_20261014j.json'
    path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'checks':checks,'tested':tested,'all_best_failures':len(fail_all),
                      'some_best_failures':len(fail_some),'least_balance':minimum['balance']}))


if __name__=='__main__':main()
