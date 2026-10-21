"""Exact mass-deficit-k zero-heavy boundary; no fitted-polynomial inference."""
from collections import Counter, defaultdict
from itertools import permutations, product
from math import comb
from pathlib import Path
import argparse
import json
from hashlib import sha256

HERE=Path(__file__).resolve().parent
OUT=HERE/'receipt.json'
CHECKS=0


def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)


def integer(x,minimum=0):
    if type(x) is not int or x<minimum: raise ValueError('exact integer outside domain')
    return x


def matchings(xs):
    if not xs:
        yield ();return
    a=xs[0]
    for b in xs[1:]:
        for rest in matchings(tuple(x for x in xs[1:] if x!=b)):
            yield ((a,b),)+rest


def compositions(total,parts):
    if parts==1:
        yield (total,);return
    for a in range(total+1):
        for rest in compositions(total-a,parts-1):
            yield (a,)+rest


def df(k):
    integer(k)
    ans=1
    for x in range(1,2*k,2): ans*=x
    return ans


def indecomp(m,k):
    return all(any(a<c<=b for a,b in m) for c in range(1,2*k))


def indecomp_counts(n):
    ans=[1]
    for k in range(1,n+1):
        ans.append(df(k)-sum(ans[j]*df(k-j) for j in range(1,k)))
    return ans


def record(k,matching,gaps,strong=False):
    integer(k,1)
    if type(strong) is not bool or type(matching) is not tuple or type(gaps) is not tuple:
        raise ValueError('literal matching/gap record required')
    if len(gaps)!=(2*k-1 if strong else 2*k+1): raise ValueError('wrong gap arity')
    if any(type(x) is not int or not 0<=x<=k for x in gaps): raise ValueError('invalid gap weight')
    flat=[]
    for e in matching:
        if type(e) is not tuple or len(e)!=2 or any(type(x) is not int for x in e) or not e[0]<e[1]:
            raise ValueError('literal ordered matching pair required')
        flat.extend(e)
    if len(matching)!=k or sorted(flat)!=list(range(2*k)): raise ValueError('not a perfect matching')
    if tuple(sorted(matching))!=matching: raise ValueError('canonical pair-list order required')
    if sum(k-g for g in gaps)!=k: raise ValueError('not the mass-deficit-k boundary')
    slots={}; weights=[]; full=(0,)+gaps+(0,) if strong else gaps
    for i,g in enumerate(full):
        if g: weights.append(g)
        if i<2*k:
            slots[i]=len(weights);weights.append(1)
    n=len(weights); rows=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    for a,b in matching:
        a,b=slots[a],slots[b];rows[a]^=1<<b;rows[b]^=1<<a
    return tuple(rows),tuple(weights)


def audit_tournament(rows,weights):
    if type(rows) is not tuple or type(weights) is not tuple or len(rows)!=len(weights):
        raise ValueError('literal complete tournament and weight tuples required')
    n=len(rows)
    for i,r in enumerate(rows):
        if type(r) is not int or not 0<=r<(1<<n) or r>>i&1: raise ValueError('invalid tournament row')
        if type(weights[i]) is not int or weights[i]<1: raise ValueError('invalid positive weight')
        for j in range(i):
            if ((r>>j)&1)+((rows[j]>>i)&1)!=1: raise ValueError('pair orientation missing or doubled')


def clone_pairs(rows):
    n=len(rows);mask=(1<<n)-1
    return tuple((i,j) for i in range(n) for j in range(i+1,n)
                 if ((rows[i]^rows[j])&(mask^(1<<i)^(1<<j)))==0)


def transform(rows,weights,order):
    return (tuple(sum(((rows[i]>>j)&1)<<p for p,j in enumerate(order)) for i in order),
            tuple(weights[i] for i in order))


def canonical(rows,weights):
    audit_tournament(rows,weights)
    big=tuple(i for i,w in enumerate(weights) if w>1)
    ordered=tuple(sorted(big,key=lambda i:-sum((rows[i]>>j)&1 for j in big)))
    if any(not rows[a]>>b&1 for i,a in enumerate(ordered) for b in ordered[i+1:]):
        raise ValueError('large-vertex induced order is not transitive')
    cells=defaultdict(list)
    for i,w in enumerate(weights):
        if w==1:
            cells[tuple((rows[i]>>b)&1 for b in ordered)].append(i)
    best=None
    for blocks in product(*(permutations(cells[s]) for s in sorted(cells))):
        order=ordered+tuple(i for block in blocks for i in block)
        rr,ww=transform(rows,weights,order)
        key=(ww,rr)
        if best is None or key<best:best=key
    return best,max(map(len,cells.values()),default=0)


def independent_canonical(rows,weights):
    # Unrelated invariant-cell order: weight and full outdegree, then exhaustive
    # permutations in each cell. Valid for arbitrary weighted tournaments.
    audit_tournament(rows,weights);cells=defaultdict(list)
    for i,w in enumerate(weights):cells[w,rows[i].bit_count()].append(i)
    best=None
    for blocks in product(*(permutations(cells[s]) for s in sorted(cells))):
        order=tuple(i for block in blocks for i in block)
        rr,ww=transform(rows,weights,order);key=(ww,rr)
        if best is None or key<best:best=key
    return best


def strong(rows):
    n=len(rows)
    for reverse in (False,True):
        adj=tuple(sum(((rows[j]>>i)&1)<<j for j in range(n)) for i in range(n)) if reverse else rows
        seen=1
        while True:
            nxt=seen
            for i in range(n):
                if seen>>i&1:nxt|=adj[i]
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<n)-1:return False
    return True


def optimum(rows,weights):
    n=len(rows);dp=[0]+[10**30]*((1<<n)-1)
    for mask in range(1,1<<n):
        best=10**30
        for i in range(n):
            if mask>>i&1:
                rest=mask^(1<<i)
                cost=weights[i]*sum(weights[j] for j in range(n) if rest>>j&1 and rows[i]>>j&1)
                best=min(best,dp[rest]+cost)
        dp[mask]=best
    return dp[-1]


def bad_gap(k,m,gaps,isstrong):
    full=(0,)+gaps+(0,) if isstrong else gaps
    return next((i for i in range(1,2*k) if full[i]==0 and (i-1,i) in m),None)


def formula(k,isstrong):
    integer(k,1)
    if type(isstrong) is not bool:raise ValueError('literal strong flag required')
    if not isstrong:return df(k)*(comb(3*k,k)-1)
    if k==1:return 0
    I=indecomp_counts(k)
    return I[k]*comb(3*k-2,k)-(2*k-3)*I[k-1]


def reject(fn,msg):
    try:fn()
    except ValueError:need(True,msg)
    else:need(False,msg)


def verify():
    global CHECKS
    CHECKS=0;layers=[];hostiles=[];small_keys={};covariance=0;dp_checks=0
    for k in range(1,5):
        for isstrong in (False,True):
            g=2*k-1 if isstrong else 2*k+1
            seen=set();secondary=set();raw=invalid=0;sectors=Counter();cellhist=Counter()
            for m in matchings(tuple(range(2*k))):
                if isstrong and not indecomp(m,k):continue
                for ds in compositions(k,g):
                    gaps=tuple(k-d for d in ds)
                    if min(gaps)<0:continue
                    rows,w=record(k,m,gaps,isstrong);raw+=1
                    bad=bad_gap(k,m,gaps,isstrong);clones=clone_pairs(rows)
                    need(bool(clones)==(bad is not None),'exact missing-gap clone iff')
                    if clones:
                        invalid+=1
                        if k<=2 and len(hostiles)<3:
                            hostiles.append({'k':k,'strong':isstrong,'matching':m,'gaps':gaps,'clone_pairs':clones,'rows':rows,'weights':w})
                        continue
                    need(sum(w)==2*k*k+(k if isstrong else 3*k)-k,'exact weighted boundary mass')
                    key,cell=canonical(rows,w)
                    need(key not in seen,'no distinct numerical-record isomorphism collision')
                    seen.add(key);cellhist[cell]+=1
                    sectors['missing' if 0 in gaps else 'full']+=1
                    if isstrong:need(strong(rows),'indecomposable matching remains strong')
                    if k<=2:
                        other=independent_canonical(rows,w)
                        need(other not in secondary,'independent full isomorphism quotient')
                        secondary.add(other)
                        need(optimum(rows,w)==k,'independent exact optimum at small boundary')
                        dp_checks+=1
                    if len(seen)<=12 or (k==3 and len(seen)%100==0):
                        for order in (tuple(reversed(range(len(rows)))),tuple(range(1,len(rows)))+(0,)):
                            rr,ww=transform(rows,w,order)
                            need(canonical(rr,ww)[0]==key,'full relabelling covariance')
                            covariance+=1
                    if k==3 and len(seen)<=4:
                        need(optimum(rows,w)==k,'independent weighted optimum k3')
                        dp_checks+=1
            I=indecomp_counts(k)
            excluded=(2*k-3)*I[k-1] if isstrong and k>=2 else (1 if isstrong else df(k))
            need(invalid==excluded,'marked adjacent matching-edge count')
            need(len(seen)==formula(k,isstrong),'proved formula matches full quotient')
            if k<=2:small_keys[k,isstrong]=secondary
            layers.append({'k':k,'strong':isstrong,'raw_records':raw,'clone_rejections':invalid,
                           'classes':len(seen),'sectors':dict(sectors),'unit_cell_sizes':dict(cellhist)})
    # The previously independent full-budget compiler is consumed through its
    # published exact mass coefficients, not assumed for the structural proof.
    old=json.loads((HERE.parent/'charge_palettes_boundary_20261020g/zero-heavy-mass.json').read_text())
    for row in old['complete_budget_compiler']:
        k=row['k']
        for isstrong in (False,True):
            key='strong_zero_heavy_by_mass' if isstrong else 'all_zero_heavy_by_mass'
            mass=2*k*k+(k if isstrong else 3*k)-k
            need(int(row[key].get(str(mass),0))==formula(k,isstrong),'independent complete old compiler agreement')
    # Exact endpoint-deficit bounds: at delta=k, r<k or shared endpoints fail.
    for k in range(1,13):
        for r in range(k+1):
            for t in range(2*r+1):
                lower=(2*k+1)*(k-r) if r<k else (k+1)*(2*k-t)
                if r<k or t<2*r:need(lower>k,'next mass layer excludes slack edges and shared endpoints')
    need(formula(2,False)==42 and df(2)*comb(6,2)==45,'literal extrapolation hostile42 versus45')
    # Reject aliases before any inferred combinatorial interpretation.
    reject(lambda:formula(True,False),'Boolean budget rejected')
    reject(lambda:formula(2,1),'integer strong flag rejected')
    reject(lambda:record(2,((0,1),(2,3)),(2,2,2,2,1)),'wrong deficit rejected')
    reject(lambda:record(2,((False,1),(2,3)),(0,2,2,2,2)),'Boolean endpoint rejected')
    rr,ww=record(2,((0,1),(2,3)),(0,2,2,2,2))
    reject(lambda:canonical(rr,tuple(float(x) for x in ww)),'float weights rejected')
    need(formula(3,False)==1245 and formula(3,True)==344,'new k3 values')
    need(formula(4,False)==51870 and formula(4,True)==15490,'new k4 values')
    # At the next layer shared endpoints really occur: the two reversed
    # endpoint edges share the middle endpoint, with every gap retained.
    next_rows=[((1<<7)-1)^((1<<(i+1))-1) for i in range(7)]
    for a,b in ((1,3),(3,5)):
        next_rows[a]^=1<<b;next_rows[b]^=1<<a
    next_rows=tuple(next_rows);next_weights=(2,1,2,1,2,1,2)
    need(not clone_pairs(next_rows),'shared-endpoint next-layer witness pair free')
    need(optimum(next_rows,next_weights)==2 and sum(next_weights)==11,'shared endpoint first occurs at deficit k+1 in witness')
    return {'status':'PROVED first missing-gap layer; FINITE-EXACT complete quotient controls',
            'checks':CHECKS,'layers':layers,'clone_hostiles':hostiles,
            'relabel_controls':covariance,'weighted_DP_controls':dp_checks,
            'next_layer_shared_endpoint':{'k':2,'deficit':3,'rows':next_rows,'weights':next_weights,'backward_edges':((1,3),(3,5))},
            'formula_values':[{'k':k,'all':formula(k,False),'strong':formula(k,True)} for k in range(1,11)]}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    payload=verify();text=json.dumps(payload,sort_keys=True,indent=2)+'\n'
    if args.verify:
        if OUT.read_text(encoding='utf-8')!=text:raise RuntimeError('saved receipt differs')
    else:OUT.write_text(text,encoding='utf-8',newline='\n')
    print(json.dumps({k:payload[k] for k in ('checks','layers','relabel_controls','weighted_DP_controls')},indent=2))
    print('LF_SHA256',sha256(text.encode()).hexdigest())


if __name__=='__main__':main()
