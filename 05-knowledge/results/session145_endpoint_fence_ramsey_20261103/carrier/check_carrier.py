#!/usr/bin/env python3
"""Exact self-contained checks. No imports from inherited or audit programs."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
import argparse, json

CHECKS=0

def require(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:
        raise RuntimeError(msg)

def inspect(n,edges):
    adj=[[] for _ in range(n)]
    for u,v in edges:
        if not (0<=u<v<n): raise ValueError('edges must be canonical')
        adj[u].append(v);adj[v].append(u)
    if any(len(a)>2 for a in adj): return None
    seen=set();isolates=[];paths=[];sign=1
    for v in range(n):
        if v in seen:continue
        if not adj[v]:isolates.append(v);seen.add(v);continue
        stack=[v];component={v}
        while stack:
            z=stack.pop()
            for w in adj[z]:
                if w not in component:component.add(w);stack.append(w)
        ends=sorted(z for z in component if len(adj[z])==1)
        if len(ends)!=2:return None
        u,w=ends;prev=-1;cur=u;walk=[u]
        while cur!=w:
            nxt=next(z for z in adj[cur] if z!=prev)
            sign*=1 if cur<nxt else -1
            prev,cur=cur,nxt;walk.append(cur)
        seen.update(component);paths.append(tuple(walk))
    paths.sort(key=lambda p:(p[0],p[-1]))
    signature=(tuple(isolates),tuple((p[0],p[-1],(len(p)-2)%2) for p in paths))
    blocks=(tuple(isolates),tuple(sorted((tuple(sorted(p)),p[0],p[-1]) for p in paths)))
    return {'signature':signature,'blocks':blocks,'sign':sign,'paths':paths,'degrees':tuple(map(len,adj)), 'isolates':isolates}

def coefficient(n,edges):
    if n==0:return Fraction(1)
    d=inspect(n,edges)
    if d is None or any(len(p)%2==0 for p in d['paths']):return Fraction(0)
    r=len(d['paths']);c=r+len(d['isolates'])
    return Fraction(d['sign']*2**r*factorial(c),2**(c-1))

def old_count(n):
    a=[1]
    for m in range(1,n+1):
        a.append(a[m-1]+sum(comb(m-1,k-1)*comb(k,2)*a[m-k] for k in range(2,m+1)))
    return a[-1]

def new_count(n):
    result=1
    for s in range(n+1):
        for r in range(1,(n-s)//2+1):
            t=n-s-2*r
            ways=sum(comb(r,j) for j in range(min(r,t)+1) if (t-j)%2==0)
            result+=factorial(n)//(factorial(s)*factorial(t)*2**r*factorial(r))*ways
    return result

def catalog(n):
    pairs=list(combinations(range(n),2));groups=defaultdict(list);old=set();forest_count=0
    for mask in range(1<<len(pairs)):
        edges=tuple(e for j,e in enumerate(pairs) if mask>>j&1)
        d=inspect(n,edges)
        if d is None:continue
        forest_count+=1;old.add(d['blocks'])
        groups[d['signature']].append({'mask':mask,'edges':edges,'sign':d['sign'],'details':d})
    return groups,old,forest_count

def direct_hamilton_coefficients(n):
    """Multiply each word's factors directly; coefficients scaled by 2^(n-1)."""
    pairs=list(combinations(range(n),2));index={e:i for i,e in enumerate(pairs)};out=Counter()
    for word in permutations(range(n)):
        factors=[]
        for u,v in zip(word,word[1:]):
            factors.append((1<<index[tuple(sorted((u,v)))],2 if u<v else -2))
        poly={0:1}
        for bit,sgn in factors:
            nxt=Counter(poly)
            for mask,c in poly.items():nxt[mask|bit]+=sgn*c
            poly=dict(nxt)
        out.update(poly)
    return out

def separator(n,A,B):
    """Construct exterior G with zero/nonzero split, without search."""
    da=A['details'];db=B['details'];G=[];N=n
    def add_edge(u,v):G.append(tuple(sorted((u,v))))
    def leaf(v):
        nonlocal N
        add_edge(v,N);N+=1
    difference=next((v for v in range(n) if da['degrees'][v]!=db['degrees'][v]),None)
    if difference is not None:
        v=difference;target=A if da['degrees'][v]<db['degrees'][v] else B
        for _ in range(2-target['details']['degrees'][v]):leaf(v)
        mechanism='degree'
    elif {(p[0],p[-1]) for p in da['paths']} != {(p[0],p[-1]) for p in db['paths']}:
        u,v=next(iter({(p[0],p[-1]) for p in da['paths']}-{(p[0],p[-1]) for p in db['paths']}))
        add_edge(u,N);add_edge(v,N);N+=1;target=B;mechanism='matching'
    else:
        target=A;mechanism='parity'
    merged=inspect(N,target['edges']+tuple(G))
    require(merged is not None,'separator target must remain a forest')
    for p in merged['paths']:
        if len(p)%2==0:leaf(p[0])
    a=coefficient(N,A['edges']+tuple(G));b=coefficient(N,B['edges']+tuple(G))
    require(bool(a)!=bool(b),'separator did not split coefficient support')
    require(all(u>=n or v>=n for u,v in G),'separator edited an old edge')
    return mechanism,N-n

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output');ap.add_argument('--reps');args=ap.parse_args()
    counts=[];catalogs={};formula_cases=0
    for n in range(7):
        groups,old,forests=catalog(n);catalogs[n]=groups
        require(len(groups)==new_count(n),f'new count {n}')
        require(len(old)==old_count(n),f'old count {n}')
        old_group_counts=Counter(len({r['details']['blocks'] for r in rows}) for rows in groups.values())
        counts.append({'n':n,'forests':forests,'old_coordinates':len(old),'new_coordinates':len(groups),'old_class_size_histogram':dict(sorted(old_group_counts.items()))})
        if n:
            direct=direct_hamilton_coefficients(n);pairs=list(combinations(range(n),2))
            for mask in range(1<<len(pairs)):
                edges=tuple(e for i,e in enumerate(pairs) if mask>>i&1)
                require(coefficient(n,edges)*2**(n-1)==direct.get(mask,0),'full Hamilton coefficient formula')
                formula_cases+=1
    reps=[rows[0] for _,rows in sorted(catalogs[6].items())]
    if args.reps:
        Path(args.reps).write_text(json.dumps({'n':6,'edge_order':list(combinations(range(6),2)),'representatives':[{'mask':r['mask'],'edges':r['edges'],'sign':r['sign'],'signature':r['details']['signature']} for r in reps]},indent=2)+'\n')
    mechanisms=Counter();max_new=0;separator_pairs=0
    for n in range(2,7):
        rs=[rows[0] for rows in catalogs[n].values()]
        for A,B in combinations(rs,2):
            mechanism,k=separator(n,A,B);mechanisms[mechanism]+=1;max_new=max(max_new,k);separator_pairs+=1
    # Exhaustively inspect all two-new-vertex exterior monomials for n<=4.
    congruence_cases=0
    for n in range(1,5):
        N=n+2;exterior=[e for e in combinations(range(N),2) if e[1]>=n]
        for mask in range(1<<len(exterior)):
            G=tuple(e for j,e in enumerate(exterior) if mask>>j&1)
            for rows in catalogs[n].values():
                ref=coefficient(N,rows[0]['edges']+G)*rows[0]['sign']
                for r in rows:
                    require(coefficient(N,r['edges']+G)*r['sign']==ref,'signed signature congruence')
                    congruence_cases+=1
    # First strict collision and the exact old-coordinate expectation.
    bymask={r['mask']:r for rows in catalogs[6].values() for r in rows}
    A,B=bymask[10376],bymask[5392]
    require(A['details']['signature']==B['details']['signature'],'first-collision signature')
    require(A['details']['blocks']!=B['details']['blocks'],'first-collision old block distinction')
    require(A['sign']==B['sign']==1,'first-collision signs')
    total_plus=Fraction(0);total_minus=Fraction(0);moment_plus=Fraction(0);moment_minus=Fraction(0);weights=set()
    for bits in range(32768):
        x=(-1)**((bits&10376).bit_count());z=(-1)**((bits&5392).bit_count())
        wp=Fraction(4+x-z,4*32768);wm=Fraction(4-x+z,4*32768)
        require(wp>0 and wm>0,'strict positivity')
        weights.update((str(wp),str(wm)));total_plus+=wp;total_minus+=wm
        moment_plus+=wp*Fraction(x,16);moment_minus+=wm*Fraction(x,16)
    require(total_plus==total_minus==1,'law normalization')
    require(moment_plus==Fraction(1,64) and moment_minus==-Fraction(1,64),'old-coordinate mismatch')
    # Hostiles: a star is acyclic but forbidden, and one odd component is fatal.
    require(coefficient(4,((0,1),(0,2),(0,3)))==0,'star hostile')
    require(coefficient(3,((0,1),(0,2),(1,2)))==0,'cycle hostile')
    require(coefficient(2,((0,1),))==0,'odd path hostile')
    require(coefficient(3,((0,1),(1,2)))==2,'even path positive control')
    require(new_count(0)==1 and new_count(1)==1,'r0 interior boundary')
    require([new_count(n) for n in range(3,7)]==[7,28,126,571],'requested dimension boundary')
    receipt={'status':'PASS','checks':CHECKS,'universe':'all labelled old graph edge sets through n=6; all Hamilton word coefficient expansions through n=6; all distinct signature pairs through n=6; all exterior subsets with two new vertices through old n=4','counts':counts,'full_hamilton_coefficients_checked':formula_cases,'separator_pairs':separator_pairs,'separator_mechanisms':dict(sorted(mechanisms.items())),'maximum_new_vertices_for_separator':max_new,'exhaustive_exterior_congruence_cases':congruence_cases,'first_collision':{'masks':[10376,5392],'law_weights':sorted(weights),'old_coordinate_expectations':[str(moment_plus),str(moment_minus)]},'dimension_counts_through_12':[{'n':n,'b':old_count(n),'c':new_count(n)} for n in range(13)],'scope':'Coordinate count and closed upper bound only; no all-order minimum rank claimed.'}
    out=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(out)
    else:print(out,end='')
if __name__=='__main__':main()
