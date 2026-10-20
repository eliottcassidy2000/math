#!/usr/bin/env python3
"""Independent finite controls; imports no producer code.
Reads existing class tables, but reconstructs quotient edges from literal paths.
"""
from pathlib import Path
from itertools import product, permutations
from collections import defaultdict
import json
ROOT=Path(__file__).resolve().parent
BASE=ROOT

def check(b,why):
    if not b: raise RuntimeError(why)
def orient(rows,a,b): rows[a] |= 1<<b

def triangle(rows):
    return sum(1 for a in range(len(rows)) for b in range(a+1,len(rows)) for c in range(b+1,len(rows))
               if (((rows[a]>>b)&1) and ((rows[b]>>c)&1) and ((rows[c]>>a)&1))
               or (((rows[b]>>a)&1) and ((rows[c]>>b)&1) and ((rows[a]>>c)&1)))
def paths(rows):
    n=len(rows)
    def extend(p,unused):
        if not unused:
            yield p;return
        available=unused if not p else rows[p[-1]] & unused
        while available:
            bit=available & -available;available-=bit
            yield from extend(p+(bit.bit_length()-1,),unused^bit)
    yield from extend((),(1<<n)-1)
def pairs(n):return [(i,j) for j in range(n-1,1,-1) for i in range(j-1)]
def decode(t,n):
    out=[0]*n
    for i in range(n-1): orient(out,i,i+1)
    for k,(i,j) in enumerate(pairs(n)):
        orient(out,j,i) if (t>>k)&1 else orient(out,i,j)
    return out

def encode(rows,p):
    return sum((1-((rows[p[i]]>>p[j])&1))<<k for k,(i,j) in enumerate(pairs(len(rows))))
def reflect(t,n):
    ps=pairs(n);rev={v:k for k,v in enumerate(ps)}
    return sum(((t>>k)&1)<<rev[(n-1-j,n-1-i)] for k,(i,j) in enumerate(ps))
def edge_key(t,n):
    mask=(1<<len(pairs(n)))-1;s=reflect(t,n)
    return min(t,s,t^mask,s^mask)
def connected(rows):
    for start in range(len(rows)):
        seen=1<<start;old=-1
        while old!=seen:
            old=seen
            for i in range(len(rows)):
                if (seen>>i)&1:seen|=rows[i]
        if seen!=(1<<len(rows))-1:return False
    return True

def literal_quotient():
    ans=[]
    for n in range(3,8):
        data=json.loads((ROOT/'inputs'/f'census_n{n}_representatives.json').read_text())
        stats={'n':n,'strong_nonself_nodes':0,'all_literal_paths':0,'descending_path_orbits':0,'minimum_descending_Q_edges':None,'cycle_cut_pairs_checked':0}
        for node in data['nodes']:
            rows=decode(node['representative'],n)
            if node['SC'] or not connected(rows):continue
            stats['strong_nonself_nodes']+=1
            deg=[r.bit_count() for r in rows]
            downhill={};ordinary=set();cycles={}
            for p in paths(rows):
                stats['all_literal_paths']+=1;t=encode(rows,p)
                check(t!=reflect(t,n),('unexpected blue',n,node['node'],p))
                if deg[p[0]]<=deg[p[-1]]:
                    ordinary.add(t);downhill.setdefault(edge_key(t,n),set()).add(t)
                    out=rows[:]
                    for a,b in zip(p,p[1:]):out[a]^=1<<b;out[b]^=1<<a
                    check(triangle(out)-triangle(rows)==deg[p[0]]-deg[p[-1]]-1,('literal triangle',n,p))
                if (rows[p[-1]]>>p[0])&1:
                    # Deduplicate unrooted cycles by canonical rotation only.
                    z=p.index(0);cy=p[z:]+p[:z];cycles[cy]=True
            check(len(downhill)==len(ordinary),('Q versus ordinary orbit',n,node['node']))
            check(all(len(v)==1 for v in downhill.values()),('folded downhill edge',n,node['node']))
            check(len(downhill)>=2,('double descent',n,node['node']))
            stats['descending_path_orbits']+=len(ordinary)
            stats['minimum_descending_Q_edges']=min(stats['minimum_descending_Q_edges'] or 10**9,len(downhill))
            for cy in cycles:
                cuts=[cy[j:]+cy[:j] for j in range(n)]
                rotations=[]
                for r in range(n):
                    mp={cy[i]:cy[(i+r)%n] for i in range(n)}
                    if all(((rows[a]>>b)&1)==((rows[mp[a]]>>mp[b])&1) for a in range(n) for b in range(n)):
                        rotations.append(r)
                good=[j for j,p in enumerate(cuts) if deg[p[0]]<=deg[p[-1]]]
                for j in good:
                    for k in good:
                        check((edge_key(encode(rows,cuts[j]),n)==edge_key(encode(rows,cuts[k]),n))==((k-j)%n in rotations),('cycle coset',n,node['node'],cy,j,k))
                        stats['cycle_cut_pairs_checked']+=1
        ans.append(stats);print('literal quotient',stats,flush=True)
    return ans

def two_orbit():
    results=[]
    for p in (3,5,7):
        skew=[]
        for signs in product((1,-1),repeat=(p-1)//2):skew.append({(i+1)*s%p for i,s in enumerate(signs)})
        count=construct=exceptions=0
        for A,B in product(skew,repeat=2):
            for mid in range(1<<(p-2)):
                S={0}|{i for i in range(1,p-1) if (mid>>(i-1))&1}
                rows=[0]*(2*p)
                for i in range(p):
                    for a in A:orient(rows,2*i,2*((i+a)%p))
                    for b in B:orient(rows,2*i+1,2*((i+b)%p)+1)
                    for j in range(p):orient(rows,2*i,2*j+1) if (j-i)%p in S else orient(rows,2*j+1,2*i)
                check(all((rows[i]>>((i+1)%(2*p)))&1 for i in range(2*p)),('original cycle',p,A,B,S))
                choice=next(((a,b) for a in A for b in B if (a+b)%p),None)
                if choice:
                    a,b=choice;x=next(x for x in S if (x-a-b)%p not in S);t=(x-a)%p
                    cy=tuple(2*((k*a)%p) for k in range(p))+tuple(2*((t+k*b)%p)+1 for k in range(p))
                    check(len(set(cy))==2*p and all((rows[cy[i]]>>cy[(i+1)%(2*p)])&1 for i in range(2*p)),('two block cycle',p,A,B,S))
                    cuts=[cy[1:]+cy[:1],cy[p+1:]+cy[:p+1]]
                    check(rows[cuts[0][0]].bit_count()!=rows[cuts[1][0]].bit_count(),('different orbit scores',p,A,B,S))
                    check(edge_key(encode(rows,cuts[0]),2*p)!=edge_key(encode(rows,cuts[1]),2*p),('two distinct Q edges',p,A,B,S))
                    for path in cuts:
                        check(rows[path[0]].bit_count()==rows[path[-1]].bit_count(),('equal score cut',p,path))
                        out=rows[:]
                        for u,v in zip(path,path[1:]):out[u]^=1<<v;out[v]^=1<<u
                        check(triangle(out)==triangle(rows)-1,('one triangle descent',p,path))
                    construct+=1
                else:
                    check(p==3,('exception prime',p,A,B))
                    c=next(c for c in range(p) if {(2*c-d)%p for d in S}==S)
                    mp={**{2*i:2*((i+c)%p)+1 for i in range(p)},**{2*j+1:2*((j-c)%p) for j in range(p)}}
                    check(all(((rows[u]>>v)&1)==((rows[mp[v]]>>mp[u])&1) for u in range(2*p) for v in range(2*p) if u!=v),('exception anti',A,B,S,c))
                    exceptions+=1
                count+=1
        row={'p':p,'total_systems':count,'literal_block_cycles':construct,'verified_self_converse_exceptions':exceptions}
        results.append(row);print('two orbits',row,flush=True)
    return results

def block_sumsets():
    results=[]
    for p,l in ((3,3),(3,4),(5,3)):
        subsets=[tuple(i for i in range(p) if (mask>>i)&1) for mask in range(1,1<<p)]
        total=signed=0
        for ss in product(subsets,repeat=l):
            chosen=list(ss);reverse=False
            if all(len(s)<p for s in ss) and sum(map(len,ss))*2<p*l:
                # reverse order 0,l-1,...,1; incoming sets are -(complement S).
                order=[0]+list(range(l-1,0,-1));chosen=[tuple(-d%p for d in range(p) if d not in ss[(order[k]-1)%l]) for k in range(l)];reverse=True
            else:order=list(range(l))
            dp={0:()}
            for s in chosen:dp={(v+d)%p:word+(d,) for v,word in dp.items() for d in s}
            check(len(dp)==p,('sumset coverage',p,l,ss,reverse))
            for eps in product((1,-1),repeat=l):
                ordered_eps=[eps[i] for i in order];ds=dp[sum(eps)%p];starts=[0]
                for k in range(l-1):starts.append((starts[-1]+ds[k]-ordered_eps[k])%p)
                check((starts[-1]+ds[-1]-ordered_eps[-1])%p==0,('cyclic closure',p,l,ss,eps))
                for k in range(l):
                    difference=(starts[(k+1)%l]-(starts[k]-ordered_eps[k]))%p
                    check(difference in chosen[k],('boundary',p,l,ss,eps,k))
                signed+=1
            total+=1
        row={'p':p,'orbits':l,'all_nonempty_cross_systems':total,'signed_boundary_certificates':signed}
        results.append(row);print('prime block sumsets',row,flush=True)
    return results

if __name__=='__main__':
    result={'status':'FINITE-EXACT; all checks passed','quotient':literal_quotient(),'two_orbit':two_orbit(),'prime_block':block_sumsets()}
    check(result==json.loads((ROOT/'results'/'independent_controls.json').read_text()),'Bundled result mismatch')
    print('Matched bundled independent control result.',flush=True)
