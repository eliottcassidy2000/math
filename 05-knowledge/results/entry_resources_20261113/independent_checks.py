"""Separate signed-cut automaton, subset partitions and contact checks."""
from pathlib import Path
from itertools import combinations,permutations
from functools import lru_cache
import numpy as np
from numba import njit
import json,gzip,random,time,sys
HERE=Path(__file__).resolve().parent

def need(p,m):
    if not p:raise RuntimeError(m)
def read(name):return json.loads(gzip.decompress((HERE/(name+'.json.gz')).read_bytes()))

@njit(cache=True)
def signed_cut(rows,O):
    n=len(rows);N=1<<n;a=np.zeros((N,n),dtype=np.int64);b=np.zeros((N,n),dtype=np.int64)
    for mask in range(1,N):
        for v in range(n):
            if not(mask>>v&1):continue
            rest=mask^(1<<v)
            if not rest:a[mask,v]=1;continue
            for u in range(n):
                if rows[u]>>v&1:a[mask,v]+=a[rest,u];b[mask,v]+=b[rest,u]
                if (O>>u&1) and (O>>v&1):b[mask,v]+=a[rest,u]
    return np.sum(a[N-1])-np.sum(b[N-1])

def build(rr,O,gg,ts,s):
    # Construct by an explicit Boolean adjacency matrix, separate from attach.
    n=len(rr);q=len(gg);a=[[False]*(n+q) for _ in range(n+q)]
    for i in range(n):
        for j in range(n):a[i][j]=bool(rr[i]>>j&1)
    for i in range(q):
        for j in range(q):a[n+i][n+j]=bool(gg[i]>>j&1)
        for j in range(n):a[n+i][j]=bool(ts[i]>>(1 if O>>j&1 else 0)&1);a[j][n+i]=not a[n+i][j]
    return np.array([sum(1<<j for j,x in enumerate(row) if x) for row in a],dtype=np.int64),O|(s<<n)

def gadgets():
    d=read('uniform_gadgets');seen=set();counts={};total=0
    for group in d['rejection_groups']:
        rr,O=group['rows'],group['selected'];need(signed_cut(np.array(rr,dtype=np.int64),O)==1,'negative witness starts at a unit seed')
        for gg,ts,s,want in group['rejected']:
            key=(tuple(gg),tuple(ts),s);need(key not in seen,'gadget partition unique');seen.add(key)
            rows,mark=build(rr,O,gg,ts,s);actual=int(signed_cut(rows,mark))
            need(actual==want and actual!=1,'independent two-phase automaton rejects gadget')
            counts[len(gg)]=counts.get(len(gg),0)+1;total+=1
    for r in d['survivors']:
        key=(tuple(r['core']),tuple(r['types']),r['selected_new']);need(key not in seen,'survivor is new');seen.add(key)
        need(r['ordinary_transitive_padding'],'all survivors have a uniform all-order proof')
        counts[len(r['core'])]=counts.get(len(r['core']),0)+1
    need(counts=={int(k):v for k,v in d['universe_by_order'].items()},'complete finite gadget universe partition')
    return total

def partition_oracle(adj):
    n=len(adj);N=1<<n;ends=[0]*N
    def visit(path,mask):
        ends[mask]|=(1<<path[0])|(1<<path[-1])
        for v in range(n):
            if not(mask>>v&1) and adj[path[-1]]>>v&1:visit(path+[v],mask|(1<<v))
    for v in range(n):visit([v],1<<v)
    K=[0]*N;E=[0]*N
    for mask in range(1,N):
        anchor=mask&-mask;best=n+1;possible=0;sub=mask
        while sub:
            if sub&anchor and ends[sub]:
                cost=1+K[mask^sub]
                if cost<best:best=cost;possible=ends[sub]|E[mask^sub]
                elif cost==best:possible|=ends[sub]|E[mask^sub]
            sub=(sub-1)&mask
        K[mask]=best;E[mask]=possible
    return K,E

def small_covers():
    from square_covers import cover_dp
    rng=random.Random(1311);checks=0
    for n in range(1,8):
        pairs=list(combinations(range(n),2));codes=range(1<<len(pairs)) if n<=5 else [rng.randrange(1<<len(pairs)) for _ in range(100)]
        for code in codes:
            adj=[0]*n
            for i,(a,b) in enumerate(pairs):
                if code>>i&1:adj[a]|=1<<b;adj[b]|=1<<a
            a,b=partition_oracle(adj);k,e=cover_dp(np.array(adj,dtype=np.uint32))
            need(a==k.tolist() and b==e.tolist(),'independent minimum-cover partition and endpoint support')
            checks+=1
    return checks

def ports():
    data=read('port_gadgets');count=0
    for r in data['failures']:
        rr,O,p=r['rows'],r['selected'],r['port'];fR,fO,edge,sel=r['rule'];n=len(rr)
        rows=list(rr)+[0]
        for v in range(n):
            old=(rr[p]>>v)&1
            new=edge if v==p else [0,1,old,1-old][fO if O>>v&1 else fR]
            if new:rows[n]|=1<<v
            else:rows[v]|=1<<n
        d=int(signed_cut(np.array(rows,dtype=np.int64),O|(sel<<n)))
        need(d==r['after_defect'] and d!=1,'stable-port rule counterexample')
        count+=1
    need(count==60 and len(data['survivors'])==4,'all64 port rules accounted for')
    return count

def run():
    t=time.time();n=gadgets();print('independent rejected gadgets',n,'seconds',round(time.time()-t,1),flush=True)
    p=small_covers();print('independent graphs',p,flush=True)
    # Independent signed-cut automaton checks all retained markings on the
    # moving-root chain, plus the original fixed marking at every order.
    chain=0
    for r in read('moving_root')['records']:
        rows=np.array(r['rows'],dtype=np.int64);O=r['fixed_selected']
        need(int(signed_cut(rows,O))==r['fixed_state'][1],'moving-root selected defect')
        for mark in r['all_unit_selected_masks']:need(signed_cut(rows,mark)==1,'all reported unit markings')
        chain+=1
    pc=ports()
    out={'status':'PASS','all_rejection_certificates':n,'port_rule_counterexamples':pc,'independent_path_cover_graphs':p,'moving_root_states':chain,'methods':['Two-phase spanning-word automaton with exactly one selected reset; independent of endpoint convolution.','All simple paths plus minimum set partitions; independent of append/reset minimum-cover recursion.','No assertion statements; normal and optimized receipts are compared separately.']}
    tag='' if __debug__ else '.optimized';(HERE/('independent_checks'+tag+'.json.gz')).write_bytes(gzip.compress((json.dumps(out,separators=(',',':'))+'\n').encode(),mtime=0))
if __name__=='__main__':run()
