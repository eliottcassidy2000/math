#!/usr/bin/env python3
"""Parent-only trace recurrence. No (n+1)-vertex canon or automorphisms in generation."""
from collections import Counter
import json
from itertools import permutations, product

def matrix(n, code):
    a = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1,n):
            a[i][j] = (code >> k) & 1
            a[j][i] = 1-a[i][j]
            k += 1
    return a

def encode(a, order):
    return sum(a[u][v] << k for k,(u,v) in enumerate(( (order[i],order[j]) for i in range(len(order)) for j in range(i+1,len(order)) )))

def degree_orders(a):
    d = [sum(row) for row in a]
    groups = [[i for i in range(len(a)) if d[i] == s] for s in sorted(set(d))]
    for choices in product(*(permutations(g) for g in groups)):
        yield tuple(v for block in choices for v in block)

def canon(a):
    return min(encode(a,p) for p in degree_orders(a))

def extension(a,s):
    n = len(a)
    return [a[i]+[1-((s>>i)&1)] for i in range(n)] + [[(s>>i)&1 for i in range(n)]+[0]]


def need(claim,msg):
    if not claim: raise RuntimeError(msg)


def rooted_key(a,v):
    remaining=[w for w in range(len(a)) if w!=v]
    parent=[[a[i][j] for j in remaining] for i in remaining]
    # Lexicographic minimum first canonizes the deletion parent, then the mask.
    return min((encode(parent,p),sum(a[v][remaining[p[i]]] << i for i in range(len(p)))) for p in degree_orders(parent))


def produce(parents,n):
    # Canonicalize attachment masks using only n-vertex deletion parents.
    states={rooted_key(extension(matrix(n,t),s),n) for t in parents for s in range(1<<n)}
    states=sorted(states)
    idx={x:i for i,x in enumerate(states)}
    rows=[]
    for t,s in states:
        child=extension(matrix(n,t),s)
        rows.append(Counter(idx[rooted_key(child,w)] for w in range(n+1)))
    return states,rows


def rank_mod(rows,p):
    basis={}
    for row in rows:
        r={j:v%p for j,v in row.items() if v%p}
        while r:
            j=min(r)
            if j not in basis:
                inv=pow(r[j],-1,p)
                basis[j]={k:(v*inv)%p for k,v in r.items()}
                break
            scale=r[j]
            for k,v in basis[j].items():
                r[k]=(r.get(k,0)-scale*v)%p
                if not r[k]: del r[k]
    return len(basis)


def main():
    parents={0}; out=[]; example=None
    for n in range(6):
        states,rows=produce(parents,n)
        N=n+1
        trace=sum(r.get(i,0) for i,r in enumerate(rows))
        need(trace%N==0,'Nonintegral trace count')
        count=trace//N
        for i,row in enumerate(rows):
            need(sum(row.values())==N,'Bad row sum')
            square=Counter()
            for j,w in row.items():
                for k,x in rows[j].items(): square[k]+=w*x
            need(square==Counter({k:N*x for k,x in row.items()}),'Q square law')
            need(all(rows[j]==row for j in row),'Rank-one block law')
        # Independent full-child oracle, NOT generation or the next parent bank:
        audit_codes=[canon(extension(matrix(n,t),s)) for t,s in states]
        children=set(audit_codes)
        need(count==len(children),'Child-canonicalization oracle mismatch')
        for i,row in enumerate(rows):
            need(set(row)=={j for j,c in enumerate(audit_codes) if c==audit_codes[i]},'Wrong fibres')
        primes={p:rank_mod(rows,p) for p in (2,3,5,7)}
        need(all(rank==count for p,rank in primes.items() if N%p),'Good-characteristic rank mismatch')
        out.append({'child_order':N,'pointed_states':len(states),'trace':trace,'unpointed_count':count,'ranks_mod_primes':primes,'checked':'row sum, Q^2=NQ, rank-one blocks, independent child-canon oracle, good-characteristic rank'})
        if n==2:
            example=[{'parent_code':t,'mask':s,'root_degree':s.bit_count(),'root_orbit_size':rows[i][i]} for i,(t,s) in enumerate(states) if s.bit_count()==1]
        # Next bank uses one raw extension per distinct row, NOT child canon.
        representatives={}
        for (t,s),row in zip(states,rows):
            signature=tuple(sorted(row.items()))
            representatives.setdefault(signature, encode(extension(matrix(n,t),s), tuple(range(N))))
        parents=set(representatives.values())
        need(len(parents)==count,'Representative bank size mismatch')
    need([x['unpointed_count'] for x in out]==[1,1,2,4,12,56],'Known count mismatch')
    return {'status':'FINITE-EXACT','results':out,'same_parent_same_attachment_weight_hostile':example}

if __name__=='__main__': print(json.dumps(main(),indent=2))
