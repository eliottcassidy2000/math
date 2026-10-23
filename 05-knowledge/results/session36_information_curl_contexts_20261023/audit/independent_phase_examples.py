#!/usr/bin/env python3
"""Small phase examples using relation-graph BFS rather than DSU."""
from independent_audit import need,mat
from independent_realization_phases import c3mask,positive_rank
import itertools,json
from collections import Counter
from pathlib import Path
OUT=Path(__file__).resolve().parent

def circulant(n,res):return [[0 if i==j else (1 if (j-i)%n in res else -1) for j in range(n)] for i in range(n)]

def solve(z):
    n=len(z);pairs=list(itertools.combinations(range(n),2));ix={p:i for i,p in enumerate(pairs)};triples=list(itertools.combinations(range(n),3));target=c3mask(z);m=len(pairs);adj=[[] for _ in pairs]
    for b,(i,j,k) in enumerate(triples):
        if target>>b&1:
            for e,f,w in ((ix[i,j],ix[j,k],0),(ix[i,j],ix[i,k],1)):
                adj[e].append((f,w));adj[f].append((e,w))
    labels={};p=0
    for e in range(m):
        if e in labels:continue
        labels[e]=(p,0);todo=[e]
        while todo:
            v=todo.pop()
            for w,bit in adj[v]:
                value=(p,labels[v][1]^bit)
                if w not in labels:labels[w]=value;todo.append(w)
                else:need(labels[w]==value,'parity inconsistency')
        p+=1
    bad,free=positive_rank(n,target);need(not bad and free==p,'phase rank mismatch')
    clauses=Counter();solutions=0
    for b,(i,j,k) in enumerate(triples):
        if target>>b&1:continue
        ls=[labels[ix[i,j]],labels[ix[j,k]],(labels[ix[i,k]][0],labels[ix[i,k]][1]^1)]
        values=set(ls)
        if any((v,e^1) in values for v,e in ls):kind='tautology'
        elif len(values)==1:kind='contradiction'
        elif len(values)==2:kind='binary_xor'
        else:kind='ternary_NAE'
        clauses[kind]+=1
    for phase in range(1<<p):
        code=sum((((phase>>v)&1)^s)<<e for e,(v,s) in labels.items())
        solutions+=c3mask(mat(n,code))==target
    prime=not any(1<len(S)<n and all(len({z[v][w] for v in S})==1 for w in range(n) if w not in S) for k in range(2,n) for S in itertools.combinations(range(n),k))
    return dict(order=n,edge_variables=m,positive_triangles=target.bit_count(),positive_phase_components=p,realizations=solutions,indecomposable=prime,clause_types=dict(clauses))

def main():
    W=mat(5,(1<<10)-1)
    for i in range(4):W[4][i]=1 if i%2==0 else -1;W[i][4]=-W[4][i]
    examples={'C5':solve(circulant(5,{1,2})),'Paley7':solve(circulant(7,{1,2,4})),'W5':solve(W)}
    producer=json.loads(Path('context/REALIZATION_PHASES_EXAMPLES.json').read_text())
    need(examples==producer,('example receipt mismatch',examples,producer))
    report=dict(status='PASS',examples=examples,definitions={'C5':'cyclic difference set {1,2} modulo5','Paley7':'cyclic difference set {1,2,4} modulo7','W5':'vertices0..3 transitive; vertex4 beats0,2 and loses to1,3'})
    (OUT/'phase_examples_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
