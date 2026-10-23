#!/usr/bin/env python3
"""Standalone exhaustive parity-DSU/NAE realization referee through n=5."""
import itertools,json,math
from collections import Counter,defaultdict
from pathlib import Path

OUT=Path(__file__).resolve().parent

def check(x,label):
    if not x:raise RuntimeError(label)

class DSU:
    def __init__(self,m):self.p=list(range(m));self.x=[0]*m
    def find(self,i):
        if self.p[i]!=i:
            r,x=self.find(self.p[i]);self.x[i]^=x;self.p[i]=r
        return self.p[i],self.x[i]
    def join(self,a,b,w):
        ra,xa=self.find(a);rb,xb=self.find(b)
        if ra==rb:return (xa^xb)==w
        self.p[ra]=rb;self.x[ra]=xa^xb^w;return True

def triangle_data(n):
    edges=list(itertools.combinations(range(n),2));ix={e:i for i,e in enumerate(edges)}
    return edges,[(ix[i,j],ix[j,k],ix[i,k]) for i,j,k in itertools.combinations(range(n),3)]

def direct_key(mask,tris):
    key=0
    for b,(a,c,d) in enumerate(tris):
        if ((mask>>a)&1)==((mask>>c)&1) and ((mask>>a)&1)!=((mask>>d)&1):key|=1<<b
    return key

def system(n,h):
    edges,tris=triangle_data(n);dsu=DSU(len(edges))
    for b,(a,c,d) in enumerate(tris):
        if h>>b&1:
            if not dsu.join(a,c,0) or not dsu.join(a,d,1):return None
    roots=sorted({dsu.find(e)[0] for e in range(len(edges))});idx={r:i for i,r in enumerate(roots)}
    literal=[(idx[r],s) for r,s in map(dsu.find,range(len(edges)))];clauses=[]
    for b,(a,c,d) in enumerate(tris):
        if not(h>>b&1):clauses.append((literal[a],literal[c],(literal[d][0],literal[d][1]^1)))
    solutions=[]
    for phases in range(1<<len(roots)):
        if all(len({(phases>>r&1)^s for r,s in cl})>1 for cl in clauses):
            mask=sum(((phases>>r&1)^s)<<e for e,(r,s) in enumerate(literal));solutions.append(mask)
    types=Counter()
    for cl in clauses:
        if any((r,s^1) in cl for r,s in cl):typ='tautology'
        elif len(set(cl))==1:typ='forced_false'
        elif len({r for r,s in cl})==2:typ='binary_xor'
        else:typ='ternary_NAE'
        types[typ]+=1
    return {'phases':len(roots),'solutions':solutions,'literal_map':literal,'clauses':clauses,'clause_types':dict(types)}

def indecomposable(n,mask):
    def arc(i,j):
        if i==j:return False
        a,b=sorted((i,j));k=edgeindex[a,b]
        return bool(mask>>k&1)==(i<j)
    edgeindex={e:i for i,e in enumerate(itertools.combinations(range(n),2))}
    for size in range(2,n):
        for vertices in itertools.combinations(range(n),size):
            s=set(vertices)
            if all(len({arc(v,w) for v in s})==1 for w in range(n) if w not in s):return False
    return True

def main():
    report={};first_unsat=None;first_prime_multiphase=None
    for n in range(1,6):
        edges,tris=triangle_data(n);expected=defaultdict(list)
        for mask in range(1<<len(edges)):expected[direct_key(mask,tris)].append(mask)
        counts=Counter();phase_all=Counter();phase_valid=Counter();savings=Counter();prime_phases=Counter();unsat_examples=[]
        for h in range(1<<len(tris)):
            z=system(n,h)
            if z is None:
                check(h not in expected,('inconsistent but realized',n,h));counts['xor_inconsistent']+=1;continue
            check(z['solutions']==expected.get(h,[]) or sorted(z['solutions'])==expected.get(h,[]),('exact realization mismatch',n,h))
            phase_all[z['phases']]+=1
            if not z['solutions']:
                counts['xor_consistent_nae_unsat']+=1
                unsat_examples.append((h.bit_count(),h,z))
            else:
                counts['realizable']+=1;phase_valid[z['phases']]+=1;savings[len(edges)-z['phases']]+=1
                if n>=2:
                    full=(1<<len(edges))-1;ss=set(z['solutions'])
                    check(all((x^full) in ss for x in ss),'complement pairing')
                    check(len(ss)%2==0,'even fibre')
                if indecomposable(n,z['solutions'][0]):
                    prime_phases[z['phases']]+=1
                    if n>=2:check(len(z['solutions'])==2,('indecomposable dual fibre',n,h))
                    if z['phases']>1 and first_prime_multiphase is None:first_prime_multiphase={'n':n,'hypergraph_mask':h,**z}
        if unsat_examples and first_unsat is None:
            _,h,z=min(unsat_examples,key=lambda t:(t[0],t[1]));first_unsat={'n':n,'hypergraph_mask':h,'positive_triples':[list(tr) for b,tr in enumerate(itertools.combinations(range(n),3)) if h>>b&1],**z}
        # Empty hypergraph has one labelled tournament per linear ordering.
        empty=system(n,0);check(len(empty['solutions'])==math.factorial(n),'transitive fibre n!')
        check(empty['phases']==len(edges),'transitive no xor compression')
        report[n]={'hypergraphs':1<<len(tris),'counts':dict(counts),'phase_histogram_consistent':dict(sorted(phase_all.items())),'phase_histogram_realizable':dict(sorted(phase_valid.items())),'saved_edge_bits_realizable':dict(sorted(savings.items())),'indecomposable_phase_histogram':dict(sorted(prime_phases.items()))}
        print(n,report[n],flush=True)
    report['first_consistent_unsat']=first_unsat
    report['first_indecomposable_multiphase']=first_prime_multiphase
    (OUT/'realization_phase_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('first consistent unsat',first_unsat)
    print('first indecomposable multiphase',first_prime_multiphase)
    print('ALL REALIZATION PHASE CHECKS PASSED')

if __name__=='__main__':main()
