"""Cheap probes for actual selected strong-atom substitutions."""
from pathlib import Path
import importlib.util
import json
from itertools import combinations, permutations

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def load(name,relative):
    spec=importlib.util.spec_from_file_location(name,ROOT/relative)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
oracle=load("root_switch","05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py")
PAIR=((66,12,41,17,39,75,30),(66,20,9,51,37,71,30))
def edge(t,i,j):return (t[i]>>j)&1
def tau(t):
    return tuple(sum((int(j==i+1) if abs(i-j)==1 else 1-edge(t,i,j))<<j
                     for j in range(len(t)) if i!=j) for i in range(len(t)))
def blow(t,p,h=(2,4,1)):
    labels=[(i,j) for i in range(len(t)) for j in range(len(h) if i==p else 1)]
    return tuple(sum((edge(h,x,y) if a==b else edge(t,a,b))<<k
                     for k,(b,y) in enumerate(labels) if (a,x)!=(b,y)) for a,x in labels)
def dent(t,p):
    n=len(t);return tuple(row|((i!=p)<<n) for i,row in enumerate(t))+(1<<p,)
def sink(t):
    return tuple(row|(1<<len(t)) for row in t)+(0,)
def induced(t,vertices):
    return tuple(sum(edge(t,i,j)<<k for k,j in enumerate(vertices)) for i in vertices)
def relabel(t,order):return induced(t,order)
def presentation(n,mask):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            bit=1 if j==i+1 else ((mask>>k)&1)
            if j>i+1:k+=1
            a,b=(i,j) if bit else (j,i);rows[a]|=1<<b
    return tuple(rows)
def register(t,ports):
    out=sink(t)
    for p in ports:out=dent(out,p)
    return out
def decode_register(t,r):
    """Recover T and ordered old ports from an unlabelled r+1 vertex register.

    Blocks report the original input vertices representing each recovered vertex;
    the terminal C3 ambiguity is collapsed, not resolved by an arbitrary label.
    """
    if r<1:raise ValueError('positive register arity')
    current=tuple(t);labels=list(range(len(t)));reverse_ports=[]
    for _ in range(r,1,-1):
        low=[i for i,row in enumerate(current) if row.bit_count()==1]
        if len(low)==1:z=low[0]
        elif len(low)==2:z=next(i for i in low if edge(current,i,low[1-low.index(i)]))
        else:raise RuntimeError('invalid register tail')
        q=next(j for j in range(len(current)) if edge(current,z,j))
        reverse_ports.append(labels[q])
        keep=[i for i in range(len(current)) if i!=z]
        current=induced(current,keep);labels=[labels[i] for i in keep]
    low=[i for i,row in enumerate(current) if row.bit_count()==1]
    if len(low)==2:
        z=next(i for i in low if not edge(current,i,low[1-low.index(i)]))
        q=next(j for j in range(len(current)) if edge(current,z,j))
        first=labels[q]
        keep=[i for i in range(len(current)) if i!=z]
        reduced=induced(current,keep);newlabels=[labels[i] for i in keep]
        ends=[i for i,row in enumerate(reduced) if row==0]
        if len(ends)!=1:raise RuntimeError('missing unique temporary sink')
        keep=[i for i in range(len(reduced)) if i!=ends[0]]
        recovered=induced(reduced,keep);blocks=[(newlabels[i],) for i in keep]
    elif len(low)==3:
        if any(sum(edge(current,i,j) for j in low)!=1 for i in low):raise RuntimeError('not a terminal triangle')
        others=[i for i in range(len(current)) if i not in low]
        if not all(edge(current,i,j) for i in others for j in low):raise RuntimeError('triangle not terminal')
        blocks=[(labels[i],) for i in others]+[tuple(labels[i] for i in low)]
        representatives=others+[low[0]];recovered=induced(current,representatives)
        first=labels[low[0]]
    else:raise RuntimeError('invalid initial register')
    owner={v:i for i,block in enumerate(blocks) for v in block}
    ports=(owner[first],)+tuple(owner[v] for v in reversed(reverse_ports))
    return recovered,ports,blocks
def modules(t):
    n=len(t);out=[]
    for k in range(2,n):
        for s in combinations(range(n),k):
            if all(len({edge(t,v,u) for u in s})==1 for v in range(n) if v not in s):out.append(s)
    return out

def main():
    records=[]
    for rev in (False,True):
        a,b=(tuple(map(tau,PAIR)) if rev else PAIR)
        records.append({'tau':rev,'aut':[len(oracle.isomorphisms(t,t,True)) for t in (a,b)],
                        'modules':[modules(t) for t in (a,b)],
                        'C3':[{'port':p,'isomorphic':bool(oracle.isomorphisms(blow(a,p),blow(b,p))),
                               'H':[oracle.hamilton(blow(t,p)) for t in (a,b)]} for p in range(7)],
                        'dent':[{'port':p,'isomorphic':bool(oracle.isomorphisms(dent(a,p),dent(b,p))),
                                 'H':[oracle.hamilton(dent(t,p)) for t in (a,b)]} for p in range(6)]})
    print(json.dumps(records,indent=2))
if __name__=='__main__':main()
