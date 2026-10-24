#!/usr/bin/env python3
"""Source-index-preserving balanced charge decomposition; no native-supply inference."""
from itertools import product
from random import Random
from pathlib import Path
import json,subprocess

CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def add(a,b):return a[0]^b[0],(a[1]+b[1])%128
def neg(a):return a[0],(-a[1])%128
def total(xs):
    z=(0,0)
    for x in xs:z=add(z,x)
    return z
def extract(seq,D=132):
    """A nonempty zero-sum subsequence among the first D, if one exists."""
    bank={}
    for j,g in enumerate(seq[:D]):
        if g==(0,0):return (j,)
        if neg(g) in bank:return bank[neg(g)]+(j,)
        updates={add(h,g):w+(j,) for h,w in bank.items()}
        updates.setdefault(g,(j,))
        for h,w in updates.items():bank.setdefault(h,w)
    raise RuntimeError('zero-sum promise or Davenport bound failed')
def split(A,B):
    need(len(A)==len(B) and total(A)==total(B),'equal cardinality and exact group charge')
    todo=list(range(len(A)));out=[]
    while todo:
        d=[add(A[i],neg(B[i])) for i in todo]
        chosen=extract(d);block=[todo[j] for j in chosen]
        need(0<len(block)<=132,'bounded atomic paired support')
        need(total(A[i] for i in block)==total(B[i] for i in block),'each block is charge neutral')
        out.append(block);used=set(block);todo=[i for i in todo if i not in used]
    need(sorted(i for c in out for i in c)==list(range(len(A))),'every indexed pair used exactly once')
    return out
def main():
    rng=Random(264132);runs=[]
    for k in (1,2,8,32,132,133,300,1000):
        for trial in range(3):
            d=[(rng.randrange(16),rng.randrange(128)) for _ in range(k-1)]
            d.append(neg(total(d)));B=[(rng.randrange(16),rng.randrange(128)) for _ in range(k)]
            A=[add(b,g) for b,g in zip(B,d)];blocks=split(A,B)
            runs.append({'pairs':k,'blocks':len(blocks),'max_slots':2*max(map(len,blocks))})
    # The standard extremal p-group sequence. Exactly all 132 terms sum to zero.
    extreme=[(1<<j,0) for j in range(4)]+[(0,1)]*127+[(15,1)]
    need(len(extreme)==132 and total(extreme)==(0,0),'extremal full relation')
    witness=extract(extreme)
    need(len(witness)==132,'no proper zero-sum relation in extremal sequence')
    need(len(split(extreme,[(0,0)]*132)[0])==132,'sharp abstract 264-slot packet')
    # Independent subset multiplicity DP, distinct from extract's witness DP.
    counts={(0,0):1}
    for g in extreme:
        nxt=counts.copy()
        for h,n in counts.items():
            q=add(h,g);nxt[q]=nxt.get(q,0)+n
        counts=nxt
    need(counts[(0,0)]==2,'only empty and full indexed subsets are zero')
    # Actual fixed-source labels: all endpoint classes in a ten-option subinventory.
    p='0f9dec9c3e1d:05-knowledge/results/session45_injective_multisource_moves_20261023/balanced-circuits/pool.json'
    pool=json.loads(subprocess.check_output(['git','show',p],cwd=Path(__file__).resolve().parents[3]))
    labels=[tuple(v) for s,v in sorted(pool['labels'].items(),key=lambda z:int(z[0]))][:10]
    buckets={};native_pairs=0
    for bits in range(1<<len(labels)):
        inds=[i for i in range(len(labels)) if bits>>i&1];key=(len(inds),total(labels[i] for i in inds))
        for old in buckets.get(key,[]):
            plus=[labels[i] for i in range(len(labels)) if bits>>i&1 and not old>>i&1]
            minus=[labels[i] for i in range(len(labels)) if old>>i&1 and not bits>>i&1]
            split(plus,minus);native_pairs+=1
        buckets.setdefault(key,[]).append(bits)
    return {'status':'FINITE-EXACT PASS; uniform theorem proved separately','group':[2,2,2,2,128],
      'Davenport':132,'abstract_sharp_atomic_slots':264,'extremal_zero_subsets':counts[(0,0)],
      'actual_ten_option_endpoint_pairs':native_pairs,'seeded_decompositions':runs,'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
