#!/usr/bin/env python3
"""Exact algebra, range-query, and CRT certificate; run with --verify to replay."""
from collections import Counter
from dataclasses import dataclass
from itertools import product, combinations
from functools import reduce
from pathlib import Path
import json, math, sys
sys.path.insert(0,str(Path(__file__).parent/'sources'))
import session55_verify as old
ALL=old.ALL

def meet(F,G): return frozenset(x&y for x in F for y in G)
def minimal(F): return frozenset(x for x in F if not any(y!=x and y&x==y for y in F))
def anti(F,G): return minimal(meet(F,G))
def convolution(F,G):
    out=Counter()
    for x,c in F.items():
        for y,d in G.items(): out[x&y]+=c*d
    return out

@dataclass
class Signature:
    ids: tuple
    hist: Counter
    witnesses: dict

def merge(F,G):
    if set(F.ids)&set(G.ids): raise ValueError('repeated component identity')
    H=convolution(F.hist,G.hist); W={}
    for x,a in F.witnesses.items():
        for y,b in G.witnesses.items(): W.setdefault(x&y,a+b)
    return Signature(F.ids+G.ids,H,W)
IDENTITY=Signature((),Counter({ALL:1}),{ALL:()})

def local(q,identity):
    Q=old.roots(q); H=Counter();W={}
    for c,a,b in product(Q,repeat=3):
        m=old.magic(q,c,a,b)
        if all(x in Q for x in m):
            mask=old.mask(m);H[mask]+=1;W.setdefault(mask,(m,))
    return Signature((identity,),H,W)

class DisjointSparseTable:
    def __init__(self,A,op):
        self.A=A;self.op=op;self.T=[]
        for k in range((len(A)-1).bit_length()):
            half=1<<k;row=[None]*len(A)
            for start in range(0,len(A),2*half):
                mid=min(start+half,len(A));end=min(start+2*half,len(A))
                if mid>start:
                    row[mid-1]=A[mid-1]
                    for j in range(mid-2,start-1,-1): row[j]=op(A[j],row[j+1])
                if mid<end:
                    row[mid]=A[mid]
                    for j in range(mid+1,end):row[j]=op(row[j-1],A[j])
            self.T.append(row)
    def query(self,l,r):
        assert 0<=l<=r<len(self.A)
        if l==r:return self.A[l]
        k=(l^r).bit_length()-1
        return self.op(self.T[k][l],self.T[k][r])

class SegmentTree:
    def __init__(self,A,op,identity):
        self.op=op;self.identity=identity;self.n=1<<(len(A)-1).bit_length()
        self.T=[identity for _ in range(2*self.n)]
        self.T[self.n:self.n+len(A)]=A
        for j in range(self.n-1,0,-1):self.T[j]=op(self.T[2*j],self.T[2*j+1])
    def query(self,l,r):
        l+=self.n;r+=self.n+1;left=right=self.identity
        while l<r:
            if l&1:left=self.op(left,self.T[l]);l+=1
            if r&1:r-=1;right=self.op(self.T[r],right)
            l//=2;r//=2
        return self.op(left,right)

class OverlapSparseTable:
    def __init__(self,A,op):
        self.op=op;self.T=[A]
        for k in range(1,len(A).bit_length()):
            h=1<<(k-1);prev=self.T[-1]
            self.T.append([op(prev[i],prev[i+h]) for i in range(len(A)-(1<<k)+1)])
    def query(self,l,r):
        k=(r-l+1).bit_length()-1
        return self.op(self.T[k][l],self.T[k][r-(1<<k)+1])

def fold(A,op,identity):return reduce(op,A,identity)
def hist_json(H):return {hex(k):v for k,v in sorted(H.items())}
def main():
    if not __debug__:
        raise RuntimeError('verification requires assertions; run without -O')
    checks={};families=[frozenset(x for x in range(8) if flag>>x&1) for flag in range(256)]
    for F in families:
        assert (meet(F,F)==F)==all(x&y in F for x in F for y in F)
        A=minimal(F);assert (anti(A,A)==A)==(len(A)<=1)
    checks['families_on_three_bits']=256
    for weights in product(range(3),repeat=4):
        h=Counter({i:w for i,w in enumerate(weights) if w})
        assert (convolution(h,h)==h)==(sum(weights)<=1)
    checks['integer_histograms_on_two_bits']=81
    ants=sorted(set(map(minimal,families)),key=lambda F:tuple(sorted(F)))
    for A,B,C in product(ants,repeat=3):assert anti(anti(A,B),C)==anti(A,anti(B,C))
    checks['antichain_associativity_triples']=len(ants)**3
    for F,G in product(families,repeat=2):assert minimal(meet(F,G))==anti(minimal(F),minimal(G))
    checks['compression_pairs']=len(families)**2
    for n in range(1,65):
        A=[(i,) for i in range(n)];op=lambda a,b:a+b
        dst=DisjointSparseTable(A,op);seg=SegmentTree(A,op,())
        for l in range(n):
            for r in range(l,n):assert dst.query(l,r)==seg.query(l,r)==tuple(range(l,r+1))
    checks['ordered_interval_queries']=sum(n*(n+1)//2 for n in range(1,65))
    rows=[]
    for mods in [(2,5,3),(5,7),(2,3,5,7),(4,9,5),(3,8,7)]:
        assert all(math.gcd(a,b)==1 for a,b in combinations(mods,2))
        A=[local(q,(i,q)) for i,q in enumerate(mods)]
        dst=DisjointSparseTable(A,merge);seg=SegmentTree(A,merge,IDENTITY)
        antitable=DisjointSparseTable([minimal(x.hist) for x in A],anti)
        bad=OverlapSparseTable([x.hist for x in A],convolution)
        for l in range(len(A)):
            for r in range(l,len(A)):
                S=dst.query(l,r);T=seg.query(l,r);naive=fold(A[l:r+1],merge,IDENTITY)
                assert S.ids==T.ids==naive.ids==tuple((i,mods[i]) for i in range(l,r+1))
                assert S.hist==T.hist==naive.hist
                assert antitable.query(l,r)==minimal(S.hist)
                N=math.prod(mods[l:r+1]);direct=local(N,'direct')
                assert S.hist==direct.hist
                assert set(S.hist)==set(old.independent_masks(N))
                for mask,locals_ in S.witnesses.items():
                    M=tuple(old.crt(v,mods[l:r+1]) for v in zip(*locals_))
                    assert old.valid(N,M) and old.mask(M)==mask
                b=bad.query(l,r)
                rows.append({'mods':mods,'interval':[l,r],'modulus':N,'histogram':hist_json(S.hist),'total_arrays':sum(S.hist.values()),'minimal_masks':list(map(hex,sorted(minimal(S.hist)))),'overlap_total':sum(b.values()),'overlap_false_distinct':0 in b and 0 not in S.hist,'distinct_local_witness':S.witnesses.get(0)})
    hostile=next(x for x in rows if x['mods']==(2,5,3) and x['interval']==[0,2])
    assert hostile['overlap_false_distinct']
    f5=local(5,'F5');f7=local(7,'F7');assert minimal(f5.hist)==frozenset((0x61089888,0x8a054250))
    assert 0 not in f5.hist and 0 in convolution(f5.hist,f5.hist)
    guard=False
    try:merge(f5,f5)
    except ValueError:guard=True
    assert guard
    overlap_guard=False
    try:
        OverlapSparseTable([local(q,(i,q)) for i,q in enumerate((2,5,3))],merge).query(0,2)
    except ValueError:overlap_guard=True
    assert overlap_guard
    assert 0 in merge(f5,f7).hist and 0 not in local(25,'Z/25Z').hist
    fixed=[0x61089888,ALL,0x8a054250,7,3,ALL,6]
    good=OverlapSparseTable(fixed,lambda x,y:x&y)
    for l in range(len(fixed)):
        for r in range(l,len(fixed)):assert good.query(l,r)==fold(fixed[l:r+1],lambda x,y:x&y,ALL)
    checks.update(crt_interval_queries=len(rows),duplicate_identity_rejected=guard,fixed_mask_queries=28,overlap_query_identity_rejected=overlap_guard)
    result={'schema':'overlap-safe-signature-v1','checks':checks,'hostile':hostile,'rows':rows,'pairs':old.PAIRS,'false_overlap_local_arrays':[local(2,'left').witnesses[ALL][0],f5.witnesses[0x61089888][0],f5.witnesses[0x8a054250][0],local(3,'right').witnesses[ALL][0]],'false_overlap_component_ids':['index0_mod2','index1_mod5','index1_mod5','index2_mod3']}
    p=Path(__file__).with_name('certificate.json')
    if '--verify' in sys.argv:assert json.loads(p.read_text())==json.loads(json.dumps(result))
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(checks,indent=2));print('PASS')
if __name__=='__main__':main()
