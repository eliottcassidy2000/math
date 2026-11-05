#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations, permutations
from collections import deque
from fractions import Fraction
from math import prod
import json, hashlib, argparse

HERE=Path(__file__).resolve().parent

def check(test,message):
    if not test:raise RuntimeError(message)

def bits(m):
    while m:
        t=m&-m;yield t.bit_length()-1;m-=t

def graph(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def tables(a):
    n=len(a);E=[[0]*n for _ in range(1<<n)];S=[[0]*n for _ in E];H=[0]*len(E);H[0]=1
    for m in range(1,len(E)):
        for i in bits(m):
            prev=m^(1<<i)
            if not prev:E[m][i]=S[m][i]=1
            else:
                E[m][i]=sum(E[prev][j] for j in bits(prev) if a[j]>>i&1)
                S[m][i]=sum(S[prev][j] for j in bits(prev) if a[i]>>j&1)
        H[m]=sum(E[m]);check(H[m]==sum(S[m]),'endpoint dual mismatch')
    return H,E,S

def signature(a,r,tab=None):
    H,E,S=tab or tables(a);full=len(H)-1;deleted=full^(1<<r)
    L=sum(E[m][r]*H[full^m] for m in range(1,full+1) if m>>r&1)
    Q=sum(S[m][r]*H[full^m] for m in range(1,full+1) if m>>r&1)
    B=sum(H[m]*H[deleted^m] for m in range(full+1) if m&deleted==m)
    return H[full],S[full][r],L,E[full][r],Q,B,int(len(a)==1),H[deleted]

def state(a,r,tab=None):
    g,s,L,E,Q,B,delta,d=signature(a,r,tab);return g,E,Q,d

def unary(a,r,op):
    n=len(a)
    if op=='D':return [x|(1<<n) for x in a]+[0],r
    b=[x|((1<<n) if i!=r else 0) for i,x in enumerate(a)]+[1<<r]
    return b,n if op=='A' else r

def glue(a,r,b,t):
    n=len(a);m=len(b);high=((1<<m)-1)<<n
    c=[x|high for x in a]+[x<<n for x in b]
    c[r]^=1<<(n+t);c[n+t]|=1<<r
    return c,n+t

def step(q,op):
    h,e,R,d=q
    if op=='A':return h+R-e,h-e,h+R,h
    if op=='K':return h+R-e,d,2*R-d,d
    if op=='D':return h,0,2*R,d
    g,s,L,E,Q,B,delta,dg=op
    return g*h-s*e+L*R,E*h-delta*e,Q*h+B*R,dg*h

def literal_state(a,r):
    n=len(a);h=e=R=0
    for word in permutations(range(n)):
        valid=all(a[word[i]]>>word[i+1]&1 for i in range(n-1))
        h+=valid;e+=valid and word[-1]==r
        cut=word.index(r)
        R+=all(a[word[i]]>>word[i+1]&1 for i in range(n-1) if i!=cut-1)
    remaining=[i for i in range(n) if i!=r]
    d=sum(all(a[word[i]]>>word[i+1]&1 for i in range(len(word)-1)) for word in permutations(remaining))
    return h,e,R,d

def saturate(bank,bound):
    seed=(1,1,1,1);seen={seed};queue=deque([seed]);support={1};ops=['A','K','D']+list(bank)
    while queue:
        q=queue.popleft()
        for op in ops:
            child=step(q,op);h,e,R,d=child
            check(0<=e<=d<=h<=R and R+e>=h+d,'physical cone failed in saturation')
            if h<=bound:support.add(h)
            if R<=bound and child not in seen:seen.add(child);queue.append(child)
    return seen,support

def determinant(a):
    a=[[Fraction(x) for x in row] for row in a];d=Fraction(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None:return 0
        if i!=j:a[i],a[j]=a[j],a[i];d=-d
        pivot=a[i][i];d*=pivot
        for j in range(i+1,len(a)):
            f=a[j][i]/pivot
            for k in range(i,len(a)):a[j][k]-=f*a[i][k]
    check(d.denominator==1,'noninteger determinant');return int(d)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',default=str(HERE.parent/'combined-growth'));parser.add_argument('--bound',type=int,default=1000);args=parser.parse_args()
    source=Path(args.source);records=json.loads((source/'inventory.json').read_text())
    signatures=[tuple(x['signature'])+(x['H_deleted_root'],) for x in records]
    check(len(signatures)==len(set(signatures))==65,'inventory cardinality')
    for rec,sig in zip(records,signatures):check(signature(rec['rows'],rec['root'])==sig,'inventory physical mismatch')
    regenerated=set();roots=graphs=unary_checks=gadget_checks=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=graph(n,code);tab=tables(a);graphs+=1
            for r in range(n):
                sig=signature(a,r,tab);regenerated.add(sig);roots+=1;q=(sig[0],sig[3],sig[4],sig[7])
                if n<=4:
                    for op in ('A','K','D'):
                        child,t=unary(a,r,op);check(state(child,t)==step(q,op),'unary physical mismatch');unary_checks+=1
                if n<=3:
                    for rec,op in zip(records,signatures):
                        child,t=glue(a,r,rec['rows'],rec['root']);check(state(child,t)==step(q,op),'gadget physical mismatch');gadget_checks+=1
    check(set(signatures)==regenerated,'inventory not complete through five')
    explicit=[]
    for word,expected_rows,expected_state in (
        ([53,'K','A'],[246,240,194,7,204,220,136,8],(163,70,168,75)),
        ([41,'K','K','A'],[502,484,448,7,462,476,392,264,8],(247,112,250,115))):
        a=[0];r=0;q=(1,1,1,1);history=[q]
        for op in word:
            if isinstance(op,int):
                rec=records[op];a,r=glue(a,r,rec['rows'],rec['root']);q=step(q,signatures[op])
            else:a,r=unary(a,r,op);q=step(q,op)
            check(state(a,r)==q,'explicit intermediate state');history.append(q)
        check(a==expected_rows and q==expected_state,'explicit witness output')
        check(literal_state(a,r)==q,'literal witness count')
        explicit.append(dict(word=word,rows=a,root=r,state=q,history=history))
    collision_graph=[2,12,9,17,7]
    collision=[]
    for root in (2,3):
        a=collision_graph;r=root;q=state(a,r)
        check(q[:3]==(13,2,18),'triple collision')
        for op in ('K','A'):
            a,r=unary(a,r,op);q=step(q,op);check(state(a,r)==q,'collision continuation')
        collision.append(q[0])
    check(collision==[55,59],'collision final counts')
    check(state([28,17,26,2,8],3)==(9,3,9,3),'gcd hostile state')
    inputs=[(1,1,1,1)]
    for word in ('A','AA','AAD'):
        q=(1,1,1,1)
        for op in word:q=step(q,op)
        inputs.append(q)
    matrix=[]
    for q in inputs:
        row=[]
        for word in ('','A','AA','KA'):
            t=q
            for op in word:t=step(t,op)
            row.append(t[0])
        matrix.append(row)
    check(determinant(matrix)==-16,'rank minor')
    seen,support=saturate(signatures,args.bound)
    misses=[x for x in range(1,args.bound+1,2) if x not in support]
    check(misses==[x for x in (7,21,67,173) if x<=args.bound],'combined bounded omissions')
    receipt=dict(status='PASS',bound=args.bound,graphs=graphs,roots=roots,signatures=len(signatures),
                 unary_physical_checks=unary_checks,gadget_physical_checks=gadget_checks,
                 retained_states=len(seen),odd_misses=misses,triple_collision_continuations=collision,
                 literal_verified_witnesses=explicit,input_by_continuation_matrix=matrix,determinant=-16,
                 inventory_sha256=hashlib.sha256((source/'inventory.json').read_bytes()).hexdigest(),
                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    rendered=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(HERE/'receipt.json').write_text(rendered);print(rendered,end='')
if __name__=='__main__':main()
