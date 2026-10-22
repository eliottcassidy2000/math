"""Exact Fourier isolation of each labelled core arc using k-1 base ports.
No imported slot realization or path-cover formula: direct endpoint DP.
"""
from common import core,response,H
from itertools import combinations

def experiment(Q,u,v,fixed=None,reverse_base=False):
    k=len(Q)
    if k<2 or u==v or not (0<=u<k and 0<=v<k):raise ValueError('distinct core endpoints required')
    m=k-1
    fixed=[0]*m if fixed is None else list(fixed)
    if len(fixed)!=m or any(not 0<=a<1<<k for a in fixed):raise ValueError('bad fixed masks')
    remaining=[w for w in range(k) if w not in (u,v)]
    R=[v,k]
    for i,w in enumerate(remaining):R += [w,k+i+1]
    edges=[]
    for a,b in zip(R,R[1:]):
        if a<k:edges.append((b-k,a,1)) # positive sign means core->base
        else:edges.append((a-k,b,0)) # positive sign means base->core
    d=len(edges)
    if d!=2*k-3 or len({(i,w) for i,w,_ in edges})!=d:raise RuntimeError('path malformed')
    total=0
    for assignment in range(1<<(d+1)):
        word=fixed[:];chi=1
        for j,(i,w,positive_bit) in enumerate(edges):
            positive=assignment>>j&1
            bit=positive_bit if positive else 1-positive_bit
            word[i]=(word[i]&~(1<<w))|(bit<<w)
            if not positive:chi=-chi
        word[-1]=(word[-1]&~(1<<u))|(((assignment>>d)&1)<<u)
        if not reverse_base:z=response(Q,word)
        else:
            out=[sum(x<<j for j,x in enumerate(row)) for row in Q]+[0]*m
            for i,a in enumerate(word):
                for w in range(k):
                    if a>>w&1:out[w]|=1<<(k+i)
                    else:out[k+i]|=1<<w
                for j in range(i):out[k+i]|=1<<(k+j)
            z=H(out)
        total+=chi*z
    return total,1<<(d+1)

if __name__=='__main__':
    tests=0;queries=0
    for k in [2,3,4]:
        for mask in range(1<<(k*(k-1)//2)):
            Q=core(k,mask)
            for u,v in combinations(range(k),2):
                fixed=[(mask+7*i+3*u+v)%(1<<k) for i in range(k-1)]
                total,n=experiment(Q,u,v,fixed)
                expected=2*(2*Q[u][v]-1)
                if total!=expected:raise RuntimeError(('faithfulness failed',k,mask,u,v,total,expected))
                tests+=1;queries+=n
        print('all labelled cores through',k,'PASS; cumulative arc tests',tests,'DP queries',queries,flush=True)
    for mask in [0,2,41,341,1023]:
        Q=core(5,mask)
        for u,v in [(0,4),(3,1)]:
            total,n=experiment(Q,u,v,[3,17,9,30])
            if total!=2*(2*Q[u][v]-1):raise RuntimeError(('k5 sample failed',mask,u,v,total))
            tests+=1;queries+=n
    # Base orientation need not be transitive in the specified order.
    for mask in range(8):
        Q=core(3,mask);total,n=experiment(Q,2,0,[7,5],reverse_base=True)
        if total!=2*(2*Q[2][0]-1):raise RuntimeError('reversed base failed')
        tests+=1;queries+=n
    print('TOTAL PASS:',tests,'arc tests,',queries,'direct DP queries; includes k5 samples and reversed-base controls.')
