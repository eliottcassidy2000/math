#!/usr/bin/env python3
"""Bounded-output exact consumers of unbounded authenticated migration words."""
from collections import deque
import json, hashlib
from pathlib import Path
from migration import compiler,tournament,multiply,extend

S0=((2,1,1,0),(0,1,0,0),(0,0,0,0),(0,0,1,1))
S1=((1,0,0,0),(1,2,1,0),(0,0,1,1),(0,0,0,0))
R1=((1,1,1,0),(0,0,0,0),(0,0,0,1),(2,0,1,0))
Z=((0,0,0,0),(0,0,0,0),(0,0,0,0),(2,0,1,1))
def require(c,msg):
    if not c:raise RuntimeError(msg)
def to_native(v):
    C,L,g,t=v
    return g+t,2*C+g+t,2*C+2*L+g,g

def positive_matrix(M):
    # Conversion uses exact integer arithmetic, not symbolic software.
    P=((0,0,1,1),(2,0,1,1),(2,2,1,0),(0,0,1,0))
    columns=[multiply(M,tuple(P[i][j] for i in range(4))) for j in range(4)]
    cols=[]
    for h,F,q,g in columns:
        require((F-h)%2==0 and (q-F+h-g)%2==0,'nonintegral conjugate')
        cols.append(((F-h)//2,(q-F+h-g)//2,g,h-g))
    return tuple(tuple(cols[j][i] for j in range(4)) for i in range(4))

def strip(r):
    if r==1:return R1
    a=2**(r-1);b=2*3**(r-1)-2**r+1
    return ((a,a,a//2,0),(b-a,b-2*a+2,(a+b-3)//2,a-1),
            (2*(a-1),2*(a-2),a,1),(2,2,1,0))

def chain_bank(max_local=4):
    bank={S0:'S0',S1:'S1'}
    for n in range(2,max_local+1):
        for code in range(1<<(n*(n-1)//2)):
            local=tournament(n,code)
            # Fresh vertices induce a transitive tournament, degree test exact.
            if sorted(sum(local[i]>>j&1 for j in range(1,n)) for i in range(1,n))!=list(range(n-1)):continue
            M=positive_matrix(compiler(local,n-1))
            require(min(x for row in M for x in row)>=0,'negative chosen-bank matrix')
            bank.setdefault(M,[n,code,n-1])
    return bank

def closure(bank,N):
    start=(0,0,0,1);seen={start:None};queue=deque([start]);values={1:start}
    while queue:
        v=queue.popleft()
        for M,recipe in bank.items():
            w=tuple(min(N+1,x) for x in multiply(M,v))
            F=2*w[0]+w[2]+w[3]
            if F>N:continue
            if F not in values:values[F]=w
            if w in seen:continue
            seen[w]=(v,recipe);queue.append(w)
    return seen,values

def trace(seen,w):
    ans=[]
    while seen[w] is not None:
        v,recipe=seen[w];ans.append(recipe);w=v
    return list(reversed(ans))

def source_marked(A,cut):
    n=len(A);size=1<<n;end=[[0]*n for _ in range(size)];start=[[0]*n for _ in range(size)]
    H=[0]*size;H[0]=1
    for mask in range(1,size):
        for v in range(n):
            if not(mask>>v&1):continue
            rest=mask^(1<<v)
            if not rest:end[mask][v]=start[mask][v]=1
            else:
                end[mask][v]=sum(end[rest][u] for u in range(n) if rest>>u&1 and A[u]>>v&1)
                start[mask][v]=sum(start[rest][u] for u in range(n) if rest>>u&1 and A[v]>>u&1)
        H[mask]=sum(end[mask])
    q=sum(H[(size-1)^Q]*sum(start[Q][v] for v in range(n) if cut>>v&1) for Q in range(1,size))
    g=sum(end[-1][v] for v in range(n) if cut>>v&1)
    apex=list(A)+[cut]
    for v in range(n):
        if not cut>>v&1:apex[v]|=1<<n
    # Separate direct end-DP on the full adjacency.
    m=n+1;D=[[0]*m for _ in range(1<<m)]
    for mask in range(1,1<<m):
        for v in range(m):
            if mask>>v&1:
                rest=mask^(1<<v)
                D[mask][v]=1 if not rest else sum(D[rest][u] for u in range(m) if rest>>u&1 and apex[u]>>v&1)
    return H[-1],sum(D[-1]),q,g

def replay(recipes):
    A=(0,);cut=0;v=(0,0,0,1)
    stages=[]
    for recipe in recipes:
        if isinstance(recipe,str):
            bit=int(recipe[-1]);n=len(A);A=tuple(x|(1<<n) for x in A)+(0,);cut|=bit<<n
            M=S0 if not bit else S1
        else:
            n,code,z=recipe;local=tournament(n,code);M=positive_matrix(compiler(local,z));A,cut=extend(A,cut,local,z)
        v=multiply(M,v);actual=source_marked(A,cut)
        require(actual==to_native(v),'actual word replay mismatch')
        stages.append({'recipe':recipe,'state':actual,'parent_order':len(A)})
    return {'parent_adjacency_masks':A,'cut':cut,'stages':stages,'state':to_native(v)}

def main():
    N=173;r0=(N).bit_length()+2 # ceil(log2(N+1))+2
    cap=lambda M:tuple(tuple(min(N+1,x) for x in row) for row in M)
    require(cap(strip(r0))==cap(strip(r0+1)),'strip tail not stabilized')
    # The proof establishes stability at all larger r; these extra controls catch coding mistakes.
    for r in range(r0,r0+50):require(cap(strip(r))==cap(strip(r0)),'tail coefficient mismatch')
    zero_bank={S0:'S0',S1:'S1',Z:'Z'}
    for r in range(1,r0+1):zero_bank.setdefault(strip(r),f'R{r}')
    seen,vals=closure(zero_bank,N)
    missing=[x for x in range(1,N+1,2) if x not in vals]
    require(missing==[7,21,67,173],'zero strip spectrum changed')
    all_local_tests=0
    for r in range(1,7):
        local=tuple([0]+[1+sum(1<<j for j in range(i+1,r+1)) for i in range(1,r+1)])
        require(positive_matrix(compiler(local,1))==strip(r),'strip/general compiler mismatch');all_local_tests+=1
    # Retained hostile against the discarded b=3^(r-1) extrapolation.
    bad_native=(1,1,1,1)
    r3=to_native(multiply(strip(3),(0,0,1,0)))
    require(r3[2]==20 and r3[2] != 18,'r3 missing cover hostile')
    bank=chain_bank();seen2,vals2=closure(bank,N)
    require([x for x in range(1,N+1,2) if x not in vals2]==[7,21],'mixed-chain result changed')
    witnesses={str(t):replay(trace(seen2,vals2[t])) for t in [67,173]}
    report={'status':'PASS','scope':'Unbounded words; strip length arbitrary; mixed-chain bank has at most three fresh vertices per step. No universal acceptance claim.',
            'strip':{'target_cap':N,'saturation_limit':N+1,'tail_representative':r0,'matrices':len(zero_bank),'states':len(seen),'missing_odds':missing,'strip_formula_compiler_crosschecks':all_local_tests},
            'mixed_chain':{'target_cap':N,'matrices':len(bank),'states':len(seen2),'missing_odds':[7,21],'repairs':witnesses},
            'hostile':{'input_native':bad_native,'r':3,'actual_output_native':r3,'discarded_q':18,'correct_q':20}}
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
