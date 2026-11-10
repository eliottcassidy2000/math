#!/usr/bin/env python3
"""Independent permutation census; imports no producer module. Explicit checks."""
import itertools as it
import json
from fractions import Fraction
from pathlib import Path

def ck(p,msg):
    if not p: raise RuntimeError(msg)

def profiles(n,bits):
    adj=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(it.combinations(range(n),2)):
        adj[i][j]=bool(bits>>k&1);adj[j][i]=not adj[i][j]
    h=[0]*(1<<n);h[0]=1
    start=[[0]*n for _ in h];end=[[0]*n for _ in h]
    for k in range(1,n+1):
        for p in it.permutations(range(n),k):
            if all(adj[u][v] for u,v in zip(p,p[1:])):
                s=sum(1<<v for v in p);h[s]+=1;start[s][p[0]]+=1;end[s][p[-1]]+=1
    full=(1<<n)-1
    reachmat=[[i==j or adj[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):reachmat[i][j] |= reachmat[i][k] and reachmat[k][j]
    out=[]
    for r in range(n):
        reach={r}
        while True:
            new=reach|{j for i in reach for j in range(n) if adj[i][j]}
            if new==reach:break
            reach=new
        U=sum(start[s][r]*h[full^s] for s in range(1,full+1))
        state=(h[full],end[full][r],U,h[full^(1<<r)])
        if len(reach)==n:out.append(state)
        # D injection: cut a Hamilton path just before r into rooted suffix + prefix.
        ck(U>=h[full],('cover bound',n,bits,r))
        outs=[a for a in range(n) if adj[r][a]]
        extra=sum(end[full^(1<<r)][a] for a in outs)
        ck(U-h[full]>=len(outs)*(len(outs)-1)+extra,('paircover deficit',n,bits,r))
        if len(outs)==1:ck(U-h[full]==extra,('outdegree-one equality',n,bits,r))
        C=sum(1<<v for v in range(n) if reachmat[r][v] and reachmat[v][r])
        P=sum(1<<v for v in range(n) if not reachmat[r][v])
        B=full^(C|P)
        UC=sum(start[x][r]*h[C^x] for x in range(1,full+1) if x&C==x)
        QB=sum(h[x]*h[B^x] for x in range(full+1) if x&B==x)
        ck(h[full]==h[P]*h[C]*h[B] and U==h[P]*UC*QB,('arbitrary SCC formula',n,bits,r))
    return out

def move(s,z):
    D,E,U,d=s
    return (D+U-E,d,2*U-d,d) if z else (D,0,2*U,d)

def image(bank,limit):
    seen=set();todo=list(bank);values={s[0] for s in bank if s[0]<=limit}
    while todo:
        s=todo.pop()
        if s in seen or s[2]>limit:continue
        seen.add(s)
        t=move(s,True)
        if t[0]<=limit:values.add(t[0])
        todo.extend([t,move(s,False)])
    return values,len(seen)

def canonical(n):
    D=1;terms=[]
    while n>2:
        m=(n+1)//2
        if n%2:terms.append(m-m%2)
        n=m
    return D+sum(terms)

bank=set();counts={};total=0
for n in range(2,6):
    for bits in range(1<<(n*(n-1)//2)):
        bank.update(profiles(n,bits));total+=1
    ck(all(0<=E<=D<=U and U>=2*d for D,E,U,d in bank),'seed bounds')
    if n>=3:
        vals,ns=image(bank,2001)
        counts[str(n)]={'seed_profiles':len(bank),'expanded_profiles':ns,'covered_odds_through_2001':len(vals),'misses_through_201':[v for v in range(1,202,2) if v not in vals]}
expected=json.loads((Path(__file__).parents[1]/'selector/RESULT.json').read_text())
ck(counts==expected['seed_banks'],'independent bank census')
J=32
lo=sum((Fraction(1,2**(j+2)+1) for j in range(J)),Fraction())
hi=lo+Fraction(1,2**(J+1))
ck(hi<Fraction(432,1000),'density rational')
for X in (100,1000,10000,100000):
    inputs=[n for n in range(1,(X-1)//2+1) if canonical(n)+2*n<=X]
    upper=sum((X//(2**(j+2)+1)+1) for j in range(X.bit_length()))
    ck(len(inputs)<=upper,('valuation bound',X))
# Verify compact CRT witnesses independently by searching the full interval,
# including the explicit small-case corrections actually given in the proof.
def prime(n):return n>=2 and all(n%d for d in range(2,__import__('math').isqrt(n)+1))
def member(n):return n>0 and n%2 and n not in (7,21) and not(n>=151 and prime(n))
for N in range(1,10001):
    if N<30:a={3:9,5:15,10:23}.get(N,2*N+1)
    else:a=next(a for a in range(2*N+1,3*N+1) if a%2 and a%3==0 and (a+2*N)%5==0)
    ck(2*N+1<=a<a+2*N<=5*N and member(a) and member(a+2*N),('model',N))
for N in range(1,61):
    for m in range(51,351,2):
        a=next(a for a in range(m-28,m+1,2) if a%3==0 and (a+2*N)%5==0)
        ck(member(a) and member(a+2*N),'anchor')
# Seven-letter pure-count compiler: literal independent permutation checks.
for N in (1,2):
    for word,coefficient in [('OJJJJOZ',40),('OJOZJOZ',42)]:
        a=[0] if N==1 else [2,0];root=0
        for letter in word:
            n=len(a);a.append(0 if letter=='O' else 1<<root)
            for i in range(n):
                if letter=='O' or i!=root:a[i]|=1<<n
            if letter=='J':root=n
        bits=sum(1<<k for k,(i,j) in enumerate(it.combinations(range(len(a)),2)) if a[i]>>j&1)
        states=profiles(len(a),bits)
        expected=(17+coefficient*N,3+6*N,21+50*N,3+6*N)
        ck(len(states)==len(a) and states[root]==expected,'pure-count compiler')
# Efficient endpoint-source profiles via independent literal permutations.
for N in range(1,17):
    a=[2,0];m=N;letters=[]
    while m>2:
        letters.append('Z' if m%2 else 'O');m=(m+1)//2
    word='Z' if N==1 else ''.join(reversed(letters))+'OZ'
    for letter in word:
        n=len(a);a.append(0 if letter=='O' else 1)
        for i in range(n):
            if letter=='O' or i!=0:a[i]|=1<<n
    root=0 if N==1 else len(a)-2
    bits=sum(1<<k for k,(i,j) in enumerate(it.combinations(range(len(a)),2)) if a[i]>>j&1)
    states=profiles(len(a),bits);g=canonical(N);eps=N%2
    ck(len(states)==len(a) and states[root]==(g+2*N,N,2*g+2*N-eps,g+N-eps),'efficient endpoint profile')
# Regression for the corrected Fminus inverse: do not check only forward maps.
for h in range(1,200,2):
    for e in range(1,200,2):
        hp,ep=20*h-3*e,6*h-e
        ck((hp-3*ep)//2==h and 3*hp-10*ep==e,'Fminus inverse sign regression')
# Literal H=173 certificate, independently reconstructed and counted by permutations.
a=[2,0]
for root,letter in [(0,c) for c in 'ZOOZ']+[(2,'Z'),(2,'Z')]:
    n=len(a);a.append(0 if letter=='O' else 1<<root)
    for i in range(n):
        if letter=='O' or i!=root:a[i]|=1<<n
bits=sum(1<<k for k,(i,j) in enumerate(it.combinations(range(len(a)),2)) if a[i]>>j&1)
literal=profiles(len(a),bits)
ck(len(literal)==8 and literal[2]==(173,9,209,9),'independent literal 173')
print(json.dumps({'status':'PASS','method':'all directed permutations on subsets, no producer imports','labelled_tournaments_orders_2_to_5':total,'seed_banks':counts,'density_constant_lower':str(lo),'density_constant_upper':str(hi),'literal_pure_count_compiler_outputs':4,'literal_efficient_endpoint_sources':16,'literal_173_profile':[173,9,209,9],'compact_model_pairs':10000,'constant_error_anchor_pairs':9000},indent=2,sort_keys=True))
