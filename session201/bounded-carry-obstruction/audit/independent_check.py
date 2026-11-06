"""Independent actual-index controls, no author imports; checks survive -O."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import json, hashlib

def need(c,why):
    if not c: raise RuntimeError(why)
@lru_cache(None)
def g(u):
    if u==1: return 1
    if u%2==0: return g(u//2)
    h=(u+1)//2
    return g(h)+h-h%2

def word(u): return ''.join('Z' if x=='0' else 'O' for x in bin(u-1)[3:])
def sample(w):
    u=2;t=F(1,2);a=0;run=[];o=0;i=0
    while i<len(w):
        if w[i]=='O':
            u*=2;t/=2;o+=1;i+=1;continue
        s=t;before=g(u)-u*t;r=0
        while i<len(w) and w[i]=='Z':
            u=2*u-1;t=(t+1)/2;r+=1;i+=1
        eta=g(u)-u*t-before
        need(eta==s+(1-s)/2**r,'independent run increment')
        need(F(1,2**(o+1))<=s<=F(1,2),'run start')
        need(eta>=F(1,2**(o+1))+F(1,2**(r+1)),'additive lower')
        need(eta**2*2**(o+r)>=1,'AMGM lower')
        run.append((o,r,eta));o=0
    return u,t,g(u)-u*t,run

def mm(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a): return list(map(list,zip(*a)))
P=[[1,0,2,1],[0,0,1,0],[1,1,2,2],[0,0,1,1]]
Pinv=[[1,-1,0,-1],[-1,1,1,-1],[0,1,0,0],[0,-1,0,1]]
Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]]
M={'O':[[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]],'Z':[[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]],'J':[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]]}
I=[[int(i==j) for j in range(4)] for i in range(4)]
need(mm(Pinv,P)==I,'inverse cone coordinates')
z=[1,1,2,2];base=mm(mm(tr(P),Q),P)
need(all(2*base[i][j]>=z[i]*z[j] for i in range(4) for j in range(4)),'base product bound')
cone={c:mm(mm(Pinv,m),P) for c,m in M.items()}
for c,T in cone.items():
    need(all(x>=0 for row in T for x in row),'cone preservation '+c)
    need(all(sum(z[i]*T[i][j] for i in range(4))>=z[j] for j in range(4)),'cover monotonicity '+c)
menu=['','J','JJ','JJJ','JOZJ','JJZJ','JJJJ','JOOZJ','JJZZJ','JOZJJ','JJZJJ','JJJJJ']
for w in menu:
    T=I
    for c in w:T=mm(M[c],T)
    C=mm(mm(mm(tr(P),Q),T),P)
    need(all(2*C[i][j]>=z[i]*z[j] for i in range(4) for j in range(4)),'menu bound '+w)
thresholds=[F(0),F(1,2),F(1),F(2)]; details=[];total=0;all_errors={}
for k in range(15):
    B=2**k; counts=[0]*((k+1)//2+1);good={E:0 for E in thresholds}
    for u in range(B+1,2*B+1):
        w=word(u);need(len(w)==k,'unique length')
        up,t,e,runs=sample(w);need(up==u,'inverse word index')
        A=sum((c=='Z')*2**i for i,c in enumerate(w))
        rev=int(format(A,'0%db'%k)[::-1],2) if k else 0
        need(u==2*B-rev and t==F(2*A+1,2*B),'bit reversal label')
        need(all(v>=0 for v in mm(Pinv,[[g(u)],[u%2],[u],[1]])[0]),'first cone coordinate')
        need(all(row[0]>=0 for row in mm(Pinv,[[g(u)],[u%2],[u],[1]])),'actual cone')
        R=len(runs);counts[R]+=1
        need(sum(o+r for o,r,_ in runs)<=k,'disjoint blocks')
        need(e==sum((x[2] for x in runs),F(0)),'actual recurrence equals whole-run sum')
        if R: need(e**(2*R)*2**k>=R**(2*R),'Jensen exact power form')
        else: need(e==0,'zero run')
        need((e==0)==(w=='O'*k),'zero error iff all O')
        all_errors[u]=e
        for E in thresholds:
            if e<=E: good[E]+=1
        total+=1
    need(counts==[comb(k+1,2*R) for R in range(len(counts))],'all run counts')
    for E in thresholds:
        bound=sum(comb(k+1,2*R) for R in range(len(counts)) if R==0 or E**(2*R)*2**k>=R**(2*R))
        need(good[E]<=bound,'necessary run upper count')
        details.append({'k':k,'E':str(E),'actual_words':good[E],'necessary_run_bound':bound})
for N in range(2,4097):
    need(sum(e==0 for u,e in all_errors.items() if u<=N)==N.bit_length()-1,'exact zero cumulative count')
for m in range(101):
    u,t,e,r=sample('ZO'*m)
    need(e==F(2*m,3)+F(1,9)*(1-F(1,4**m)),'alternating exact formula including empty')
root=Path(__file__).resolve().parent
for row in json.loads((root/'INPUTS.json').read_text()):
    need(hashlib.sha256((root/row['frozen']).read_bytes()).hexdigest()==row['sha256'],'frozen input integrity')
author=json.loads((root/'AUTHOR-normal.json').read_text())
need(details==author['exact_bounded_error_counts'],'independent author count match')
print(json.dumps({'status':'PASS','actual_indices':total,'max_word_length':14,'thresholds':list(map(str,thresholds)),'exact_bounded_error_counts':details,'zero_error_cumulative_N':4096,'alternating_family_m':100,'cone_matrices':cone,'twelve_templates_product_bound':'PASS','all_order_claim':'Established by written audit, not finite controls'},indent=2,sort_keys=True))
