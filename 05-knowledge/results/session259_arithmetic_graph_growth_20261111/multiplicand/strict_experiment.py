from math import isqrt

def prime(a): return a>1 and all(a%d for d in range(2,int(a**.5)+1))
def barcode(N,strict):
    faces=[]
    def visit(s):
        faces.append(s)
        for b in range(s[-1]+1,N+1):
            if all(b%a==0 and not(strict and prime(a) and b==a*a) for a in s): visit(s+(b,))
    for a in range(2,N+1): visit((a,))
    faces.sort(key=lambda s:(s[-1],len(s),s))
    ix={s:i for i,s in enumerate(faces)};cols={};zeros=[];pairs={}
    for j,s in enumerate(faces):
        col=0
        if len(s)>1:
            for k in range(len(s)): col^=1<<ix[s[:k]+s[k+1:]]
        while col and col.bit_length()-1 in cols: col^=cols[col.bit_length()-1]
        if col: p=col.bit_length()-1;cols[p]=col;pairs[p]=j
        else: zeros.append(j)
    return [(len(faces[i])-1,faces[i][-1],faces[pairs[i]][-1] if i in pairs else None) for i in zeros if i not in pairs or faces[i][-1]!=faces[pairs[i]][-1]]
def squarefree_factors(m):
    ps=[]
    for p in range(2,m+1):
        if m%p==0 and prime(p):
            if m%(p*p)==0: return None
            ps.append(p)
    return ps
from collections import Counter
import json
receipts=[]
for N in [30,100,300]:
    w=barcode(N,False);s=barcode(N,True)
    predicted=[]
    for p in range(2,isqrt(N)+1):
        if not prime(p): continue
        predicted.append((0,p*p,2*p*p if 2*p*p<=N else None))
        for m in range(3,N//(p*p)+1,2):
            ps=squarefree_factors(m)
            if ps is not None:
                b=p*p*m
                predicted.append((len(ps),b,2*b if 2*b<=N else None))
    if Counter(s)!=Counter(w)+Counter(predicted):
        raise RuntimeError(('strict barcode mismatch',N))
    receipts.append({'cutoff':N,'weak_bars':len(w),'strict_bars':len(s),
                     'correction_bars':len(predicted),'status':'FINITE-EXACT PASS'})
# The last full matrix reduction is at N=300; its bars determine all prefix Euler values.
mu={1:1}
for n in range(2,301):
    mu[n]=-sum(mu[d] for d in range(1,n) if n%d==0)
M={n:sum(mu[j] for j in range(1,n+1)) for n in range(1,301)}
reduced=[bar for bar in s if bar!=(0,2,None)]
events=Counter()
for k,b,d in reduced:
    events[b]+=1
    if d is not None: events[d]+=1
if any(count!=1 for count in events.values()):
    raise RuntimeError('event collision')
for n in range(2,301):
    chi=sum((-1)**k for k,b,d in reduced if b<=n and (d is None or n<d))
    predicted=-M[n]+sum(M[n//(p*p)] for p in range(2,isqrt(n)+1) if prime(p))
    if chi!=predicted: raise RuntimeError(('strict Euler',n,chi,predicted))
print(json.dumps({'event_uniqueness_and_Euler_cutoff':300,'checks':receipts,
                  'first_higher_correction':[1,12,24],
                  'claim':'weak bars plus prime-square-scaled suspended weak bars and prime-square H0 bars'},indent=2))
