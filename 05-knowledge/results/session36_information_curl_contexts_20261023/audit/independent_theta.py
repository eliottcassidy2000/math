#!/usr/bin/env python3
"""Direct coordinate-vector root enumeration, independent of producer code."""
import itertools,collections,json,math
from pathlib import Path
OUT=Path(__file__).resolve().parent

def need(x,why):
    if not x:raise RuntimeError(why)

def dot(x,y):return sum(a*b for a,b in zip(x,y))
def roots(n):
    return [tuple(int(v==i)-int(v==j) for v in range(n)) for i in range(n) for j in range(n) if i!=j]

def h(M):
    n=len(M);d={(1<<v,v):1 for v in range(n)}
    for k in range(1,n):
        e=collections.defaultdict(int)
        for (s,v),a in d.items():
            for w in range(n):
                if not s>>w&1 and M[v][w]==1:e[s|1<<w,w]+=a
        d=e
    return sum(d.values())

def main():
    bases=[]
    for n in range(2,7):
        rr=roots(n);G=[[dot(x,y) for y in rr] for x in rr];found=[]
        def rec(word):
            if len(word)==n-1:found.append(tuple(word));return
            for j in range(len(rr)):
                if all(G[i][j]==1 for i in word):rec(word+[j])
        rec([])
        expected=2 if n==2 else 2*math.factorial(n)
        need(len(found)==expected,('root basis count',n,len(found)))
        bases.append(dict(n=n,bases=len(found)))
    n=8;C=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            C[i][j]=1 if i==0 or (j-i)%7 in (1,2,4) else -1;C[j][i]=-C[i][j]
    need(all(sum(C[i][k]*C[j][k] for k in range(n))==(7 if i==j else 0) for i in range(n) for j in range(n)),'base conference identity')
    rr=roots(n);G=[[dot(x,y) for y in rr] for x in rr];N=len(rr)
    triples=[(a,b,c) for a in range(N) for b in range(N) if G[a][b]==-1 for c in range(N) if G[a][c]==G[b][c]==0]
    ref=None;classes=collections.Counter();records=[]
    for mask in range(128):
        signs=[1]+[1 if not mask>>(i-1)&1 else -1 for i in range(1,n)]
        M=[[signs[i]*C[i][j]*signs[j] for j in range(n)] for i in range(n)]
        Mr=[[sum(M[i][j]*y[j] for j in range(n)) for i in range(n)] for y in rr]
        W=[[dot(x,my) for my in Mr] for x in rr]
        pair=collections.Counter((G[a][b],W[a][b]) for a in range(N) for b in range(N))
        if ref is None:ref=pair
        need(pair==ref,('genus two differs',mask))
        count=sum(W[a][b]==-3 and W[a][c]==-4 and W[b][c]==2 for a,b,c in triples)
        H=h(M);universal=any(abs(sum(row))==7 for row in M)
        need((count,H)==((84,189) if universal else (36,333)),('genus three/H',mask,count,H,universal))
        need(H==441-3*count,('finite affine formula',mask))
        classes[count,H,universal]+=1;records.append(dict(mask=mask,coefficient=count,H=H,universal=universal))
    report=dict(status='PASS',root_basis_counts=bases,conference_switches=len(records),genus_two_support_size=len(ref),gram_filtered_root_triples=len(triples),classes=[dict(coefficient=k[0],H=k[1],universal=k[2],switches=v) for k,v in sorted(classes.items())],records=records)
    (OUT/'theta_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
if __name__=='__main__':main()
