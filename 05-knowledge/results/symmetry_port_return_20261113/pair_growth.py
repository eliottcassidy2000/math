"""Exact paired vertex-birth probe for the nontrivial determinant pair."""
from pathlib import Path
from math import gcd
import json

ROOT=Path(__file__).resolve().parent
LEFT=(4,9,82,21,35,143,187,31)
RIGHT=(6,12,48,21,131,219,159,15)

def require(ok,why):
    if not ok: raise RuntimeError(why)

def inventory(rows):
    n=len(rows); full=(1<<n)-1
    require(all(not(row>>i&1) and 0<=row<=full for i,row in enumerate(rows)), 'loop/range')
    require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)), 'pair')
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,full+1):
        for v,val in enumerate(dp[mask]):
            if not val: continue
            rem=rows[v]&~mask
            while rem:
                b=rem&-rem;rem-=b
                dp[mask|b][b.bit_length()-1]+=val
    hh=[1]+[sum(row) for row in dp[1:]]
    return hh[-1],sum(hh[s]*hh[full^s] for s in range(full+1))

def add(rows,mask):
    n=len(rows); full=(1<<n)-1
    return tuple(row|((mask>>i&1)<<n) for i,row in enumerate(rows))+(full^mask,)

def strong(rows):
    n=len(rows)
    for start in range(n):
        seen=1<<start; old=-1
        while seen!=old:
            old=seen
            for v,row in enumerate(rows):
                if seen>>v&1: seen|=row
        if seen!=(1<<n)-1: return False
    return True

def extensions(rows):
    profiles={}
    for mask in range(1<<len(rows)):
        graph=add(rows,mask);h,p=inventory(graph)
        profiles.setdefault((h,p),[]).append(dict(mask=mask,rows=graph,strong=strong(graph)))
    return profiles

def analyze(left,right):
    L=extensions(left);R=extensions(right);hits=[]
    for (h,p),lm in L.items():
        if gcd(h,p)!=1 or h%3==0:continue
        for (j,q),rn in R.items():
            e=j*(6*p+h)-6*h*q
            if abs(e)!=1: continue
            A=12*p*q+6*h*q+2*p*j-h*j
            hits.append(dict(h=h,p=p,j=j,q=q,e=e,source=A,target=3*A+2*e,
                             left=lm,right=rn))
    return dict(left_profiles=len(L),right_profiles=len(R),hits=hits)

def main():
    require(inventory(LEFT)==(181,3424),'left control')
    require(inventory(RIGHT)==(179,3416),'right control')
    out=dict(status='FINITE-EXACT',universe='Every single-vertex extension on each supplied labelled order-eight module; all ordered pairs of resulting response profiles',base_left=LEFT,base_right=RIGHT,result=analyze(LEFT,RIGHT))
    (ROOT/'pair_growth.json').write_text(json.dumps(out,indent=2)+'\n')
    print('profiles',out['result']['left_profiles'],out['result']['right_profiles'],'hits',len(out['result']['hits']))
    for hit in out['result']['hits']:print({k:hit[k] for k in ['h','p','j','q','e','source','target']})

if __name__=='__main__':main()
