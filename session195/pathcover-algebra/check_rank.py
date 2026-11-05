#!/usr/bin/env python3
"""Independent context-matrix rank control: vertex DP versus cyclic quotient words."""
from fractions import Fraction
from math import factorial,prod
import hashlib,json,pathlib
from verify import hamilton,require
HERE=pathlib.Path(__file__).resolve().parent

def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def context(n,m):
    if m==0:return transitive(n)
    sizes=(n,m,m);offsets=(0,n,n+m);out=[0]*(n+2*m)
    for b,size in enumerate(sizes):
        for x in range(size):
            for y in range(x+1,size):out[offsets[b]+x]|=1<<(offsets[b]+y)
            nxt=(b+1)%3
            for y in range(sizes[nxt]):out[offsets[b]+x]|=1<<(offsets[nxt]+y)
    return tuple(out)

def ordered_covers(n):
    s=[[0]*(n+1) for _ in range(n+1)];s[0][0]=1
    for i in range(1,n+1):
        for r in range(1,i+1):s[i][r]=s[i-1][r-1]+r*s[i-1][r]
    return [factorial(r)*s[n][r] for r in range(n+1)]

def word_coefficient(r,s,t):
    target=(r,s,t);length=r+s+t;answer=0
    for start in range(3):
        counts=[0,0,0]
        for j in range(length):counts[(start+j)%3]+=1
        if tuple(counts)==target:answer+=1
    return answer

def coefficients(m):
    if m==0:return [0,1]
    c=ordered_covers(m)
    return [0]+[sum(word_coefficient(r,s,t)*c[s]*c[t] for s in range(1,m+1) for t in range(1,m+1)) for r in range(1,m+2)]

def determinant(a):
    a=[[Fraction(x) for x in row] for row in a];result=Fraction(1);n=len(a)
    for i in range(n):
        pivot=next((j for j in range(i,n) if a[j][i]),None)
        if pivot is None:return 0
        if pivot!=i:a[i],a[pivot]=a[pivot],a[i];result=-result
        value=a[i][i];result*=value
        for j in range(i+1,n):
            factor=a[j][i]/value
            for k in range(i,n):a[j][k]-=factor*a[i][k]
    require(result.denominator==1,'integer determinant')
    return result.numerator

def main():
    N=5;matrix=[];ks=[]
    for m in range(N):
        k=coefficients(m);ks.append(k)
        require(k[-1]==factorial(m)**2,'highest run coefficient')
    for n in range(1,N+1):
        row=[];c=ordered_covers(n)
        for m in range(N):
            actual=hamilton(context(n,m))
            transfer=sum(c[r]*ks[m][r] for r in range(1,min(n,m+1)+1))
            require(actual==transfer,f'context mismatch n={n},m={m}')
            row.append(actual)
        matrix.append(row)
    determinants=[]
    for j in range(1,N+1):
        actual=determinant([row[:j] for row in matrix[:j]])
        expected=prod(factorial(n)*factorial(n-1)**2 for n in range(1,j+1))
        require(actual==expected,'rank determinant')
        determinants.append(actual)
    receipt={'status':'PASS','matrix':matrix,'context_coefficients':ks,'leading_principal_determinants':determinants,
       'largest_order':3*N-2,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
    (HERE/'rank_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
