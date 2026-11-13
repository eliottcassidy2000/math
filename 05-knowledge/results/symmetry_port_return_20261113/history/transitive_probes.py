#!/usr/bin/env python3
"""Exact controls for the bounded-run physical probe theorems."""
from itertools import product
from math import factorial, comb
from fractions import Fraction
from pathlib import Path
import json
from hypergraph_probe import inventory, blowup, require

def stirling(n,r):
    a=[[0]*(r+1)for _ in range(n+1)];a[0][0]=1
    for i in range(1,n+1):
        for j in range(1,r+1):a[i][j]=a[i-1][j-1]+j*a[i-1][j]
    return a[n][r]

def TT(n):return tuple(((1<<n)-1)^((1<<(i+1))-1)for i in range(n))

def c(n,r):return factorial(r)*stirling(n,r)

def main():
    controls=0
    for n in range(1,7):
        for r in range(1,7):
            direct=sum(len(set(w))==r for w in product(range(r),repeat=n))
            exponential=sum((-1)**(r-k)*comb(r,k)*k**n for k in range(1,r+1))
            require(direct==c(n,r)==exponential,'ordered-cover basis mismatch')
            controls+=1
    C3=(2,4,1);Q4=blowup(C3,(TT(2),TT(1),TT(1)));require(inventory(Q4)[0][-1]==5,'Q4 source count mismatch')
    rows=[]
    for n in range(1,5):
        M=TT(n)
        W=blowup(C3,(M,TT(1),TT(1)))
        L3=blowup((2,0),(C3,M));L5=blowup((2,0),(Q4,M))
        hs=[inventory(G)[0][-1]for G in (W,L3,L5)]
        require(hs==[2**n+1,3,5],'literal separator context mismatch')
        rows.append(dict(n=n,counts=hs,differences=[hs[0]-hs[1],hs[0]-hs[2]]))
    require(all(0 in z['differences']for z in rows[:2]),'simultaneous-separator hostile mismatch')
    require(all(x!=0 for x in rows[2]['differences']),'common separator failed')
    source=[[2,-1,3],[0,4,-2]];R,S=2,3
    values=[[sum(source[r-1][s-1]*c(n,r)*c(m,s)for r in range(1,R+1)for s in range(1,S+1))for m in range(1,S+1)]for n in range(1,R+1)]
    recovered=[[Fraction(0)]*S for _ in range(R)]
    for n in range(1,R+1):
        for m in range(1,S+1):
            rest=sum(recovered[r-1][s-1]*c(n,r)*c(m,s)for r in range(1,n+1)for s in range(1,m+1)if(r,s)!=(n,m))
            recovered[n-1][m-1]=Fraction(values[n-1][m-1]-rest,c(n,n)*c(m,m))
    require(recovered==source,'signed tensor reconstruction mismatch')
    frame=list(TT(9))
    for a,b in [(0,4),(0,5),(2,7),(2,8)]:frame[a]^=1<<b;frame[b]^=1<<a
    frame_samples=[]
    for n in range(1,4):
        row=[]
        for m in range(1,3):
            blocks=[TT(1)]*9;blocks[3]=TT(n);blocks[1]=TT(m)
            value=inventory(blowup(frame,blocks))[0][-1]
            require(value==60*3**n*2**m+1,'literal Pi-frame count mismatch')
            row.append(value)
        frame_samples.append(row)
    frame_coeff=[[Fraction(0)]*2 for _ in range(3)]
    for n in range(1,4):
        for m in range(1,3):
            rest=sum(frame_coeff[r-1][s-1]*c(n,r)*c(m,s)for r in range(1,n+1)for s in range(1,m+1)if(r,s)!=(n,m))
            frame_coeff[n-1][m-1]=Fraction(frame_samples[n-1][m-1]-rest,c(n,n)*c(m,m))
    require(frame_coeff==[[361,180],[360,180],[120,60]],'Pi-frame tensor reconstruction mismatch')
    result=dict(status='FINITE-EXACT',universe='onto-map basis for n,r1..6; literal single-module tournament contexts n1..4; one signed2x3 tensor; literal nine-site Pi frame on3x2 transitive grid',basis_controls=controls,diagonal=[factorial(n)for n in range(1,7)],literal_contexts=rows,product_bases=[1,2,4],first_common_separator=3,tensor_source=source,tensor_samples=values,tensor_recovered=[[int(x)for x in row]for row in recovered],pi_frame=dict(rows=frame,M_site=3,N_site=1,transitive_samples=frame_samples,run_coefficients=[[int(x)for x in row]for row in frame_coeff]))
    Path(__file__).with_name('transitive_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
