#!/usr/bin/env python3
"""Exact positive-integer termination certificate for the cheaper ternary queue."""
from fractions import Fraction as F
import json

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def step(n):
    t=3*n+1
    a=(t & -t).bit_length()-1
    return t >> a, a

def potential(n,q):
    return (n+1)*F(3,2)**q

def branch_bound(n):
    m=0
    left=n+1
    right=2
    while left*27 >= right*28:
        left*=27
        right*=28
        m+=1
    return m

def run(n):
    source=n; q=1; j=d=c=m=0; word=[]; prior=potential(n,q)
    limit=3*branch_bound(n)+1
    while q:
        y,a=step(n)
        branch=a%3==0
        nextq=q+(2 if branch else -1)
        now=potential(y,nextq)
        require(now<=prior, 'potential increase')
        if branch:
            require(n>=13 and now*28<=prior*27,'branch contraction')
            m+=1
        c=3*c+(1<<d); d+=a; j+=1; word.append(a)
        require((3**j*source+c)==(1<<d)*y,'carry')
        modulus=1<<(d+1)
        residue=((1<<d)-c)*pow(3**j,-1,modulus)%modulus
        require(source%modulus==residue,'source guard')
        if nextq:
            require(y<=source,'active orbit bound')
        n,q,prior=y,nextq,now
        require(j<=limit,'termination clock certificate')
    require(j==3*m+1,'tree size')
    require(d<=j*((3*source+1).bit_length()-1),'division clock bound')
    require(n+1<=F(3,2)*(source+1)*F(27,28)**m,'terminal budget')
    return {'source':source,'target':n,'word':word,'branches':m,'J':j,'D':d,'C':c,'source_residue':residue,'source_modulus':modulus}

def main():
    # Every native edge of every positive odd n through 200001.
    equality=[]
    for n in range(1,200002,2):
        y,a=step(n)
        ratio=F(y+1,n+1)*F(3,2)**(2 if a%3==0 else -1)
        require(ratio<=1,'edge')
        if a%3==0:
            require(ratio<=F(27,28),'branch edge')
            if ratio==F(27,28): equality.append(n)
    require(equality==[13],'sharp branch constant')
    records=[run(n) for n in range(1,20002,2)]
    controls={str(n):run(n) for n in [1,3,13,41,27,165,129,2**100-1]}
    require(controls['1']['target']==1 and controls['3']['target']==5,'hostiles')
    require(controls['41']['target']==31,'41 changed grammar')
    # Rational non-integer fixed source 1/5 has exact valuation 3 forever.
    require((3*F(1,5)+1)/8==F(1,5),'rational hostile')
    require(F(1,5)<13,'arithmetic guard essential')
    # Twelve branches suffice for first-stop native descent.
    require(F(3,2)*F(27,28)**12<1,'descent threshold')
    require(F(3,2)*F(27,28)**11>1,'threshold sharp for this budget')
    print(json.dumps({'edge_inputs':100001,'full_sources':10001,'sharp_branch_input':13,
       'max_test_J':max(r['J'] for r in records),'max_test_branches':max(r['branches'] for r in records),
       'controls':controls},sort_keys=True))
if __name__=='__main__': main()
