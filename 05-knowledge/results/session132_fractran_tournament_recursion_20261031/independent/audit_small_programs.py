#!/usr/bin/env python3
"""Read the two short programs, then run literal ordered fraction semantics.
Does not import either generating/verifying implementation. Standard library.
"""
import ast
from math import comb,gcd,prod
from pathlib import Path
import json
import sys
from audit_fractions import factor, require

ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent


def step(v,fs):
    for a,b in fs:
        if (v*a)%b==0:
            return v*a//b
    return None


def direct16():
    source=(ROOT/'fractran/literal16.py').read_text()
    tree=ast.parse(source)
    fs=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
            and any(isinstance(t,ast.Name) and t.id=='FRACTIONS' for t in n.targets))
    require(len(fs)==16 and all(gcd(a,b)==1 for a,b in fs),'Invalid16fraction artifact')
    result=[]
    for n in range(3,8):
        val=306*13**(n-3);steps=0
        while True:
            nxt=step(val,fs)
            if nxt is None:break
            val=nxt;steps+=1
            require(steps<200000,'16fraction progress bound')
        expected=2**((n-1)*(n-2)//2-1)
        require(val==59*2**expected*3**(n-1),'16fraction terminal integer')
        t=n-3
        require(steps==6*expected-5+(7*t*t+29*t)//2,'16fraction exact step count')
        result.append(dict(n=n,output=expected,steps=steps))
    return result


def stream22():
    fs=[]
    for line in (ROOT/'pure_blue_fractran/pure_blue.fractran').read_text().splitlines():
        line=line.split('#')[0].strip()
        if line:fs.append(tuple(map(int,line.split('/'))))
    require(len(fs)==22 and all(gcd(a,b)==1 for a,b in fs),'Invalid22fraction artifact')
    data=(2,3,5,7,11,13,17,19)
    controls=(23,29,31,37,41,43,47,53,59,61,67,71,73,79,83)
    require(not set(data)&set(controls),'Overlap')
    for a,b in fs:
        an,bn=factor(a),factor(b)
        qa=[(p,e) for p,e in an.items() if p in controls]
        qb=[(p,e) for p,e in bn.items() if p in controls]
        require(len(qa)==len(qb)==1 and qa[0][1]==qb[0][1]==1 and qa[0][0]!=qb[0][0], 'Lost one-hot guard')
    def narayana(k):
        return sum(comb(k-2*j,j) for j in range(k//3+1))
    def expected(k):
        return sum(narayana((k-c)//2) for c in [0,1,3,5] if k>=c and (k-c)%2==0)
    val=817086270;steps=0;result=[]
    for n in range(25):
        target=23*19**n*prod(p**expected(n+j) for j,p in enumerate(data[:6]))
        require(val==target,'Checkpoint integer mismatch')
        result.append(dict(n=n,output=expected(n),steps=steps))
        if n==24:break
        start=steps
        while True:
            val=step(val,fs);steps+=1
            require(val is not None,'Unexpected stream halt')
            if val%23==0:break
        a,b,c,d,e,f=[expected(n+j) for j in range(6)]
        require(steps-start==8+2*(2*a+b+c+d+2*e+f),'Checkpoint budget mismatch')
    # Exhaust all 0/1 windows, including zero-register branches.
    for mask in range(64):
        window=tuple((mask>>j)&1 for j in range(6))
        val=23*prod(p**k for p,k in zip(data,window))
        while True:
            val=step(val,fs)
            require(val is not None,'Unexpected zero-window halt')
            if val%23==0:break
        target=window[1:]+(window[0]+window[4],)
        require(val==23*19*prod(p**k for p,k in zip(data,target)),'Incorrect transfer')
    return dict(checkpoints=result,zero_one_windows=64)


print(json.dumps(dict(status='PASS',direct16=direct16(),stream22=stream22()),indent=2))
