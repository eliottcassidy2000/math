#!/usr/bin/env python3
"""Independent integer first-eligible PRIMEGAME and square-recognizer audit.
No producer executor is imported; producer outputs are compared only at end.
"""
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
P=(2,3,5,7,11,13,17,19,23,29)
F=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
def check(ok,message):
    if not ok: raise RuntimeError(message)
def fac(n):
    es=[]
    for p in P:
        e=0
        while n%p==0:n//=p;e+=1
        es.append(e)
    check(n==1,('factor',n))
    return tuple(es)
V=tuple((fac(a),fac(b)) for a,b in F)
def transition(x):
    for row,(a,b) in enumerate(V):
        if all(z>=w for z,w in zip(x,b)):
            return tuple(z+u-v for z,u,v in zip(x,a,b)),row
    raise RuntimeError('no fraction')
def raw(n):
    for i,(a,b) in enumerate(F):
        if n%b==0:return n//b*a,i
    raise RuntimeError('no fraction')
def budgets(x):return x[0]+x[2]+x[8],x[1]+x[3]+x[6]+x[7]+x[8]+x[9],sum(x[4:])
def canonical(x):
    a,b,c,d,*q=x
    return a==d==0 and b>=1 and c==b+1 and q==[1,0,0,0,0,0]
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
memo={}
max_path=0
actual_periods=set()
def independent_basin(x):
    global max_path
    seen={};path=[]
    while True:
        if x in memo: label=memo[x];break
        if canonical(x):label='prime';break
        if x in seen:
            period=len(path)-seen[x];actual_periods.add(period)
            check(period in (2,6),('unexpected cycle',period,x))
            label='cycle'+str(period);break
        seen[x]=len(path);path.append(x)
        old=budgets(x);nxt,row=transition(x);new=budgets(nxt)
        if old[2]>0:
            check(row<11,'positive control enabled L/M/N')
            check(new[0]==old[0] and new[1]<=old[1] and new[2]<=old[2],('budgets',x,nxt))
        x=nxt
    max_path=max(max_path,len(path))
    for state in path:memo[state]=label
    return label

def main():
    producer=load(ROOT/'primegame_machine/primegame_basin.py','producer_prime')
    onehot=0
    for data in itertools.product(range(5),repeat=4):
        for q in range(7):
            ctrl=tuple(int(j==q-1) for j in range(6)) if q else (0,)*6
            x=data+ctrl
            expected=independent_basin(x)
            got=producer.classify_vector(x)['basin']
            check(expected==got,('onehot',x,expected,got));onehot+=1
    capped=0
    for x in itertools.product(range(3),repeat=10):
        expected=independent_basin(x)
        got=producer.classify_vector(x)
        check(expected==got['basin'],('capped',x,expected,got))
        check(got['transitions']<=got['bound'],('bound',x));capped+=1
    rawcases=0
    for n in range(1,2001):
        # Strip spectators independently before vector comparison.
        residual=n;es=[]
        for p in P:
            e=0
            while residual%p==0:residual//=p;e+=1
            es.append(e)
        y,row=transition(tuple(es));v,i=raw(n)
        check(i==row and v==residual*math.prod(p**e for p,e in zip(P,y)),('raw',n))
        rawcases+=1
    # Square machine: parse the actual emitted fixed fraction list and use only
    # whole positive integers with a global ordered divisibility scan.
    square=[]
    for token in (ROOT/'square_machines/square.fractions.txt').read_text().strip().split(','):
        a,b=map(int,token.strip().split('/'));check(math.gcd(a,b)==1,'unreduced');square.append((a,b))
    squarecases=0;maxsteps=0
    for n in range(1025):
        value=7*(1<<n);steps=0
        while True:
            for a,b in square:
                if value%b==0:value=value//b*a;steps+=1;break
            else:break
        k=math.isqrt(n);r=n-k*k
        expected=31 if n==0 else (31*5**(2*k-1) if r==0 else 37*3**(2*k-r)*5**(r+1))
        cost=1 if n==0 else 5*n+2 if r==0 else 3*n+2*k*k+4*k+4
        check(value==expected and steps==cost,('square',n,value,expected,steps,cost))
        squarecases+=1;maxsteps=max(maxsteps,steps)
    print(json.dumps({'status':'PASS','primegame_onehot':onehot,'primegame_cap_0_2':capped,'primegame_independent_states':len(memo),'primegame_max_unmemoized_path':max_path,'observed_periods':sorted(actual_periods),'integer_transition_cases':rawcases,'square_literal_cases':squarecases,'square_maxsteps':maxsteps},indent=2))
if __name__=='__main__':main()
