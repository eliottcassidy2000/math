#!/usr/bin/env python3
"""Independent Fraction/polynomial audit of the frozen weight22 certificate.
Reads only immutable coefficient data, never imports the producer or SymPy.
"""
import ast
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

checks=0

def need(ok,why):
    global checks
    checks+=1
    if not ok: raise RuntimeError(why)

def clean(p):return {k:v for k,v in p.items() if v}
def add(a,b):
    c=dict(a)
    for k,v in b.items():c[k]=c.get(k,F(0))+v
    return clean(c)
def scale(a,v):return clean({k:c*v for k,c in a.items()})
def mul(a,b):
    c={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            if j+l<=16:c[i+k,j+l]=c.get((i+k,j+l),F(0))+v*w
    return clean(c)
def power(a,n):
    p={(0,0):F(1)}
    for _ in range(n):p=mul(p,a)
    return p
def derivative(a,axis):
    out={}
    for key,v in a.items():
        n=key[axis]
        if n:
            k=list(key);k[axis]-=1;out[tuple(k)]=v*n
    return out
def bracket(a,b):return add(mul(derivative(a,0),derivative(b,1)),scale(mul(derivative(a,1),derivative(b,0)),-1))
def trunc(a,n):return {k:v for k,v in a.items() if k[1]<=n}
def literal(v):return {(0,0):F(str(v))}
def parse(expr):
    def go(n):
        if isinstance(n,ast.Name):return {'x':{(1,0):F(1)},'t':{(0,1):F(1)}}[n.id]
        if isinstance(n,ast.Constant):return literal(n.value)
        if isinstance(n,ast.UnaryOp):return scale(go(n.operand),-1 if isinstance(n.op,ast.USub) else 1)
        if isinstance(n,ast.BinOp):
            a,b=go(n.left),go(n.right)
            if isinstance(n.op,ast.Add):return add(a,b)
            if isinstance(n.op,ast.Sub):return add(a,scale(b,-1))
            if isinstance(n.op,ast.Mult):return mul(a,b)
            if isinstance(n.op,ast.Div):return scale(a,1/b[(0,0)])
            if isinstance(n.op,ast.Pow):return power(a,int(b[(0,0)]))
        raise TypeError(ast.dump(n))
    return go(ast.parse(expr,mode='eval').body)

def main():
    root=Path(__file__).resolve().parents[1]
    raw=(root/'05-knowledge/results/planar_jc48_sep06_weight14.out').read_bytes()
    need(hashlib.sha256(raw).hexdigest()=='74d12cf42c23d9402fdb892faab12fd4ce7309d43561b3dc63ee92e9c6b2e132','exact coefficient certificate pinned')
    data=json.loads(next(s for s in raw.decode().splitlines() if s.startswith('{')))
    # Independently transcribed from THM4308 equation4 after Phi0,K=-32/5.
    A=parse('1+x**2/4+(4/3+2*x**2)*t+(-32/9-4*x**2/5)*t**2+(2176/135+(1088/315+128/35)*x**2-32*x**4/15)*t**3')
    C=parse('-3*x/4-x**3/8+(-4*x-3*x**3/2)*t+(88*x/15-12*x**3/5)*t**2+((-8128/315-192/35)*x+(736/105-96/35)*x**3+8*x**5/5)*t**3')
    dA={};dC={}
    for rows,dest in [(data['delta_A'],dA),(data['delta_C'],dC)]:
        for n,expr in rows.items():
            term=mul(parse(expr),parse('t**'+n));dest.update(term)
    p=parse('t*(1+x**2*t)');y=mul(parse('x*t'),p);u=parse('x**2*t')
    K=add(add(mul(power(p,5),power(y,4)),scale(mul(p,power(y,6)),-F(7,10))),scale(mul(power(p,2),power(y,6)),-F(508,135)))
    T=add(K,scale(mul(power(p,3),power(y,6)),F(235202,27945)))
    j=add(bracket(A,dC),bracket(dA,C))
    variation=add(add(scale(mul(C,dC),2),scale(mul(power(A,2),dA),-3)),scale(dA,F(3,4)))
    need(not trunc(j,14),'every bracket coefficient through t14')
    need(not trunc(add(variation,scale(T,-1)),15),'every defect coefficient through t15')
    need(not trunc(bracket(dA,dC),14),'quadratic finite action bracket')
    need(not trunc(add(power(dC,2),scale(mul(A,power(dA,2)),-3)),15),'quadratic finite action defect')
    need(not trunc(power(dA,3),15),'cubic finite action defect')
    for depth,terms,target in [(2,data['P2_witness'],dA),(3,data['P3_witness'],dC)]:
        witness={}
        for (a,b,c,e),co in terms:
            need(a+b<=depth,'actual P_depth generator')
            need(b+c+2*e>=12,'no earlier witness row')
            term=mul(mul(power(parse('x'),a),power(u,b)),mul(power(p,c),power(y,e)))
            witness=add(witness,scale(term,F(co)))
        need(not trunc(add(witness,scale(target,-1)),15),'complete depth representative equality')
    for k in dA:need(k[1]>=12 and k[0]<=k[1]+1,'A row and degree support')
    for k in dC:need(k[1]>=12 and k[0]<=k[1]+2,'C row and degree support')
    # The first background row not retained really enters the next equations.
    need(bool({k:v for k,v in bracket(parse('t**4'),dC).items() if k[1]==15}),'next bracket row sees unretained background')
    need(bool({k:v for k,v in mul(mul(A,parse('t**4')),dA).items() if k[1]==16}),'next defect row sees unretained background')
    residual=trunc(add(variation,scale(K,-1)),15)
    need(residual==parse('235202*x**6*t**15/27945'),'exact lower-packet response and sign')
    need(F(-27945,235202)*F(235202,27945)==-1,'cancels old weight24 term')
    need(F(27945,235202)*F(145,24)==F(1350675,1881616),'parameter section constant')
    print('PASS independent Fraction coefficient convolution, full depth representatives and finite action')
    print('SCOPE fixed low jet; first unretained background enters next row; no termination claim')
    print('EXPLICIT_GATES',checks)

if __name__=='__main__':main()
