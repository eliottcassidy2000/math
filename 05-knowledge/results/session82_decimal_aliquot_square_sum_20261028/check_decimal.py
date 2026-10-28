#!/usr/bin/env python3
"""Exact, dependency-free Session82 decimal-family certificate."""
from math import gcd,isqrt,prod
import json

def check(p,msg):
    if not p: raise RuntimeError(msg)
def factor(n):
    out={}; p=2
    while p*p<=n:
        while n%p==0: out[p]=out.get(p,0)+1; n//=p
        p+=1
    if n>1: out[n]=out.get(n,0)+1
    return out

def sigma(n): return prod((p**(e+1)-1)//(p-1) for p,e in factor(n).items())
def s(n): return sigma(n)-n
def paired(n):
    out=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            out+=d
            if d*d!=n: out+=n//d
    return out-n
def square(n): return n>=0 and isqrt(n)**2==n
def edges(v): return [(i,j) for i in range(len(v)) for j in range(i+1,len(v)) if square(v[i]+v[j])]
def family(t):return t*t+1,t+1,t*(t+1)
def petal(t):
    h=(t*t+1)**2; d=8*t*(t*t-1)
    return 2*(h-d),2*h,2*(h+d)

def main():
    for n in range(1,5001): check(s(n)==paired(n),('sigma',n))
    pell=[]
    for t in range(2,10001):
        a,b,c=family(t)
        check(a+b==c+2 and b+c==b*b,('identities',t))
        check(gcd(a,b)==gcd(a,c)==gcd(2,t+1) and gcd(b,c)==b,('gcd',t))
        check(not square(a+b),('AB exclusion',t))
        r=isqrt(a+c)
        check(square(a+c)==((4*t+1)**2-8*r*r==-7),('Pell',t))
        if square(a+c): pell.append([t,r])
        check((a-2)**2+(2*t)**2==a*a,('Pythagorean',t))
        if t<=1000:
            check(s(c)==t*s(b)+b*s(t)+s(t)*s(b),('aliquot transport',t))
        if t>=8:
            x,y,z=petal(t)
            check(x>0 and x+z==2*y and edges((x,y,z))==[(0,1),(0,2),(1,2)],('petal',t))
        if 2<=t<=7: check(petal(t)[0]<0,('positivity hostile',t))
    check(pell==[[16,23],[45,64],[552,781],[1537,2174]],'finite Pell list')
    base=[(x,isqrt((x*x+7)//8)) for x in range(1,45,4) if square((x*x+7)//8)]
    check(base==[(1,1),(5,2)],'Pell descent base')
    for t,r in [(0,1),(1,2)]:
        for _ in range(30):
            check(2*t*t+t+1==r*r,'Pell recurrence')
            old=(4*t+1,r)
            t,r=17*t+12*r+4,24*t+17*r+6
            x=4*t+1
            check((17*x-48*r,17*r-6*x)==old,'Pell inverse')
            check(0<17*x-48*r<x and 17*r-6*x>0,'Pell descent')
    for t in [2,4,16,256]:
        check(len(factor(t+1))==1 and list(factor(t+1).values())==[1], 'prime guard')
        check(s(t+1)+s(t*(t+1))==(t+1)**2-2, 'guarded edge loss')
    vals=family(1000); images=tuple(s(v) for v in vals)
    check(vals==(1000001,1001,1001000),'anchor')
    check(images==(10003,343,2143960),'anchor aliquot')
    for n in vals: check(s(n)==paired(n),('anchor independent',n))
    check(edges(vals)==[(1,2)] and edges(images)==[],'anchor edge loss')
    check(edges((3,6))==[(0,1)] and edges((s(3),s(6)))==[],'small hostile')
    check(prod(factor(prod(vals)))==10010010010,'radical')
    for b in range(2,101):
        check(b**6+1==(b*b+1)*(b**4-b*b+1),'base factor A')
        check(b**3+1==(b+1)*(b*b-b+1),'base factor B')
    check(pow(10,136,289)==288 and pow(10,272,289)==1 and pow(10,16,289)==137,'order witnesses')
    check([k for k in range(136) if (pow(10,2*k,289)+1)%289==0]==[68],'complete residue window')
    for k in range(1,7):
        a=10**(2*k)+1; f=factor(a); rad=prod(f)
        check(all(p%4==1 and p>=13 for p in f),'prime factor restriction')
        check((10*rad<a)==any(e>=2 for e in f.values()),'radical iff')
    for n in range(2,1001):
        a,b,c=family(n)
        check(a-2+b-2==c-2,'translated additive')
        check(gcd(gcd(a-2,b-2),c-2)==n-1,'translation gcd')
        check(tuple(v//(n-1) for v in (a-2,b-2,c-2))==(n+1,1,n+2),'translation primitive')
        check((b//b,c//b,b*b//b)==(1,n,n+1),'edge primitive')
    print(json.dumps({'status':'PASS','family_range':[2,10000],'sigma_crosscheck':[1,5000],
      'anchor':[{'n':n,'factorization':factor(n),'aliquot':s(n)} for n in vals],
      'anchor_edges':edges(vals),'aliquot_image_edges':edges(images),
      'anchor_AP_triangle':petal(1000),'pell_t_le_10000':pell},indent=2))
if __name__=='__main__':main()
