from itertools import permutations, combinations, product
from math import comb, isqrt
from pathlib import Path
from collections import Counter
import json

from math import isqrt,gcd

def F(x):
 a=1
 while x:
  if not x&1:
   u=x>>1;a+=u+(u&1)
  x>>=1
 return a

def H(x,y):return x*y+F(x)+F(y)-1+((x&1)&(y&1))

def selector(N):
 for x in range(max(1,isqrt(N)//2-2),isqrt(N)+1):
  s=F(x); lo=max(x,(N-s+x-1)//(x+1)); hi=min(4*x,(N-s)//x)
  for y in range(lo,hi+1):
   if N==x*y+s+F(y)-1+((x&1)&(y&1)):
    d=gcd(x,y);a,b=x//d,y//d
    if d!=isqrt((N-1)//(a*b)):raise ArithmeticError((x,y))
    yield x,y,a,b,d

from math import isqrt
from collections import Counter

def paths(n):
 f=(1<<n)-1
 m=[sum(1<<j for j in range(n) if j!=i and isqrt(i+j+2)**2==i+j+2) for i in range(n)]
 E=Counter()
 def dfs(S,v):
  if not S:
   E[v]+=1;return
  B=m[v]&S
  while B:
   b=B&-B;B^=b;dfs(S^b,b.bit_length()-1)
 C={}
 for i in range(n):
  E.clear()
  dfs(f^(1<<i),i)
  for j,k in E.items():
   if j>=i:C[i+1,j+1]=k
 return C

def need(x):
    if not x: raise ArithmeticError("certificate failed")

def hamilton(T):
    n=len(T); f=(1<<n)-1
    need(all(0<=r<=f and not r>>i&1 for i,r in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(n) for j in range(i)))
    a=[[0]*n for s in range(f+1)]
    for i in range(n): a[1<<i][i]=1
    for s in range(1,f+1):
        for j,v in enumerate(a[s]):
            if not v: continue
            e=T[j]&~s
            while e:
                b=e&-e; e^=b
                a[s^b][b.bit_length()-1]+=v
    return sum(a[f])

def reverse(T,i,j):
    out=list(T); out[i]^=1<<j;out[j]^=1<<i
    return tuple(out)

def inflate(T,sizes):
    blocks=[];z=0
    for i in range(len(T)):
        m=sizes.get(i,1);blocks.append(tuple(range(z,z+m)));z+=m
    out=[]
    for i,I in enumerate(blocks):
        v=sum(1<<j for k,J in enumerate(blocks) if T[i]>>k&1 for j in J)
        for a in range(len(I)):out.append(v+sum(1<<j for j in I[a+1:]))
    return tuple(out)

def slot_coeffs(T,w):
    others=[i for i in range(len(T)) if i!=w]; ans=Counter()
    for P in permutations(others):
        req=0; opt=(T[w]>>P[0]&1)+(T[P[-1]]>>w&1)
        for i,j in zip(P,P[1:]):
            can=(T[i]>>w&1)&(T[w]>>j&1)
            if T[i]>>j&1: opt+=can
            elif can:req+=1
            else:break
        else:
            for k in range(opt+1): ans[req+k]+=comb(opt,k)
    ans.pop(0,None)
    return dict(ans)

def both_slot_coeffs(T,w,z):
    out=Counter()
    for i,j in ((1,1),(1,2),(2,1),(2,2),(3,1),(1,3),(2,3),(3,2),(3,3)):
        # For this seven-vertex scaffold both clone labels each have a
        # one-vertex boundary that forces at most two uses. Independently
        # invert the literal small transitive inflation in increasing order.
        a=hamilton(inflate(T,{w:i,z:j}))
        for (r,s),c in out.items():
            if r<=i and s<=j:
                a-=c*surj(i,r)*surj(j,s)
        need(a%(surj(i,i)*surj(j,j))==0)
        out[i,j]=a//(surj(i,i)*surj(j,j))
    return tuple(out[i,j] for i,j in ((1,1),(1,2),(2,1),(2,2))),[v for k,v in out.items() if 3 in k]

def surj(n,r):return sum((-1)**i*comb(r,i)*(r-i)**n for i in range(r+1))

def apex_rows(word):
    return (sum((b=='1')<<(i+1) for i,b in enumerate(word)),
            *(sum(1<<(j+1) for j in range(i+1,len(word)))+(b=='0')
              for i,b in enumerate(word)))

def rc(s):return ''.join('1' if b=='0' else '0' for b in s[::-1])

def fermat_apex(N):
    u=isqrt(N//8)
    return [(x,y) for x in (2*u,2*u-1) if x>=0 for y in (2*u+2,2*u+1)
            if 2*x*y+2*u+F(x)+F(y)-1==N]

def prime(n):return n>1 and all(n%d for d in range(2,isqrt(n)+1))

def canonical(rows):
    n=len(rows)
    cells=[tuple(i for i,r in enumerate(rows) if r.bit_count()==k) for k in range(n)]
    return min(tuple(sum((rows[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        for choice in product(*(permutations(c) for c in cells if c))
        for p in [sum(choice,())])

def exhaustive_census():
    layer={(0,)}
    target={9:0,25:0,81:0,289:0}
    for n in range(1,9):
        need(len(layer)==(1,1,2,4,12,56,456,6880)[n-1])
        for T in layer:
            h=hamilton(T)
            if h in target and any(hamilton(reverse(T,i,j))==h+2
                    for i,j in combinations(range(n),2)):
                target[h]+=1
        if n<8:
            layer={canonical(tuple(v+((1-(m>>i&1))<<n) for i,v in enumerate(T))+(m,))
                   for T in layer for m in range(1<<n)}
    need(list(target.values())==[20,3,2,6])
    print("Complete isomorphism census through order eight: PASS",flush=True)

def run(large_square_sums):
    root=Path(__file__).parent
    data=json.loads((root/'flip_edges.json').read_text())
    lows={p:set() for p in range(1,5)}
    for p,A,B,i,j in data:
        need(tuple(B)==reverse(A,i,j))
        a=(2**p+1)**2
        need(hamilton(A)==a and hamilton(B)==a+2)
        lows[p].add(tuple(A))
    need([len(lows[p]) for p in range(1,5)]==[20,3,2,6])
    A=(8,13,25,176,99,71,143,55);B=reverse(A,4,7)
    need(slot_coeffs(A,0)=={1:289,2:151} and slot_coeffs(B,0)=={1:291,2:151})
    for m in range(1,5):
        need(hamilton(inflate(A,{0:m}))==(151<<m)-13)
        need(hamilton(inflate(B,{0:m}))==(151<<m)-11)
    need(((151<<18)-11)==7*5654819 and prime(5654819))
    U=(118,124,120,49,96,64,8);V=reverse(U,2,3)
    need(both_slot_coeffs(U,1,4)==((25,10,10,4),[0]*5))
    need(both_slot_coeffs(V,1,4)==((27,12,9,4),[0]*5))
    for i in range(1,4):
        for j in range(1,4):
            q=(1+(2<<i))*(1+(2<<j))
            need(hamilton(inflate(U,{1:i,4:j}))==q)
            need(hamilton(inflate(V,{1:i,4:j}))==q+(2<<j)-(1<<i))
    D=[]
    for p in range(1,501):
        for x,y in fermat_apex((2**p+1)**2):D.append((p,(x,y)))
    need(D==[(1,(1,3)),(2,(2,4))])
    N=7*30013; need(prime(30013))
    need(next(selector(N))==(275,762,275,762,1))
    w=bin(275)[2:]+rc(bin(762)[2:])
    need(w=='1000100111010000010' and H(275,762)==210091)
    need(hamilton(apex_rows(w))==N)
    for p in range(1,501):
        n=(2**p+1)**2
        A=isqrt(2*n-1)**2==2*n-1 or isqrt(2*n-2)**2==2*n-2
        need(A==(p in (1,2,4)))
    P=paths(27)
    need(sum(P.values())==35 and len(P)==16 and all(18 in s for s in P))
    need(Counter(P.values())=={1:4,2:6,3:5,4:1})
    if large_square_sums:
        result=[sum(paths(n).values()) for n in range(1,36)]
        need(result==[1]+[0]*13+[1,1,1]+[0]*5+[3,0,10,12,35,52,19,20,349,392,669,4041,17175])
    print("Fermat flip witnesses, both modular frames, one-apex inverse, and square-sum controls: PASS")

if __name__=="__main__":
    import sys
    run("--square-sum" in sys.argv)
    if "--census" in sys.argv: exhaustive_census()
