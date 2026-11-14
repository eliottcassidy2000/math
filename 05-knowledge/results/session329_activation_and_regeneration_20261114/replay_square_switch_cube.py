from collections import Counter
from math import isqrt

P=(1,8,28,21,4,32,17,19,6,30,34,15,10,26,23,13,12,24,25,11,5,20,29,35,14,2,7,18,31,33,16,9,27,22,3)
Q=(1,8,28,21,4,32,17,19,6,30,34,15,10,26,23,2,7,18,31,33,16,9,27,22,14,35,29,20,5,11,25,24,12,13,3)
# Integers are singleton vertices. Tuples give (offset, reverse).
G=(1,(-1,0),(1,0),(-7,1),(6,0),(-6,0),(0,1),11,(-5,1),5,4,12,
   (-12,0),(12,0),(7,1),(-8,0),(2,1),(-3,0),9,7,(4,0),(-4,0),10,6,
   (5,0),(-11,1),2,(-2,0),8,(3,0),(-9,1),(9,1),(-10,0),(10,0),(11,0),(8,1),3)
def edges(path):
    return {tuple(sorted(e)) for e in zip(path,path[1:])}
def signature(path):
    return frozenset(path[::2])
def check(path):
    if path[0]!=1 or path[-1]!=3 or set(path)!=set(range(1,len(path)+1)):
        raise RuntimeError("not a fixed-port spanning path")
    if any(isqrt(x+y)**2!=x+y for x,y in zip(path,path[1:])):
        raise RuntimeError("not a square-sum path")
def glue(d):
    p=[]
    for x in G:
        if isinstance(x,int): p.append(x); continue
        c,rev=x
        a=[25*v+(-1)**i*c for i,v in enumerate(d[c])]
        p.extend(a[::-1] if rev else a)
    return tuple(p)
def census(n):
    A=[sum(1<<(y-1) for y in range(1,n+1) if x!=y and isqrt(x+y)**2==x+y) for x in range(1,n+1)]
    C=Counter()
    def rec(x,rem,even):
        if not rem:
            C[even]+=1
            return
        b=A[x]&rem
        while b:
            r=b&-b;b-=r
            y=r.bit_length()-1
            if (y==2)!=(rem==r):continue
            rec(y,rem^r,even|r if rem.bit_count()%2==1 else even)
    rec(0,((1<<n)-1)^1,1)
    return C

def main():
    nine_ports=set(range(1,14))|set(range(23,32))
    if {x for x in nine_ports if x!=27 and isqrt(x+27)**2==x+27}!={9}:
        raise RuntimeError("9-block obstruction failed")
    print("9: zero-offset port 27 only meets its opposite port 9")
    check(P);check(Q)
    if signature(P)!=signature(Q) or len(edges(P)-edges(Q))!=3:
       raise RuntimeError('seed pair incorrect')
    for n in range(3,36,2):
        counts=census(n)
        if n<35 and max(counts.values(),default=0)>1:
           raise RuntimeError('earlier switchable same-parity pair')
        if counts: print(n,sorted(counts.values()))
    copies={c:P for c in range(-12,13)}
    base=glue(copies);check(base)
    scars=[]
    for c in copies:
        changed=glue({**copies,c:Q});check(changed)
        if signature(changed)!=signature(base) or len(edges(changed)-edges(base))!=3:
            raise RuntimeError('failed switched block')
        scars.append(edges(changed)^edges(base))
    if len(set.union(*scars))!=25*6: raise RuntimeError("switches not independent")
    print("35: disjoint 3-edge repairs:", edges(P)-edges(Q), edges(Q)-edges(P))
    print("887: 25 independent switches, all fixed 1→3 and same-signature")
    depth2=glue({c:glue({j:Q if ((j+c)%3) else P for j in copies}) for c in copies})
    check(depth2)
    if signature(depth2)!=signature(glue({c:base for c in copies})):
       raise RuntimeError("parity memory lost in composition")
    print("22187: mixed two-level literal replay passed")

if __name__=="__main__":main()
