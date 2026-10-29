"""Independent reciprocal-orbit count and named Golay generators audit."""
from collections import Counter
from math import gcd,comb
import json

def need(b,m):
    if not b:raise ArithmeticError(m)

def ord2(d):
    t=1;k=0
    while True:
        t=2*t%d;k+=1
        if t==1:return k

def divisors(n):return [d for d in range(3,n+1,2) if n%d==0]
def prime(n):return n>1 and all(n%d for d in range(2,int(n**0.5)+1))

def cosets(n):
    left=set(range(1,n));out=[]
    while left:
        t=min(left);c=set()
        while t not in c:c.add(t);t=2*t%n
        left-=c;out.append(frozenset(c))
    return out

def code(n,g):
    basis=[g<<j for j in range((n-1)//2)]
    need(all(v.bit_length()<=n for v in basis),'generator degree')
    need(all((u&v).bit_count()%2==0 for u in basis for v in basis),'orthogonal basis')
    words={0}
    for v in basis:words|={w^v for w in words}
    need(len(words)==2**len(basis),'dimension')
    mask=(1<<n)-1
    need(all(((w<<1)&mask)|(w>>(n-1)) in words for w in words),'cyclic')
    extended=words|{w^((1<<(n+1))-1) for w in words}
    punctured={w&mask for w in extended}
    need(len(punctured)==len(extended),'puncture injective')
    return Counter(w.bit_count() for w in extended),min(w.bit_count() for w in punctured if w)

def main():
    samples={}
    for n in range(3,512,2):
        cs=cosets(n);obstruction=any(frozenset(-x%n for x in C)==C for C in cs)
        prime_obstruction=any(ord2(p)%2==0 for p in divisors(n) if prime(p))
        need(obstruction==prime_obstruction,'prime obstruction equivalence')
        phi=lambda d:sum(gcd(i,d)==1 for i in range(d))
        factors=sum(phi(d)//ord2(d) for d in divisors(n))
        need(factors==len(cs),'factor count')
        r=0 if obstruction else 2**(factors//2)
        if r:need(factors%2==0,'pairs')
        if n in [7,21,23,31,49,71,73,89,127]:samples[n]=r
    hamming={0:1,4:14,8:1};golay={0:1,8:759,12:2576,16:759,24:1}
    for n,gs,expected,d in [(7,[29,23],hamming,3),(23,[7973,5279],golay,7)]:
        for g in gs:
            enum,mind=code(n,g)
            need(dict(enum)==expected and mind==d,'enumerator and punctured distance')
            k=(n+1)//2;r=(d-1)//2
            need((1<<k)*sum(comb(n,j) for j in range(r+1))==1<<n,'perfect balls')
    print(json.dumps({'odd_n_checked':255,'upper_n':511,'counts':samples,
        'Golay_weight_enumerator':golay,'independent_generators':[7973,5279]},sort_keys=True))
if __name__=='__main__':main()
