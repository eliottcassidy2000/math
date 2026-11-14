from itertools import combinations, permutations, product
from math import comb, factorial, gcd
from fractions import Fraction
from collections import Counter

def check(x):
    if not x: raise RuntimeError("control failed")

def divisors(n):
    return [d for d in range(1,n+1,2) if n%d==0]

def rot_exp(n,d):
    return n*(n-3)//(2*d)

def ref_exp(n,d):
    m=n//d
    return (n-3)*(m+1)//4 if n%2 else (n//2-2)*(m//2)+n//4-1+(m+3)//4

def sectors(n):
    D=divisors(n)
    ex={}; blue={}
    for d in reversed(D):
        ex[d]=2**rot_exp(n,d)-sum(ex[e] for e in D if e>d and not e%d)
        blue[d]=2**ref_exp(n,d)-sum(blue[e] for e in D if e>d and not e%d)
    a={d:Fraction(d*ex[d],n) for d in D}
    b={d:Fraction(blue[d],1+(n%2==0)) for d in D}
    return a,b

def literal(n):
    pairs=list(combinations(range(n),2))
    free=[(a,b) for a,b in pairs if b-a not in (1,n-1)]
    forced={(a,a+1) for a in range(n-1)}
    forced.add((n-1,0))
    out=Counter(); outb=Counter()
    all_rot={k:[((i+k)%n,(j+k)%n) for i,j in free] for k in range(n)}
    anti={s:[((s-j)%n,(s-i)%n) for i,j in free] for s in range(n)}
    for z in product((0,1),repeat=len(free)):
        D=forced|{p if v else p[::-1] for p,v in zip(free,z)}
        d=sum(all((q in D)==v for q,v in zip(L,z)) for L in all_rot.values())
        b=sum(all((q in D)==v for q,v in zip(L,z)) for L in anti.values())
        check(b in (0,d))
        out[d]+=1; outb[d]+=bool(b)
    return ({d:Fraction(d*out[d],n) for d in divisors(n)},
            {d:Fraction(d*outb[d],n) for d in divisors(n)})

def multi(n,k):
    out=factorial(n-1)
    for t in range(1,(k-1)//2+1):
        out+=factorial(n-2*t-1)*sum((1<<r)*comb(t-1,r-1)*comb(k-2*t,r)
                                       for r in range(1,min(t,k-2*t)+1))
    return out

def literal_multi(n,k):
    out=0
    for p in permutations(range(1,n)):
        r=1
        for a,b in zip((0,)+p,p+(0,)):
            if a==b+1 and a<k:
                r=0;break
            if b==a+1 and b<k:r*=2
        out+=r
    return out

for n in range(3,8):
    a,b=sectors(n)
    check((a,b)==literal(n))
    check(sum((n//d)*a[d] for d in a)==2**rot_exp(n,1))
    check(sum(b.values())*(1 if n%2 else 2)==2**(((n-1)**2)//4-1))
    print(n,sum(a.values()),sum(b.values()),multi(n,n))
for n in range(3,10):
    for k in range(1,n+1):
        check(multi(n,k)==literal_multi(n,k))
for n in range(3,101):
    a,b=sectors(n)
    check(all(0<=b[d]<=a[d] and a[d].denominator==b[d].denominator==1 for d in a))
    check(sum(b.values())*(1 if n%2 else 2)==2**(((n-1)**2)//4-1))
    check(sum(a.values())*n==sum(sum(gcd(j,d)==1 for j in range(1,d+1))*2**rot_exp(n,d)
                                  for d in divisors(n)))
print("Exact sectors and blue closure counts; literal chord orientations through n=7; independent mixed cycle moments through n=9.")

def delta(n):
    z=Fraction(2*(n%2==0))
    for t in range(1,(n-1)//2+1):
        z+=factorial(n-2*t-1)*sum(Fraction((1<<r)*2*t,r)*comb(t-1,r-1)*comb(n-2*t-1,r-1)
                                  for r in range(1,min(t,n-2*t)+1))
    return z/2**(n-1)
def literal_delta(n):
    z=0
    for p in permutations(range(1,n)):
        t=0; r=1
        for a,b in zip((0,)+p,p+(0,)):
            if a==n-1 and b==0: t=1
            if b==n-1 and a==0: t=-1
            if a==b+1:
                r=0; break
            if b==a+1:r*=2
        z+=t*r
    return Fraction(z,2**(n-1))
for n in range(3,10):
    check(delta(n)==literal_delta(n) and delta(n)>0)
A=(22,60,8,17,36,13); B=(42,4,57,50,35,2)
for i,j in combinations(range(6),2):
    bit=(A[i]>>j)&1; bit2=(B[i]>>j)&1
    check(bit==bit2 if j==i+1 else bit!=bit2)
def literal_cycles(a):
    return sum(all((a[i]>>j)&1 for i,j in zip((0,)+P,P+(0,)))
               for P in permutations(range(1,len(a))))
check((literal_cycles(A),literal_cycles(B))==(2,3))
print("Literal metagraph mean-bias and opposite-sign line checked.")
from functools import cache

def spec(rows):
    n=len(rows)
    ph=tuple(sum(all((rows[i]>>j)&1 for i,j in zip(P,P[1:]))
                 for P in permutations(range(n),k)) for k in range(1,n+1))
    cy=tuple(sum(all((rows[i]>>j)&1 for i,j in zip((S[0],)+P,P+(S[0],)))
                 for S in combinations(range(n),k) for P in permutations(S[1:]))
             for k in range(3,n+1))
    aut=sum(all(((rows[i]>>j)^(rows[P[i]]>>P[j]))&1==0
                for i,j in combinations(range(n),2)) for P in permutations(range(n)))
    return ph,cy,aut

def hp(rows):
    n=len(rows)
    @cache
    def p(s,v):
        t=s^(1<<v)
        return sum(p(t,i) for i in range(n) if (rows[v]&t)>>i&1) if t else 1
    return sum(p((1<<n)-1,v) for v in range(n))

def min_neighbor(rows):
    n=len(rows); all_=(1<<n)-1
    H=[]; closed=[]
    for p in permutations(range(n)):
        if not all((rows[i]>>j)&1 for i,j in zip(p,p[1:])):
            continue
        Q=[all_^x^(1<<i) for i,x in enumerate(rows)]
        for i,j in zip(p,p[1:]):
            Q[i]^=1<<j; Q[j]^=1<<i
        h=hp(tuple(Q))
        H.append(h)
        if (rows[p[-1]]>>p[0])&1: closed.append(h)
    return min(H),min(closed),len(H),len(closed)

A=(34,68,33,87,71,90,5); B=(70,64,98,39,47,3,56)
check(spec(A)==spec(B)==((7,21,53,107,161,151,65),(9,15,16,11,3),1))
sa,sb=[sorted(x.bit_count() for x in Z) for Z in (A,B)]
check(sa!=sb and sorted(6-x for x in sa)!=sb)
for T in (A,B):
    check(not any(all(((T[i]>>j)^(T[P[i]]>>P[j]))&1
                      for i,j in combinations(range(7),2))
                  for P in permutations(range(7))))
check(min_neighbor(A)==(33,33,65,21) and min_neighbor(B)==(23,23,65,21))
print("All-length path/cycle spectral collision survives converse and has distinct metagraph incidence spectra.")
