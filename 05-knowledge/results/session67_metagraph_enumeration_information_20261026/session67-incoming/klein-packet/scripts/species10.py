"""n=10 species check: primes of size <=9 with full cycle indices (p9.txt) + the size-10 prime count (each contributes 1)."""
import os
DATA=os.environ.get("DATA",".")
import re, collections
from fractions import Fraction
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'species.py')).read().split("print(f'distinct primes by k:'")[0].replace("NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 7","NMAX=10").replace("files=sys.argv[2:] if len(sys.argv)>2 else [f'p{NMAX}.txt']","files=[os.path.join(DATA,'p9.txt')]"))
P10=7906161; P10SC=5525
def polya(cts, series, n):
    tot=Fraction(0)
    for ct in cts:
        poly=[Fraction(0)]*(n+1); poly[0]=Fraction(1)
        for l,c in ct.items():
            base=[Fraction(0)]*(n+1)
            for j in range(1,n+1):
                if l*j<=n: base[l*j]=Fraction(series[j])
            for _ in range(c):
                new=[Fraction(0)]*(n+1)
                for a in range(n+1):
                    if poly[a]==0: continue
                    for b in range(n+1-a):
                        if base[b]: new[a+b]+=poly[a]*base[b]
                poly=new
        tot+=poly[n]
    return tot/len(cts)
a=[0]*11; s=[0]*11; a[1]=1; s[1]=1
def seqpow(s,n):
    f=[Fraction(0)]*(n+1); f[0]=Fraction(1)
    for m in range(1,n+1): f[m]=sum(Fraction(s[j])*f[m-j] for j in range(1,m+1))
    return f[n]
for n in range(2,11):
    sn=sum(polya(primes[p][0],a,n) for p in primes if p[0]<=n) + (P10 if n==10 else 0)
    s[n]=sn; a[n]=seqpow(s,n)
print('A000568 from species (n=1..10):', [a[n] for n in range(1,11)])
print('strong (n=1..10):', [s[n] for n in range(1,11)])
c=[0]*11; cs=[0]*11; c[1]=1; cs[1]=1
def euler_mset(cs,n):
    f=[Fraction(0)]*(n+1); f[0]=Fraction(1)
    for j in range(1,n+1):
        for _ in range(int(cs[j])):
            for m in range(j,n+1): f[m]+=f[m-j]
    return f[n]
for n in range(2,11):
    tot=Fraction(0)
    for (k,key),(A,Apm) in primes.items():
        if k>n: continue
        w=Fraction(1) if len(Apm)>len(A) else Fraction(1,2)
        tot+=w*polya(Apm,c,n)
    if n==10: tot+=Fraction(P10+P10SC,2)
    cs[n]=tot; c[n]=euler_mset(cs,n)
print('C3-structure classes c(n), n=1..10:', [c[n] for n in range(1,11)])
print('strong C3-structure classes, n=1..10:', [cs[n] for n in range(1,11)])
