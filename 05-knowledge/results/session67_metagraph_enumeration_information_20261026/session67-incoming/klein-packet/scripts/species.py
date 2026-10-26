"""Species check: A000568 from the prime census via Polya substitution, and the C3-structure (intransitivity)
class counts via Aut+- cycle indices.  Reads p{n}.txt files written by tourn_enum (--primes)."""
import os
DATA=os.environ.get("DATA",".")
import sys, re, math, csv, collections
from fractions import Fraction
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 7
files=sys.argv[2:] if len(sys.argv)>2 else [f'p{NMAX}.txt']
# parse primes: key -> (k, list of cycle types for Aut, list for Aut+- )
primes={}
for fn in files:
    for line in open(os.path.join(DATA,fn)):
        if not line.startswith('PRIME'): continue
        m=re.match(r'PRIME k=(\d+) key=(\d+) aut=(\d+) AUT:(.*?) SC=(\d)(?: ANTI:(.*?))? ARCS:',line)
        k=int(m.group(1)); key=int(m.group(2)); aut=int(m.group(3)); sc=int(m.group(5))
        def parse(s):
            cts=[]
            for tok in s.split():
                ct={}
                for part in tok.strip(',').split(','):
                    if not part: continue
                    l,e=part.split('^'); ct[int(l)]=int(e)
                cts.append(ct)
            return cts
        A=parse(m.group(4)); assert len(A)==aut
        anti=parse(m.group(6)) if sc else []
        primes[(k,key)]=(A,A+anti)
print(f'distinct primes by k:',collections.Counter(k for k,_ in primes))
def polya(cts, series, n):
    """[x^n] of (1/|G|) sum_g prod_l series(x^l)^{c_l(g)} where series is list s[0..n], s[0]=0"""
    tot=Fraction(0)
    for ct in cts:
        # product over cycle lengths l of (sum_j s_j x^{lj})^{c_l}
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
# A000568 recursion: a[n]; s[n] strong; a = sum_{k>=1} s^k  (ordered sequences);  s[n] = [n==1] + sum_primes polya(Aut, a)
a=[0]*(NMAX+1); s=[0]*(NMAX+1); a[1]=1; s[1]=1
def seqpow(s,n):
    # number of ordered sequences of strong classes with total size n (>=1 part)
    f=[Fraction(0)]*(n+1); f[0]=Fraction(1)
    for m in range(1,n+1): f[m]=sum(Fraction(s[j])*f[m-j] for j in range(1,m+1))
    return f[n]
for n in range(2,NMAX+1):
    sn=sum(polya(primes[p][0],a,n) for p in primes if p[0]<=n)
    s[n]=sn; a[n]=seqpow(s,n)
    print(f'n={n}: strong={s[n]} A000568={a[n]}')
# C3-structure (intransitivity hypergraph) iso classes: c[n]; strong C3 classes cs[n]
# c = MSET_{>=1}(cs): Euler transform; cs[n] = [n==1] + sum_{primes up to duality} polya(Aut+-, c)
# primes up to duality: a non-self-converse prime and its converse give the same C3 class -> count each pair once.
# tourn_enum lists each prime key once; for non-SC primes the converse is a different key also listed. Pair them by using
# the ANTI data: SC=1 -> one class with Aut+-; SC=0 -> group with converse (same Aut cycle types; Aut+- = Aut); divide by 2.
c=[0]*(NMAX+1); cs=[0]*(NMAX+1); c[1]=1; cs[1]=1
def euler_mset(cs,n):
    # number of multisets of strong C3 classes with total size n (>=1 element): coefficients of prod (1-x^j)^{-cs_j} - 1
    f=[Fraction(0)]*(n+1); f[0]=Fraction(1)
    for j in range(1,n+1):
        for _ in range(int(cs[j])):
            for m in range(j,n+1): f[m]+=f[m-j]
    return f[n]
for n in range(2,NMAX+1):
    tot=Fraction(0)
    for (k,key),(A,Apm) in primes.items():
        if k>n: continue
        w=Fraction(1) if len(Apm)>len(A) else Fraction(1,2)   # non-self-converse primes come in converse pairs
        tot+=w*polya(Apm,c,n)
    cs[n]=tot; c[n]=euler_mset(cs,n)
    print(f'n={n}: strong C3-classes={cs[n]} C3-classes={c[n]}')
