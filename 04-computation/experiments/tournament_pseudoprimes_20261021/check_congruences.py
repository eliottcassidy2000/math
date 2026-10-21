# check_congruences.py -- THM-4747 checks on circulant data (circ_counts outputs).  mac-mini-2026-10-21-fanopsp
import re, glob
from math import gcd
def isprime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
rows=[]; seen=set()
for f in ['circ_small.out','circ_13_17.out','circ_19.out','circ_21.out','circ_23.out','paley19.out','paley23.out']:
    for l in open(f):
        m=re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+) H/n=(\d+) C=(\d+) rho0=(\d+) A=(\d+)',l)
        if m:
            n,S,mu,H,h0,C,rho,A=m.groups()
            nn=int(n); Sl=[int(x) for x in S.split(',')]
            key=(nn,min(tuple(sorted((a*x)%nn for x in Sl)) for a in range(1,nn) if gcd(a,nn)==1))
            if key in seen: continue
            seen.add(key); rows.append((int(n),[int(x) for x in S.split(',')],int(mu),int(H),int(C),int(rho),int(A)))
bad=0; checked=0; comp={}
for n,S,mu,H,C,rho,A in rows:
    checked+=1
    if H%n: bad+=1; print('H not divisible by n',n,S)
    if rho%2==0: bad+=1; print('rho0 even',n,S)
    if isprime(n):
        if (C-(n-1)//2)%(n*mu): bad+=1; print('C congruence fails',n,S,C,mu)
        if (rho-mu)%(2*mu): bad+=1; print('rho congruence fails',n,S,rho,mu)
        if H%(n*mu): bad+=1; print('H not divisible by n*mu',n,S)
    else:
        comp.setdefault(n,[0,0]); comp[n][1]+=1
        if (C-A)%n==0: comp[n][0]+=1
print('distinct circulant classes checked',checked,'violations',bad)
print('composite n: classes passing naive C = A (mod n):',{n:f'{a}/{b}' for n,(a,b) in sorted(comp.items())})
