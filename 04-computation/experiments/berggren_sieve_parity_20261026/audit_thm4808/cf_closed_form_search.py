# Audit search: XOR of simple continued-fraction / parity features that equals the first letter or chi (m < 400).
# Finding: chi odd == (m odd) XOR (sum of RCF partial quotients odd); proof in REPORT.md; checked to m <= 1500 in tree_checks.py.
from math import gcd
def word(m,n):
    w=[]
    while (m,n)!=(2,1):
        if m>3*n: w.append('D'); m-=2*n
        elif m>2*n: w.append('A'); m,n=n,m-2*n
        else: w.append('U'); m,n=n,2*n-m
    return ''.join(reversed(w))
def cf(m,n):
    d=[]
    while n: d.append(m//n); m,n=n,m%n
    return d
rows=[]
for m in range(3,400):
    for n in range(1,m):
        if (m-n)%2==0 or gcd(m,n)!=1: continue
        w=word(m,n); d=cf(m,n)
        rows.append((m,n,w,d))
feats={
 'cL odd': lambda m,n,d: d[-1]%2==1,
 'c0 odd': lambda m,n,d: d[0]%2==1,
 'sum odd': lambda m,n,d: sum(d)%2==1,
 'L odd': lambda m,n,d: (len(d)-1)%2==1,
 'm odd': lambda m,n,d: m%2==1,
 'sum even-index digits odd': lambda m,n,d: sum(d[0::2])%2==1,
 'sum odd-index digits odd': lambda m,n,d: sum(d[1::2])%2==1,
 '#odd digits odd': lambda m,n,d: sum(1 for x in d if x%2)%2==1,
}
import itertools
tgt=[w[0]=='U' for m,n,w,d in rows]
chi=[w.count('A')%2==1 for m,n,w,d in rows]
names=list(feats)
vals={k:[feats[k](m,n,d) for m,n,w,d in rows] for k in names}
for r in (1,2,3):
    for combo in itertools.combinations(names,r):
        x=[False]*len(rows)
        for k in combo: x=[a^b for a,b in zip(x,vals[k])]
        for t,name in ((tgt,'firstU'),(chi,'chi odd')):
            if x==t or [not a for a in x]==t: print(name,'== XOR of',combo, '(complemented)' if x!=t else '')
print(len(rows))
