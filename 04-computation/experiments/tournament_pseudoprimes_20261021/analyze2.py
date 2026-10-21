import gzip, os
def _open_census(n):
    p=f'census_n{n}.txt'
    return open(p) if os.path.exists(p) else gzip.open(p+'.gz','rt')
import math
from collections import Counter, defaultdict
def load(n): return [tuple(map(int,l.split())) for l in _open_census(n)]
first=defaultdict(lambda:None); firststrong=defaultdict(lambda:None)
nonstrongvals={}
for n in range(3,10):
    for H,C,s,A,AA,c3 in load(n):
        if first[H] is None: first[H]=n
        if s and firststrong[H] is None: firststrong[H]=n
# values attained (n<=9) only via strong tournaments at their first order?
R9=load(9)
strong9=Counter(r[0] for r in R9 if r[2]); all9=Counter(r[0] for r in R9)
nonstrong9=Counter(r[0] for r in R9 if not r[2])
print('first appearance order m(h) for odd h<=200 (None = not by n=9):')
row=[]
for h in range(1,201,2):
    row.append(f'{h}:{first[h]}/{firststrong[h]}')
print(' '.join(row))
print()
print('values <=700 present at n=9 but NEVER by a non-strong 9-tournament (strong-only values):')
so=[h for h in sorted(all9) if h<=700 and nonstrong9[h]==0]
print(so)
comp=[h for h in so if not all(h%p for p in range(3,int(h**.5)+1,2)) and h>1]
print('composite strong-only (H-pseudoprime candidates):', comp)
def fac(h):
    f=[];d=3;x=h
    while d*d<=x:
        while x%d==0: f.append(d); x//=d
        d+=2
    if x>1: f.append(x)
    return f
print([ (h,fac(h)) for h in comp])
# divisibility frequencies at n=9 (classes and labeled)
fact=math.factorial(9)
for p in (3,5,7,11,13):
    cl=sum(1 for r in R9 if r[0]%p==0)/len(R9)
    lab=sum(fact//r[3] for r in R9 if r[0]%p==0)/2**36
    print(f'p={p}: frac classes with p|H {cl:.4f}  labeled {lab:.4f}  (1/p={1/p:.4f})')
