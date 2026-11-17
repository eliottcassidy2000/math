from pathlib import Path
from collections import Counter
from itertools import permutations,combinations
P=Path(__file__).with_name('witnesses.tsv')
def ck(v):
    if not v:raise ArithmeticError('control failed')
def paths(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s,row in enumerate(D):
        for i,x in enumerate(row):
            if not x:continue
            for j in range(n):
                if (a[i]>>j&1) and not (s>>j&1):D[s|1<<j][j]+=x
    return sum(D[-1])
def flip(a,i,j):
    v=list(a);v[i]^=1<<j;v[j]^=1<<i;return tuple(v)
def de(a,k):
    return tuple((r&((1<<k)-1))|((r>>(k+1))<<k) for i,r in enumerate(a) if i!=k)
def valid(a):
    n=len(a)
    ck(all(0<=r<2**n and not(r>>i&1) for i,r in enumerate(a)))
    ck(all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i,j in combinations(range(n),2)))
V={};O=Counter()
for line in P.read_text().splitlines()[1:]:
    n,h,z=line.split('\t');n=int(n);h=int(h);a=tuple(map(int,z.split(',')))
    ck(len(a)==n);valid(a);ck(paths(a)==h and h not in V);V[h]=n;O[n]+=1
ck(O=={1:1,3:1,4:1,5:4,6:12,7:58,8:243,9:1200})
for a,ed,v in [((8,9,3,148,7,31,63,119),((0,7),(0,6)),1),((8,1,3,198,15,31,55,119),((0,7),(1,7)),6)]:
    ck(paths(a)==63 and V[63]==8)
    ck(set(paths(de(a,j)) for j in range(8))=={1,11,25,33,37})
    ck(all(V[paths(de(a,j))]<7 for j in range(8)))
    for e,h,g in zip(ed,(163,211),(49,59)):
        b=flip(a,*e);ck(paths(b)==h and V[h]==8 and v not in e)
        ck(paths(de(b,v))==g and V[g]==7)
ca=(2,24,3,37,13,23);cb=(4,9,10,49,7,23)
for a in (ca,cb):
    ck(paths(a)==25 and V[25]==V[27]==6)
    ck(sorted(paths(flip(a,i,j)) for i,j in combinations(range(6),2) if a[i].bit_count()==a[j].bit_count() and a[i].bit_count() in (2,3))==[19,23])
print('PASS all 1520 literal birth witnesses, relay squares, smallest missing unit link')
