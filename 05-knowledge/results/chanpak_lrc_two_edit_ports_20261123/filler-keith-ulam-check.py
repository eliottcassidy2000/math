"""Exact finite inputs for the Keith/Ulam/Kaprekar filler note."""
from collections import Counter
from math import prod
from pathlib import Path
from runpy import run_path

bank = run_path(str(Path(__file__).with_name("third-slot-bank.py")))["data"]
entries = {f"{i}{tag}":(d,q,r[0])
           for i,(d,q,u,v) in enumerate(bank) for tag,r in zip("ab",(u,v))}
coeffs=(11,42,20)
P=(3,5,7,11,13)
S1=frozenset("0a 0b 1a 2a 2b 3a 3b 6a 6b 7a 7b 11a 11b".split())
S2=frozenset("0a 1a 2a 2b 3a 3b 6a 6b 7a 7b 11a 13a 13b".split())
families={frozenset(entries)}
for p in P:
    families={s & frozenset(name for name,t in entries.items() if t[j]%p)
              for s in families for j,c in enumerate(coeffs) if c%p}
assert max(map(len,families))==13
assert {s for s in families if len(s)==13}=={S1,S2}
for s,tail in ((S1,1),(S2,2)):
    for p,survivor in zip(P,(0,1,2,1,tail)):
        assert [j for j,c in enumerate(coeffs) if c%p and
                all(entries[name][j]%p for name in s)]==[survivor]

G20={name for name,(d,q,r) in entries.items() if d%3}
G18={name for name in G20 if entries[name][1]%5}
assert len(G20)==20 and G20-G18=={"5a","5b"}
assert all(entries[name][0]%5==0 for name in G20-G18)
assert all(q%5 and r%5 and q%4==0 for name,(d,q,r) in entries.items() if name in G18)
assert entries["0a"][0]%5==0 and pow(entries["0a"][1],4,25)==11

M=(16,16,5,5)
def seed(name):
    d,q,r=entries[name]
    return d%16,r%16,q%5,r%5
def close(gens):
    states={(1,1,1,1)}
    queue=list(states)
    for x in queue:
        for y in gens:
            z=tuple(a*b%m for a,b,m in zip(x,y,M))
            if z not in states:
                states.add(z)
                queue.append(z)
    return states
group=close(tuple(map(seed,("0a","1a","2a","11a"))))
assert len(group)==128
assert all(close(tuple(map(seed,S)))==group for S in (S1,S2,G18))
assert Counter(((11*d+20*r)%16,42*q%5) for d,r,q,_ in group)=={
    (x,y):4 for x in range(1,16,2) for y in range(1,5)}
for n,v in (("0a",(0,1)),("2a",(1,0))):
    d,_,r=entries[n]
    assert ((d*d-1)//8%2,(r*r-1)//8%2)==v

u=[1,2]
reps={}
for n in range(3,484):
    reps[n]=[(a,n-a) for a in u if 2*a<n and n-a in u]
    if len(reps[n])==1:
        u.append(n)
assert [n for n in u if n<=48]==[1,2,3,4,6,8,11,13,16,18,26,28,36,38,47,48]
assert reps[19]==[(1,18),(3,16),(6,13),(8,11)]
assert 28 in u and 47 in u and len(reps[61])==3 and 11 in u and 77 in u

def keith(n,B):
    a,b=divmod(n,B)
    assert 1<=a<B
    while b<n:
        a,b=b,a+b
    return b==n
assert [n for n in range(10,100) if keith(n,10)]==[14,19,28,47,61,75]
def K(n,B):
    a,b=divmod(n%(B*B),B)
    return (B-1)*abs(a-b)
assert [(K(n,10),K(7*n,10)) for n in [19,47,61]]==[(72,0),(27,63),(45,45)]
assert K(5,8)==35 and K(77,12)==11 and K(91,14)==13
assert keith(13,12) and keith(91,12)
assert Counter(K(n,10) for n in range(100) if n%2 and n%5)==dict(
    zip(range(0,90,9),(4,7,6,5,4,4,4,3,2,1)))
assert {69,341,483}.issubset(u) and 277 not in u
assert [reps[n] for n in (69,341,483)]==[[(16,53)],[(2,339)],[(69,414)]]
assert all(keith(n,61) for n in (69,483))
assert all(keith(n,26) for n in (55,385)) and K(55,26)==25 and K(385,26)==175

def H(rows):
    n=len(rows)
    full=(1<<n)-1
    d=[[0]*n for _ in range(full+1)]
    for j in range(n):
        d[1<<j][j]=1
    for seen in range(1,full+1):
        for j,num in enumerate(d[seen]):
            if num:
                z=rows[j]&~seen
                while z:
                    b=z&-z
                    d[seen|b][b.bit_length()-1]+=num
                    z^=b
    return sum(d[full])
def toggle(rows,arcs):
    rows=list(rows)
    for a,b in arcs:
        rows[a]^=1<<b
        rows[b]^=1<<a
    return rows

squares=[
 ((1022,508,504,496,992,448,384,768,0,366),((3,8),(4,8)),(69,341,277,483)),
 ((1022,1020,1016,752,224,960,640,512,216,272),((1,9),(2,9)),(55,283,167,385)),
 ((1022,1020,1016,496,992,960,384,512,128,328),((2,9),(3,8)),(25,63,59,175))]
for rows,arcs,expected in squares:
    assert all(((rows[i]>>j)&1)^((rows[j]>>i)&1)
               for j in range(len(rows)) for i in range(j))
    assert tuple(H(toggle(rows,aa)) for aa in ((),arcs[:1],arcs[1:],arcs))==expected

support_order="0a 0b 1a 2a 2b 3a 3b 6a 6b 7a 7b 11a 11b".split()
exponents=(9,12,7,2,8,3,1,10,12,8,8,9,6)
assert set(support_order)==S1 and sum(exponents)*9+10==865
D,Q,R=(prod(entries[n][j]**e for n,e in zip(support_order,exponents)) for j in range(3))
M=25*16**4
h=11*D+42*Q+20*R
assert (D%M,Q%M,R%M)==(1189625,1048576,830099)
assert h%M==47 and h%15015==12197
z=(h-47)//M
base=(27,82,41,124,62,31,94,47)
for _ in range(4):
    assert h%100==47
    h=h-47+K(h,10)
    for a in base[:-1]:
        assert h%2==a%2
        h=3*h+1 if h%2 else h//2
assert h==47+25*27**4*z

print("Exact support sieve, decimal lift seeds, and digit examples verified.")
