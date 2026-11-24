from itertools import combinations
from random import Random
R=Random(57)
def defects(A,E):
    n=len(A)
    t=sum((A[x]&A[y]).bit_count() for x,y in E)//3
    c=sum(u*(u-1)//2 for x,y in combinations(range(n),2) if (u:=(A[x]&A[y]).bit_count())>1)//2
    return 6*t+8*c,t,c
def area(E,A):
    for (a,b),(c,d) in combinations(E,2):
        if len({a,b,c,d})!=4:continue
        for x,y in ((c,d),(d,c)):
            if not(A[a]>>x&1) and not(A[b]>>y&1):
                yield (a,b),(c,d),tuple(sorted((a,x))),tuple(sorted((b,y)))
def delta(A,t):
    q=0
    for u,v in t:
        old=A[u]>>v&1
        a=A[u];b=A[v]
        f=0
        while a:
            p=a&-a;a^=p
            f+=(A[p.bit_length()-1]&b).bit_count()
        val=6*(A[u]&b).bit_count()+8*(f-(A[u].bit_count()+b.bit_count()-1 if old else 0))
        q+=-val if old else val
        A[u]^=1<<v;A[v]^=1<<u
    for u,v in t:A[u]^=1<<v;A[v]^=1<<u
    return q
def slide(A,E,t):
    for u,v in t:A[u]^=1<<v;A[v]^=1<<u
    E.difference_update(t[:2]);E.update(t[2:])

E=set([(0, 2), (0, 18), (0, 25), (0, 32), (0, 35), (0, 38), (0, 48), (1, 3), (1, 4), (1, 15), (1, 44), (1, 45), (1, 48), (1, 49), (2, 5), (2, 11), (2, 38), (2, 42), (2, 44), (2, 46), (3, 4), (3, 21), (3, 22), (3, 28), (3, 35), (3, 42), (4, 7), (4, 13), (4, 29), (4, 37), (4, 41), (5, 6), (5, 12), (5, 20), (5, 22), (5, 23), (5, 30), (6, 7), (6, 9), (6, 31), (6, 35), (6, 39), (6, 43), (7, 8), (7, 9), (7, 17), (7, 25), (7, 46), (8, 10), (8, 11), (8, 15), (8, 23), (8, 40), (8, 44), (9, 11), (9, 21), (9, 25), (9, 34), (9, 45), (10, 11), (10, 13), (10, 20), (10, 24), (10, 35), (10, 47), (11, 12), (11, 14), (11, 41), (12, 14), (12, 18), (12, 26), (12, 37), (12, 49), (13, 16), (13, 26), (13, 33), (13, 38), (13, 47), (14, 16), (14, 28), (14, 29), (14, 39), (14, 48), (15, 16), (15, 32), (15, 40), (15, 43), (15, 49), (16, 17), (16, 19), (16, 21), (16, 32), (17, 19), (17, 20), (17, 37), (17, 42), (17, 48), (18, 20), (18, 21), (18, 33), (18, 36), (18, 43), (19, 22), (19, 31), (19, 36), (19, 41), (19, 44), (20, 23), (20, 28), (20, 45), (21, 23), (21, 24), (21, 46), (22, 25), (22, 30), (22, 40), (22, 47), (23, 26), (23, 41), (23, 48), (24, 27), (24, 30), (24, 38), (24, 39), (24, 49), (25, 26), (25, 27), (25, 28), (26, 29), (26, 42), (26, 49), (27, 28), (27, 30), (27, 37), (27, 43), (27, 44), (28, 31), (28, 46), (29, 30), (29, 32), (29, 39), (29, 42), (30, 32), (30, 33), (31, 33), (31, 34), (31, 41), (31, 49), (32, 35), (32, 45), (33, 34), (33, 36), (33, 48), (34, 37), (34, 40), (34, 42), (34, 45), (35, 36), (35, 37), (36, 39), (36, 44), (36, 46), (37, 40), (38, 39), (38, 40), (38, 41), (39, 40), (41, 43), (42, 43), (43, 47), (44, 45), (45, 47), (46, 47), (46, 49), (47, 48)])

A=[0]*50
for u,v in E:A[u]|=1<<v;A[v]|=1<<u
assert all(a.bit_count()==7 for a in A)
assert defects(A,E)==(518,29,43)
q=[]
for t in area(E,A):
    B=A[:];F=E.copy();slide(B,F,t)
    DD=defects(B,F)[0]-518
    assert DD==delta(A,t)
    if DD==0:q.append(t)
    if DD<0:raise AssertionError('not local maximum')
assert sum(1 for _ in area(E,A))==22223
print('verified every local neighbor; neutral count',len(q),flush=True)
for t in q:
    B=A[:];F=E.copy();slide(B,F,t)
    for u in area(F,B):
        if (dd:=delta(B,u))<0:
            print('escape via neutral move then',dd, t,u,flush=True)
            quit()
print('no two-step plateau escape')
