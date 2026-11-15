"""Finite closure certificate for arbitrary two-sided padding around H=5."""
from collections import deque
from functools import cache
from itertools import permutations, product

C = (6, 12, 8, 1)
BITS = (0, 1)

def require(ok):
    if not ok:
        raise AssertionError("five-core closure control failed")

@cache
def paths(rows):
    return tuple(p for p in permutations(range(len(rows)))
                 if all(rows[u] >> v & 1 for u, v in zip(p, p[1:])))

def induce(R, J):
    return tuple(sum(((R[i] >> j) & 1) << k for k, j in enumerate(J))
                 for i in J)

def rev(rows):
    return tuple(sum((r >> i & 1) << j for j, r in enumerate(rows))
                 for i in range(len(rows)))

def apex(rows, mask):
    k = len(rows)
    return tuple(r | ((1 - (mask >> i & 1)) << k)
                 for i, r in enumerate(rows)) + (mask,)

def local(rows, u):
    return (len(paths(apex(rows, u))),
            sum(u >> p[-1] & 1 for p in paths(rows)))

def starts(R):
    ans = [0] * len(R)
    for v in range(1, 1 << len(R)):
        A = [i for i in range(len(R)) if v >> i & 1]
        B = [i for i in range(len(R)) if not(v >> i & 1)]
        for p in paths(induce(R,A)):
            ans[A[p[0]]] += len(paths(induce(R,B)))
    return tuple(ans)

F = [local(C,u) for u in range(16)]
B = [local(rev(C),15-u) for u in range(16)]
require(all(f[0] == b[0] for f,b in zip(F,B)))
require(len(paths(C)) == 5)
S,T = starts(C),starts(rev(C))
require((S,T) == ((7,9,6,5),(5,6,9,7)))
require(sum(S)+5 == sum(T)+5 == 32)
marks = {(sum(S[i] for i in range(4) if u >> i & 1),
          sum(T[i] for i in range(4) if u >> i & 1)):u
         for u in range(16)}
require(len({s for (s,r) in marks}) == len(marks) == 16)

G = {(s,r):[] for s,r in product(range(-5,33),repeat=2)}
for s,r in G:
    for z,w in product(BITS,repeat=2):
        sp, q = 2*s+5*z-32*w, r+32*z-5*w
        if q%2==0 and (sp,q//2) in G:
            G[s,r].append(((sp,q//2),z,w))
R = {u:[] for u in G}
for u in G:
    for v,z,w in G[u]:
        R[v].append((u,z,w))

def closure(A, graph):
    seen, q = set(A), deque(A)
    while q:
        u = q.popleft()
        for v,_,_ in graph[u]:
            if v not in seen:
                seen.add(v)
                q.append(v)
    return seen

live = closure(marks,G) & closure(marks,R)
out, visits = set(), 0
for b,M in marks.items():
    for v,z,w in G[b]:
        if v not in live:
            continue
        V0 = (1-z)*b[0] - w*(27-v[1]) + F[M][1]*(1-z)
        seen, q = set(), deque([(v,z,w,V0)])
        while q:
            v,pz,pw,V = key = q.popleft()
            if key in seen:
                continue
            seen.add(key)
            visits += 1
            require(-100 <= V <= 100)
            if v in marks:
                N = marks[v]
                dH = F[M][0] - F[N][0] + V - B[N][1]*pw
                ds = 5*pz - F[N][1]
                de = B[M][1] - 5*(1-w)
                for f,e in product(BITS,repeat=2):
                    out.add(dH + f*de + e*ds)
            for v2,z2,w2 in G[v]:
                if v2 in live:
                    V2 = V+(1-z2)*v[0]-w2*(27-v2[1])+5*pz*(1-z2)-5*pw*(1-w2)
                    q.append((v2,z2,w2,V2))
require(len(live) == 34 and out == {-2,0,2})
print("All live states:",len(live),"closure visits:",visits,
      "full-context count differences:",sorted(out))
