from collections import deque
from itertools import combinations

m = 11

def count(x):
    answer = 1
    while x:
        x, t = divmod(x, 2)
        if not t:
            answer += 2*((x+1)//2)
    return answer

def literal(word):
    n = m+1
    out = [sum(1<<j for j in range(i+1, m)) for i in range(m)] + [0]
    for i in range(m):
        if word >> (m-1-i) & 1:
            out[m] |= 1<<i
        else:
            out[i] |= 1<<m
    p = [[0]*n for _ in range(1<<n)]
    for i in range(n):
        p[1<<i][i] = 1
    for mask in range(1, 1<<n):
        for i, v in enumerate(p[mask]):
            if v:
                nxt = out[i] & ~mask
                while nxt:
                    j = nxt & -nxt
                    p[mask|j][j.bit_length()-1] += v
                    nxt ^= j
    return sum(p[-1])

route = [1158,1238,1114,1146]
if [literal(x) for x in route] != [721,723,725,727]:
    raise RuntimeError("literal HP mismatch")
if [x^y for x,y in zip(route,route[1:])] != [80,140,32]:
    raise RuntimeError("route mismatch")

H = [count(x) for x in range(1<<m)]
steps = [sum(1<<i for i in I) for k in (1,2) for I in combinations(range(m),k)]
seen = {route[0]}
todo = deque(seen)
while todo:
    x = todo.popleft()
    for z in steps:
        y = x ^ z
        if y not in seen and abs(H[y]-H[x]) <= 2:
            seen.add(y)
            todo.append(y)
if seen != {1158,1190,1270,1238}:
    raise RuntimeError(seen)
print("Four literal HP counts, 2/3/1-mark crossing, and exact two-mark component")
