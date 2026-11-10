# Exact O/Z/J orbit in nonnegative coordinates (h, E, Op, M): Op = O - M.
#   J: (h+2E, M+Op, E, h)   Z: (h+2E, E+Op, E+Op, M)   O: (h, E+Op+M, E+Op, M)
# start (1,0,0,1) = single vertex.  value = h.  Next J/Z value = h + 2E.
import sys
from collections import defaultdict
N = int(sys.argv[1])
start = (1,0,0,1)
seen = {start: 0}   # state -> min number of J's
frontier = [start]
minJ = defaultdict(lambda: 99)   # value -> min number of J letters needed
cnt = defaultdict(int)
minJ[1] = 0
while frontier:
    nf = []
    for st in frontier:
        h,E,Op,M = st; j = seen[st]
        if h + 2*E <= N:
            for ns, dj in (((h+2*E, M+Op, E, h), 1), ((h+2*E, E+Op, E+Op, M), 0)):
                if ns not in seen or seen[ns] > j+dj:
                    seen[ns] = j+dj; nf.append(ns)
                if minJ[ns[0]] > j+dj: minJ[ns[0]] = j+dj
        ns = (h, E+Op+M, E+Op, M)
        if ns[0] + 2*ns[1] <= N and (ns not in seen or seen[ns] > j):
            seen[ns] = j; nf.append(ns)
    frontier = nf
vals = sorted(v for v in minJ if v <= N)
miss = [n for n in range(1, N+1, 2) if minJ[n] == 99]
print("states:", len(seen), " values<=N:", len(vals), " missing odd:", miss[:60], "count", len(miss))
from collections import Counter
c = Counter(minJ[n] for n in range(1, N+1, 2))
print("min #J distribution:", sorted(c.items()))
big = [n for n in range(1, N+1, 2) if minJ[n] >= 3 and minJ[n] < 99]
print("need >=3 J:", big[:80], len(big))
