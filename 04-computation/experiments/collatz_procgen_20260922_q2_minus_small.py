# minus sheet E_- (n -> 3n-1, n -> n/2 for even n): BFS from 1, report reachability of all m<=5000, 3 !| m (values capped)
from collections import deque
CAP=10**7; seen=bytearray(CAP+1); q=deque([1]); seen[1]=1
while q:
    x=q.popleft()
    for y in (3*x-1, x//2 if x%2==0 else 0):
        if 1<=y<=CAP and not seen[y]: seen[y]=1; q.append(y)
bad=[m for m in range(1,5001) if m%3 and not seen[m]]
print("minus sheet: m<=5000, 3!|m, not reached from 1 (cap 1e7):", bad)
print("m=2 reached:",bool(seen[2]),"m=4 reached:",bool(seen[4]))
