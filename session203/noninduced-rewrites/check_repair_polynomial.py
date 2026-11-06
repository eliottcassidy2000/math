import json,math
from itertools import permutations
import contextlib,io
with contextlib.redirect_stdout(io.StringIO()):
    from verify775 import graph,count
A=graph(4);n=len(A);b=3+(1<<7)+(1<<8)
C=[0]*(n+2)
for P in permutations(range(n)):
    allowed=int(bool(b>>P[0]&1))+int(not(b>>P[-1]&1));mandatory=0
    for x,y in zip(P,P[1:]):
        gap=not(b>>x&1) and bool(b>>y&1)
        allowed+=int(gap)
        if not(A[x]>>y&1):
            if not gap:break
            mandatory+=1
    else:
        for j in range(mandatory+1):C[allowed-j]+=(-1)**j*math.comb(mandatory,j)
if C!=[0,-7,10,26]+[0]*(n-2):raise RuntimeError(C)
for m in range(6):
    if count(graph(m+4))!=26*3**m+10*2**m-7:raise RuntimeError(m)
print(json.dumps({'core_rows':A,'tail_outgoing_mask':b,'exact_exponential_coefficients':C,'permutations_checked':math.factorial(n),'identity':'H(k)=26*3^(k-4)+10*2^(k-4)-7 for every k>=4','target_hit':'H(k)=6*(2^k+1)+1 iff k=7'},indent=2))
