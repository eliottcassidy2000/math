from profiles_near_30030 import profile, canon, extend, strong
from inverse_carries import logical
from math import isqrt, gcd
from pathlib import Path
from json import loads
D=loads(Path(__file__).with_name('phase_table.json').read_text())
for a,b,q,A,B,pad in (
    (5,19,256,(2,4,9,3),(2,20,9,3,45,15),3),
    (79,91,1184,(6,12,24,97,11,23,55),(8,5,17,102,43,71,23),0),
):
    logical(strong(A) and strong(B) and profile(A)==(a,q>>pad) and profile(B)==(b,q))
    logical(gcd(a,b)==1 and gcd(99*q,a*b)==1)
logical([D[str(c)] for c in (19,5,91,79)]==[0,0,0,2])
for n in (25344+95,117216*64+7189):
    logical(all(n%d for d in range(2,isqrt(n)+1)))
T={(0,)};hits=[]
for n in range(2,7):
    T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
    for A in T:
        if strong(A) and profile(A)[1] in (100,200):
            hits.append(A)
logical(hits==[(2,20,9,3,13)])
A=hits[0];logical(canon(tuple(31^r^(1<<i) for i,r in enumerate(A)))==A)
print('two exact matched-Q target switches; both prime controls; round-clock class certified')
