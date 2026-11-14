from math import gcd
from collections import Counter
from json import loads
from pathlib import Path
from profiles_near_30030 import profile, strong
from inverse_carries import logical, pair
D={int(c):r for c,r in loads(Path(__file__).with_name('phase_table.json').read_text()).items() if gcd(int(c),99)==1}
logical([Counter(D.values())[i] for i in range(5)]==[20,18,10,9,3])
for n in (1,30013,30017,30029,30047):
    logical(gcd(n,198)==1)
    first=0; both=0
    by={c%198:c for c in D}
    for h in range(19008):
        if gcd(h,198)!=1: continue
        c=by[n*pow(h,-1,198)%198]
        z=n-c*h
        m=99*2**(D[c]+2)
        if z%m==0:
            first+=1;both+=int(z%(3*m)!=0)
    logical((first,both)==(1575,1050))
A=(2,12,9,17,39,15)
logical(strong(A) and profile(A)==(23,272))
logical(pair((1,1,0,0,0))==(23,296))
logical(99*296+31*23==30017 and (30017-99*272)%23!=0 and 2*99*272>30017)
print('actual phase gates: 1050 classes among 5760 units; fixed-H non-greedy control passed')
