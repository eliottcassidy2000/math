from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from counter36 import H
from itertools import product
from functools import lru_cache
H=lru_cache(None)(H)
zero=(1,0,0); two=(0,1,0);three=(0,0,1)
def block(p,q,r):
    A,B=2**(p+r),3**(p+r)
    Q,R=2**q,3**q
    S=3**p+3**r; U=2**p*3**r+3**p*2**r
    return 3-2*Q+2*Q*S+(4+8*Q+2*R)*A+(4*Q+2*R)*U+2*R*B
assert all(H('21'+'0'*p+'2'*q+'3'*r+'12')==block(p,q,r)
           for p,q,r in product((1,2,3),repeat=3))
def matched(a,b):
    return (a&1 and b&2) or (a&2 and b&1)
ends=[''.join(map(str,p)) for p in product(range(4),repeat=2) if matched(*p)]
rends=[''.join(str(3-int(d)) for d in e) for e in ends]
C=[(lo+b+hi) for lo in ends for hi in rends for b in map(''.join,product('0123',repeat=3))]
assert min(map(H,C))==281
def strong(w):
    return w[0]!='0' and w[-1]!='3' and ('1' in w or ('3' in w and '0' in w[w.index('3')+1:]))
for N in (385,573):
    todo=[(w,(0,1,5,6)) for w in C if H(w)<=N]; seen=set(todo)
    print('target',N,'viable anchored roots',len(todo))
    hit=set()
    for w,P in todo:
        if H(w)==N and strong(w):hit.add(w)
        if len(w)>=N.bit_length()-1+3:continue
        a,b,c,d=P
        # Inserting between b and c would add a triple site; both targets are below the next floor.
        for j in list(range(0,b+1))+list(range(c+1,len(w)+1)):
            Q=tuple(k+(k>=j) for k in P)
            for digit in '0123':
                W=w[:j]+digit+w[j:]
                if (W,Q) not in seen and H(W)<=N:
                    seen.add((W,Q));todo.append((W,Q))
    assert (not bool(hit)) == (N==385)
    if N==573:assert '20302212' in hit
    print('anchored states',len(seen),'strong target hits',len(hit))
