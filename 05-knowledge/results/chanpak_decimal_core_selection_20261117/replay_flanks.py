from pathlib import Path
import importlib.util
source=Path(__file__).resolve().with_name('base_flanks.py')
spec=importlib.util.spec_from_file_location('apex_interval_previous',source)
ap=importlib.util.module_from_spec(spec);spec.loader.exec_module(ap)
from generic import F,select

def require(x):
    if not x:raise ArithmeticError('control failed')
def sig(a,q):
    h,P,L,R,d=ap.local(a,q)
    return h,P,sum(L),sum(R),L[0]-L[1],R[0]-R[1],h+d
def answer(u,w,S):
    h,P,L,R,ls,rs,k=S
    return k+h*(F(u)+F(w)-2)+R*u+L*w+P*u*w-rs*(u%2)-ls*(w%2)
def blocks(a,q,u,w):
    b=[]
    if u:
        s=list(map(int,bin(u)[2:]));b.append((ap.TT(len(s)),sum(x<<i for i,x in enumerate(s))))
    b.append((a,q))
    if w:
        s=[1-int(x) for x in reversed(bin(w)[2:])];b.append((ap.TT(len(s)),sum(x<<i for i,x in enumerate(s))))
    return b
D=(56,33,35,6,14,24);C=(4,1,2);E=(18,12,17,5,10);Q=(2,4,9,3)
require(sig(D,6)==(45,432,129,258,-15,-30,137))
require(sig(D,9)==(45,432,129,258,-15,-30,139))
for n in (257,631,8191,65537,458759,57337,4417,1799):
    require(not any(select(n,*sig(D,k)) for k in (6,9)))
fixtures=[(C,1,13,1,257),(C,1,11,4,631),(C,1,1219,0,8191),(C,1,65,83,65537),(C,1,0,94342,458759),(C,1,1,3370,57337),(C,1,230,0,1799),(Q,1,82,1,4417),(E,3,0,144,8191),(E,18,0,144,8189),(E,5,0,8064,458759),(E,24,0,8064,458761)]
for a,q,u,w,H in fixtures:
    S=sig(a,q)
    require(answer(u,w,S)==H)
    require(ap.law(blocks(a,q,u,w))==H)
    literal=ap.make(blocks(a,q,u,w))
    if len(literal)<=17:require(ap.H(literal)==H)
    print(len(literal),H,*literal)
# exhaustive small source formula vs direct DP for odd and even arguments
from itertools import product
for a,q in [(C,1),(E,3),(E,18),(E,5),(E,24),(D,6),(D,9)]:
    S=sig(a,q)
    for u,w in product(range(6),repeat=2):
        n=answer(u,w,S)
        require(ap.H(ap.make(blocks(a,q,u,w)))==n)
        z=select(n,*S);require(z is not None and answer(*z,S)==n)
print('PASS alternate-core and unrestricted parity checks')
