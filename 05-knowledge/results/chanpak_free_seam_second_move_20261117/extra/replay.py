from itertools import product
from math import isqrt
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
root=Path(__file__).resolve().parents[2]
p=next(root/name/'base_flanks.py' for name in ('decimal_target','chanpak_decimal_core_selection_20261117') if (root/name).is_dir())
spec=spec_from_file_location('marked',p);a=module_from_spec(spec);spec.loader.exec_module(a)
from try_target import F,read,select

D=(54,32,2,119,6,20,55);tag=26
def require(x):
    if not x:raise ArithmeticError('control failed')
def rm(m):
    t=list(a.TT(m+2));t[0]^=1<<(m+1);t[-1]=1;return tuple(t)
def tr(w):
    s=bin(w)[2:] if w else ''
    return a.TT(len(s)),sum(int(b)<<j for j,b in enumerate(s))
def dual(X):
    x,t=X;return a.swap(x),((1<<len(x))-1)^t
def pre(w,X):
    x,t=X;y,s=tr(w);return a.join(y,x),s|(t<<len(y))
def sig(X):
    if not X[0]:return 1,1,[0,0],[0,0],0
    return a.local(*X)
def val(X,Y):
    h,p,L,R,d=sig(X);hh,pp,l,r,dd=sig(Y)
    return 91*h*hh+5*(d*hh+dd*h)+142*(sum(L)*hh+sum(l)*h)+256*sum(L)*sum(l)-2*(l[0]-l[1])*h
def append(z,X,Y):
    return a.make([tr(z),X,(D,tag),dual(Y)])
def all_literal(X,Y):
    return a.H(a.make([X,(D,tag),dual(Y)]))
def shifted(z,X,Y):
    s=a.law([tr(z),X,(D,tag),dual(Y)]) if z else val(X,Y)
    t=a.law([Y,(D,tag),dual(X),dual(tr(z))]) if z else val(Y,X)
    return s,t
local=sig((D,tag))
require((local[0],local[1],sum(local[2]),sum(local[3]),local[2][0]-local[2][1],local[3][0]-local[3][1],local[-1])==(5,256,142,142,-2,0,86))
for m in range(1,5):
    T=rm(m)
    expect=(1,)+tuple(1<<j for j in range(m))+(1,)
    require(tuple(a.htable(T)[-1])==expect)
    for z in range(1<<len(T)):
        h,P,L,R,d=sig((T,z));require(P==4*2**m+2*(3**m-1))
        require(L[1]-L[0]==sum(x for j,x in enumerate(expect) if z>>j&1))
for m,u,w in [(1,0,2),(1,1,2),(2,3,3),(2,0,3),(1,2,1),(1,3,5),(3,0,4)]:
    A=rm(m);X,Y=pre(u,(A,1 if w%2==0 else (1<<len(A))-2)),tr(w)
    p,q=all_literal(X,Y),all_literal(Y,X)
    require((p,q)==(val(X,Y),val(Y,X)))
    require(abs(p-q)==2 and min(p,q)==read(m,u,w))
    for z in (1,2,5):
        U,V=shifted(z,X,Y)
        require(U-p==V-q and U>p)
for m,u,w in [(1,0,0),(1,1,2),(2,0,1)]:
    X=pre(u,(rm(m),1));T=rm(m+1);Y=pre(w,(T,1|(1<<(len(T)-1))))
    p,q=all_literal(X,Y),all_literal(Y,X)
    require((p,q)==(val(X,Y),val(Y,X)) and p-q==2)
for m,u,w,p in [(1,0,21,42461),(2,0,10,26699),(3,2,25,1525031)]:
    require(read(m,u,w)==p)
    require(all(all(n%d for d in range(2,isqrt(n)+1)) for n in (p,p+2)))
for N in (8191,65537,458759):require(not select(N))
print('PASS Farey profiles, independent full DP controls, guard commutation and literal prime pairs')
