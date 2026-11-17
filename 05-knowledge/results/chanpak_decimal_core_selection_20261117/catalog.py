from generic import select
from replay_flanks import sig,answer,D,C,E,Q
def choose(N,unit=False):
    for a in (D,C,Q,E):
        fs=[(q,sig(a,q)) for q in range(1,(1<<len(a))-1)]
        for q,s in fs:
            other=[(r,t[-1]-s[-1]) for r,t in fs if t[:-1]==s[:-1] and abs(t[-1]-s[-1])==2]
            if unit and not other:continue
            u=select(N,*s)
            if u is None:continue
            if answer(*u,s)!=N:raise ValueError(N)
            return a,q,u,other
    return None
for n,unit in ((10007,False),(10009,False),(30011,False),(30013,False),(65537,False),
               (458759,True),(10**5+3,False),(10**6+3,False),(10**6+7,False),(10**9+7,True)):
    t=choose(n,unit)
    if t is None:raise ValueError((n,unit))
    print(n,t)
