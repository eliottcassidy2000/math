"""Recover slope-matched carrier windows without enumerating all suffix words."""
from collections import defaultdict
from functools import lru_cache
from sys import argv
from scan_carrier_windows import rows, signature

@lru_cache(None)
def suf(d,r):
    if d==0:
        return (1,0,0,0,0)
    D,R,A,B,C=suf(d-1,r//2)
    if r%2:
        return (2*D,r,A,B,C)
    return (2*D,r,A+D,B+(d==1),C+R+(R%2))

def mirrored(a,q,c,d,j):
    """c -> d, transferring j low bits from Y to X."""
    _,_,s,e,t,f,n=c
    _,_,S,E,T,F,N=d
    L=(T+a)*(1<<j)-(t+a)
    R=(s+a)*(1<<j)-(S+a)
    if not(0<=L and 0<=R):
        return ()
    maxcarry=q+a-1
    if max(L>>j,R>>j)>maxcarry:
        return ()
    @lru_cache(None)
    def paths(i,r,l):
        if i==j:
            return ((0,0),) if r==R>>j and l==0 else ()
        k=j-i-1
        ans=[]
        for u in (0,1):
            for v in (0,1):
                r0=r+q*u+a*v-(R>>i&1)
                l0=(L>>k&1)+2*l-q*v-a*u
                if (r0%2==0 and 0<=r0//2<=maxcarry and 0<=l0<=maxcarry):
                    for U,V in paths(i+1,r0//2,l0):
                        ans.append((U+(u<<i),V+(v<<k)))
        return tuple(ans)
    return paths(0,0,L>>j)

def main(K):
    data={}
    for (a,q),group in rows().items():
        expanded={}
        for (a,q,s,e,t,f,n),loc in group:
            expanded.setdefault((a,q,s,e,t,f,n),(loc,0))
            expanded.setdefault((a,q,t,f,s,e,n),(loc,1))
        data[a,q]=expanded
    best=0
    for j in range(K+1):
        total=0
        for (a,q),group in sorted(data.items()):
            seen=set()
            for c in group:
                for d in group:
                    for u,v in mirrored(a,q,c,d,j):
                        S=signature(c,suf(0,0),suf(j,v))[0]
                        T=signature(d,suf(j,u),suf(0,0))[0]
                        if S[:4]!=T[:4]:
                            if j==0:
                                continue
                            raise RuntimeError("slope")
                        key=(S[:4],u,v)
                        if key in seen:
                            continue
                        seen.add(key)
                        total+=1
                        vals=[{} for _ in range(4)]
                        for i in range(j+1):
                            x=suf(i,u>>(j-i)); y=suf(j-i,v>>i)
                            for g,desc in group.items():
                                sg,z=signature(g,x,y)
                                if sg[:4]!=S[:4]:
                                    continue
                                for parity in range(4):
                                    z0=z+sg[4]*(parity&1)+sg[5]*(parity>>1)
                                    vals[parity].setdefault(z0,(*desc,i,x[1],y[1]))
                        for parity,V in enumerate(vals):
                            for t in sorted(V):
                                if t-2 in V:
                                    continue
                                stop=t
                                while stop in V:
                                    stop+=2
                                L=(stop-t)//2
                                if L>best:
                                    best=L
                                    print('record span',j,'length',L,'signature',S[:4],'parity',parity)
                                    for m in range(t,stop,2):
                                        print(m,V[m])
        print('span',j,'matched extremal suffix pairs',total,'record',best,flush=True)
if __name__=='__main__':
    main(int(argv[1]) if len(argv)>1 else 16)
