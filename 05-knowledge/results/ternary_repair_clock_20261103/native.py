"""Literal execution of the inherited one-temporary-vertex native compiler."""
from functools import lru_cache
from itertools import combinations


def need(ok,message):
    if not ok:raise RuntimeError(message)


def flip(T,a,b):
    U=list(T);U[a]^=1<<b;U[b]^=1<<a;return tuple(U)


def legal(T,p):
    return sorted(p)==list(range(len(T))) and all(T[a]>>b&1 for a,b in zip(p,p[1:]))


def native(T,p):
    need(legal(T,p),'native path guard')
    n=len(T);U=[((1<<n)-1)^r^(1<<v) for v,r in enumerate(T)]
    for a,b in zip(p,p[1:]):U[a]|=1<<b;U[b]&=~(1<<a)
    U=tuple(U);need(legal(U,p),'native retained path');return U


def extend(T,out):
    n=len(T)
    return tuple(r|(0 if out>>j&1 else 1<<n) for j,r in enumerate(T))+(out,)


def arc_path(T,u,v):
    """Contract the required directed arc into a two-ended unit; exact subset DFS."""
    need(T[u]>>v&1,'arc orientation')
    units=tuple((j,) for j in range(len(T)) if j not in (u,v))+((u,v),)
    full=(1<<len(units))-1
    @lru_cache(None)
    def visit(mask,last):
        if mask==full:return ()
        for j,w in enumerate(units):
            if mask>>j&1 or (last>=0 and not T[units[last][-1]]>>w[0]&1):continue
            tail=visit(mask|1<<j,j)
            if tail is not None:return w+tail
        return None
    return visit(0,-1)


def traced(T,u,v,p):
    n=len(T);k=p.index(u)
    if n==2:cut=0;p1=p+(n,);q=(v,n,u)
    elif v!=p[-1]:cut=1<<u;p1=p+(n,);q=p[:k+1]+(n,)+p[k+1:]
    else:cut=((1<<n)-1)^(1<<v);p1=(n,)+p;q=p[:k+1]+(n,)+p[k+1:]
    A=extend(T,cut);B=native(A,p1);C=native(B,q)
    need(tuple(r&((1<<n)-1) for r in C[:n])==flip(T,u,v),'traced old relation')
    return A,B,C,p1,q


def arc_command(T,a,b):
    u,v=(a,b) if T[a]>>b&1 else(b,a);p=arc_path(T,u,v);reverse=p is None
    if reverse:
        S=flip(T,u,v);p=arc_path(S,v,u);need(p is not None,'one-guard separation')
        C,B,A,q,p1=traced(S,v,u,p)
    else:A,B,C,p1,q=traced(T,u,v,p)
    need(tuple(r&((1<<len(T))-1) for r in A[:len(T)])==T,'birth retains source')
    need(native(A,p1)==B and native(B,q)==C,'both legal operations')
    out=tuple(r&((1<<len(T))-1) for r in C[:len(T)])
    need(out==flip(T,u,v),'exact arc command')
    return out,{'arc':[a,b],'birth_out':A[-1],'paths':[p1,q],'inverse_case':reverse}


def compile_repair(source,target):
    """Pad with sinks, execute named arc commands, then delete named surplus."""
    source=tuple(source);target=tuple(target);n=max(len(source),len(target))
    T=source;U=target
    while len(T)<n:T=extend(T,0)
    while len(U)<n:U=extend(U,0)
    commands=[]
    for a,b in combinations(range(n),2):
        if ((T[a]>>b)^(U[a]>>b))&1:T,c=arc_command(T,a,b);commands.append(c)
    need(T==U,'complete target before named deletion')
    final=tuple(r&((1<<len(target))-1) for r in T[:len(target)])
    need(final==target,'complete literal target')
    return {'source_rows':source,'target_rows':target,'padding_births':n-len(source),
            'terminal_deletions':list(range(len(target),n)), 'commands':commands,
            'temporary_births':len(commands),'native_moves':2*len(commands),
            'temporary_deletions':len(commands),'peak_order':n+1 if commands else n}
