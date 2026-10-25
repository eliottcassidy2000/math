"""Exact labelled endpoint-four checks; no producer imports or source filters.

Exhausts the I={2,n-3,n-2}, c->b, c->x slice at n=7..11.
Only Case A is claimed below n11; smaller Case B outcomes are diagnostics.
Random controls are deterministic and explicitly bounded, not an all-order proof.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from random import Random
import json

def need(ok,why):
    if not ok: raise RuntimeError(why)

def edge(u,v):return (min(u,v),max(u,v))
def edges(q):return {edge(u,v) for u,v in zip(q,q[1:])}
def directed(T,q):return all((u,v) in T for u,v in zip(q,q[1:]))
def blue(T,q):
    r=dict(zip(q,q[::-1]))
    return all((r[v],r[u]) in T for u,v in T)
def move(T,q):
    E=edges(q)
    return {(u,v) if edge(u,v) in E else (v,u) for u,v in T}
def as_rows(T,n):
    out=[0]*n
    for u,v in T:out[u]|=1<<v
    return out
def row_blue(rows,q):
    n=len(rows);r=[0]*n
    for i,v in enumerate(q):r[v]=q[n-1-i]
    return all(bool(rows[u]>>v&1)==bool(rows[r[v]]>>r[u]&1)
               for u in range(n) for v in range(n) if u!=v)
def row_move(rows,q):
    n=len(rows);out=[((1<<n)-1)^row^(1<<u) for u,row in enumerate(rows)]
    for u,v in zip(q,q[1:]):out[u]^=1<<v;out[v]^=1<<u
    return out

def orbits(n):
    seen=set();fixed=[];free=[]
    for u,v in combinations(range(1,n-1),2):
        if (u,v) in seen:continue
        orb={(u,v),(n-1-v,n-1-u)};seen|=orb
        (fixed if v==u+1 else free).append(orb)
    return fixed,free
def assemble(n,h,fixed,free,bits):
    T=set().union(*fixed)
    for orb,bit in zip(free,bits):T.update(orb if bit else {(v,u) for u,v in orb})
    I={h,n-3,n-2};J={n-1-v for v in I}
    for v in range(1,n-1):
        T.add((v,0) if v in I else (0,v))
        T.add((n-1,v) if v in J else (v,n-1))
    T.add((n-1,0));return T
def random_source(n,h,rng,guards):
    fixed,free=orbits(n);bits=[]
    for orb in free:
        f=bool(orb&guards);b=bool({(v,u) for u,v in orb}&guards)
        need(not(f and b),'consistent actual guards')
        bits.append(1 if f else 0 if b else rng.getrandbits(1))
    T=assemble(n,h,fixed,free,bits);need(guards<=T,'all source guards');return T

def ledger_a(n):
    D=tuple(range(3,n-3));a,b,x,y,c,z=0,1,2,n-3,n-2,n-1
    return ((z,a,b,x)+D+(y,c),(z,a,x)+D+(b,y,c),(a,x)+D+(b,y,c,z))
def ledger_b(n,h=2):
    a,b,y,c,z=0,1,n-3,n-2,n-1;L=tuple(range(2,h+1));d=h+1;U=tuple(range(h+2,n-3))
    return ((b,)+L+(a,d)+U+(y,c,z),(b,)+L+(a,d,c,z)+U+(y,),
            (a,b)+L+(d,c,z)+U+(y,))

def check_source(T,n,h):
    need(len(T)==n*(n-1)//2 and all((v,u) not in T for u,v in T),'tournament')
    P=tuple(range(n));need(directed(T,P) and blue(T,P),'closed blue source')
    need({u for u,v in T if v==0}=={n-1,h,n-3,n-2},'exact incoming ports')
    need({v for u,v in T if u==n-1}=={0,1,2,n-1-h},'exact outgoing ports')
    need((n-2,1) in T and (n-2,2) in T,'closed middle and c->x')

def witness(T,n,h,which,qs):
    a,b,x,d,y,c,z=0,1,2,h+1,n-3,n-2,n-1
    S=move(T,qs[0]);B=move(T,tuple(range(n)))
    def mismatch(A,q,pair,other):
        r=dict(zip(q,q[::-1]));u,v=pair
        need((r[v],r[u])==other,'actual reflected ordered pair')
        need((pair in A)!=(other in A),'literal blackness witness')
    if which=='A':
        mismatch(T,qs[0],(c,a),(y,z))
        mismatch(S,qs[1],(z,x),(b,c))
        mismatch(B,qs[2],(a,y),(3,z))
    else:
        r=dict(zip(qs[0],qs[0][::-1]));v=n-h-2
        g=lambda u:r[n-1-u]
        need(g(a)==b and g(v)==v,'prefix reflection action')
        need((a,v) in T and (b,v) not in T,'prefix fixed-vertex obstruction')
        if h==2 and n==11:mismatch(S,qs[1],(b,z),(z,y))
        else:mismatch(S,qs[1],(a,n-h-6),(z,n-h-3))
        mismatch(B,qs[2],(a,n-4),(b,y))

def replay(T,n,h,which,strict=True):
    qs=ledger_a(n) if which=='A' else ledger_b(n,h)
    check_source(T,n,h);now=T;rows=as_rows(T,n);xor=set();ok=True
    for q in qs:
        need(len(q)==n and set(q)==set(range(n)),'full labelled word')
        if not directed(now,q):
            need(not strict,'directed actual word')
            return None
        flag=blue(now,q);need(flag==row_blue(rows,q),'independent anti predicate')
        ok &= not flag
        now=move(now,q);rows=row_move(rows,q);xor^=edges(q)
        need(rows==as_rows(now,n),'independent K evaluator')
    need(xor==edges(tuple(range(n))),'path-edge XOR')
    need(now==move(T,tuple(range(n))),'literal K partner')
    if strict:
        need(ok,'all three words black');witness(T,n,h,which,qs)
    return ok

def main():
    counts=Counter();digest=sha256();wrong={}
    for n in range(7,12):
        fixed,free=orbits(n)
        for bits in product((0,1),repeat=len(free)):
            T=assemble(n,2,fixed,free,bits)
            if (n-2,1) not in T or (n-2,2) not in T:continue
            which='A' if (3,n-2) in T else 'B'
            strict=which=='A' or n>=11
            good=replay(T,n,2,which,strict)
            label='black' if good else 'small_illegal_step' if good is None else 'small_blue_step'
            counts[f'exhaustive_n{n}_{which}_'+label]+=1
            digest.update(bytes([n,which=='A',2 if good is None else good]));digest.update(str(as_rows(T,n)).encode())
            if which not in wrong:
                qs=ledger_b(n) if which=='A' else ledger_a(n);now=T;failure=None
                for i,q in enumerate(qs):
                    if not directed(now,q):failure=i+1;break
                    now=move(now,q)
                need(failure is not None,'wrong complementary guard hostile')
                wrong[which]={'n':n,'rows':as_rows(T,n),'wrong_ledger_first_illegal_step':failure}
    rng=Random(202610250041)
    for n in range(11,81):
        for which in ('A','B'):
            guards={(n-2,1),(n-2,2),(3,n-2) if which=='A' else (n-2,3)}
            for _ in range(8):
                T=random_source(n,2,rng,guards);replay(T,n,2,which)
                counts['sampled_second_port_'+which]+=1
        for h in range(3,(n-8)//2+1):
            T=random_source(n,h,rng,{(n-2,1),(n-2,2),(n-2,h+1)})
            replay(T,n,h,'B');counts['sampled_prefix']+=1
            if n==2*h+8:counts['prefix_exact_threshold']+=1
    old=(6,4,8,51,39,7);T={(u,v) for u,row in enumerate(old) for v in range(6) if row>>v&1}
    B=move(T,tuple(range(6)))
    from itertools import permutations
    paths=[q for q in permutations(range(6)) if directed(B,q)]
    need(paths and all(blue(B,q) for q in paths),'pure-blue hostile remains')
    print(json.dumps({'status':'PASS; finite exact and deterministic bounded controls',
        'counts':dict(sorted(counts.items())),'wrong_guard_hostiles':wrong,
        'pure_target_order6_paths':len(paths),'exhaustive_digest':digest.hexdigest(),
        'filters':'only stated source arcs; no strongness, mixedness, isomorphism or clone filter',
        'independence':'local arc-set and bit-row K and anti evaluators; no producer imports'},sort_keys=True,indent=2))
if __name__=='__main__':main()
