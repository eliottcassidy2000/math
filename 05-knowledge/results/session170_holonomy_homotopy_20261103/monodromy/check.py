#!/usr/bin/env python3
"""Exact controls for signed-unit quotient rigidity and literal flat torus covers.
No external packages, assertions, or source imports.
"""
from functools import lru_cache
from itertools import permutations, product
from collections import deque
import json

def need(test, label):
    if not test:
        raise RuntimeError(label)

@lru_cache(None)
def cd(n,a,b):
    if n == 0:
        return 0,0
    half=1<<(n-1)
    A,u=divmod(a,half); B,v=divmod(b,half)
    if not A and not B:
        sign,index=cd(n-1,u,v)
    elif not A and B:
        sign,index=cd(n-1,v,u)
    elif A and not B:
        sign,index=cd(n-1,u,v); sign ^= bool(v)
    else:
        sign,index=cd(n-1,v,u); sign ^= 1 ^ bool(v)
    return int(sign),index+half*(A^B)

def left(n,a):
    d=1<<n
    out=[]
    for s,x in product(range(2),range(d)):
        t,y=cd(n,a,x)
        out.append((s^t)*d+y)
    return tuple(out)

def compose(a,b):
    return tuple(a[x] for x in b)

def inv(a):
    b=[0]*len(a)
    for x,y in enumerate(a): b[y]=x
    return tuple(b)

def closure_partition(size, generators, pairs, partition=None):
    root=list(range(size))
    def find(x):
        while root[x]!=x:
            root[x]=root[root[x]]; x=root[x]
        return x
    todo=deque(pairs)
    if partition is not None:
        todo.extend((x,partition[x]) for x in range(size))
    while todo:
        a,b=todo.popleft(); ra,rb=find(a),find(b)
        if ra==rb: continue
        if ra>rb: ra,rb=rb,ra
        root[rb]=ra
        todo.extend((g[a],g[b]) for g in generators)
    return tuple(find(x) for x in range(size))

def all_congruences(n):
    size=1<<(n+1); gens=[left(n,a) for a in range(1,1<<n)]
    initial=tuple(range(size)); found={initial}; todo=[initial]
    while todo:
        part=todo.pop()
        for v in range(1,size):
            if part[0]==part[v]: continue
            q=closure_partition(size,gens,[(0,v)],part)
            if q not in found: found.add(q); todo.append(q)
    for part in found:
        if part==initial: continue
        d=1<<n
        need(all(part[x]==part[x+d] for x in range(d)), 'proper quotient kept a sign')
        W={x for x in range(d) if part[x]==part[0]}
        need(all((x^y) in W for x in W for y in W),'index block not a subspace')
        need(all((part[x]==part[y]) == ((x^y) in W) for x in range(d) for y in range(d)), 'not a coset partition')
    return len(found)

def group(gens):
    identity=tuple(range(len(gens[0])))
    todo=[identity]; seen={identity}
    while todo:
        g=todo.pop()
        for h in gens:
            v=compose(h,g)
            if v not in seen: seen.add(v); todo.append(v)
    return seen

SEEDS=(153,8345,8337,8369,177,185)
P=(4,2,1,5,0,3)

def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]

def adjacency(n,mask):
    a=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs(n)):
        a[i][j]=bool((mask>>k)&1); a[j][i]=not a[i][j]
    return a

def relabel(n,mask,p):
    a=adjacency(n,mask); result=0
    pairindex={pair:k for k,pair in enumerate(pairs(n))}
    for i,j in pairs(n):
        u,v=p[i],p[j]
        if u<v:
            bit=a[i][j]
        else:
            u,v=v,u; bit=a[j][i]
        if bit: result |= 1<<pairindex[u,v]
    return result

def strong(a):
    n=len(a)
    for s in range(n):
        seen={s}; todo=[s]
        while todo:
            u=todo.pop()
            for v in range(n):
                if a[u][v] and v not in seen: seen.add(v); todo.append(v)
        if len(seen)!=n: return False
    return True

def profile_dp(a):
    n=len(a); size=1<<n
    dp=[[0]*n for _ in range(size)]
    for j in range(n): dp[1<<j][j]=1
    for mask in range(1,size):
        for j in range(n):
            if not ((mask>>j)&1): continue
            old=mask^(1<<j)
            if old:
                dp[mask][j]=sum(dp[old][i] for i in range(n) if a[i][j])
    h=[sum(row) for row in dp]
    covers=[[0]*(n+1) for _ in range(size)]; covers[0][0]=1
    for mask in range(1,size):
        bit=mask&-mask; sub=mask
        while sub:
            if sub&bit:
                rest=mask^sub
                for k in range(1,n+1): covers[mask][k]+=h[sub]*covers[rest][k-1]
            sub=(sub-1)&mask
    return covers[-1]

def profile_permcuts(a):
    n=len(a); ans=[0]*(n+1)
    for perm in permutations(range(n)):
        for cuts in range(1<<(n-1)):
            paths=[]; start=0; valid=True
            for k in range(n-1):
                if (cuts>>k)&1:
                    paths.append(tuple(perm[start:k+1])); start=k+1
                elif not a[perm[k]][perm[k+1]]:
                    valid=False; break
            if not valid: continue
            paths.append(tuple(perm[start:]))
            # Each unordered path cover has precisely one ordering by smallest vertex.
            minima=[min(path) for path in paths]
            if minima==sorted(minima): ans[len(paths)]+=1
    return ans

def ordinal(h):
    N=len(h); n=6*N; a=[[False]*n for _ in range(n)]
    for block,z in enumerate(h):
        b=adjacency(6,SEEDS[z])
        for i in range(6):
            for j in range(6): a[6*block+i][6*block+j]=b[i][j]
    for i in range(n):
        for j in range(i+1,n):
            if i//6 != j//6: a[i][j]=True
    mask=0
    for k,(i,j) in enumerate(pairs(n)):
        if a[i][j]: mask |= 1<<k
    return mask

def signed_coordinate(r):
    q=r[0]; y=sum(bit<<j for j,bit in enumerate(r[1:-1])); t=r[-1]
    return q^(t*bool(y)),y,t

def e_action(state,j,n):
    s,y,t=state
    if j==0: return s^1,y,t
    if j==n: return s^bool(y),y,t^1
    u=1<<(j-1)
    return s^(t*(bool(y)^bool(y^u))),y^u,t

def main():
    result={}
    # Independent subgroup enumeration for the exact finite sizes.
    sizes={}
    for n in (2,3,4):
        G=group([left(n,a) for a in range(1,1<<n)])
        sizes[str(n)]=len(G)
    need(sizes=={'2':8,'3':128,'4':4096}, 'group orders')
    result['transport_group_orders']=sizes
    counts={str(n):all_congruences(n) for n in (2,3,4)}
    need(counts=={'2':6,'3':17,'4':68}, 'equivariant quotient counts')
    result['all_equivariant_quotient_partitions_including_discrete']=counts
    stabilizer_checks=0
    for n in range(3,8):
        d=1<<n; low=d//2
        for x in range(1,d):
            y=x%low; t=x//low
            if t:
                f=lambda z: z//low
            else:
                z=next(k for k in range(1,low) if k!=y)
                f=lambda w,y=y,z=z: int(w%low in (y,z))
            need(f(0)==0 and f(x)==1,'stabilizer separator')
            stabilizer_checks+=1
    result['stabilizer_separators_n3_to_n7']=stabilizer_checks
    # Full monodromy spreads one associator defect to the whole sign layer.
    n=3; d=1<<n; gens=[left(n,a) for a in range(1,d)]
    defect=compose(inv(left(n,3)),compose(left(n,1),left(n,2)))
    need(sum(x==defect[x] for x in range(2*d))==8,'defect fixed sheets')
    need(defect==tuple(((x//d)^bool((x%d)&4))*d+x%d for x in range(2*d)), 'defect formula')
    quotient=closure_partition(2*d,gens,[(x,defect[x]) for x in range(2*d)])
    need(all((quotient[x]==quotient[y])==((x%d)==(y%d)) for x in range(2*d) for y in range(2*d)), 'relator quotient')
    result['associator_relator']={'fixed_sheets':8,'original_sheets':16,'maximal_extending_quotient_sheets':len(set(quotient))}
    # Literal seed controls independent of inherited verifiers.
    expected=[0,29,131,162,77,15,1]
    for h,mask in enumerate(SEEDS):
        a=adjacency(6,mask)
        need(strong(a),'strong seed')
        need(profile_dp(a)==expected,'DP profile')
        need(profile_permcuts(a)==expected,'permutation/cut profile')
        need((mask^SEEDS[(h+1)%6]).bit_count()==1,'hexagon edge')
        need(relabel(6,mask,P)==SEEDS[(h+3)%6],'physical deck relabelling')
    for a,b in product(range(6),repeat=2):
        need(((SEEDS[a]^SEEDS[b]).bit_count()==1)==((a-b)%6 in (1,5)), 'induced hexagon')
    canon=[min(relabel(6,mask,p) for p in permutations(range(6))) for mask in SEEDS]
    need(all((canon[i]==canon[j])==((i-j)%3==0) for i,j in product(range(6),repeat=2)), 'exact seed isomorphism classes')
    # Expand the falling factorials to verify the complete labelled-slot polynomial.
    slot=[0]*7; falling=[1]
    for k in range(1,7):
        nxt=[0]*(len(falling)+1)
        for degree,c in enumerate(falling):
            nxt[degree+1]+=c; nxt[degree]-=(k-1)*c
        falling=nxt
        for degree,c in enumerate(falling): slot[degree]+=expected[k]*c
    need(slot==[0,0,16,0,12,0,1], 'slot polynomial')
    result['seed_canonical_masks']=canon
    result['seed_slot_coefficients']=slot
    result['seed_profile']=expected
    result['seed_masks']=list(SEEDS)
    product_checks={}
    for N in range(1,4):
        states=list(product(range(6),repeat=N)); masks={h:ordinal(h) for h in states}
        for ai,a in enumerate(states):
            for b in states[ai+1:]:
                dif=[i for i in range(N) if a[i]!=b[i]]
                adjacent=len(dif)==1 and (a[dif[0]]-b[dif[0]])%6 in (1,5)
                need(((masks[a]^masks[b]).bit_count()==1)==adjacent,'product induced graph')
            logical=tuple(x%3 for x in a)
            projected=[]
            for i in range(N):
                for sign in (-1,1):
                    b=list(a); b[i]=(b[i]+sign)%6
                    projected.append(tuple(x%3 for x in b))
            need(len(set(projected))==2*N,'local cover bijection')
        product_checks[str(N)]={'physical_vertices':6**N,'logical_vertices':3**N,'sheets':2**N,'order':6*N}
    result['literal_product_controls']=product_checks
    action_checks=0
    for n in range(3,8):
        for r in product(range(2),repeat=n+1):
            s=signed_coordinate(r)
            for j in range(n+1):
                new=list(r); new[j]^=1
                need(signed_coordinate(new)==e_action(s,j,n),'E_n sheet chart')
                action_checks+=1
    result['E_n_coordinate_action_checks_n3_to_n7']=action_checks
    # Every E_3 generator is a three-literal-arc loop on actual order-24 tournaments.
    literal_edges=0
    for r in product(range(2),repeat=4):
        initial=tuple(3*z for z in r)
        for j in range(4):
            h=list(initial)
            for step in range(3):
                old=ordinal(h); h[j]=(h[j]+1)%6; new=ordinal(h)
                need((old^new).bit_count()==1,'E3 actual arc loop')
                literal_edges+=1
            need(all(h[k]%3==0 for k in range(4)),'logical return')
            decoded=signed_coordinate([h[k]//3 for k in range(4)])
            need(decoded==e_action(signed_coordinate(r),j,3),'physical E3 monodromy')
    result['literal_E3']={'tournament_order':24,'signed_sheets_per_logical_vertex':16,'checked_arc_flips':literal_edges,'Hamilton_count':29**4}
    # Flat square returns and triangular generator returns.
    for h in product(range(6),repeat=3):
        for i,j in ((0,1),(0,2),(1,2)):
            v=list(h)
            for k,s in ((i,1),(j,1),(i,-1),(j,-1)): v[k]=(v[k]+s)%6
            need(tuple(v)==h,'mixed square flatness')
        for i in range(3):
            v=list(h); v[i]=(v[i]+3)%6
            need(tuple(x%3 for x in v)==tuple(x%3 for x in h) and tuple(v)!=h,'triangle deck twist')
    result['mixed_square_checks']=6**3*3
    result['coordinate_triangle_checks']=6**3*3
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__': main()
