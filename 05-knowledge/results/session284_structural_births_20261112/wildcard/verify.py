#!/usr/bin/env python3
"""Standalone exact ternary-host and marked-symmetry controls."""
import itertools, json, math, collections, functools

def check(ok, msg):
    if not ok: raise RuntimeError(msg)

def mul(a,b):
    out=collections.Counter()
    for (s,v),c in a.items():
        for (t,w),d in b.items(): out[s+t,v+w]+=c*d
    return out

def orbit_step(a):
    out=mul(mul(a,a),a)
    for (s,v),c in a.items():
        out[3*s,3*v]-=c
        out[3*s,3*v+1]+=3*c
    check(all(v%3==0 for v in out.values()), 'orbit integrality')
    return {k:v//3 for k,v in out.items() if v}

def canon_col(word):
    if len(word)==1: return word[0],0
    n=len(word)//3
    parts=[canon_col(word[n*i:n*(i+1)]) for i in range(3)]
    keys=tuple(x[0] for x in parts)
    rot=min(keys[i:]+keys[:i] for i in range(3))
    return rot,sum(x[1] for x in parts)+(keys[0]==keys[1]==keys[2])

X=('x',)
def seq(items):
    parts=[]
    for p in items:
        if p is None: continue
        parts.extend(p[1] if p[0]=='s' else (p,))
    if not parts: return None
    if len(parts)==1: return parts[0]
    return ('s',tuple(parts))

def cyc(items):
    keys=tuple(items)
    return ('c',min(keys[i:]+keys[:i] for i in range(3)))

def induced(mask,e):
    if not mask: return None
    if e==0: return X
    n=3**(e-1); z=(1<<n)-1
    kids=[induced((mask>>(i*n))&z,e-1) for i in range(3)]
    have=[i for i,x in enumerate(kids) if x is not None]
    if len(have)==3: return cyc(kids)
    if len(have)==1: return kids[have[0]]
    i,j=have
    if (j-i)%3!=1: i,j=j,i
    return seq([kids[i],kids[j]])

@functools.lru_cache(None)
def birth(t):
    if t[0]=='x': return 0
    if t[0]=='c': return 1+max(map(birth,t[1]))
    a=t[1]; m=len(a)
    dp={(i,i+1):birth(a[i]) for i in range(m)}
    for w in range(2,m+1):
        for i in range(m-w+1):
            j=i+w
            dp[i,j]=1+min(max(dp[i,k],dp[k,j]) for k in range(i+1,j))
    return dp[0,m]

def greedy_birth(costs):
    p=0
    for c in costs:
        w=1<<c;p=((p+w-1)//w+1)*w
    return (p-1).bit_length()

def embed(t,e):
    """Selected leaves in T_e; requires e >= exact birth(t)."""
    check(e>=birth(t),'insufficient host')
    if t[0]=='x': return 1
    n=3**(e-1)
    if t[0]=='c': return sum(embed(c,e-1)<<(i*n) for i,c in enumerate(t[1]))
    a=t[1]
    for k in range(1,len(a)):
        l=seq(a[:k]);r=seq(a[k:])
        if max(birth(l),birth(r))<=e-1:
            return embed(l,e-1)|(embed(r,e-1)<<n)
    raise RuntimeError('birth split absent')

def tree_pointwise(mask,e):
    if e==0: return 0
    if not mask: return (3**e-1)//2
    n=3**(e-1);z=(1<<n)-1
    return sum(tree_pointwise((mask>>(i*n))&z,e-1) for i in range(3))

def E(e,s): return sum(max(3**d-s,0) for d in range(e))

def tower_rows(e):
    n=3**e; out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            for d in range(e-1,-1,-1):
                a=(i//3**d)%3;b=(j//3**d)%3
                if a!=b:
                    if (b-a)%3==1: out[i]|=1<<j
                    else: out[j]|=1<<i
                    break
    return out

def H(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            z=dp[s][i]
            if not z:continue
            todo=rows[i]&~s
            while todo:
                bit=todo&-todo; todo-=bit;j=bit.bit_length()-1
                dp[s|bit][j]+=z
    return sum(dp[-1])

def val(n):
    s=0
    while n%3==0:n//=3;s+=1
    return s

def run():
    result={}
    for r in range(1,8):
        for costs in itertools.product(range(4),repeat=r):
            dp={(i,i+1):costs[i] for i in range(r)}
            for width in range(2,r+1):
                for i in range(r-width+1):
                    j=i+width
                    dp[i,j]=1+min(max(dp[i,k],dp[k,j]) for k in range(i+1,j))
            check(greedy_birth(costs)==dp[0,r],'greedy alphabetic code')
    result['greedy_vs_minimax_cost_lists']=sum(4**r for r in range(1,8))
    P={(0,0):1,(1,0):1};poly=[]
    for e in range(5):
        n=3**e;I=(n-1)//2
        check(set(a for _,a in P)==set(range(I+1)), 'missing stabilizer exponent')
        check(sum(c*3**(I-a) for (s,a),c in P.items())==2**n,'orbit-size mass')
        check(sum(c for (s,a),c in P.items() if a==0)==2,'free orbits')
        poly.append({'depth':e,'vertices':n,'rooted_extension_types':sum(P.values()),'stabilizer_exponents':[0,I], 'coefficients':len(P)})
        if e<=2:
            seen={}
            for w in itertools.product((0,1),repeat=n):
                key,a=canon_col(w);seen[key]=(sum(w),a)
            brute=collections.Counter(seen.values())
            check(dict(brute)==P,'literal orbit classification')
        if e<4:P=orbit_step(P)
    result['marked_extension_orbits']=poly
    birthsets=[]
    for e in range(3):
        seen={induced(m,e) for m in range(1,1<<(3**e))}
        check(all(birth(t)<=e for t in seen),'observed birth lower bound')
        for t in seen:
            h=birth(t);check(induced(embed(t,h),h)==t,'canonical embedding')
        births=collections.Counter(birth(t) for t in seen)
        size=collections.Counter(m.bit_count() for m in range(1,1<<(3**e)))
        birthsets.append({'depth':e,'nonempty_induced_types':len(seen),'birth_depth_distribution':dict(sorted(births.items()))})
    result['induced_host_types']=birthsets
    # Generate complete types from child types, not subsets, at the next level.
    types={None,X}
    for e in range(1,4):
        nxt={None}
        for kids in itertools.product(types,repeat=3):
            have=[i for i,t in enumerate(kids) if t is not None]
            if not have:continue
            if len(have)==3:t=cyc(kids)
            elif len(have)==1:t=kids[have[0]]
            else:
                i,j=have
                if (j-i)%3!=1:i,j=j,i
                t=seq([kids[i],kids[j]])
            nxt.add(t)
        types=nxt
        check(all(birth(t)<=e for t in types if t is not None),'grammar birth')
    nonempty=types-{None}
    for t in nonempty:
        h=birth(t);check(induced(embed(t,h),h)==t,'all-type embedding')
    result['depth3_complete_induced_grammar']={
        'types':len(nonempty),
        'birth_depth_distribution':dict(sorted(collections.Counter(birth(t) for t in nonempty).items()))}

    for e in range(3):
        n=3**e
        mins={s:100000 for s in range(n+1)}
        for mask in range(1<<n): mins[mask.bit_count()]=min(mins[mask.bit_count()],tree_pointwise(mask,e))
        check(all(mins[s]==E(e,s) for s in mins),'exact stabilizer floor')
    # Literal Hamilton controls for every one-vertex extension of the nine-tower.
    base=tower_rows(2);hc=H(base);check(hc==3159,'inherited H(T2)')
    # Independent full adjacency-based automorphism enumeration (not tree recursion).
    pairs=[(i,j,(base[i]>>j)&1) for i in range(9) for j in range(i+1,9)]
    autos=[p for p in itertools.permutations(range(9))
           if all(((base[p[i]]>>p[j])&1)==b for i,j,b in pairs)]
    check(len(autos)==81,'full automorphism group')
    seen_masks=set();actual=collections.Counter()
    for mask in range(512):
        if mask in seen_masks:continue
        orb={sum(((mask>>i)&1)<<p[i] for i in range(9)) for p in autos}
        seen_masks.update(orb)
        stabilizer=81//len(orb);actual[mask.bit_count(),val(stabilizer)]+=1
    expected=orbit_step(orbit_step({(0,0):1,(1,0):1}))
    check(dict(actual)==expected,'full adjacency orbit audit')
    result['independent_adjacency_automorphisms']={'order9_automorphisms':len(autos),'rooted_orbits':sum(actual.values())}

    gcd=0; histogram=collections.Counter()
    for mask in range(512):
        rows=[r|(((mask>>i)&1)<<9) for i,r in enumerate(base)]
        rows.append(511^mask)
        count=H(rows);gcd=math.gcd(gcd,count)
        word=tuple((mask>>i)&1 for i in range(9));a=canon_col(word)[1]
        check(count%3**max(a,E(2,2))==0,'physical divisibility')
        histogram[val(count)]+=1
    result['all_512_one_vertex_extensions']={'base_H':hc,'gcd':gcd,'valuation_histogram':dict(sorted(histogram.items()))}
    result['endpoint_protection']=[{'depth':e,'vertices':3**e,'exponents_k0_to_5':[E(e,2*k) for k in range(6)]} for e in range(1,9)]
    # Ordered-cost hostile: same multiset, different exact host birth.
    one=X;tri=cyc([X,X,X]);a=seq([tri,X,X]);b=seq([X,tri,X])
    check(birth(a)==2 and birth(b)==3,'alphabetic order hostile')
    result['ordered_cost_hostile']={'costs_100':birth(a),'costs_010':birth(b),'orders':[5,5]}
    # Many genuinely larger instances constructed without listing host tournaments.
    examples=[]
    for r in [1,2,3,4,7,8,9,50,200]:
        t=seq([tri]*r);h=1+(r-1).bit_length()
        check(birth(t)==h,'homogeneous alphabetic law')
        examples.append({'triangles':r,'order':3*r,'birth_depth':h,'smallest_balanced_host_order':3**h})
    result['all_order_family_examples']=examples
    return result

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
