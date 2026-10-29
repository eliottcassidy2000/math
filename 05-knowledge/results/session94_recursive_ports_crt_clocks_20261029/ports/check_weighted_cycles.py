#!/usr/bin/env python3
"""Exact cyclic-tree palette inventories; stdlib, always-active hostile tests."""
from fractions import Fraction
from itertools import product, combinations
from collections import Counter
from math import comb, gcd

CHECKS=0
def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

def mul(a,b,R):
    c=[0]*(R+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:R+1-i]):
                c[i+j]+=x*y
    return c

def power(a,n,R):
    b=[1]+[0]*R
    while n:
        if n&1: b=mul(b,a,R)
        a=mul(a,a,R); n//=2
    return b

def subst(a,l,R):
    b=[0]*(R+1)
    for i,x in enumerate(a):
        if i*l<=R: b[i*l]=x
    return b

def phi(n):
    return sum(gcd(i,n)==1 for i in range(1,n+1))

def step(a,d,R):
    """Cycle-index wreath substitution, not a fixed-point-only recursion."""
    out=[0]*(R+1)
    for l in range(1,d+1):
        if d%l==0:
            term=power(subst(a,l,R),d//l,R)
            for i,t in enumerate(term): out[i]+=phi(l)*t
    if any(x%d for x in out): raise ArithmeticError('nonintegral inventory')
    return [x//d for x in out]

def inventory(d,k,Q,R,centralizer=False):
    """Q is a finite signed-weight palette inventory. Truncation exact to R."""
    if d<3 or not d&1: raise ValueError('odd cyclic alphabet >=3 required')
    q=(list(Q)+[0]*(R+1))[:R+1]
    if not centralizer:
        for _ in range(k): q=step(q,d,R)
        return q
    out=q[:]; a=mul(q,q,R)
    for _ in range(k):
        out=mul(out,power(a,(d-1)//2,R),R)
        a=step(a,d,R)
    return out

def tree_group(d,k):
    """Literal source group; only small depths used."""
    if k==0: return [(0,)]
    child=tree_group(d,k-1); n=d**(k-1)
    return [tuple(((i//n+r)%d)*n+gs[i//n][i%n] for i in range(d*n))
            for r in range(d) for gs in product(child,repeat=d)]

def neg_vertex(v,d,k):
    out=0; place=1
    for _ in range(k):
        out+=((-v)%d)*place; v//=d; place*=d
    return out

def central_group(d,k):
    if k==0: return [(0,)]
    n=d**(k-1); prev=central_group(d,k-1); child=tree_group(d,k-1)
    out=[]
    for z in prev:
        for gs in product(child,repeat=(d-1)//2):
            blocks=[None]*d; blocks[0]=z
            for b,g in enumerate(gs,1):
                blocks[b]=g
                blocks[d-b]=tuple(neg_vertex(g[neg_vertex(i,d,k-1)],d,k-1) for i in range(n))
            out.append(tuple(b*n+blocks[b][i] for b in range(d) for i in range(n)))
    return out

def cycles(g):
    seen=set(); out=[]
    for i in range(len(g)):
        if i in seen: continue
        j=i; l=0
        while j not in seen:
            seen.add(j); l+=1; j=g[j]
        out.append(l)
    return sorted(out)

def literal_inventory(group,Q,R):
    total=[0]*(R+1); q=(list(Q)+[0]*(R+1))[:R+1]
    for g in group:
        a=[1]+[0]*R
        for l in cycles(g): a=mul(a,subst(q,l,R),R)
        for i,x in enumerate(a): total[i]+=x
    check(all(x%len(group)==0 for x in total),'literal integral Burnside')
    return [x//len(group) for x in total]

def cycle_dist(group):
    return Counter(tuple(cycles(g)) for g in group)

def cycle_mul(a,b):
    out=Counter()
    for x,p in a.items():
        for y,q in b.items():out[tuple(sorted(x+y))]+=p*q
    return out

def cycle_step(a,d):
    out=Counter()
    for l in range(1,d+1):
        if d%l: continue
        scaled={tuple(l*x for x in key):p for key,p in a.items()}
        term={():Fraction(1)}
        for _ in range(d//l):term=cycle_mul(term,scaled)
        for key,p in term.items():out[key]+=Fraction(phi(l),d)*p
    return out

def cycle_law(d,k,centralizer=False):
    a={(1,):Fraction(1)}
    if not centralizer:
        for _ in range(k):a=cycle_step(a,d)
        return a
    out={(1,):Fraction(1)}
    for _ in range(k):
        paired={tuple(sorted(key+key)):p for key,p in a.items()}
        for _ in range((d-1)//2):out=cycle_mul(out,paired)
        a=cycle_step(a,d)
    return out

def binomial_table(d,R):
    """Exact Newton coefficients; theorem bounds each column's degree by r."""
    q=[1,1]+[0]*max(0,R-1); q=q[:R+1]
    a=mul(q,q,R); current=q[:]; rows=[current[:]]
    for _ in range(R):
        current=mul(current,power(a,(d-1)//2,R),R)
        rows.append(current[:]); a=step(a,d,R)
    out=[]
    for r in range(R+1):
        row=[x[r] for x in rows[:r+1]]; coeff=[]
        while row:
            coeff.append(row[0]); row=[b-a for a,b in zip(row,row[1:])]
        out.append(coeff)
    return out

def fast_depth(k,table):
    return [sum(b*comb(k,l) for l,b in enumerate(row) if l<=k) for row in table]

def active_levels(subset,k):
    """Intrinsic ternary centralizer port depth set, from folded pair LCAs."""
    groups={}; active=set()
    for v in subset:
        ds=[]; z=v
        for _ in range(k): ds.append(z%3); z//=3
        ds=ds[::-1]
        if not any(ds):continue
        a=next(i for i,x in enumerate(ds) if x)
        active.add(a)
        suffix=tuple((x if ds[a]==1 else -x)%3 for x in ds[a+1:])
        groups.setdefault(a,set()).add(suffix)
    for a,words in groups.items():
        for x,y in combinations(words,2):
            l=next(i for i,(u,v) in enumerate(zip(x,y)) if u!=v)
            active.add(a+1+l)
    return active

def main():
    groups=[]
    for d, depths, cdepths in [(3,range(3),range(4)),(5,range(2),range(3)),(9,range(2),range(3))]:
        for central,ks in [(False,depths),(True,cdepths)]:
            for k in ks:
                G=central_group(d,k) if central else tree_group(d,k)
                check(len(set(G))==len(G),'group duplicates')
                literal=cycle_dist(G)
                expected=cycle_law(d,k,central)
                check({a:Fraction(v,len(G)) for a,v in literal.items()}==expected,'full cycle law')
                for Q in ([1,1],[2,1],[1,2,1],[1,-1]):
                    R=min(8,d**k)
                    check(inventory(d,k,Q,R,central)==literal_inventory(G,Q,R),'weighted palette')
                fixed=[sum(i==g[i] for i in range(len(g))) for g in G]
                mean=Fraction(sum(fixed),len(G)); second=Fraction(sum(x*x for x in fixed),len(G))
                if central:
                    check(mean==1+(d-1)*k,'central mean')
                    check(second-mean*mean==(d-1)**2*k*(k-1),'central variance')
                else:
                    check(mean==1,'critical mean')
                    check(second-1==(d-1)*k,'critical variance')
                groups.append((d,k,central,len(G)))
    # Literal fixed subsets, independently of cycle products, on all 512 subsets.
    G=central_group(3,2); counts=[0]*10; seen=set()
    for mask in range(1<<9):
        if mask in seen:continue
        orb={sum(1<<g[i] for i in range(9) if mask>>i&1) for g in G}
        seen.update(orb); counts[mask.bit_count()]+=1
    check(counts==inventory(3,2,[1,1],9,True),'explicit subset orbits')
    # Fixed-point-only tuple moments miss invariant cycles of selected vertices.
    naive=sum(Fraction(comb(sum(i==g[i] for i in range(9)),3),len(G)) for g in G)
    check(naive!=counts[3],'fixed-point-only subset hostile')
    # Paired C3 action: 4 colors on 3 vertices, not 2 independent necklaces.
    actual=sum(inventory(3,2,[1,1],9,True))
    independent=8*4*4  # three fixed vertices, two independent binary C3 necklaces
    doubled=8*4        # erroneous replacement of paired 3-cycles by 6-cycles
    check(actual==192 and independent==128 and doubled==32,'paired palette hostile')
    # Cross-side ordered one-port profiles: three joint types, one marginal pair.
    cross={(min(((a+t)%3,(b+t)%3) for t in range(3))) for a in range(3) for b in range(3)}
    check(len(cross)==3,'diagonal shared ports')
    # A composite cyclic alphabet needs ALL divisors, not merely identity/full cycle.
    z9=cycle_law(9,1)
    check(z9[(3,3,3)]==Fraction(2,9),'composite rotation clock')
    check(inventory(9,1,[1,1],9)[3]==10,'C9 triple necklaces')
    # Same entire fixed-point law, different invariant subset inventories.
    cyclic=tree_group(9,1)
    elementary=[tuple(((i//3+a)%3)*3+(i%3+b)%3 for i in range(9))
                for a in range(3) for b in range(3)]
    f=lambda group:Counter(sum(i==g[i] for i in range(9)) for g in group)
    check(f(cyclic)==f(elementary),'same all-order moments hostile')
    check(literal_inventory(cyclic,[1,1],3)[3]==10,'cyclic nine subsets')
    check(literal_inventory(elementary,[1,1],3)[3]==12,'elementary nine subsets')
    # Exact degree-truncated recurrence at depths impossible for group enumeration.
    for d in (3,5,9):
        low=inventory(d,11,[1,1],12,True)
        high=inventory(d,11,[1,1],15,True)
        check(low==high[:13],'truncation compatibility')
        check(low[0]==1 and low[1]==1+(d-1)*11,'deep low coefficients')
    # All-order polynomial theorem: degree r, leading (d-1)^r.
    for d in (3,5,9):
        samples=[inventory(d,k,[1,1],8,True) for k in range(12)]
        for r in range(1,9):
            row=[x[r] for x in samples]
            for _ in range(r): row=[b-a for a,b in zip(row,row[1:])]
            import math
            check(all(x==math.factorial(r)*(d-1)**r for x in row),'subset polynomial leading term')
    for k in range(14):
        inv=inventory(3,k,[1,1],3,True)
        check(inv[2]==k*(4*k-1),'quadratic subsets')
        check(inv[3]==k*(2*k-1)*(4*k-3),'cubic subsets')
    # Active-level compression, verified against actual subset orbits.
    table=binomial_table(3,12)
    for k in (2,3):
        G=central_group(3,k)
        for r in range(1,5):
            seen=set(); hist=Counter()
            for sub in combinations(range(3**k),r):
                if sub in seen:continue
                orb={tuple(sorted(g[v] for v in sub)) for g in G}
                seen.update(orb)
                levels=active_levels(sub,k)
                check(len(levels)<=r,'active depth bound')
                check(all(active_levels(x,k)==levels for x in orb),'intrinsic active depths')
                hist[len(levels)]+=1
            check(dict(hist)=={l:b*comb(k,l) for l,b in enumerate(table[r]) if l<=k and b},'saturated template counts')
    for d in (3,5,9):
        tab=binomial_table(d,12)
        check(all(b>=0 for row in tab for b in row),'nonnegative templates')
        check(fast_depth(50,tab)==inventory(d,50,[1,1],12,True),'fast depth exact')
        huge=fast_depth(10**6,tab)
        check(huge[1]==1+(d-1)*10**6,'huge depth one port')
    # Tag-weighted named ports retain signed costs and paired source multiplicities.
    k=3; u=[2,-3,7]; v=[-1,4,11]; root_weight=5
    weights=[]
    for vertex in range(27):
        ds=[]; z=vertex
        for _ in range(k):ds.append(z%3); z//=3
        ds=ds[::-1]
        if not any(ds):weights.append(root_weight);continue
        a=next(i for i,x in enumerate(ds) if x); j=k-a-1
        weights.append((u if ds[a]==1 else v)[j])
    C=central_group(3,3)
    actual=Counter(sum(weights[i] for i in range(27) if g[i]==i) for g in C)
    expected=Counter({root_weight:1})
    for j in range(3):
        A=tree_group(3,j)
        dist=Counter(sum(i==g[i] for i in range(len(g))) for g in A)
        nxt=Counter()
        for old,mass in expected.items():
            for x,count in dist.items():nxt[old+(u[j]+v[j])*x]+=mass*count
        expected=nxt
    check(actual==expected,'weighted tag fixed-source law')
    check(all(weights[i]==weights[g[i]] for g in C for i in range(27)),'invariant tag weights')
    dist=Counter(sum(i==g[i] for i in range(9)) for g in tree_group(3,2))
    check(Fraction(sum(dist),len(dist))==Fraction(9,2),'support uniform hostile')
    check(Fraction(sum(x*c for x,c in dist.items()),sum(dist.values()))==1,'group uniform law')
    print('groups:',groups)
    print('C2 subset inventory:',counts)
    print('fixed-point-only triples:',naive,'actual subsets:',counts[3])
    print('C11 subsets through 12:',inventory(3,11,[1,1],12,True))
    print('C12 subsets through 12:',inventory(3,12,[1,1],12,True))
    print('always-active checks:',CHECKS)

if __name__=='__main__':main()
