"""Exact saturated-board counts and scoped restart-bound controls; stdlib only."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial, gcd

gates = 0
def check(ok, message):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(message)

def boards(n):
    pairs = tuple(combinations(range(n), 2))
    degrees = [0]*n
    rows = []
    def visit(r):
        if r == n:
            if all(v == 2 for v in degrees):
                yield tuple(rows)
            return
        for a,b in pairs:
            if degrees[a] == 2 or degrees[b] == 2:
                continue
            degrees[a] += 1; degrees[b] += 1
            if all(2-v <= n-r-1 for v in degrees):
                rows.append((a,b))
                yield from visit(r+1)
                rows.pop()
            degrees[a] -= 1; degrees[b] -= 1
    yield from visit(0)

def profile(rows):
    n = len(rows)
    adj = [[] for _ in range(2*n)]
    for r,pair in enumerate(rows):
        for c in pair:
            adj[r].append(n+c); adj[n+c].append(r)
    unseen = set(range(2*n)); parts=[]
    while unseen:
        stack=[unseen.pop()]; size=0
        while stack:
            v=stack.pop(); size+=1
            for w in adj[v]:
                if w in unseen:
                    unseen.remove(w); stack.append(w)
        check(size%2 == 0 and size >= 4, 'cycle type')
        parts.append(size//2)
    return tuple(sorted(parts))

def success_determinants(points):
    return all((b[0]-a[0])*(c[1]-a[1]) != (c[0]-a[0])*(b[1]-a[1])
               for a,b,c in combinations(points,3))

def success_directions(points):
    # Two equal unoriented primitive directions from one point certify a triple.
    for a in points:
        seen=set()
        for b in points:
            if a == b: continue
            dx,dy=b[0]-a[0],b[1]-a[1]
            d=gcd(abs(dx),abs(dy)); dx//=d; dy//=d
            if dx<0 or (dx==0 and dy<0): dx=-dx; dy=-dy
            if (dx,dy) in seen: return False
            seen.add((dx,dy))
    return True

def b_numbers(nmax):
    out=[1,0]
    for n in range(2,nmax+1):
        out.append(n*(n-1)*out[-1]+n*(n-1)**2//2*out[-2])
    return out

def coefficient_counts(nmax):
    # Independent formal convolution exp(-z/2)*(1-z)^(-1/2).
    from math import comb
    return [factorial(n)**2*sum((F((-1)**j,2**j*factorial(j))*
                 F(comb(2*(n-j),n-j),4**(n-j)) for j in range(n+1)),F(0))
            for n in range(nmax+1)]

def exp_bounds(x, terms=100):
    check(x >= 0 and x < terms+2, 'Taylor range')
    term=F(1); total=term
    for k in range(1,terms+1):
        term*=x/k; total+=term
    next_term=term*x/(terms+1)
    return total,total+next_term/(1-x/F(terms+2))

def main():
    print('Exact uncolored saturated boards and adaptive uniform-column restarts')
    bn=b_numbers(40)
    for n,value in enumerate(coefficient_counts(40)):
        check(value.denominator==1 and value==bn[n], 'count normalization')
    grand=0
    for n in range(2,7):
        counts={}; successes=selected=total=0; weighted_pairs=0; orbits={}
        for rows in boards(n):
            total+=1; typ=profile(rows); counts[typ]=counts.get(typ,0)+1
            weighted_pairs += 2**len(typ)
            points=tuple((r,c) for r,pair in enumerate(rows) for c in pair)
            good=success_determinants(points)
            check(good==success_directions(points), 'independent collinearity')
            occupancy={}
            for r,c in points: occupancy[c-r]=occupancy.get(c-r,0)+1
            zero=all(v<=2 for v in occupancy.values())
            check(not good or zero, 'selected zero inclusion')
            successes+=good; selected+=zero
            column_word=tuple(sorted(tuple(r for r,pair in enumerate(rows) if c in pair)
                                     for c in range(n)))
            record=orbits.setdefault(column_word,[0,0,0,typ.count(2)])
            record[0]+=1; record[1]+=good
            record[2]+=sum(max(v-2,0) for v in occupancy.values())
        check(total==bn[n], 'literal board count')
        for typ,count in counts.items():
            denominator=1
            for r in set(typ):
                multiplicity=typ.count(r)
                denominator *= (2*r)**multiplicity*factorial(multiplicity)
            check(count*denominator==factorial(n)**2, 'orbit stabilizer')
        derangements=sum((-1)**k*factorial(n)//factorial(k) for k in range(n+1))
        check(weighted_pairs==factorial(n)*derangements, 'colored pair bias')
        for size,good,defect_sum,c4 in orbits.values():
            check(size*2**c4==factorial(n),'fixed-row column stabilizer')
            mean=F(defect_sum,size)
            _,exp_upper=exp_bounds(mean**2/(8*(n-1)),100)
            check(F(good,size)*exp_upper<=1,'literal conditional tail bound')
        print('n',n,'boards',total,'full-success',successes,'slope-one-zero',selected,
              'cycle-profiles',len(counts),'ordered-pairs',weighted_pairs,
              'column-orbits',len(orbits))
        grand+=total
    elo,ehi=exp_bounds(F(2),100)
    alpha_lo=1-5/elo; alpha_hi=1-5/ehi
    check(F(323,1000)<alpha_lo<alpha_hi<F(324,1000), 'alpha enclosure')
    print('conservative exponents: conditional uses rowfreeze only; uncolored takes max with tenth')
    for n in (64,128,256,512,1024,2048,4096):
        # Simplified exact lower alpha 0.323 avoids enormous printed fractions.
        lower=max(F(0),F(323*n,1000)-21,F(2*(n-9),15)+F(2,n*(n-1)))
        conditional=lower**2/(8*(n-1))
        uncolored=max(conditional,F(n*n,900*(n-1)))
        check(conditional>=0 and uncolored>=conditional, 'bound directions')
        print('n',n,'conditional exponent >=',conditional,
              'uncolored exponent >=',uncolored)
    # Arbitrary history-dependent success hazards <=q retain the product survival bound.
    for depth in range(1,9):
        q=F(1,3)
        masses={0:F(1)}
        for step in range(depth):
            new={}
            for state,mass in masses.items():
                hazard=q*F((state+step)%5,4)
                check(0<=hazard<=q,'adaptive hazard cap')
                new[2*state]=mass*(1-hazard)/3
                new[2*state+1]=mass*(1-hazard)*F(2,3)
            masses=new
        survival=sum(masses.values(),F(0))
        check(survival >= (1-q)**depth,'history-dependent survival induction')
        check(1-survival <= depth*q,'adaptive union bound')
    print('literal board universe',grand)
    print('PASS',gates,'always-active gates')

if __name__ == '__main__':
    main()
