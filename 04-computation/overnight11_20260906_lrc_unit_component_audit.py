"""Independent exact referee of the primitive-unit actual LRC entry theorem.

Standard-library arithmetic only; no producer imports. Bounded toy universes,
all primitive pairs under the stated sum cap, explicit rational clock lifts.
"""
from fractions import Fraction as F
from math import gcd, prod, isqrt
from functools import reduce
from itertools import combinations
import sys

sys.stdout.reconfigure(newline="\n")
Q = 91**6
gates = 0


def need(ok, message):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(message)


def dist(w, x):
    y = w*x
    residue = y - y.numerator//y.denominator
    return min(residue, 1-residue)


def ceil(x):
    return -((-x.numerator)//x.denominator)


def pair_center(p, q):
    s = p+q
    j = pow(p, -1, s)*(s//2) % s
    z = F(j, s)
    need(dist(p,z) == dist(q,z) == F(s//2,s), "exact simultaneous pair center")
    need(dist(p,z) >= F(1,3), "elementary pair one-third clearance")
    return z


def transport(fixed_clock, multiplier, modulus, left, right):
    # Physical lifts x=(fixed_clock+j)/modulus; select the integer k first.
    need(gcd(multiplier,modulus)==1, "complete lifted clock requires coprimality")
    need(right-left >= F(1,modulus), "closed interval covers grid spacing")
    k = ceil(modulus*left-multiplier*fixed_clock)
    target = (multiplier*fixed_clock+k)/modulus
    j = pow(multiplier,-1,modulus)*k % modulus
    x = (fixed_clock+j)/modulus
    need(left <= target <= right, "ceil construction lands in closed interval")
    delta = multiplier*x-target
    need(delta.denominator==1, "clock equality modulo integers")
    return x


def main():
    divisions = 0
    for bound in range(1,13):
        for K in range(1,bound+1):
            represented = {a*K+r for a in range(bound+1) for r in range(bound+1)}
            need(represented == set(range(bound*(K+1)+1)), "independent Minkowski sum is full integer interval")
            for x in range(1,bound*(K+1)+1):
                a = min(x//K,bound)
                r = x-a*K
                need(x==a*K+r and 0<=a<=bound and 0<=r<=bound, "clamped quotient covers full box")
                divisions += 1
            need(bound*(K+1)+1 not in represented, "first positive exterior point excluded")
            need(max(abs(a*K+r) for a in (-bound,bound) for r in (-bound,bound))==bound*(K+1), "signed box cannot enlarge the radius")
    need(3 not in {6*a+r for a in range(-2,3) for r in range(-2,3)}, "dropping K<=Q creates a hole")

    crosses = 0
    for bound in (5,8,13):
        for K in range(2,bound+1):
            for t in range(1,bound+1):
                for p in (1,4,6,11):
                    for g in range(1,20):
                        if gcd(t,g)!=1:
                            continue
                        delta = gcd(t,p)
                        x = g*p//delta
                        if x > bound*(K+1):
                            continue
                        a, r = min(x//K,bound), x-min(x//K,bound)*K
                        c = t//delta
                        need(c*g*p == a*t*K+r*t, "literal cleared relation")
                        need(1<=c<=bound and 0<=a<=bound and 0<=r<=bound and a+r>0, "height and nonzero component restrictions")
                        need(c*g*p>0, "nonzero pair restriction excludes decoder span")
                        crosses += 1

    pairs = [(p,s-p) for s in range(3,357) for p in range(1,(s+1)//2) if gcd(p,s)==1]
    need(max(p for p,q in pairs)==177 and max(q for p,q in pairs)==355, "all primitive sum-cap pair extrema")
    max_t = 0
    for p,q in pairs:
        need(p<q and p+q<=356 and gcd(p,q)==1, "pair universe directly enumerated by sum")
        tmax = (21*q-1)//11
        max_t = max(max_t,tmax)
        need(tmax<=677<Q and 11*tmax<21*q<=11*(tmax+1), "strict scale branch integer cutoff")
        need(Q>42*p, "all-parameter scale forcing constant")
        pair_center(p,q)
    need(max_t==677, "core-scale cutoff attained")
    need(Q==567869252041, "literal inherited box height")
    need(2*(F(1,12)-F(1,14))==F(1,42), "eleven-core safe arc constant")
    need(2*(F(1,3)-F(1,14))==F(11,21), "pair safe arc constant")

    clock_cases = []
    core = tuple(range(1,12))
    y = F(1,12)
    for p,q in ((1,2),(1,3),(2,11),(1,355),(177,178)):
        z = pair_center(p,q)
        # Both equality (p,q,t)=(2,11,21) and strict grid containment occur.
        for t in sorted({ceil(F(21*q,11)), ceil(F(21*q,11))+1}):
            for g in (1,t+1,3*t+1):
                x = transport(y,g,t,z-F(11,42*q),z+F(11,42*q))
                speeds = tuple(t*v for v in core)+(g*p,g*q)
                need(min(dist(v,x) for v in speeds)>=F(1,14), "large core-scale literal physical safety")
                clock_cases.append(("large",t,g,p,q))
        for t in (1,2,6,17,677):
            if 11*t>=21*q:
                continue
            delta = gcd(t,p)
            g = 2*Q*t*max(core)+1
            need(g*p>delta*Q*(max(core)+1) and gcd(t,g)==1, "forced-scale positive control")
            need(F(delta*Q*(max(core)+1),p)>42*max(core), "exact general lower-bound comparison")
            x = transport(z,t,g,y-F(1,84*max(core)),y+F(1,84*max(core)))
            speeds = tuple(t*v for v in core)+(g*p,g*q)
            need(min(dist(v,x) for v in speeds)>=F(1,14), "small core-scale literal physical safety")
            clock_cases.append(("small",t,g,p,q))
    need(any(t==21 and q==11 for _,t,g,p,q in clock_cases), "closed-arc equality represented")

    grid_cases = 0
    for modulus in range(2,20):
        for multiplier in range(1,modulus+1):
            if gcd(multiplier,modulus)!=1:
                continue
            for fixed in (F(0),F(2,7)):
                for left in (F(-3,11),F(0),F(7,13)):
                    x = transport(fixed,multiplier,modulus,left,left+F(1,modulus))
                    need(0<=x<1+F(1,modulus), "literal lift representative finite")
                    grid_cases += 1
    residues = {F(2*j,4)%1 for j in range(4)}
    need(not any(F(1,8)<=x<=F(3,8) for x in residues), "noncoprime clock misses a nominal-spacing arc")
    need(F(3,8)-F(1,8)==F(1,4), "hostile retains claimed grid-width comparison")

    U = (1,4,6,8,10,12,14,15,16,18,22)
    g = 2**45
    speeds = U+(g,3*g)
    # Reconstruct allowed pair sums by multiplicative generation, not trial
    # factorization of each ratio as in the producer.
    inert_primes = [p for p in range(2,357) if p%3==2 and all(p%d for d in range(2,isqrt(p)+1))]
    allowed = {1}
    for p in inert_primes:
        allowed |= {s*p**e for s in tuple(allowed) for e in (1,2) if s*p**e<=356}
    need(sum(p+q in allowed for p,q in pairs)==5855, "independently generated inherited atlas cardinality")
    adjacency = [set() for _ in speeds]
    for i,j in combinations(range(13),2):
        d = gcd(speeds[i],speeds[j])
        if (speeds[i]+speeds[j])//d in allowed:
            adjacency[i].add(j)
            adjacency[j].add(i)
    unseen = set(range(13))
    components = []
    while unseen:
        stack = [min(unseen)]
        component = set()
        while stack:
            i = stack.pop()
            if i in component:
                continue
            component.add(i)
            stack.extend(adjacency[i]-component)
        unseen -= component
        components.append(component)
    need(components==[set(range(11)),{11,12}], "independent decoder BFS of nonvacuous control")
    need(sum(speeds)<=Q**2 and g>2*Q*max(U), "physical finite box and crossing exclusion by dominance")
    need(reduce(gcd,speeds)==1 and len(set(speeds))==13, "physical primitive distinct nonvacuity control")
    phase = F(9,19)
    need(min(dist(v,phase) for v in U)>=F(1,12), "literal canonical core supplier phase")
    z = pair_center(1,3)
    x = transport(z,1,g,phase-F(1,84*max(U)),phase+F(1,84*max(U)))
    need(min(dist(v,x) for v in speeds)>=F(1,14), "actual-entry positive control literal witness")

    primes = (37,43,61,67,73,79,97,103,127)
    P = 15*prod(primes)
    V = tuple(2*P//r for r in primes)+(P//3,P//5)
    need(reduce(gcd,V)==1 and min(V)>1, "primitive nonunit canonical hostile")
    need(max(V)==237907127334685115>Q, "nonunit maximum is not bounded by Q")
    need(max(max(a,b)//gcd(a,b) for a,b in combinations(V,2))==127, "all internal heights still small")
    for K in (11,Q,42*95):
        for p in (1,17,177):
            need(F(14*Q*(K+1),3*p)>29*K, "d3 consumer scale threshold")
    print("STATUS: PASS; independent actual-entry primitive-unit theorem referee")
    print("BOUND CONTROLS:",divisions,"full interval points;",crosses,"cleared crossing identities")
    print("PAIR UNIVERSE:",len(pairs),"all coprime p<q, p+q<=356; maximum residual t",max_t)
    print("CLOCK CONTROLS:",len(clock_cases),"physical branch witnesses;",grid_cases,"closed-grid transports")
    print("NONVACUITY: independently generated inert sums + decoder BFS + crossing dominance + literal witness")
    print("HOSTILES: omitted K bound; noncoprime clock; primitive nonunit with K>Q")
    print("SCOPE: inherited actual-box equality; cited eleven-speed LRC; weak clearance >=1/14")
    print("ACTIVE GATES:",gates)


if __name__=="__main__":
    main()
