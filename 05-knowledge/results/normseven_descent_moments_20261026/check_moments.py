"""Exact finite universes; no assert statements and no external dependencies."""
import itertools
import json
from collections import Counter


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def span(xs):
    s = {0}
    for x in xs:
        s |= {v ^ x for v in tuple(s)}
    return frozenset(s)


def add(a, b):
    return frozenset(x ^ y for x in a for y in b)


def subspaces(n):
    found = {frozenset({0})}
    todo = list(found)
    while todo:
        s = todo.pop()
        for x in range(1 << n):
            if x not in s:
                t = add(s, (0, x))
                if t not in found:
                    found.add(t)
                    todo.append(t)
    return sorted(found, key=lambda s: (len(s), sorted(s)))


def cosets(n, h):
    return sorted({frozenset(x ^ y for y in h) for x in range(1 << n)},
                  key=lambda s: min(s))


def finite_orbits():
    cases = 0
    # Every inertia subspace, every translation, canonical complementary
    # H/K in every dimension split. Arbitrary complements are conjugate by GL.
    for n in range(1, 5):
        size = 1 << n
        for hdim in range(n + 1):
            h = span(1 << i for i in range(hdim))
            k = span(1 << i for i in range(hdim, n))
            for inertia in subspaces(n):
                hi, ki = add(h, inertia), add(k, inertia)
                j = hi & ki
                source = cosets(n, inertia)
                left, right = cosets(n, hi), cosets(n, ki)
                target = list(itertools.product(left, right))
                # An independent literal orbit representation, including map.
                map_to = []
                for s in source:
                    x = min(s)
                    map_to.append((next(a for a in left if x in a),
                                   next(b for b in right if x in b)))
                need(set(map_to) == set(target), 'diagonal transitivity')
                need(set(Counter(map_to).values()) == {len(j)//len(inertia)}, 'sheet degree')
                d, e = len(source), len(target)
                for a in range(size):
                    values = []
                    for g in range(size):
                        f = a ^ g
                        x = sum(frozenset(v ^ f for v in s) == s for s in source)
                        y = sum(frozenset(v ^ f for v in l) == l and
                                frozenset(v ^ f for v in r) == r for l, r in target)
                        need(x == (d if f in inertia else 0), 'source fixed points')
                        need(y == (e if f in j else 0), 'coarse fixed points')
                        values.append((x, y))
                    need(sum(x for x,y in values) == size, 'source mean')
                    need(sum(y for x,y in values) == size, 'target mean')
                    need(sum(x*x for x,y in values) == size*d, 'source second moment')
                    need(sum(y*y for x,y in values) == size*e, 'target second moment')
                    need(sum(x*y for x,y in values) == size*e, 'mixed moment')
                    need(sum((x-y)**2 for x,y in values) == size*(d-e), 'energy')
                    cases += 1
    return {'dimensions': [1,2,3,4], 'inertia_translation_split_cases': cases,
            'scope': 'all inertia subspaces and translations; coordinate complements in every split'}


class Field:
    def __init__(self, p, degree):
        self.p, self.degree, self.q = p, degree, p**degree
        self.nu = next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
        self.roots = Counter(self.mul(x,x) for x in range(self.q))
        self.ns = next(x for x in range(1,self.q) if not self.roots[x])
    def plus(self, a, b):
        p = self.p
        return ((a%p+b%p)%p) + p*((a//p+b//p)%p)
    def neg(self,a):
        return self.plus(0,(-a%self.p)+self.p*((- (a//self.p))%self.p))
    def mul(self, a,b):
        p = self.p
        if self.degree == 1:
            return a*b%p
        x,y,z,w = a%p,a//p,b%p,b//p
        return (x*z+self.nu*y*w)%p + p*((x*w+y*z)%p)
    def inv(self,a):
        need(a != 0,'zero inverse')
        return next(b for b in range(1,self.q) if self.mul(a,b)==1)


def character_code(labels):
    return [sum((1 << i) for i,x in enumerate(labels) if (x & c).bit_count()%2)
            for c in range(8)]


PAIRS = [([1,2,4,5,6,3,7],[1,2,4,6,3,7,5]),
         ([1,2,4,5,7,3,6],[1,2,4,7,3,6,5])]


def field_checks(p, degree):
    f = Field(p,degree)
    q = f.q
    points = [(1,x,z) for x in range(q) for z in range(q)]
    points += [(0,1,z) for z in range(q)] + [(0,0,1)]
    twists = [[1]+[f.ns if mask>>i&1 else 1 for i in range(6)] for mask in range(64)]
    rows = []
    for pair in PAIRS:
        codes = list(map(character_code,pair))
        # Projective sign gauge: bit 6 is zero, inertia_6 is all six bits.
        h,k = [frozenset(g for g in range(64) if all((g&c).bit_count()%2==0 for c in code)) for code in codes]
        need(len(h)==len(k)==8 and h&k=={0}, 'complementary kernels')
        strata, sums, triples = Counter(), [0]*5, []
        for E,x,z in points:
            L = [E,f.plus(E,x),f.plus(E,f.neg(x)),f.plus(E,z),f.plus(E,f.neg(z)),
                 f.plus(f.plus(E,x),z),f.plus(f.plus(E,f.neg(x)),f.neg(z))]
            zero = [i for i,v in enumerate(L) if v==0]
            zmask = sum(1<<i for i in zero)
            inertia = span((1<<i) if i<6 else 63 for i in zero)
            j = add(h,inertia)&add(k,inertia)
            r = len(zero)
            e = (len(j)//len(inertia)).bit_length()-1
            strata[(r,e)] += 1
            if r==3:
                triples.append({'lines':zero,'h_intersection':len(h&inertia),
                                'k_intersection':len(k&inertia),'missing_bits':e})
            active = [[c for c in code if not c&zmask] for code in codes]
            pivot = next(i for i,v in enumerate(L) if v)
            local = []
            for delta in twists:
                vals = [f.mul(a,b) for a,b in zip(delta,L)]
                inv = f.inv(vals[pivot])
                # Actual projective root fibre, by enumerating squares in F_q.
                count_x = 1
                for i,v in enumerate(vals):
                    if i != pivot:
                        count_x *= f.roots[f.mul(v,inv)]
                # Full invariant algebra: all surviving unit characters.
                signmask = sum(1<<i for i,v in enumerate(vals) if v and not f.roots[v])
                count_y = 1
                for code in active:
                    count_y *= len(code) if all((signmask&c).bit_count()%2==0 for c in code) else 0
                g = (signmask ^ (127 if signmask>>6&1 else 0)) & 63
                need(count_x == (64//len(inertia) if g in inertia else 0),'root/coset agreement')
                need(count_y == (64//len(j) if g in j else 0),'character/coset agreement')
                local.append((count_x,count_y))
            need(sum(x for x,y in local)==64 and sum(y for x,y in local)==64,'local means')
            terms = [sum(x*x for x,y in local),sum(y*y for x,y in local),
                     sum(x*y for x,y in local),sum((x-y)**2 for x,y in local),sum(x-y for x,y in local)]
            for i,t in enumerate(terms):
                sums[i] += t
        expected = [64*q*q-160*q+(136 if p==3 else 152),
                    64*q*q-272*q+(351 if p==3 else 403)]
        need(sums[0]==64*expected[0] and sums[1]==64*expected[1],'global local-moment polynomials')
        need(sums[2]==sums[1] and sums[3]==sums[0]-sums[1] and sums[4]==0,'global identities')
        rows.append({'strata':{f'r{r}_e{e}':v for (r,e),v in sorted(strata.items())},
                     'triples':sorted(triples,key=lambda d:d['lines']),
                     'source_second':sums[0]//64,'coarse_second':sums[1]//64,
                     'mixed':sums[2]//64,'energy':sums[3]//64})
    return {'p':p,'degree':degree,'q':q,'base_points':len(points),'twists':64,'pairs':rows}


def rational_hostile():
    # delta=5/6, t=4/5 gives residual class 5. Direct finite reductions,
    # excluding 2,3,5 (bad models/denominators); no density claim.
    primes = [p for p in range(7,301) if all(p%d for d in range(2,int(p**.5)+1))]
    rows = []
    for p in primes:
        roots = Counter(x*x%p for x in range(p))
        inv5 = pow(5,-1,p)
        vals = [1,0,1,9*inv5%p,inv5,4*inv5%p,1]
        xcount = 1
        for v in vals[1:]: xcount *= roots[v]
        ycount = 1
        for code in map(character_code,PAIRS[0]):
            active = [c for c in code if not c&2]
            def sq(c):
                prod = 1
                for i,v in enumerate(vals):
                    if c>>i&1: prod=prod*v%p
                return roots[prod]>0
            ycount *= len(active) if all(sq(c) for c in active) else 0
        need(ycount==16, 'fixed rational quotient pair persists')
        need(xcount==(32 if roots[5] else 0),'residual five criterion')
        rows.append([p,xcount,ycount])
    return {'excluded_primes':[2,3,5],'range':'7 <= p <= 300','rows':rows}


if __name__ == '__main__':
    result = {'finite_orbits':finite_orbits(),
              'fields':[field_checks(p,d) for p,d in [(3,1),(5,1),(7,1),(11,1),(13,1),(3,2),(5,2)]],
              'rational_hostile':rational_hostile()}
    print(json.dumps(result,indent=2,sort_keys=True))
