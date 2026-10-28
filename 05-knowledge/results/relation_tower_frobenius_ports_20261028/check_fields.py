"""Exact full form spaces, normal-port recovery, and hostile controls; stdlib only."""
import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path

CHECKS = 0
def need(b, message):
    global CHECKS
    CHECKS += 1
    if not b:
        raise RuntimeError(message)

def rank(rows):
    basis = {}
    for x in rows:
        while x:
            k = x.bit_length()-1
            if k in basis:
                x ^= basis[k]
            else:
                basis[k] = x
                break
    return len(basis)

def kernel(rows, n):
    piv = {}
    for x in rows:
        while x:
            k = x.bit_length()-1
            if k in piv:
                x ^= piv[k]
            else:
                piv[k] = x
                break
    ans = []
    for j in range(n):
        if j in piv:
            continue
        x = 1 << j
        for k in sorted(piv):
            if (piv[k] & x).bit_count() % 2:
                x ^= 1 << k
        ans.append(x)
    need(all(not ((r & x).bit_count() % 2) for r in rows for x in ans), 'kernel residual')
    return ans

def span(basis):
    out = [0]
    for x in basis:
        out += [y ^ x for y in out]
    return out

class Field:
    def __init__(self, r, polynomial):
        self.r, self.poly, self.N = r, polynomial, 1 << r
        self.table = [[self.multiply(a,b) for b in range(self.N)] for a in range(self.N)]
        self.sq = [self.table[x][x] for x in range(self.N)]
        # In a finite commutative quotient ring, invertibility of every nonzero
        # element proves that the chosen polynomial is irreducible.
        need(all(1 in self.table[x] for x in range(1,self.N)), 'not a field')
        self.tr = []
        for x in range(self.N):
            y, t = x, 0
            for _ in range(r):
                t ^= y
                y = self.sq[y]
            need(y == x and t in (0,1), 'trace/Frobenius')
            self.tr.append(t)
        self.qpow = [self.frob(x,r//2) for x in range(self.N)]

    def multiply(self, a, b):
        z = 0
        while b:
            if b & 1: z ^= a
            b >>= 1
            a <<= 1
            if a & self.N: a ^= self.poly
        return z

    def frob(self, x, k=1):
        for _ in range(k % self.r): x = self.sq[x]
        return x

    def orbit(self, x):
        out, y = [x], self.sq[x]
        while y != x:
            out.append(y)
            y = self.sq[y]
        return out

    def B(self, x, y):
        return self.tr[self.table[x][self.qpow[y]]]

    def profile(self,x,y):
        return sum(self.B(self.frob(x,j),y) << j for j in range(self.r))

def audit(r, polynomial):
    f = Field(r,polynomial)
    pairs = list(combinations(range(r),2))
    def evaluation(x,y):
        return sum((((x>>i&1)*(y>>j&1)) ^ ((x>>j&1)*(y>>i&1))) << k for k,(i,j) in enumerate(pairs))
    def form_rank(mask):
        return rank([sum((((mask & evaluation(1<<i,1<<j)).bit_count())%2)<<j for j in range(r)) for i in range(r)])
    fr = [evaluation(f.sq[1<<i], f.sq[1<<j]) ^ evaluation(1<<i,1<<j) for i,j in pairs]
    subfield = [x for x in range(f.N) if f.qpow[x] == x]
    need(len(subfield) == 1 << (r//2), 'subfield')
    kb = []
    for x in subfield:
        if rank(kb+[x]) > len(kb): kb.append(x)
    spread = set(tuple(sorted(f.table[x][a] for a in subfield)) for x in range(1,f.N))
    sr = []
    for line in spread:
        x = next(a for a in line if a)
        sr += [evaluation(f.table[x][a],f.table[x][b]) for a,b in combinations(kb,2)]
    spaces = {}
    full = {}
    for name, constraints in [('frobenius',fr),('spread',sr),('joint',fr+sr)]:
        forms = span(kernel(constraints,len(pairs)))
        full[name] = forms
        spaces[name] = {'dimension':len(forms).bit_length()-1,'forms':len(forms), 'nondegenerate':sum(form_rank(x)==r for x in forms)}
    canonical = sum(f.B(1<<i,1<<j) << k for k,(i,j) in enumerate(pairs))
    need(set(full['joint']) == {0,canonical}, 'joint form not unique')
    need(spaces['spread']['forms']==len(subfield), 'spread classification')
    if r <= 6:
        # Independent exhaustive filtering of all alternating forms, including
        # singular forms, rather than just a proposed trace-form family.
        direct = [mask for mask in range(1<<len(pairs)) if all(not ((mask & row).bit_count()%2) for row in fr+sr)]
        need(direct == sorted(full['joint']), 'full ambient enumeration')

    orbits = sorted({tuple(sorted(f.orbit(x))) for x in range(1,f.N)})
    normal = [o for o in orbits if rank(o)==r]
    full_period = [o for o in orbits if len(o)==r]
    profile_cost = Counter()
    example = None
    for ox,oy in product(normal, repeat=2):
        x,y = ox[0],oy[0]
        words = [f.profile(x,f.frob(y,k)) for k in range(r)]
        need(len(set(words))==r, 'normal relative phase collision')
        need(all(w.bit_count()%2==1 for w in words), 'normal correlation parity')
        cost = next(k for k in range(r+1) if any(len({w & mask for w in words})==r for mask in (sum(1<<j for j in js) for js in combinations(range(r),k))))
        profile_cost[cost] += 1
        if cost == r-1 and example is None:
            example = {'x':x,'y':y,'profiles':words}
    need(max(profile_cost)==r-1,'sharp probe cost')
    # Every normal anchor distinguishes *all* field values, not just one orbit.
    for ox in normal:
        need(len({f.profile(ox[0],y) for y in range(f.N)})==f.N,'normal decoder')
    hostile = None
    for ox in full_period:
        x = ox[0]
        if rank(ox)!=r-1: continue
        for oy in full_period:
            for y in oy:
                if f.qpow[y] == y^1 and f.profile(x,y)==f.profile(x,y^1):
                    hostile = {'x':x,'x_orbit':f.orbit(x),'x_span_rank':rank(ox),'y':y,'other_y':y^1,'y_orbit':f.orbit(y),'profile':f.profile(x,y)}
                    break
            if hostile: break
        if hostile: break
    if r in (6,8): need(hostile is not None, 'missing full-period hostile')
    # Three normal named ports; exhaustive r^3 tuples reduce to r^2 edge
    # phase states. The missing third profile is recovered along the path.
    roots = [normal[i%len(normal)][0] for i in range(3)]
    carrier_states = set()
    for phases in product(range(r),repeat=3):
        vals = [f.frob(x,a) for x,a in zip(roots,phases)]
        deltas = ((phases[1]-phases[0])%r,(phases[2]-phases[1])%r)
        carrier_states.add(deltas)
        delta02 = sum(deltas)%r
        need(f.profile(vals[0],vals[2])==f.profile(roots[0],f.frob(roots[2],delta02)), 'LCA transport')
    need(len(carrier_states)==r*r,'minimal state count')
    return {'degree':r,'polynomial_mask':polynomial,'form_spaces':spaces,'canonical_form_mask':canonical,'spread_lines':len(spread),'normal_elements':sum(map(len,normal)),'full_period_elements':sum(map(len,full_period)),'normal_orbit_pairs_minimum_bit_queries':dict(sorted(profile_cost.items())),'sharp_query_witness':example,'full_period_nonnormal_hostile':hostile,'three_port_relative_states':len(carrier_states)}

def tower():
    # H_4 has the nonzero four-bit vectors as its columns.
    free = [j for j in range(1,16) if j not in (1,2,4,8)]
    basis = [(1 << (j-1)) ^ sum(1 << ((1<<k)-1) for k in range(4) if j>>k&1) for j in free]
    words = span(basis)
    need(len(set(words))==2048, 'relation basis')
    need(all(not __import__('functools').reduce(int.__xor__,(j for j in range(1,16) if w>>(j-1)&1),0) for w in words),'Hamming syndrome')
    count = 0
    for a in range(1,2048):
        for b in range(a+1,2048):
            c = a^b
            if b < c:
                need(words[a]^words[b]^words[c]==0,'relation of relations')
                count += 1
    need(count==2047*2046//6,'line count')
    return {'point_counts':[7,15,2047],'ranks':[3,4,11],'next_rank':2036,'next_points_expression':'2^2036 - 1','hamming15_kernel_basis_support_masks':basis,'new_triple_relations':count,'hamming15_weight_enumerator':dict(sorted(Counter(x.bit_count() for x in words).items()))}

if __name__=='__main__':
    result = {'status':'FINITE-EXACT','fields':[audit(r,p) for r,p in [(2,7),(4,19),(6,67),(8,283)]],'relation_tower':tower()}
    result['active_checks'] = CHECKS
    path = Path(__file__).with_name('fields.json')
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
