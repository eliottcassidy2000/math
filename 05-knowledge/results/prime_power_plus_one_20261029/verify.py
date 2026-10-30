"""Exact prime-power-plus-one controls; all truth gates survive python -O."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from functools import lru_cache
import argparse
import hashlib
import importlib.util
import json
import math

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent/'regular_parent_symmetry_20261029'/'verify.py'
spec = importlib.util.spec_from_file_location('regular_parent_helpers', PREVIOUS)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
GATES = 0


def need(value, why):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(why)


def factor(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def poly_remainder(f, g, p):
    a = f[:]
    while len(a) >= len(g):
        lead = a[-1]
        off = len(a)-len(g)
        for j, c in enumerate(g):
            a[off+j] = (a[off+j]-lead*c) % p
        while a and a[-1] == 0:
            a.pop()
    return a


def irreducible(p, k):
    if k == 1:
        return [0, 1]
    for coeff in product(range(p), repeat=k):
        if not coeff[0]:
            continue
        f = list(coeff)+[1]
        if all(poly_remainder(f, list(g)+[1], p)
               for d in range(1, k//2+1) for g in product(range(p), repeat=d)):
            return f
    raise RuntimeError('No irreducible polynomial')


class Field:
    def __init__(self, p, k):
        self.p, self.k, self.q = p, k, p**k
        self.modulus = irreducible(p, k)
        self.digits = [tuple((a//p**j) % p for j in range(k)) for a in range(self.q)]
        self.add = [[sum(((x+y) % p)*p**j for j, (x, y) in enumerate(zip(da, db)))
                     for db in self.digits] for da in self.digits]
        self.neg = [sum((-x % p)*p**j for j, x in enumerate(da)) for da in self.digits]
        self.mul = [[self.raw_mul(a, b) for b in range(self.q)] for a in range(self.q)]

    def raw_mul(self, a, b):
        coeff = [0]*(2*self.k-1)
        for i, x in enumerate(self.digits[a]):
            for j, y in enumerate(self.digits[b]):
                coeff[i+j] = (coeff[i+j]+x*y) % self.p
        rem = poly_remainder(coeff, self.modulus, self.p)
        return sum(x*self.p**j for j, x in enumerate(rem))

    def sub(self, a, b):
        return self.add[a][self.neg[b]]

    def power(self, a, n):
        out = 1
        while n:
            if n & 1:
                out = self.mul[out][a]
            a = self.mul[a][a]
            n //= 2
        return out

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        return self.power(a, self.q-2)

    def div(self, a, b):
        return self.mul[a][self.inv(b)]

    def chi(self, a):
        if not a:
            return 0
        return 1 if self.power(a, (self.q-1)//2) == 1 else -1


def paley(field):
    q = field.q
    need(q % 4 == 3, 'quadratic signs must be antisymmetric')
    return [sum(1 << j for j in range(q) if field.chi(field.sub(j, i)) == 1)
            for i in range(q)]


def transitive_max(rows):
    @lru_cache(None)
    def visit(mask):
        if not mask:
            return 0
        best = 1
        left = mask
        while left:
            bit = left & -left; left ^= bit
            v = bit.bit_length()-1
            tail = mask & rows[v]
            if tail.bit_count()+1 > best:
                best = max(best, 1+visit(tail))
        return best
    return visit((1 << len(rows))-1)


def separating_mask(rows, u, v):
    return sum(1 << w for w in range(len(rows)) if w not in [u, v]
               and ((rows[w] >> u) & 1) != ((rows[w] >> v) & 1))


def deletion_checks(field, supplied_rows=None):
    rows = paley(field) if supplied_rows is None else supplied_rows
    q = len(rows); full = (1 << q)-1; half = (q-1)//2
    need(all(r.bit_count() == half for r in rows), 'Paley regular')
    need(all((rows[u] & rows[v]).bit_count() == (q-3)//4
             for u, v in combinations(range(q), 2)), 'doubly regular premise')
    # Every potential proper module has an outside uniform witness, hence is
    # contained in one out- or in-neighborhood. This is the complete universe.
    candidates = set()
    cones = set()
    for v, out in enumerate(rows):
        inc = full ^ out ^ (1 << v)
        cones.update([out, inc])
        for bank in [out, inc]:
            sub = bank
            while sub:
                if sub.bit_count() >= 2:
                    candidates.add(sub)
                sub = (sub-1) & bank
    minima = {}
    minimum_cuts = {}
    for mask in candidates:
        m = mask.bit_count()
        uniform = sum(1 << v for v, row in enumerate(rows)
                      if not(mask >> v & 1) and (row & mask) in [0, mask])
        need(bool(uniform), 'candidate has outside witness')
        deleted = full ^ mask ^ uniform
        r = deleted.bit_count()
        lower = Fraction(m*(q-m), m+1) if m % 2 else Fraction((m-1)*(q-m)+1, m)
        need(r >= lower >= half, 'parity-sensitive spectral deletion bound')
        minima[m] = min(minima.get(m, q), r)
        if r == half:
            need(m in [2, half], 'sharp modules only pair or half-size')
            minimum_cuts.setdefault(deleted, []).append(mask)
    pair_cuts = {separating_mask(rows, u, v) for u, v in combinations(range(q), 2)}
    need(len(pair_cuts) == q*(q-1)//2, 'distinct pair deletion cuts')
    need(len(cones) == 2*q and not(pair_cuts & cones), 'cone and pair cuts distinct')
    need(set(minimum_cuts) == pair_cuts | cones, 'complete sharp deletion-cut classification')
    need(len(minimum_cuts) == q*(q+3)//2, 'closed sharp cut count')
    need(all(len(mods) == 1 for mods in minimum_cuts.values()), 'unique proper nontrivial module at threshold')
    # Independent direct subset-module audit for small orders.
    direct = None
    if q <= 11:
        size = (q+1)//2
        direct = sum(bool(old.nontrivial_modules(old.induced(rows, S)))
                     for S in combinations(range(q), size))
        need(direct == q*(q+3)//2, 'independent threshold induced-subtournament count')
        need(all(not old.nontrivial_modules(old.induced(rows, S))
                 for S in combinations(range(q), size+1)), 'independent just-above-threshold primality')
    return {'q': q, 'field_modulus': field.modulus if supplied_rows is None else None,
            'constructor': 'Paley' if supplied_rows is None else 'skew-Hadamard doubling',
            'candidate_modules': len(candidates),
            'deletion_minimum_by_module_size': dict(sorted(minima.items())),
            'sharp_deletion_count': len(minimum_cuts), 'pair_cuts': len(pair_cuts),
            'cone_cuts': len(cones), 'direct_decomposable_threshold_count': direct}


def carousel_checks():
    data = []
    for p in [3, 5, 7, 11, 13, 17, 19]:
        base = old.carousel(p)
        group = old.iso_maps(base, base)
        need(len(group) == p, 'native carousel group')
        need(transitive_max(base) == (p+1)//2, 'carousel transitive maximum')
        need(all(old.permute_mask(1, g) != 1 for g in group if g != tuple(range(p))),
             'singleton is distinguishing')
        for k in range(1, 6):
            n = p**k; delta = (p-2)**k
            a = (n-delta)//2
            need(2*a+delta == n, 'parity weight identity')
            data.append({'p': p, 'depth': k, 'child_n': n+1, 'cut_weight': a,
                         'imbalance': delta, 'child_transitive_max': ((p+1)//2)**k+1,
                         'parent_group_exponent': (n-1)//(p-1)})
    mixed = []
    children = []
    for sizes in [(3, 5), (5, 3)]:
        rows, colors = [0], [0]
        order = 1
        for b in reversed(sizes):
            rows = old.lex_product(old.carousel(b), rows)
            colors = [1-c if v == 0 else c for v in range(b) for c in colors]
            order = b*order**b
        actual = old.iso_maps(rows, rows)
        need(len(actual) == order, 'mixed tower full native group')
        child = old.add_vertex(rows, sum(c << i for i, c in enumerate(colors)))
        need(len(old.iso_maps(child, child)) == 1, 'mixed child rigid')
        need(transitive_max(child) == 7, 'mixed child exact transitive maximum')
        children.append(child)
        mixed.append({'outer_to_inner': sizes, 'parent_group_order': order,
                      'parent_n': len(rows), 'cut_weight': sum(colors),
                      'child_transitive_max': transitive_max(child)})
    need(not old.iso_maps(children[0], children[1], True), 'same products do not determine child isotype')
    return {'carousel_families': data, 'mixed_order_hostile': mixed}


def norm_clock(field):
    F = field; q = F.q; N = q+1
    delta = next(a for a in range(1, q) if F.chi(a) == -1)
    def mul(x, y):
        a, b = x; c, d = y
        return (F.add[F.mul[a][c]][F.mul[delta][F.mul[b][d]]],
                F.add[F.mul[a][d]][F.mul[b][c]])
    def power(x, n):
        ans = (1, 0)
        while n:
            if n & 1: ans = mul(ans, x)
            x = mul(x, x); n //= 2
        return ans
    params = [(1, 0)]
    for t in range(q):
        t2 = F.mul[t][t]; den = F.sub(t2, delta)
        params.append((F.div(F.add[t2][delta], den), F.div(F.add[t][t], den)))
    unit_set = {(a, b) for a in range(q) for b in range(q)
                if F.sub(F.mul[a][a], F.mul[delta][F.mul[b][b]]) == 1}
    need(set(params) == unit_set and len(params) == N, 'complete projective norm parametrization')
    primes = list(factor(N))
    z = next(x for x in params if all(power(x, N//ell) != (1, 0) for ell in primes))
    need(power(z, N) == (1, 0), 'primitive norm clock')
    need(next(j for j in range(1, 2*F.k+1) if power(z, F.p**j) == z) == 2*F.k,
         'primitive Frobenius degree is twice exponent')
    trace = lambda x: F.add[x[0]][x[0]]
    tr = trace(z)
    need(next(j for j in range(1, F.k+1) if F.power(tr, F.p**j) == tr) == F.k,
         'trace Frobenius degree is exponent')
    need(len({power(z, 2*j) for j in range(N//2)}) == N//2, 'squaring halves primitive clock')
    for x in params:
        inv = (x[0], F.neg[x[1]])
        need(power(x, q) == inv and mul(x, inv) == (1, 0), 'q-Frobenius is inversion')
        if x[1]:
            need(trace(mul(z, x)) != trace(mul(z, inv)), 'future multiplication detects trace branch')
    # Parameter infinity is represented by None; basis delta remains explicit.
    ts = [None]+list(range(q)); z_by_t = dict(zip(ts, params))
    for t in ts:
        for s in ts:
            u = s if t is None else t if s is None else (
                None if F.add[t][s] == 0 else F.div(F.add[F.mul[t][s]][delta], F.add[t][s]))
            need(mul(z_by_t[t], z_by_t[s]) == z_by_t[u], 'transported projective addition law')
    primitive_t = ts[params.index(z)]
    need(primitive_t is not None, 'nonidentity primitive clock')
    s = primitive_t
    det = F.sub(F.mul[s][s], delta)
    need(F.chi(det) == -1, 'primitive torus generator lies outside determinant-square half')
    perm, scales = [], []
    for x in ts:
        if x is None:
            perm.append(s); scales.append(1)
        else:
            den = F.add[x][s]
            num = F.add[F.mul[s][x]][delta]
            perm.append(None if den == 0 else F.div(num, den))
            scales.append(num if den == 0 else den)
    def wedge(x, y):
        return 1 if x is None else F.neg[1] if y is None else F.sub(x, y)
    for i, x in enumerate(ts):
        for j in range(i+1, N):
            y = ts[j]
            need(F.chi(wedge(perm[i], perm[j])) ==
                 F.chi(det)*F.chi(scales[i])*F.chi(scales[j])*F.chi(wedge(x, y)),
                 'exact determinant switching/converse phase')
    need(F.chi(F.neg[1]) == (-1 if q % 4 == 3 else 1), 'quadratic orientation parity guard')
    return {'p': F.p, 'k': F.k, 'q': q, 'N': N, 'modulus': F.modulus,
            'delta': delta, 'primitive_z': z, 'primitive_parameter': primitive_t,
            'frobenius_degree': 2*F.k, 'trace_degree': F.k, 'square_clock': N//2,
            'quadratic_relation': 'tournament' if q % 4 == 3 else 'symmetric'}


def arithmetic_checks():
    data = []
    for p in [3, 5, 7, 11, 13, 17, 19]:
        for k in range(1, 9):
            N = p**k+1
            factors = factor(N)
            need(math.prod(ell**e for ell, e in factors.items()) == N, 'exact complete factorization')
            need(factors[2] == (factor(p+1)[2] if k % 2 else 1), 'two-adic plus-one law')
            orders = {ell: next(i for i in range(1, 2*k+1) if pow(p, i, ell) == 1)
                      for ell in factors if ell != 2}
            primitive = [ell for ell, order in orders.items() if order == 2*k]
            if k >= 2:
                need(bool(primitive), 'finite Zsigmondy specialization')
            need(all((ell-1) % (2*k) == 0 for ell in primitive), 'primitive clock divisibility')
            data.append({'p': p, 'k': k, 'N': N, 'factors': factors,
                         'prime_orders': orders, 'primitive_divisors': primitive})
        for a in range(1, 13):
            for b in range(1, 13):
                g = math.gcd(a, b)
                expected = p**g+1 if (a//g) % 2 and (b//g) % 2 else 2
                need(math.gcd(p**a+1, p**b+1) == expected, 'plus-one gcd with exponent phase')
    need(factor(2**3+1) == {3: 2} and pow(2, 2, 3) == 1, 'base-two primitive-prime exception')
    return data


def paley27_check(F):
    rows = paley(F)
    group = old.iso_maps(rows, rows)
    need(len(group) == 1053, 'native full Paley-27 group')
    pair_mask = 3
    pair_group = [g for g in group if old.permute_mask(pair_mask, g) == pair_mask]
    need(len(pair_group) == 3, 'two field points retain Frobenius C3')
    rigid = next(S for S in combinations(range(27), 3)
                 if sum(old.permute_mask(sum(1 << u for u in S), g) == sum(1 << u for u in S)
                        for g in group) == 1)
    child_pair = old.add_vertex(rows, pair_mask)
    child_rigid = old.add_vertex(rows, sum(1 << u for u in rigid))
    need(len(old.iso_maps(child_pair, child_pair)) == 3, 'unmarked two-point field attachment group')
    need(len(old.iso_maps(child_rigid, child_rigid)) == 1, 'unmarked three-point field attachment rigid')
    tau = transitive_max(rows)
    need(tau == 5, 'Paley27 exact transitive maximum')
    return {'field_modulus': F.modulus, 'core_group_order': len(group), 'two_point_group': len(pair_group),
            'rigid_three_set': rigid, 'core_transitive_max': tau, 'triangle_tower_same_size_group': 3**13,
            'triangle_tower_same_size_transitive_max': 8}


def local_base_checks(fields):
    result = []
    for key in [(7,1),(11,1),(19,1),(23,1),(3,3),(59,1)]:
        F = fields[key]; rows = paley(F)
        vertices = [v for v in range(F.q) if rows[0] >> v & 1]
        local = old.induced(rows, vertices)
        b = len(local); group = old.iso_maps(local, local)
        need(all(row.bit_count() == (b-1)//2 for row in local), 'local base regular')
        # The proved deletion boundary implies primeness for every q;
        # small bases get an additional exhaustive module check.
        if b <= 13:
            need(not old.nontrivial_modules(local), 'native local base prime')
        tau = transitive_max(local)
        selected = None
        for size in range(1, 4):
            for S in combinations(range(b), size):
                mask = sum(1 << u for u in S)
                if sum(old.permute_mask(mask, g) == mask for g in group) != 1:
                    continue
                complement = old.induced(local, [v for v in range(b) if v not in S])
                if transitive_max(complement) == tau:
                    selected = S; break
            if selected is not None: break
        need(selected is not None, 'constructive rigid cut leaves maximum transitive set')
        child = old.add_vertex(local, sum(1 << u for u in selected))
        need(len(old.iso_maps(child, child)) == 1, 'local-base one-port child rigid')
        need(transitive_max(child) == tau+1, 'local-base child transitive maximum')
        result.append({'ambient_q': F.q, 'base_n': b, 'base_rows': local,
                       'base_group_order': len(group), 'tau': tau, 'cut': selected,
                       'imbalance_base': b-2*len(selected),
                       'family': [{'depth': k, 'child_n': b**k+1,
                                   'cut_weight': (b**k-(b-2*len(selected))**k)//2,
                                   'child_transitive_max': tau**k+1} for k in range(1,5)]})
    return result


def necklace_checks():
    result = []
    for p in [3,5,7,11,13]:
        F, R = 2, 2
        for k in range(1,4):
            F = (F**p+(p-1)*F)//p
            R = (R**p-R)//p
            need(R <= F, 'regular orbit count within all orbits')
            result.append({'p':p,'depth':k,'parent_n':p**k,'subset_orbits':F,
                           'unmarked_children':F-1,'rigid_children':R-(1 if k==1 else 0)})
        if p==3: need(R==2,'ternary binary regular palette fixed point')
        else: need(R>2,'higher prime binary regular palette grows')
    return result


def composite_necklace_checks():
    result=[]
    for b in [9,15]:
        representatives=set(); regular=set()
        full=(1<<b)-1
        for mask in range(1<<b):
            orbit={((mask<<j)|(mask>>(b-j)))&full for j in range(b)}
            rep=min(orbit)
            representatives.add(rep)
            if len(orbit)==b: regular.add(rep)
        divisors=[d for d in range(1,b+1) if b%d==0]
        phi=lambda d:sum(math.gcd(i,d)==1 for i in range(1,d+1))
        def mu(d):
            fac=factor(d)
            return 0 if any(e>1 for e in fac.values()) else (-1)**len(fac)
        need(b*len(representatives)==sum(phi(d)*2**(b//d) for d in divisors),
             'independent composite necklace orbit enumeration')
        need(b*len(regular)==sum(mu(d)*2**(b//d) for d in divisors),
             'independent composite primitive-necklace enumeration')
        result.append({'b':b,'all_orbits':len(representatives),'regular_orbits':len(regular)})
    return result


def guarded_radius(base, forbidden, cap):
    b, q = len(base), len(forbidden)
    # A direct embedding witnesses the lower bound; prime obstructions give the upper.
    target_bank = {s: [old.induced(base, S) for S in combinations(range(b), s)]
                   for s in range(3, b+1)}
    witnesses = []
    for S in combinations(range(q), cap+1):
        obstruction = None
        for size in range(len(S), 2, -1):
            for A in combinations(S, size):
                atom = old.induced(forbidden, A)
                if old.nontrivial_modules(atom): continue
                if not any(old.iso_maps(atom, T, True) for T in target_bank.get(size, [])):
                    obstruction = A; break
            if obstruction is not None: break
        need(obstruction is not None, 'every larger old intersection has a forbidden prime atom')
        witnesses.append({'subset': S, 'prime_obstruction': obstruction})
    found = None
    for S in combinations(range(q), cap):
        A = old.induced(forbidden, S)
        for T in combinations(range(b), cap):
            maps = old.iso_maps(A, old.induced(base, T), True)
            if maps:
                found = (S,T,maps[0]); break
        if found: break
    need(found is not None, 'direct embedding attains attachment capacity')
    return {'base_n': b, 'forbidden_n': q, 'max_old_intersection': cap,
            'minimum_new_vertices': q-cap, 'embedding_witness': found,
            'upper_witnesses': witnesses}


def frobenius_seed_counts(fields, locals_data):
    local13 = next(x['base_rows'] for x in locals_data if x['base_n'] == 13)
    result = []
    for name, rows, ell in [('P7',paley(fields[(7,1)]),3),
                            ('P11',paley(fields[(11,1)]),5),('local13',local13,3)]:
        p = len(rows); group = old.iso_maps(rows, rows)
        c = 1+(p-1)//ell
        cycles = Counter(len(old.cycle_partition(g)) for g in group)
        need(cycles == {p:1,1:p-1,c:p*(ell-1)}, 'Frobenius group cycle inventory')
        regular_masks = sum(sum(old.permute_mask(mask,g)==mask for g in group)==1
                            for mask in range(1<<p))
        free = (2**p-p*2**c+(p-1)*2)//(p*ell)
        need(regular_masks == p*ell*free, 'native binary free-color polynomial')
        for Q in [2,3,4,5]:
            burnside = sum(Q**j*count for j,count in cycles.items())//len(group)
            need(burnside == (Q**p+p*(ell-1)*Q**c+(p-1)*Q)//(p*ell), 'Frobenius Burnside polynomial')
        result.append({'base':name,'p':p,'ell':ell,'fixed_cycle_count':c,
                       'minimum_index':p-c,'binary_regular_orbits':free,
                       'binary_regular_masks':regular_masks,
                       'binary_unmarked_children':sum(2**j*n for j,n in cycles.items())//len(group)-1})
    return result


def lattice_checks(fields):
    from sympy import Matrix, eye, ones, ZZ
    from sympy.matrices.normalforms import smith_normal_form
    result = []
    examples=[paley(fields[key]) for key in [(3,1),(7,1),(11,1),(19,1),(23,1),(3,3),(59,1)]]
    examples.append(doubled_drt(paley(fields[(7,1)])))
    for rows in examples:
        q=len(rows); r=(q-1)//2; C=(q+1)//4
        S=Matrix(q,q,lambda i,j:0 if i==j else 1 if rows[i]>>j&1 else -1)
        need(S*S == ones(q)-q*eye(q), 'skew Hadamard Gram identity')
        numerator=Matrix(q-1,q-1,lambda i,j:int(i==j)+S[i,j]-S[i,q-1])
        need(all(x%2==0 for x in numerator), 'integral root-lattice half operator')
        T=numerator/2; I=eye(q-1); G=I+ones(q-1)
        need(T*(I-T)==C*I, 'complementary integral norm factorization')
        need(T.T*G==G*(I-T), 'adjoint uses the root-lattice Gram matrix')
        need(G.det()==q and math.gcd(q,C)==1, 'Gram form nondegenerate modulo norm')
        diag=smith_normal_form(T,domain=ZZ)
        entries=[abs(int(diag[i,i])) for i in range(q-1)]
        need(entries==[1]*r+[C]*r, 'complete integral Smith classification')
        result.append({'q':q,'C':C,'r':r,'smith_diagonal':entries,
                       'single_cokernel_order':C**r,'paired_group_order':C**(2*r),
                       'lagrangians_if_C_prime':math.prod(C**i+1 for i in range(1,r+1))
                       if len(factor(C))==1 and list(factor(C).values())==[1] else None})
    return result


def doubled_drt(rows):
    q=len(rows); n=q+1
    K=[[0 if i==j else 1 if i==0 else -1 if j==0 else
        1 if rows[i-1]>>(j-1)&1 else -1 for j in range(n)] for i in range(n)]
    H=[[K[i][j]+int(i==j) for j in range(n)] for i in range(n)]
    H2=[H[i]+H[i] for i in range(n)]+[[-H[j][i] for j in range(n)]+[H[j][i] for j in range(n)] for i in range(n)]
    K2=[[H2[i][j]-int(i==j) for j in range(2*n)] for i in range(2*n)]
    signs=[1]+K2[0][1:]
    normalized=[[signs[i]*K2[i][j]*signs[j] for j in range(2*n)] for i in range(2*n)]
    return [sum(1<<(j-1) for j in range(1,2*n) if normalized[i][j]==1) for i in range(1,2*n)]


def lagrangian_checks():
    result=[]
    for p in [2,3,5,7]:
        all_spaces=0; isotropic=0
        for pivots in combinations(range(4),2):
            free=[(i,j) for i,pivot in enumerate(pivots)
                  for j in range(pivot+1,4) if j not in pivots]
            for values in product(range(p),repeat=len(free)):
                rows=[[0]*4 for _ in range(2)]
                for i,j in enumerate(pivots): rows[i][j]=1
                for (i,j),v in zip(free,values): rows[i][j]=v
                x,y=rows
                omega=x[0]*y[2]+x[1]*y[3]-x[2]*y[0]-x[3]*y[1]
                all_spaces+=1
                isotropic+=omega%p==0
        need(all_spaces*(p**2-1)*(p**2-p)==(p**4-1)*(p**4-p),
             'complete RREF universe of two-dimensional subspaces')
        need(isotropic==(p+1)*(p**2+1),'independent Lagrangian recurrence count')
        result.append({'p':p,'dimension':4,'all_planes':all_spaces,'lagrangians':isotropic})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'certificate.json')
    args = parser.parse_args()
    fields = {(p,k): Field(p,k) for p,k in [
        (3,1),(3,2),(3,3),(3,4),(5,1),(5,2),(5,3),(7,1),(7,2),
        (11,1),(11,2),(13,1),(13,2),(19,1),(23,1),(59,1)]}
    locals_data = local_base_checks(fields)
    hostile = old.hostile_checks()
    need([x['regular_binary_masks'] for x in hostile['tournaments']] == [1620,1728], 'retained 81/27 control')
    result = {'deletion': [deletion_checks(fields[key]) for key in [(7,1),(11,1),(19,1),(23,1),(3,3)]],
              'recursive': carousel_checks(),
              'norm_clocks': [norm_clock(F) for F in fields.values()],
              'arithmetic': arithmetic_checks(),
              'field_vs_tower_27': paley27_check(fields[(3,3)]),
              'local_bases': locals_data, 'prime_necklaces': necklace_checks(),
              'composite_necklaces':composite_necklace_checks(),
              'frobenius_seeds': frobenius_seed_counts(fields,locals_data),
              'exact_attachment_radii': [guarded_radius(old.carousel(5),paley(fields[(7,1)]),4),
                                         guarded_radius(paley(fields[(7,1)]),paley(fields[(11,1)]),5)],
              'hostile': hostile, 'integral_lattices': lattice_checks(fields)}
    result['lagrangian_rank_two']=lagrangian_checks()
    drt15=doubled_drt(paley(fields[(7,1)]))
    result['deletion'].append(deletion_checks(None,drt15))
    result['nonfield_15']={'rows':drt15,'group_order':len(old.iso_maps(drt15,drt15)),
                           'transitive_max':transitive_max(drt15)}
    result['gates'] = GATES
    raw = json.dumps(result, indent=2)+'\n'
    args.output.write_text(raw, encoding='utf-8', newline='\n')
    print(json.dumps({'gates': GATES, 'sha256': hashlib.sha256(raw.encode()).hexdigest()}))


if __name__ == '__main__':
    main()
