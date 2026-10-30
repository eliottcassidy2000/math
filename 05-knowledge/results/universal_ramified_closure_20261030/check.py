"""Literal cyclotomic colors, universal-group upper bounds and selector counts."""
from collections import Counter
from fractions import Fraction
from itertools import product, combinations
from math import comb, factorial
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import random

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


class Ring:
    def __init__(self, p, k, d=1):
        self.p, self.k, self.q, self.D = p, k, p**k, p*d
        self.n = self.q**2
        self.points = list(product(range(self.q), repeat=2))
        inv = pow(1-self.D, -1, self.q)
        self.g = ((1+self.D)*inv % self.q, 2*inv % self.q)
        self.powers = [(1, 0)]
        for _ in range(self.q-1):
            self.powers.append(self.mul(self.powers[-1], self.g))
        need(self.mul(self.powers[-1], self.g) == (1, 0)
             and len(set(self.powers)) == self.q, 'full principal norm clock')
        z = self.powers[p**(k-1)]
        self.eta = ((z[0]-1) % self.q, z[1])
        self.owner = [-1]*self.n
        self.orbits = []
        for i, v in enumerate(self.points):
            if self.owner[i] >= 0:
                continue
            orb, w = [], v
            while self.owner[self.idx(w)] < 0:
                self.owner[self.idx(w)] = len(self.orbits)
                orb.append(self.idx(w))
                w = self.mul(w, self.g)
            need(w == v, 'literal orbit closes')
            self.orbits.append(orb)
        need(len(self.orbits) == 2*self.q-1, 'all directed orbit colors')

    def idx(self, v):
        return v[0]*self.q+v[1]

    def mul(self, v, w):
        return ((v[0]*w[0]+self.D*v[1]*w[1]) % self.q,
                (v[0]*w[1]+v[1]*w[0]) % self.q)

    def color(self, i, j):
        x, y = self.points[i], self.points[j]
        return self.owner[((y[0]-x[0]) % self.q)*self.q+(y[1]-x[1]) % self.q]

    def valuation(self, v):
        def vp(a):
            if not a:
                return self.k
            z = 0
            while a % self.p == 0:
                z, a = z+1, a//self.p
            return z
        return min(2*vp(v[0]), 2*vp(v[1])+1, 2*self.k)

    def affine_fibre(self, t, b, h):
        ans = []
        for v in self.points:
            w = self.mul(self.powers[t], v)
            z = h[v[0] % self.p]
            ans.append(self.idx(tuple((w[i]+b[i]+self.eta[i]*z) % self.q for i in range(2))))
        return ans


def colored_group_certificate(ring):
    # Any automorphism fixing 0 and 1 must preserve these two color labels.
    p, n, q = ring.p, ring.n, ring.q
    landmarks = [j*q for j in range(p)]
    colors = [[ring.color(i, j) for j in range(n)] for i in range(n)]
    domains = [[v for v in range(n) if (colors[0][v], colors[q][v]) ==
                (colors[0][a], colors[q][a])] for a in landmarks[2:]]
    profiles = list(map(len, domains))
    need(profiles == [p]*(p-2), 'two-point domains have exactly p candidates')
    extras = [j*q+1 for j in range(p)]
    full_landmarks = landmarks+extras
    original_signatures = [tuple(colors[a][v] for a in full_landmarks) for v in range(n)]
    need(len(set(original_signatures)) == n, 'two-stage landmark colors separate every vertex')
    accepted, rejected, maps_digest = 0, 0, hashlib.sha256()
    for choices in product(*domains):
        images = [0, q, *choices]
        for a in extras:
            signature = tuple(colors[b][a] for b in landmarks)
            candidates = [v for v in range(n) if tuple(colors[b][v] for b in images[:p]) == signature]
            need(len(candidates) == 1, 'auxiliary landmark reconstructed without another selector')
            images.append(candidates[0])
        reverse = {tuple(colors[a][v] for a in images): v for v in range(n)}
        if len(reverse) != n or any(sig not in reverse for sig in original_signatures):
            rejected += 1
            continue
        perm = [reverse[sig] for sig in original_signatures]
        need(len(set(perm)) == n, 'candidate permutation bijective')
        valid = all(colors[x][y] == colors[perm[x]][perm[y]]
                    for x in range(n) for y in range(x+1, n))
        if not valid:
            rejected += 1
            continue
        h = {}
        for i, j in enumerate(perm):
            v, w = ring.points[i], ring.points[j]
            shift = ((w[0]-v[0]) % q, (w[1]-v[1]) % q)
            possible = [s for s in range(p) if shift == tuple(s*a % q for a in ring.eta)]
            need(len(possible) == 1, 'map lies on exact socle line')
            residue = v[0] % p
            need(residue not in h or h[residue] == possible[0], 'shift constant on entire residue fibre')
            h[residue] = possible[0]
        need(h[0] == h[1] == 0, 'root and unit selector retained')
        accepted += 1
        maps_digest.update(json.dumps(perm, separators=(',', ':')).encode())
    need(accepted == p**(p-2) and rejected == 0, 'exhaustive full two-point stabilizer')
    return {'p': p, 'k': ring.k, 'D': ring.D, 'landmark_domains': profiles,
            'accepted': accepted, 'rejected': rejected,
            'full_universal_order': n*ring.q*accepted,
            'all_maps_sha256': maps_digest.hexdigest()}


def predicted_stabilizer(ring, pins):
    if not pins:
        return ring.p**(3*ring.k+ring.p-2)
    j = len({ring.points[i][0] % ring.p for i in pins})
    if j >= 2:
        return ring.p**(ring.p-j)
    root = ring.points[pins[0]]
    tau = min(ring.valuation(tuple((ring.points[i][z]-root[z]) % ring.q for z in range(2))) for i in pins)
    return ring.p**(ring.p-2+(tau+1)//2)


def direct_stabilizer(ring, pins):
    if not pins:
        return ring.p**(3*ring.k+ring.p-2)
    root = ring.points[pins[0]]
    vecs = [tuple((ring.points[i][z]-root[z]) % ring.q for z in range(2)) for i in pins]
    shifts = {tuple(s*a % ring.q for a in ring.eta): s for s in range(ring.p)}
    total = 0
    for g in ring.powers:
        forced = {0: 0}
        for v in vecs:
            w = ring.mul(g, v)
            delta = tuple((v[z]-w[z]) % ring.q for z in range(2))
            s = shifts.get(delta)
            z = v[0] % ring.p
            if s is None or (z in forced and forced[z] != s):
                break
            forced[z] = s
        else:
            total += ring.p**(ring.p-len(forced))
    need(total % ring.p == 0, 'top clock duplication removed integrally')
    return total//ring.p


def selector_checks():
    rng, reports = random.Random(301031), []
    for p, k, d in ((5, 1, 1), (5, 2, 1), (5, 3, 1), (7, 1, 1), (7, 2, 1), (11, 1, 2)):
        ring = Ring(p, k, d)
        hist = Counter()
        for v in range(ring.n):
            exact = direct_stabilizer(ring, [0, v])
            need(exact == predicted_stabilizer(ring, [0, v]), 'every anchored pair stabilizer')
            hist[exact] += 1
        for _ in range(120):
            pins = rng.sample(range(ring.n), rng.randrange(1, min(10, ring.n)+1))
            need(direct_stabilizer(ring, pins) == predicted_stabilizer(ring, pins), 'arbitrary joint-pin law')
        for j in range(1, p+1):
            pins = [z*ring.q for z in range(j)]
            need(direct_stabilizer(ring, pins) == predicted_stabilizer(ring, pins), 'transversal stabilizer staircase')
        reports.append({'p': p, 'k': k, 'D': ring.D, 'complete_anchored_pairs': ring.n,
                        'random_joint_sets': 120, 'pair_stabilizer_histogram': sorted(hist.items())})
    return reports


def full_small_group():
    ring = Ring(5, 1)
    maps, fixed_masks, cycle_hist = set(), Counter(), Counter()
    for t, b, values in product(range(5), ring.points, product(range(5), repeat=3)):
        perm = tuple(ring.affine_fibre(t, b, (0, 0)+values))
        maps.add(perm)
        fixed_masks[sum(1 << i for i, j in enumerate(perm) if i == j)] += 1
        seen, cycles = set(), Counter()
        for v in range(25):
            if v in seen:
                continue
            x, size = v, 0
            while x not in seen:
                seen.add(x)
                x, size = perm[x], size+1
            need(x == v, 'complete small permutation cycles')
            cycles[size] += 1
        cycle_hist[tuple(sorted(cycles.items()))] += 1
    need(len(maps) == 5**6, 'complete small group has no duplicated normal forms')
    for size in (1, 2, 3):
        for pins in combinations(range(25), size):
            mask = sum(1 << i for i in pins)
            count = sum(mass for fixed, mass in fixed_masks.items() if fixed & mask == mask)
            need(count == predicted_stabilizer(ring, pins), 'all small singleton/pair/triple supports')
    hist = Counter()
    for mask, mass in fixed_masks.items():
        hist[mask.bit_count()] += mass
    return {'group_elements': len(maps), 'all_pin_sets_of_sizes_1_2_3': sum(comb(25,j) for j in (1,2,3)),
            'fixed_point_histogram': sorted(hist.items()),
            'complete_cycle_inventory': [{'cycles': list(t), 'mass': m} for t, m in sorted(cycle_hist.items())]}


def stirling(n, j):
    return sum((-1)**(j-i)*comb(j,i)*i**n for i in range(j+1))//factorial(j)


def tuple_burnside(p, k, r):
    N = p**(2*k)
    numerator = p**(p-2)*sum((p-1)*p**(3*m-2)*p**((2*k-2*m+1)*r) for m in range(2,k+1))
    numerator += (N//p)**r*sum(comb(p,z)*(p-1)**(p-z)*z**r for z in range(p+1))
    value = Fraction(numerator, p**(3*k+p-2))
    need(value.denominator == 1, 'tuple Burnside integrality')
    return int(value)


def tuple_strata(p, k, r):
    N, order = p**(2*k), p**(3*k+p-2)
    total = Fraction(1)
    for t in range(1, 2*k):
        mass = N*(p**((2*k-t)*(r-1))-p**((2*k-t-1)*(r-1)))
        total += Fraction(mass*p**(p-2+(t+1)//2), order)
    for j in range(2, min(p,r)+1):
        mass = factorial(p)//factorial(p-j)*stirling(r,j)*(N//p)**r
        total += Fraction(mass*p**(p-j), order)
    need(total.denominator == 1, 'configuration stratum orbit integrality')
    return int(total)


def counting_checks():
    reports = []
    for p, k in product((5,7,11,13), range(1,5)):
        roots = [z*p**(2*k-1) for z in range(1,p+1)]+[p**(2*k-2*m+1) for m in range(2,k+1)]
        poly = [1]
        for root in roots:
            nxt = [0]*(len(poly)+1)
            for i, coef in enumerate(poly):
                nxt[i] -= root*coef
                nxt[i+1] += coef
            poly = nxt
        values = [tuple_burnside(p,k,r) for r in range(1,2*len(roots)+5)]
        for r, value in enumerate(values, 1):
            need(value == tuple_strata(p,k,r), 'independent all-configuration orbit count')
        for start in range(len(values)-len(roots)):
            need(sum(c*values[start+i] for i,c in enumerate(poly)) == 0, 'exact joint-selector recurrence')
        need(values[:2] == [1, 2*p**k-1], 'transitivity and directed-color orbit controls')
        reports.append({'p':p,'k':k,'minimal_recurrence_order':len(roots),
                        'characteristic_coefficients_ascending':list(map(str,poly)),
                        'first_eight_ordered_tuple_orbits':list(map(str,values[:8]))})
    return reports


def bank_kernel_checks():
    reports=[]
    for Q in (2,3,4,8):
        rows=[]
        for j in range(3):
            mass=Q**4*comb(2,j)*(Q**2)**j*(Q**10-Q**2)**(2-j)
            need(mass % 5**(2-j) == 0,'bank stratum orbit integrality')
            rows.append({'fixed_blocks':j,'effective_stabilizer':5**j,'labelled_states':str(mass),
                         'orbits':str(mass//5**(2-j)),'lifted_vertex_subgroup_order':5**(9+j)})
        need(sum(int(row['labelled_states']) for row in rows)==Q**24,'complete bank mass')
        orbit_sum=sum(int(row['orbits']) for row in rows)
        need(25*orbit_sum==Q**24+8*Q**16+16*Q**8,'independent inherited bank Burnside polynomial')
        reports.append({'Q':Q,'rows':rows,'orbits':str(orbit_sum)})
    return reports


def hostiles_and_packets():
    reports = []
    for p,k in ((5,1),(5,2),(5,3),(7,1),(7,2),(11,1)):
        ring = Ring(p,k)
        M = ring.q-1
        bits = 1+(M-1).bit_length()
        codes, used = {}, set()
        for color, orbit in enumerate(ring.orbits):
            if orbit == [0] or color in codes:
                continue
            x,y = ring.points[orbit[0]]
            opp = ring.owner[((-x)%ring.q)*ring.q+(-y)%ring.q]
            label = len(codes)//2
            codes[color] = label
            codes[opp] = ((1<<bits)-1)^label
            used.add(label)
        need(len(codes) == 2*M and len(set(codes.values())) == 2*M, 'injective complementary tournament code')
        need(2**(bits-1) < 2*M <= 2**bits, 'minimal bit count for color-faithful coding')
        landmarks = [j*ring.q for j in range(p)]
        extras = [j*ring.q+1 for j in range(p)]
        signatures = [tuple(ring.color(a,v) for a in landmarks+extras) for v in range(ring.n)]
        need(len(set(signatures)) == ring.n, 'two-stage resolution throughout larger universes')
        reports.append({'p':p,'k':k,'paired_orbits':M,'named_tournament_layers':bits,
                        'two_stage_resolved_vertices':ring.n})
    ring=Ring(5,2)
    marks=[j*25 for j in range(5)]
    collision=[ring.idx((5,y)) for y in range(0,25,5)]
    need(len({tuple(ring.color(a,v) for a in marks) for v in collision}) == 1, 'base is not a direct resolving set')
    need(predicted_stabilizer(ring,marks)==1, 'collision persists despite trivial base stabilizer')
    # Excluded prime: the Cayley clock modulo9 with delta^2=6 has order3.
    g=(4,5)
    mul=lambda a,b: ((a[0]*b[0]+6*a[1]*b[1])%9,(a[0]*b[1]+a[1]*b[0])%9)
    need(mul(mul(g,g),g)==(1,0) and g!=(1,0),'p3 defective clock hostile')
    spec=importlib.util.spec_from_file_location('native_parent',HERE.parent/'regular_parent_symmetry_20261029/verify.py')
    native=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    control=native.hostile_checks()
    return {'packets':reports,'direct_signature_collision':[[5,y] for y in range(0,25,5)],
            'p3_hostile':{'p':3,'k':2,'D':6,'clock':g,'order':3},
            'control81_27':control,'inherited_gates':native.GATES}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'receipt.json')
    args = parser.parse_args()
    reports = [colored_group_certificate(Ring(p, k)) for p, k in ((5, 1), (5, 2), (7, 1))]
    result = {'status': 'FINITE-EXACT stated colored universes, not every bank group',
              'colored_groups': reports,'selectors':selector_checks(),
              'full_small_group':full_small_group(),'counts':counting_checks(),
              'bank_kernel':bank_kernel_checks(),
              'hostiles':hostiles_and_packets(),'gates': GATES}
    data = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    args.output.write_bytes(data)
    print(json.dumps({'gates': GATES, 'sha256': hashlib.sha256(data).hexdigest(), 'groups': reports}))


if __name__ == '__main__':
    main()
