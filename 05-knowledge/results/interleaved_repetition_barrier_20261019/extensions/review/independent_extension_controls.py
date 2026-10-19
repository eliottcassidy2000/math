#!/usr/bin/env python3
"""Independent exact controls; imports no source checker or proof package."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json

checks = 0

def require(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)


def affine(bits):
    """Closed odd-position sum, evaluated from right to left."""
    h = bits.count('1')
    suffix_power, b = 1, 0
    for i in range(len(bits)-1, -1, -1):
        if bits[i] == '1':
            b += suffix_power << i
            suffix_power *= 3
    return 3**h, b, 1 << len(bits)


def lift(bits):
    """Lift a lawful residue one bit at a time, without modular inversion."""
    residue, source_modulus, image_step, image = 0, 1, 1, 0
    for bit in bits:
        if image % 2 != int(bit):
            residue += source_modulus
            image += image_step
        if bit == '1':
            image = (3*image+1)//2
            image_step *= 3
        else:
            image //= 2
        source_modulus *= 2
    return residue


def step(x):
    return (3*x+1)//2 if x % 2 else x//2


def replay(x, bits):
    for b in bits:
        require(x % 2 == int(b), 'lawful chronological guard')
        x = step(x)
    return x


def log_interval(a, b, terms=96):
    """Exact positive atanh series and geometric remainder majorant."""
    z = Q(a-b, a+b)
    require(0 < z < 1, 'positive atanh argument')
    lo = 2*sum((z**(2*j+1)/Q(2*j+1) for j in range(terms)), Q(0))
    remainder = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return lo, lo+remainder


K = '110110110110110101011011101010'
V, F, H = K*4+'10', '10', '110101000'
a, b, M = affine(V)
require((a, M) == (3**77, 2**122), 'V prime exponents')
v = Q(a, M)
require(1 < v and v**9 < Q(4,3) < v**10, 'greedy reset spacing')
require(Q(1,200)-Q(1,360) == Q(1,450), 'one-percent lower slope')
L0, B = 2**48, 10**10
require(Q(L0,450) > 53*B+174, 'partial-block overlap base')
require(L0 > 450*(B+1), 'overlap increasing beyond base')
require(15*L0 < 2**52, 'short and long ranges overlap')
partial_prefixes = 0
for name, block in [('K',K),('V',V),('F',F),('H',H)]:
    require(len(block) <= 122, 'bank maximum length')
    for q in range(len(block)+1):
        aq, bq, Mq = affine(block[:q])
        require(aq < (1 << 122)*Mq and bq < (1 << 122)*Mq,
                'uniform partial-block affine envelope')
        if name in ('K','V') and q:
            require(aq > Mq, 'every literal expansion prefix expands')
        partial_prefixes += 1
for x in range(1,129):
    y = x
    for q in range(123):
        require(Q(y+1, x+1) <= Q(3,2)**q, 'stepwise augmented-height envelope')
        require(y < (1 << 122)*(x+1), '2^122 envelope on actual state')
        y = step(y)

ln2lo, ln2hi = log_interval(2,1)
ln32lo, ln32hi = log_interval(3,2)
ln3lo, ln3hi = ln2lo+ln32lo, ln2hi+ln32hi
numlo, numhi = 77*ln3lo-122*ln2hi, 77*ln3hi-122*ln2lo
denlo, denhi = 2*ln2lo-ln3hi, 2*ln2hi-ln3lo
require(0 < numlo < numhi and 0 < denlo < denhi, 'theta quotient is positive')
thlo, thhi = numlo/denhi, numhi/denlo
require(Q(1,10) < thlo < thhi < Q(1,9), 'theta run-length interval')

def floor_theta(n):
    lo, hi = (n*thlo, n*thhi) if n >= 0 else (n*thhi, n*thlo)
    p, q = lo.numerator//lo.denominator, hi.numerator//hi.denominator
    require(p == q, 'certified integer part of n theta')
    return p

max_m = 2048
floors = [floor_theta(n) for n in range(max_m+1)]
mechanical = ''.join(str(floors[n]-floors[n-1]) for n in range(1,max_m+1))
A, D, conditional = 1, 1, 0
for n in range(1,max_m+1):
    A *= 3**77
    D <<= 122
    reset = 3*A >= 4*D
    require(int(reset) == int(mechanical[n-1]), 'exact rational greedy equals mechanical coding')
    if reset:
        A *= 3
        D <<= 2
        conditional += 1
    require(conditional == floors[n], 'reset count equals floor n theta')
    require(D <= A and 3*A < 4*D, 'greedy interval invariant')

partition_cases = 0
for k in range(1,41):
    points = [(Q(0), Q(0), 0)]
    for j in range(1,k+1):
        ceil_j = floor_theta(j)+1
        points.append((ceil_j-j*thhi, ceil_j-j*thlo, j))
    points.sort(key=lambda p:p[0])
    patterns = set()
    for t, (low, high, j) in enumerate(points):
        right = points[t+1][0] if t+1 < len(points) else Q(1)
        require(high < right, 'rotation breakpoints have certified ordering')
        phase = (high+right)/2
        fs = []
        for q in range(k+1):
            f1 = (phase+q*thlo).__floor__()
            f2 = (phase+q*thhi).__floor__()
            require(f1 == f2, 'interval interior has certified coding')
            fs.append(f1)
        interior = ''.join(str(fs[q+1]-fs[q]) for q in range(k))
        # At z={-j theta}, integral constants cancel in the successive differences.
        edge_floors = [floor_theta(q-j) for q in range(k+1)]
        endpoint = ''.join(str(edge_floors[q+1]-edge_floors[q]) for q in range(k))
        require(endpoint == interior, 'floor convention selects right interval')
        patterns.add(interior)
        partition_cases += 1
    seen = {mechanical[i:i+k] for i in range(len(mechanical)-k+1)}
    require(len(patterns) <= k+1 and seen <= patterns, 'rotation partition bounds all observed factors')

endpoint_words = 0
cylinder_pairs = 0
for n in range(1,10):
    for tup in product('01',repeat=n):
        bits = ''.join(tup)
        a, b, M = affine(bits)
        residue = lift(bits)
        require((a*residue+b) % M == 0, 'lift equals exact affine residue')
        x = residue if residue else M
        x2 = x+M
        y, y2 = replay(x,bits), replay(x2,bits)
        require(M*y == a*x+b and y2-y == a, 'closed affine and cylinder gap')
        require((x2-x) % (1 << n) == 0, 'matching prefix gap is power of two')
        cylinder_pairs += 1
        if a < M:
            delta = M-a
            least = b*pow(M,-1,delta) % delta if delta > 1 else 0
            hits = []
            for displacement in range(b//M+2):
                divisibility = b > M*displacement and (b-M*displacement) % delta == 0
                residue_condition = displacement >= least and (displacement-least) % delta == 0 and M*displacement < b
                require(divisibility == residue_condition, 'endpoint and residue filters equivalent')
                if divisibility:
                    source = (b-M*displacement)//delta
                    require(replay(source,bits) == source+displacement, 'endpoint equality authenticates all guards')
                    hits.append(displacement)
            require(bool(hits) == (M*least < b), 'least residue decides nonnegative endpoint existence')
            endpoint_words += 1

pigeonhole_cases = 0
for n in range(1,10):
    for tup in product('01',repeat=n):
        bits = ''.join(tup)
        for r in range(1,n+1):
            factors = [bits[i:i+r] for i in range(n-r+1)]
            P = len(set(factors))
            if P < n-r+1:
                seen = {}
                pair = None
                for j in range(P+1):
                    if factors[j] in seen:
                        pair = (seen[factors[j]],j)
                        break
                    seen[factors[j]] = j
                require(pair is not None, 'first P+1 starts contain duplicate')
                i,j = pair
                require(i <= P-1 and 1 <= j-i <= P, 'cycle-entry and period bounds')
                pigeonhole_cases += 1
require(len({'000','001'}) == len('0001')-3+1, 'strict threshold equality can have no repeated factor')

matched_pairs = 0
for source in range(1,65):
    orbit = [source]
    for _ in range(96):
        orbit.append(step(orbit[-1]))
    parity = ''.join(str(x%2) for x in orbit[:-1])
    for r in (1,2,4,8,12,24):
        first = {}
        for j in range(len(parity)-r+1):
            word = parity[j:j+r]
            if word in first:
                i = first[word]
                require((orbit[j]-orbit[i]) % (1 << r) == 0, 'equal factors give exact dyadic divisibility')
                if orbit[j] != orbit[i]:
                    require(abs(orbit[j]-orbit[i]) >= (1 << r), 'distinct matched-state gap')
                else:
                    overlap = len(parity)-j
                    require(parity[i:i+overlap] == parity[j:j+overlap], 'return forces every available aligned continuation')
                matched_pairs += 1
            else:
                first[word] = j

# Scope hostiles: a contracting exact return is real, and a preperiod need not vanish.
require(replay(1,'10'*24) == 1, 'genuine contracting exact-return hostile')
require(affine('10') == (3,1,4), 'contracting return separator')
trajectory = [3]
for _ in range(20): trajectory.append(step(trajectory[-1]))
require(trajectory[4] == trajectory[6] == 2 and 3 not in trajectory[1:], 'internal cycle does not imply source periodic')
require(affine('10')[0::2] == affine('01')[0::2], 'equal slopes hostile')
require(replay(1,'10') == 1 and replay(2,'01') == 2 and (2-1) % 4 != 0, 'equal slopes cannot replace equal parity factors')

records = []
for m in (1,5,10,11,100,901,916,1024,2048):
    reset = mechanical[:m]
    bits = ''.join(V+(F if c == '1' else '') for c in reset)+F
    n, L = len(bits), 122*m
    a,b,M = affine(bits)
    require(a < M, 'greedy final crossing')
    require(6*b < n*M and 6*b < n*a, 'greedy strict carry tariffs')
    Delta = M-a
    d0 = b*pow(M,-1,Delta) % Delta if Delta > 1 else 0
    require(M*d0 >= b, 'finite greedy controls have no nonnegative displacement')
    r = (L+99)//100
    lengths = sorted({1,2,3,7,121,122,123,124,125,r,n,n+1})
    counts = {}
    for size in lengths:
        P = len({bits[i:i+size] for i in range(max(0,n-size+1))})
        bound = Q(62,61)*size+498
        require(P <= bound, 'literal variable-image complexity bound')
        if size == r:
            counts = {'r':r,'P':P}
            if m >= 5: require(P < n-r+1, 'one-percent pigeonhole')
    if m >= 5:
        require(Q(123,61)*r+498 < n+1, 'all-m numerical threshold control')
    runs = []
    run = 0
    for c in reset:
        run += 1
        if c == '1': runs.append(run);run = 0
    if run: runs.append(run)
    require(max(runs) <= 10, 'bounded expansion runs')
    require(all(q in (9,10) for q in runs[:-1]), 'completed greedy runs nine or ten')
    records.append({'m':m,'N':n,'L':L,'conditional_resets':reset.count('1'),
                    **counts,'max_run':max(runs),'literal_sha256':sha256(bits.encode()).hexdigest()})

require(Q(123*5,50)+Q(123,61)+498 < 122*5+1, 'complexity threshold m=5')
require(Q(123,50) < 122, 'complexity margin grows with m')
require(100*1098 < 120*916, 'parse-robust one-percent threshold')
# Keep first V^10; split the remaining 891 V letters in W_901.
require(mechanical[:10] == '0000000001', 'initial ten-V expansion run')
mixed_L = 10*122+891*120
require(mixed_L == 108140 and 100*1098 >= mixed_L, 'm>900 is not arbitrary-parse robust')

print(json.dumps({'status':'PASS: exact independent finite controls; universal claims audited separately',
                  'checks':checks,'partial_bank_prefixes':partial_prefixes,
                  'rotation_partition_intervals':partition_cases,'mechanical_terms':max_m,
                  'dyadic_cylinder_pairs':cylinder_pairs,'contracting_endpoint_words':endpoint_words,
                  'strict_pigeonhole_cases':pigeonhole_cases,'trajectory_matched_pairs':matched_pairs,
                  'parse_threshold_hostile':{'m':901,'r':1098,'L':mixed_L},
                  'greedy_controls':records},indent=2,sort_keys=True))
