"""Exact 105 bridges, with independent source maps and hostile controls.

No large open conjecture is assumed or decided. All gates run under -O.
Universe: all 105 matchings; all phases mod105; eight reduced forms;
permutation cycle types through degree16; all parallel matchings at ranks2--5;
the displayed arithmetic and weighted sign controls.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from math import comb, gcd, isqrt, lcm
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def phi_mul(x, y, modulus=None):
    a, b = x
    c, d = y
    z = (a*c+b*d, a*d+b*c+b*d)
    return tuple(t % modulus for t in z) if modulus else z


def phi_pow(n, modulus=None):
    value, base = (1, 0), (0, 1)
    while n:
        if n & 1:
            value = phi_mul(value, base, modulus)
        base = phi_mul(base, base, modulus)
        n //= 2
    return value


def matchings(labels):
    if not labels:
        yield ()
        return
    a = labels[0]
    for b in labels[1:]:
        rest = tuple(x for x in labels if x not in (a, b))
        for m in matchings(rest):
            yield tuple(sorted(((a, b),) + m))


def cycles(p):
    unseen, result = set(p), []
    while unseen:
        start = min(unseen)
        row, x = [], start
        while x not in row:
            row.append(x)
            unseen.remove(x)
            x = p[x]
        result.append(tuple(row))
    return result


def partition(n, lo=1):
    if not n:
        yield ()
    for a in range(lo, n+1):
        for tail in partition(n-a, a):
            yield (a,) + tail


def induced_orbits(objects, permutation):
    left, result = set(objects), []
    while left:
        first = min(left)
        x, orbit = first, []
        while x not in orbit:
            orbit.append(x)
            left.remove(x)
            x = tuple(sorted(permutation[a] for a in x))
        check(x == first, "orbit closes at source")
        result.append(orbit)
    return result


def arithmetic():
    # Independent exact evaluation of the inherited cyclotomic coefficient word.
    from sympy import Poly, Symbol, cyclotomic_poly
    x = Symbol("x")
    coeff = [int(a) for a in Poly(cyclotomic_poly(105, x), x).all_coeffs()]
    value = (0, 0)
    for a in coeff:
        value = phi_mul(value, (0, 1))
        value = (value[0] + a, value[1])
    A, B = value
    check(value == (5224431949, 8453308464), "cyclotomic evaluation")
    N = A*A+A*B-B*B
    check(N == 21211*767131 == 16271615641, "norm")
    check(gcd(B, N) == gcd(5, N) == 1, "quotient and conjugate guards")
    r = -A*pow(B, -1, N) % N
    check((r*r-r-1) % N == 0, "golden polynomial")
    check((A+B*r) % N == 0, "ideal killed")
    check(pow(r, 105, N) == 1, "period105")
    for p in (3, 5, 7):
        check(pow(r, 105//p, N) != 1, "exact105")
    rb = (1-r) % N
    check(pow(rb, 105, N) == N-1, "conjugate sign")
    check(pow(rb, 210, N) == 1, "conjugate210")
    for p in (2, 3, 5, 7):
        check(pow(rb, 210//p, N) != 1, "exact210")
    check(gcd(r-rb, N) == 1, "comaximal conjugate factors")
    check(phi_pow(80, 105) == (1, 0), "scalar golden80")
    for p in (2, 5):
        check(phi_pow(80//p, 105) != (1, 0), "exact scalar80")
    check(pow(2, 12, 105) == 1, "binary12")
    check(pow(2, 80, 105) == 46, "missing binary phase")
    primitive_periods = Counter()
    units = 0
    for a, b in product(range(105), repeat=2):
        units += gcd(a*a+a*b-b*b, 105) == 1
        if gcd(gcd(a, b), 105) != 1:
            continue
        cur = (a, b)
        for k in range(1, 81):
            cur = (cur[1], (cur[0]+cur[1]) % 105)
            if cur == (a, b):
                primitive_periods[k] += 1
                break
        else:
            raise RuntimeError("phase period exceeds80")
    check(units == 7680, "units")
    check(dict(primitive_periods) == {16: 1536, 80: 7680}, "primitive phases")
    return dict(alpha=[A, B], norm=N, root=r, conjugate_root=rb,
                quotient_order=105, conjugate_order=210,
                scalar_modulus_order=80, binary_order=12, joint_order=240,
                primitive_periods=dict(primitive_periods), units=units)


def class_characters():
    D, forms = -420, []
    for a in range(1, isqrt(140)+1):
        for b in range(-a, a+1):
            if (b*b-D) % (4*a):
                continue
            c = (b*b-D)//(4*a)
            if a > c or gcd(gcd(a,b),c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a,b,c))
    expected = [(1,0,105),(2,2,53),(3,0,35),(5,0,21),
                (6,6,19),(7,0,15),(10,10,13),(11,8,11)]
    check(forms == expected, "complete reduced forms")
    check(all(b == 0 or abs(b) == a or a == c for a,b,c in forms),
          "all forms ambiguous")
    primes = [109,53,47,41,19,43,13,11]
    vectors = [(-2,-1),(-1,1),(-2,-1),(-2,-1),
               (-1,1),(-2,-1),(-1,1),(-1,0)]
    signs = []
    for (a,b,c), p, (u,v) in zip(forms, primes, vectors):
        check(a*u*u+b*u*v+c*v*v == p, "represented prime")
        check(all(p % d for d in range(2,isqrt(p)+1)), "prime")
        bits = tuple(int(pow(d,(p-1)//2,p) == p-1) for d in (3,5,7))
        signs.append(bits)
    check(len(set(signs)) == 8, "genus signs separate eight classes")
    radical = [1,3,5,15,7,21,35,105]
    masks = [sum(((d & g).bit_count() % 2) << g for g in range(8))
             for d in range(8)]
    check(masks == [0,170,204,102,240,90,60,150], "character masks")
    H8 = {s ^ f for s in (0,255) for f in masks}
    C = {h ^ t for h in H8 for t in (0,768)}
    check(Counter(h.bit_count() for h in H8) == {0:1,4:14,8:1}, "H8 weights")
    golden_pairs = sorted(tuple(sorted((radical[d], radical[d^2])))
                          for d in range(8) if d < (d^2))
    check(golden_pairs == [(1,5),(3,15),(7,35),(21,105)], "golden atlas")
    return C, masks, dict(forms=forms, prime_witnesses=primes,
                          genus_signs=signs, characters=radical,
                          masks=masks, golden_pairs=golden_pairs)


def matching_dictionary(C, masks):
    ms = list(matchings(tuple(range(8))))
    check(len(ms) == len(set(ms)) == 105, "all matchings")
    E = (0,255,768,1023)
    fibres = {a: {masks[a]^e for e in E} for a in range(8)}
    atlases = set()
    by_direction = Counter()
    cross = Counter()
    passports = Counter()
    ordinary_cross = Counter()
    records = []
    # A Singer7-cycle is linear on F2^3, fixing0: x -> x*X in F8.
    singer = {x: ((x << 1) ^ (11 if x & 4 else 0)) for x in range(8)}
    check(sorted(singer.values()) == list(range(8)), "Singer permutation")
    check(sorted(map(len,cycles(singer))) == [1,7], "Singer7")
    for m in ms:
        atlas = tuple(sorted(tuple(sorted(fibres[a] | fibres[b])) for a,b in m))
        check(set().union(*(set(cell) for cell in atlas)) == C, "atlas covers")
        check(sum(map(len,atlas)) == len(C), "atlas disjoint")
        atlases.add(atlas)
        nd = len({a ^ b for a,b in m})
        by_direction[nd] += 1
        root_pair = next(p for p in m if 0 in p)
        unmatched = root_pair[1]
        partial = tuple(p for p in m if 0 not in p)
        check(tuple(sorted(partial + ((0,unmatched),))) == m, "puncture inverse")
        involution = {x:x for x in range(1,8)}
        for a,b in partial:
            involution[a],involution[b] = b,a
        beta = {x:involution[singer[x]] for x in range(1,8)}
        passport = tuple(sorted(map(len,cycles(beta)), reverse=True))
        genus = (4-len(passport))//2
        check(genus in (0,1) and 4-len(passport) == 2*genus, "Euler genus")
        # Independent genus-zero test: no two chords interlace in the circle.
        order = cycles(singer)[1]
        position = {x:i for i,x in enumerate(order)}
        chords = [sorted((position[a],position[b])) for a,b in partial]
        noncrossing = not any(
            a<c<b<d or c<a<d<b for (a,b),(c,d) in combinations(chords,2)
        )
        check(noncrossing == (genus == 0), "chord/genus independent agreement")
        cross[(nd,genus)] += 1
        if genus == 0:
            passports[passport] += 1
        records.append(dict(matching=m, directions=nd, genus=genus, passport=passport))
        ordinary = {x:involution[x%7+1] for x in range(1,8)}
        ordinary_cross[(nd,(4-len(cycles(ordinary)))//2)] += 1
    check(len(atlases) == 105, "atlas map injective")
    check(dict(by_direction) == {1:7,2:42,4:56}, "direction histogram")
    check(sum(v for (d,g),v in cross.items() if not g) == 35, "planar35")
    check(dict(passports) == {(4,1,1,1):7,(3,2,1,1):21,(2,2,2,1):7},
          "inherited heptic passports")
    check(dict(cross) == {(1,1):7,(2,0):14,(2,1):28,(4,0):21,(4,1):35},
          "Singer genus/direction intersection")
    check([ordinary_cross[d,0] for d in (1,2,4)] == [3,15,17],
          "changing cyclic sidecar changes the filter")
    # Another old105-set: a marked root and a 2-subset of the other six.
    # Compare its S7 action with the near-perfect matching action.
    transposition = {x:(2 if x == 1 else 1 if x == 2 else x) for x in range(8)}
    transform_m = lambda m: tuple(sorted(tuple(sorted((transposition[a],transposition[b])))
                                        for a,b in m))
    matching_fixed = sum(transform_m(m) == m for m in ms)
    root_pair_objects = [(a,bc) for a in range(1,8)
                         for bc in combinations([b for b in range(1,8) if b != a],2)]
    root_pair_fixed = sum(transposition[a] == a and
                          tuple(sorted(transposition[b] for b in bc)) == bc
                          for a,bc in root_pair_objects)
    check(len(root_pair_objects) == 105, "marked-root/pair105")
    check((matching_fixed,root_pair_fixed) == (15,35), "S7-equivariance hostile")
    return dict(direction_histogram=dict(by_direction),
                singer_cycle=cycles(singer),
                direction_genus=[dict(directions=d,genus=g,count=n)
                                 for (d,g),n in sorted(cross.items())],
                planar_passports={str(k):v for k,v in passports.items()},
                ordinary_cycle_planar_by_directions=[3,15,17],
                transposition_fixed=[matching_fixed,root_pair_fixed],
                records=records)


def clock_carriers():
    results = []
    for lengths in [(3,5,7),(3,5,8)]:
        blocks, perm, start = [], {}, 0
        for length in lengths:
            block = list(range(start,start+length))
            blocks.append(block)
            perm.update({x:block[(i+1)%length] for i,x in enumerate(block)})
            start += length
        order = lcm(*lengths)
        pairs = list(combinations(range(start),2))
        triples = [tuple(x) for x in product(*blocks)]
        pair_orbits = induced_orbits(pairs,perm)
        triple_orbits = induced_orbits(triples,perm)
        check(len(pairs) == len(triples) == order, "equal cardinalities")
        check(len(triple_orbits) == 1, "regular transversal clock")
        check(all(len(row) < order for row in pair_orbits), "pair clock loses phase")
        minimum = None
        for degree in range(1,start+1):
            if any(lcm(*p) == order for p in partition(degree)):
                minimum = degree
                break
        check(minimum == start, "minimum permutation degree")
        results.append(dict(blocks=lengths,order=order,minimum_degree=minimum,
                            pair_orbit_lengths=sorted(map(len,pair_orbits)),
                            triple_orbit_lengths=list(map(len,triple_orbits))))
    check(results[0]["pair_orbit_lengths"] == [3,5,5,7,7,7,15,21,35], "105 orbit list")
    orders8 = sorted({lcm(*p) for p in partition(8)})
    check(80 not in orders8 and 105 not in orders8, "eight-label clock obstruction")
    crt = {(k%3,k%5,k%7):k for k in range(105)}
    check(len(crt) == 105, "CRT reversible")
    for (a,b,c),k in crt.items():
        check((70*a+21*b+15*c)%105 == k, "CRT decoder")
    compatible = [(a,b) for a,b in product(range(15),range(35)) if (a-b)%5 == 0]
    check(len(compatible) == 105, "shared-mod5 fibre product")
    for a,b in compatible:
        k=(a+15*((b-a)%7))%105
        check((k%15,k%35) == (a,b), "coupled-pair decoder")
    return dict(carriers=results,S8_orders=orders8,
                compatible_pair_states=len(compatible),
                unrestricted_pair_states=15*35)


def small_bridges():
    # A certificate threshold, not a consequence of prime factorization.
    f = lambda H: Fraction(14,6*H+13)
    check(f(1)+12*f(104) == Fraction(1730,1729), "Kraft hostile104")
    check(f(1)+12*f(105) == Fraction(12194,12217), "Kraft positive105")
    # Fibonacci -> Gaussian square -> primitive Pythagorean triple.
    check(13*13-8*8 == 105 and 2*13*8 == 208 and 13*13+8*8 == 233, "Euclid")
    check(105*105+208*208 == 233*233, "Pythagoras")
    check(gcd(105,208) == 1, "primitive triple")
    # Independent shortcut Collatz prefixes preserve the exact source difference.
    def shortcut(n): return (3*n+1)//2 if n%2 else n//2
    a,b = 105,233
    odd=0
    for _ in range(7):
        check(a%2 == b%2, "shared parity")
        odd += a%2
        a,b = shortcut(a),shortcut(b)
    check(odd == 5 and b-a == 3**5, "literal105/233 shared prefix")
    # Inherited sibling clock retains translation, not just linear slope.
    hist=Counter()
    for n in range(105):
        v=n
        for t in range(1,7):
            v=(4*v+2)%105
            if v==n:
                hist[t]+=1
                break
        else: raise RuntimeError("sibling clock")
    check(set(hist) == {3,6}, "sibling periods")
    def oddstep(n):
        q=3*n+1
        a=0
        while q%2 == 0:
            q//=2
            a+=1
        return q,a
    check(1%105 == 211%105, "same scalar residue")
    check(oddstep(1) == (1,2) and oddstep(211) == (317,1), "native branch hostile")
    return dict(Kraft104=str(f(1)+12*f(104)),Kraft105=str(f(1)+12*f(105)),
                pythagorean=[105,208,233],collatz_after7=[a,b],
                collatz_odd_steps=odd,sibling_point_periods=dict(hist),
                same_residue_native_hostile=[[1,*oddstep(1)],[211,*oddstep(211)]])


def singer_ranks():
    records=[]
    for rank, polynomial in [(2,7),(3,11),(4,19),(5,37)]:
        Q=1 << rank
        def multiply_alpha(x):
            shifted=x << 1
            return shifted ^ polynomial if shifted & Q else shifted
        c={x:multiply_alpha(x) for x in range(1,Q)}
        check(sorted(map(len,cycles(c))) == [Q-1], "primitive Singer cycle")
        for a in range(1,Q):
            # Directly build the actual punctured involution and its response.
            tau={x:(x ^ a) if x != a else a for x in range(1,Q)}
            response={x:tau[c[x]] for x in range(1,Q)}
            lengths=sorted(map(len,cycles(response)))
            check(lengths == [1,Q-2], "all-r predicted passport")
            genus=1+((Q-1)-len(cycles(tau))-len(cycles(c))-len(cycles(response)))//2
            check(genus == Q//4-1, "all-r genus")
            full={x:multiply_alpha(x) ^ a for x in range(Q)}
            check(sorted(map(len,cycles(full))) == [1,Q-1], "independent affine passport")
            bypass={x:(full[0] if full[x] == 0 else full[x]) for x in range(1,Q)}
            check(bypass == response, "puncture retains correct successor")
        records.append(dict(rank=rank,field_size=Q,parameters=Q-1,
                            passport=[Q-2,1],genus=Q//4-1))
    return records


def weighted_controls():
    # Three two-dimensional orthogonal layers; a rank-one marked observer
    # versus the diagonal norm. Exact sign averaging is a separate path.
    blocks=((1,2),(3,4),(5,6))
    rank_one_energies=[]
    diagonal_energies=[]
    for signs in product((-1,1),repeat=3):
        values=[sign*x for sign,block in zip(signs,blocks) for x in block]
        rank_one_energies.append(sum(values)**2)
        diagonal_energies.append(sum(x*x for x in values))
    average=Fraction(sum(rank_one_energies),8)
    block_energy=sum(sum(block)**2 for block in blocks)
    cross_debt=2*sum(sum(a)*sum(b) for a,b in combinations(blocks,2))
    check(average == block_energy == 179, "sign average equals pinching")
    check(21**2-block_energy == cross_debt == 262, "exact discarded cross debt")
    check(set(diagonal_energies) == {91}, "block diagonal positive control")
    check(len(set(rank_one_energies)) > 1, "marked rank-one hostile")
    return dict(rank_one_energies=rank_one_energies,pinched_energy=block_energy,
                cross_debt=cross_debt,diagonal_energy=91)


def main():
    arithmetic_result=arithmetic()
    C,masks,characters=class_characters()
    out=dict(status="PROVED elementary maps; FINITE-EXACT enumerations",
             arithmetic=arithmetic_result,class_characters=characters,
             matching=matching_dictionary(C,masks),
             clocks=clock_carriers(),small_bridges=small_bridges(),
             singer_ranks=singer_ranks(),weighted_controls=weighted_controls())
    out["active_checks"]=CHECKS
    target=ROOT/"05-knowledge/results/bridges_105_20261018f.json"
    encoded=json.dumps(out,indent=2,sort_keys=True)+"\n"
    if "--verify" in sys.argv:
        check(target.read_text(encoding="utf-8")==encoded,"frozen receipt")
    else:
        target.write_text(encoded,encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k!="matching"},indent=2))
    print(json.dumps({k:v for k,v in out["matching"].items() if k!="records"},indent=2))


if __name__=="__main__":
    main()
