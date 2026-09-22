"""Exact odd-square triangle chart, inherited descent, and shell permutation.

Standard library only. All checks use exceptions and survive python -O.
No floating-point calculation decides an arithmetic or dynamical claim.
"""
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json


CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


def triple(u, v):
    return u*v, (u*u-v*v)//2, (u*u+v*v)//2


def primitive(t):
    a,b,c = t
    return a > 0 and b > 0 and a % 2 == 1 and b % 2 == 0 and a*a+b*b == c*c and gcd(gcd(a,b),c) == 1


def roots(t):
    a,b,c = t
    u,v = isqrt(c+b),isqrt(c-b)
    check(u*u == c+b and v*v == c-b, "inverse roots are exact squares")
    return u,v


def children(u,v):
    return {"L":(u+2*v,v),"M":(2*u+v,u),"R":(2*u-v,u)}


def parent(u,v):
    if u == 3*v:
        return None
    if u > 3*v:
        return "L",(u-2*v,v)
    if u > 2*v:
        return "M",(v,u-2*v)
    return "R",(v,2*v-u)


def shell_step(t):
    a,b,c = t
    return abs(b+c-2*a),2*(a+b-c),3*c-b-2*a


def shell_inverse(u,w):
    candidates = [(u-w)//2,(u+w)//2]
    odd = [v for v in candidates if v % 2]
    check(len(odd) == 1, "unique parity-selected inverse")
    return odd[0]


def phi(u):
    return sum(gcd(v,u) == 1 for v in range(1,u))


def signed_order_two(u):
    r = 1
    for k in range(1,phi(u)+1):
        r = 2*r % u
        if r in (1,u-1):
            return k
    raise RuntimeError("unit orbit failed to return")


def cycles(u, fibre):
    unseen = set(fibre)
    answer = []
    while unseen:
        start = min(unseen)
        orbit = []
        v = start
        while v not in orbit:
            check(v in unseen, "shell cycles disjoint")
            orbit.append(v)
            unseen.remove(v)
            v = abs(u-2*v)
        check(v == start, "shell cycle returns to own start")
        answer.append(orbit)
    return answer


def trace_identity(u,v,w):
    # Laurent traces in Z[z]/(z^u-1), an exact path independent of cosine floats.
    trace = Counter({v % u:1,(-v) % u:1})
    square = Counter()
    for i,a in trace.items():
        for j,b in trace.items():
            square[(i+j) % u] += a*b
    square[0] -= 2
    square = {i:a for i,a in square.items() if a}
    target = {w % u:1,(-w) % u:1}
    return square == target


def integer_family(k):
    a,b,c = k*k-1,2*k,k*k+1
    d = gcd(gcd(a,b),c)
    a,b,c = a//d,b//d,c//d
    if a % 2 == 0:
        a,b = b,a
    return a,b,c


def family_count(x):
    if x < 5:
        return 0
    return isqrt(x-1)//2+(isqrt(2*x-1)-1)//2-1


def accelerated(n,b):
    value = 3*n+b
    if value == 0:
        raise ValueError("singular numerator")
    return value // (abs(value) & -abs(value))


def main():
    cap = 301
    nodes = set()
    shell_records = []
    nonprimitive = 0
    for u in range(3,cap+1,2):
        fibre = [v for v in range(1,u,2) if gcd(u,v) == 1]
        check(len(fibre)*2 == phi(u), "totient fibre count")
        for v in range(1,u,2):
            t = triple(u,v)
            g = gcd(u,v)
            check(gcd(gcd(t[0],t[1]),t[2]) == g*g, "exact odd-square content")
            if g != 1:
                nonprimitive += 1
                continue
            nodes.add((u,v))
            check(primitive(t) and roots(t) == (u,v), "lossless primitive chart")
            for branch,child in children(u,v).items():
                check(gcd(*child) == 1 and primitive(triple(*child)), "Berggren child admissible")
                check(parent(*child) == (branch,(u,v)), "Berggren child-parent inverse")
            point = u,v
            depth = 0
            while parent(*point) is not None:
                branch,up = parent(*point)
                check(0 < up[1] < up[0] < point[0] and gcd(*up) == 1, "strict admissible descent")
                point = up
                depth += 1
            check(point == (3,1), "primitive descent ends at root")
            boundary = v == 1 or u-v == 2
            check((depth == (u-3)//2) == boundary, "boundary rays are maximal-depth shells")
            w = abs(u-2*v)
            dt = shell_step(t)
            check(dt == triple(u,w) and primitive(dt), "side-length shell operation")
            check(dt[1]+dt[2] == u*u and dt[1] == 4*((t[0]+t[1]-t[2])//2), "shell and inradius preserved law")
            check(shell_inverse(u,w) == v, "arithmetic shell inverse")
            check({w % u,(-w) % u} == {(2*v) % u,(-2*v) % u}, "unit-coset doubling")
            check(trace_identity(u,v,w), "exact Laurent Chebyshev conjugacy")
        orbits = cycles(u,fibre)
        order = signed_order_two(u)
        check(all(len(orbit) == order for orbit in orbits), "uniform shell period")
        check(2*order*len(orbits) == phi(u), "shell orbit count")
        shell_records.append({"u":u,"fibre_size":len(fibre),"period":order,"cycles":orbits})

    # Independent tree enumeration has no coordinate-fibre selection filter.
    generated = set()
    queue = deque([(3,1)])
    while queue:
        node = queue.popleft()
        check(node not in generated, "Berggren tree no repeated node")
        generated.add(node)
        queue.extend(child for child in children(*node).values() if child[0] <= cap)
    check(generated == nodes, "tree and totative enumeration agree")

    height = 10000
    euclid = set()
    for m in range(2,isqrt(height)+1):
        for n in range(1,m):
            if gcd(m,n) == 1 and (m-n) % 2 and m*m+n*n <= height:
                euclid.add((m*m-n*n,2*m*n,m*m+n*n))
    odd_roots = {triple(u,v) for u in range(3,isqrt(2*height)+1,2) for v in range(1,u,2)
                 if gcd(u,v) == 1 and (u*u+v*v)//2 <= height}
    check(euclid == odd_roots and len(euclid) == 1593, "independent height census")
    direct = set()
    for c in range(3,1001,2):
        for a in range(1,c,2):
            b = isqrt(c*c-a*a)
            if b and b % 2 == 0 and a*a+b*b == c*c and gcd(a,b) == 1:
                direct.add((a,b,c))
    check(direct == {t for t in euclid if t[2] <= 1000} and len(direct) == 158,
          "Euclid-free direct side census")

    integers = {}
    for k in range(2,2001):
        t = integer_family(k)
        check(primitive(t), "integer k reduction primitive")
        u,v = roots(t)
        check((u,v) == ((k+1,k-1) if k % 2 == 0 else (k,1)), "integer family exactly boundary rays")
        integers.setdefault(t,[]).append(k)
    collisions = [ks for ks in integers.values() if len(ks)>1]
    check(collisions == [[2,3]], "only repeated integer shape")
    cutoffs = [0,4,5,13,17,29,100,1000,10000,100000,1000000]
    sparse_counts = []
    for x in cutoffs:
        actual = sum(t[2] <= x for t in integers)
        check(actual == family_count(x), "exact integer-family counting function")
        sparse_counts.append({"X":x,"count":actual})
    omitted = min((t for t in euclid if not (roots(t)[1] == 1 or roots(t)[0]-roots(t)[1] == 2)),key=lambda t:t[2])
    check(omitted == (21,20,29), "first integer-family omitted triple")

    rational_controls = 0
    for u,v in nodes:
        k = Fraction(u+v,u-v)
        mirror = Fraction(u,v)
        check((k+1)/(k-1) == mirror and (mirror+1)/(mirror-1) == k, "two rational k sheets")
        for x in (k,mirror):
            p,q = x.numerator,x.denominator
            raw = (p*p-q*q,2*p*q,p*p+q*q)
            d = gcd(gcd(raw[0],raw[1]),raw[2])
            check(d == (2 if p % 2 and q % 2 else 1), "rational k content law")
            a,b,c = (z//d for z in raw)
            if a % 2 == 0:
                a,b = b,a
            check((a,b,c) == triple(u,v), "both rational k sheets reconstruct same ordered triple")
            rational_controls += 1

    special = {r['u']:r for r in shell_records if r['u'] in (3,5,7,9,13,17,21,63,65)}
    period3 = [r['u'] for r in shell_records if r['period'] == 3]
    period6 = [r['u'] for r in shell_records if r['period'] == 6]
    check(period3 == [7,9] and period6 == [13,21,63,65], "period-three and period-six conductor controls")
    check(special[17]['cycles'] == [[1,15,13,9],[3,11,5,7]], "first shell doubling-transitivity hostile")
    incomplete = [r['u'] for r in shell_records if len(r['cycles']) > 1]
    check(incomplete[0] == 17, "first incomplete closure shell")
    check(triple(17,3) == (51,140,149), "first missing node of boundary-shell closure")
    check(shell_step((35,12,37)) == (21,20,29), "Collatz edge triangle leaves edge locus")
    check(accelerated(5,-1) == 7 and accelerated(7,-1) == 5, "known Collatz two-cycle control")
    check(all(accelerated(a,b) != z for b in (-1,1) for a,z in ((7,3),(3,7))),
          "shell permutation does not preserve either Collatz edge relation")
    check(shell_step((27,36,45)) == (27,36,45) and not primitive((27,36,45)), "nonprimitive fixed-point hostile")
    power_neighbours = []
    for exponent in range(1,17):
        for sign in (-1,1):
            u = 2**exponent+sign
            if u < 3:
                continue
            expected = 1 if (exponent,sign) == (2,-1) else exponent
            period = signed_order_two(u)
            check(period == expected, "power-of-two neighbour exact signed order")
            power_neighbours.append({"exponent":exponent,"sign":sign,"u":u,
                                     "period":period,"cycles":phi(u)//(2*period)})
    matrix = [[-2,1,1],[2,2,-2],[-2,-1,3]]
    metric = [1,1,-1]
    for i in range(3):
        for j in range(3):
            check(sum(matrix[k][i]*metric[k]*matrix[k][j] for k in range(3)) == (4*metric[i] if i == j else 0),
                  "conformal Lorentz matrix identity")

    source = Path(__file__)
    output = {"status":"PASS","source_sha256_lf":sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
              "scope":"Elementary written proofs plus finite exact controls; no Collatz conjugacy or convergence claim.",
              "universe":{"odd_root_max":cap,"primitive_root_pairs":len(nodes),"nonprimitive_root_pairs":nonprimitive,
                          "independent_hypotenuse_max":height,"independent_hypotenuse_count":len(euclid),
                          "direct_side_hypotenuse_max":1000,"direct_side_count":len(direct),
                          "integer_k_range":[2,2000],"rational_k_reconstructions":rational_controls},
              "shells":shell_records,"selected_shells":special,"integer_family_counts":sparse_counts,
              "integer_family_only_collision":collisions,"first_integer_family_omission":omitted,
              "period_three_conductors":period3,"period_six_conductors":period6,
              "period_six_total_cycles":sum(special[u]['fibre_size']//6 for u in period6),
              "first_incomplete_boundary_closure":{"u":17,"omitted_triple":[51,140,149]},
              "power_of_two_neighbour_controls":power_neighbours,
              "checks_passed":CHECKS}
    source.with_suffix('.json').write_text(json.dumps(output,indent=2)+"\n",encoding='utf-8')
    print(json.dumps({"status":"PASS","checks_passed":CHECKS,"primitive_root_pairs":len(nodes),
                      "height_count":len(euclid),"output":source.with_suffix('.json').name}))


if __name__ == '__main__':
    main()
