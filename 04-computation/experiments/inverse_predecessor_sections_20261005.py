"""Exact inverse sections and the size-versus-flow-mass boundary.

No positive support assertion is inferred from inverse branching.
All checks remain active under python -O.
"""
from fractions import Fraction as F
from math import factorial
import json

CHECKS = 0
RESIDUES = (1, 2, 4, 5, 7, 8)
PALETTE = {
    1: (6, 2, 4), 2: (5, 1, 3), 4: (4, 6, 2),
    5: (1, 3, 5), 7: (2, 4, 6), 8: (3, 5, 1),
}
UNIT = dict(zip(RESIDUES, (2, 1, 2, 3, 4, 1)))
GROWING = dict(zip(RESIDUES, (2, 3, 2, 3, 4, 5)))


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def nat(n, positive=False):
    if type(n) is not int or n < int(positive):
        raise ValueError("exact natural integer required")


def odd(n):
    nat(n, True)
    if n % 2 != 1:
        raise ValueError("positive odd integer required")


def unit(n):
    odd(n)
    if n % 3 == 0:
        raise ValueError("target must be coprime to three")


def U(n):
    odd(n)
    x = 3*n+1
    return x // (x & -x)


def S(n, k):
    odd(n)
    nat(k)
    return 4**k*n+(4**k-1)//3


def inverse(n, exponent):
    unit(n)
    nat(exponent, True)
    numerator = 2**exponent*n-1
    if numerator % 3:
        raise ValueError("illegal inverse exponent")
    return numerator//3


def sibling(n, j):
    unit(n)
    nat(j)
    a = 2 if n % 3 == 1 else 1
    return inverse(n, a+2*j)


def section(n, kind):
    """Least strict predecessor of the stated type; ROOT self-loop omitted."""
    unit(n)
    if kind == "any":
        a = 2 if n % 3 == 1 else 1
    elif kind == "unit":
        a = UNIT[n % 9]
    elif kind == "divisible":
        a = PALETTE[n % 9][0]
    elif kind == "growing_unit":
        a = GROWING[n % 9]
    else:
        raise ValueError("unknown section")
    if n == 1 and kind in ("any", "unit", "growing_unit"):
        a = 4
    return a, inverse(n, a)


def weight(length, depth):
    nat(length)
    nat(depth)
    return F(2*factorial(depth)*factorial(length+1),
             factorial(length+depth+2))


def tail_fraction(length, depth, count):
    """Exact incoming tail after count minimal-exponent siblings."""
    nat(length)
    nat(depth)
    nat(count)
    result = F(1)
    for i in range(length+2):
        result *= F(depth+1+i, depth+count+1+i)
    return result


def ceil_fraction(x):
    return -((-x.numerator)//x.denominator)


def count_bound(length, depth, delta):
    nat(length)
    nat(depth)
    if type(delta) not in (int, F) or not 0 < delta < 1:
        raise ValueError("exact rational tail tolerance in (0,1) required")
    delta = F(delta)
    return ceil_fraction((1/delta-1)*F(depth+length+2, length+2))


def optimal_count(length, depth, delta):
    """Minimal first-sibling head with relative tail at most delta."""
    hi = count_bound(length, depth, delta)
    lo = 0
    while lo < hi:
        mid = (lo+hi)//2
        if tail_fraction(length, depth, mid) <= delta:
            hi = mid
        else:
            lo = mid+1
    return lo


def uniform_size_head(cap):
    """All x<cap*n lie in j<count, for every positive unit target n."""
    if type(cap) not in (int, F) or cap <= 0:
        raise ValueError("positive exact size cap required")
    count = 0
    while 2*4**count < 3*cap+1:
        count += 1
    return count


def main():
    report = {}
    targets = [n for n in range(1,4096,2) if n % 3]
    limits = {"any": F(4,3), "unit": F(16,3),
              "divisible": F(64,3), "growing_unit": F(32,3)}
    for n in targets:
        for row in range(3):
            a = PALETTE[n % 9][row]
            raw = inverse(n,a)
            need(raw % 3 == row and U(raw) == n, "complete residue palette")
            possible = [e for e in range(1,7)
                        if (2**e*n-1) % 9 == 3*row]
            need(possible == [a], "independent full exponent phase")
        for kind in limits:
            a, x = section(n,kind)
            need(x > 0 and x % 2 and U(x) == n and x != n,
                 "strict actual inverse section")
            candidates = []
            for e in range(1,a+1):
                if (2**e*n-1) % 3:
                    continue
                y = inverse(n,e)
                if y == n:
                    continue
                if kind == "unit" and y % 3 == 0:
                    continue
                if kind == "divisible" and y % 3:
                    continue
                if kind == "growing_unit" and (y % 3 == 0 or y <= n):
                    continue
                candidates.append((e,y))
            need(candidates[0] == (a,x), "literal minimality among earlier exponents")
            if n > 1 or kind != "any":
                need(F(x,n) < limits[kind], "sharp uniform ratio upper bound")
        if n > 1:
            need((section(n,"any")[1] < n) == (n % 3 == 2),
                 "iff unrestricted smaller predecessor")
            need((section(n,"unit")[1] < n) == (n % 9 in (2,8)),
                 "iff unit smaller predecessor")
            need((section(n,"divisible")[1] < n) == (n % 9 == 5),
                 "iff three-divisible smaller predecessor")
            small = [sibling(n,j) for j in range(7) if sibling(n,j) < 22*n]
            need(small == [sibling(n,j) for j in range(3)], "exact factor22 palette")
            need(sorted(x % 3 for x in small) == [0,1,2], "one child in every row")
    need(section(1,"any") == (4,5), "root self-loop excluded")
    need(inverse(1,2) == 1 and inverse(1,8) == 85, "root row-one exception")
    report["target_universe"] = {"positive_odd_units_below":4096,"count":len(targets)}
    report["palette_exponents_rows_0_1_2"] = PALETTE
    report["sharp_bounds"] = {k:str(v) for k,v in limits.items()}

    for kind, residue in (("any",1),("unit",7),("divisible",1),("growing_unit",17)):
        for t in (1,10,1000):
            n = 18*t+residue
            _, x = section(n,kind)
            need(limits[kind]-F(x,n) == F(1,3*n), "sharpness along unbounded arithmetic progression")

    rays = []
    for start in (1,5,7,11,13,17,19,25,31,37,41,43):
        n = start
        exponents = []
        for step in range(16):
            a, x = section(n,"growing_unit")
            need(n < x and x % 3 and U(x) == n, "infinite-ray rule exact finite replay")
            need(F(x,n) < F(32,3), "recursive growth bound")
            exponents.append(a)
            n = x
        rays.append({"start":start,"end_after16":n,"exponents":exponents})
    report["increasing_unit_rays"] = rays[:3]
    report["ray_controls"] = len(rays)*16

    adaptive_controls = 0
    for length in range(13):
        for depth in range(49):
            w = weight(length,depth)
            for count in (0,1,2,3,7,12):
                finite = sum(weight(length+1,depth+j) for j in range(count))
                need(finite/w+tail_fraction(length,depth,count) == 1,
                     "independent factorial versus product exact tail")
                need(tail_fraction(length,depth,count)
                     == weight(length,depth+count)/w,
                     "telescoping infinite remainder")
            for delta in (F(1,2),F(1,4),F(1,10)):
                bound = count_bound(length,depth,delta)
                optimal = optimal_count(length,depth,delta)
                need(0 < optimal <= bound, "adaptive finite bound")
                need(tail_fraction(length,depth,bound) <= delta, "Bernoulli capture guarantee")
                need(tail_fraction(length,depth,optimal) <= delta
                     < tail_fraction(length,depth,optimal-1), "exact minimal sibling count")
                adaptive_controls += 1
    report["formal_array_adaptive_controls"] = adaptive_controls

    escapes = []
    for k in (1,4,10,31,100,301):
        n = S(1,k)
        need(n % 3 and U(n) == 1, "literal rooted unit target")
        w = weight(0,k)
        head = sum(weight(1,k+j) for j in range(3))
        ratio = head/w
        need(ratio == 1-F((k+1)*(k+2),(k+4)*(k+5)), "factor22 exact mass escape")
        for j in range(3):
            x = sibling(n,j)
            need(U(x) == n and x < 22*n, "actual first3 rooted children")
        count = optimal_count(0,k,F(1,2))
        farthest = sibling(n,count-1)
        need(3*farthest < 4**count*n, "adaptive maximum size cost")
        escapes.append({"K":k,"first3_mass_fraction":str(ratio),
                        "minimal_half_mass_count":count,
                        "target_bits":n.bit_length(),"last_selected_source_bits":farthest.bit_length()})
    report["root_ray_fixed_cap_escape_and_adaptive_repair"] = escapes
    for cap in (F(1,2),F(1),F(2),F(22),F(100)):
        count = uniform_size_head(cap)
        for n in targets[:100]:
            for j in range(count+3):
                need(not (sibling(n,j) < cap*n) or j < count, "uniform fixed size cap")
    need(uniform_size_head(22) == 3, "factor22 uniform count")

    heavy = []
    for t in range(16):
        k = 9*t+1
        n = S(1,k)
        z = inverse(n,1)
        need(n % 9 == 5 and z % 3 == 0 and U(z) == n and U(n) == 1,
             "explicit three-divisible source with two-step root certificate")
        need(weight(1,k)/weight(0,k) == F(2,k+3), "near predecessor still has vanishing relative mass")
        need(weight(1,k) == F(4,(k+1)*(k+2)*(k+3)), "injection atom exact polynomial tail")
        if t in (0,1,15):
            heavy.append({"t":t,"K":k,"n":n,"z":z,"lambda_atom":str(weight(1,k))})
    report["certified_heavy_size_tail_examples"] = heavy

    bad = [
        lambda: section(True,"unit"),lambda: section(3,"unit"),
        lambda: section(7,"unknown"),lambda: inverse(7,1),
        lambda: inverse(7,False),lambda: tail_fraction(1,2,-1),
        lambda: count_bound(1,2,0.5),lambda: count_bound(1,2,F(1)),
        lambda: optimal_count(False,2,F(1,2)),lambda: uniform_size_head(22.0),
    ]
    for job in bad:
        try:
            job()
        except ValueError:
            need(True,"invalid exact input rejected")
        else:
            raise ValueError("hostile input accepted")
    report["malformed_controls"] = len(bad)
    report["exact_checks"] = CHECKS
    report["scope"] = "Sections and conditional weight transport proved; no universal positivity claim."
    print(json.dumps(report,indent=2,sort_keys=True))
    print("PASS")


if __name__ == "__main__":
    main()
