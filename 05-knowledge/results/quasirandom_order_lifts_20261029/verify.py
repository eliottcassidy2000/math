"""Exact controls for README propositions; spectral samples are labelled numerical.

Run: python3 -B 05-knowledge/results/quasirandom_order_lifts_20261029/verify.py
No random sampling or floating arithmetic is used to certify a theorem.
"""
from collections import Counter
from itertools import product
from math import gcd
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    result = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            result[d] = result.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def order(a, modulus, multiple):
    require(pow(a, multiple, modulus) == 1, "invalid order multiple")
    d = multiple
    for p in factor(multiple):
        while d % p == 0 and pow(a, d // p, modulus) == 1:
            d //= p
    return d


def subgroup(n, g):
    out, x = [], 1
    while x not in out:
        out.append(x)
        x = x * g % n
    require(x == 1, "not a unit")
    return out


def packets(n, g):
    left = set(range(1, n))
    out = []
    while left:
        x = min(left)
        orbit = {x * h % n for h in subgroup(n, g)}
        negative = {-y % n for y in orbit}
        require(not orbit & negative, "negation collision")
        packet = [(y, 1) for y in sorted(orbit)]
        packet += [(y, -1) for y in sorted(negative)]
        left.difference_update(orbit | negative)
        out.append(packet)
    return out


def covariance(n, ps):
    c = [0] * n
    for packet in ps:
        for a, s in packet:
            for b, t in packet:
                c[(a - b) % n] += s * t
    return c


def exhaustive_covariance(n, ps):
    # Independent route: enumerate actual tournaments, then all cyclic row products.
    sums = [0] * n
    for signs in product((-1, 1), repeat=len(ps)):
        symbol = [0] * n
        for sign, packet in zip(signs, ps):
            for x, v in packet:
                symbol[x] = sign * v
        require(all(symbol[x] == -symbol[-x % n] for x in range(n)), "oddness")
        for d in range(n):
            sums[d] += sum(symbol[x] * symbol[(x + d) % n] for x in range(n))
    return sums


def iterate(x, steps, modulus):
    for _ in range(steps):
        x = (x * x - x + 1) % modulus
    return x


def main():
    report = {}
    sylvester = [2]
    for _ in range(5):
        a = sylvester[-1]
        sylvester.append(a * a - a + 1)
    require(sylvester == [2, 3, 7, 43, 1807, 3263443], "Sylvester seeds")
    report["Sylvester_factors"] = {str(x): factor(x) for x in sylvester}
    require(factor(1807) == {13: 1, 139: 1}, "1807")
    require(all(gcd(x, y) == 1 for i, x in enumerate(sylvester) for y in sylvester[:i]), "coprimality")
    require(all(p >= 2*j+1 for j,n in enumerate(sylvester) if j>=2 for p in factor(n)), "prime factor first-arrival bound")
    deterministic_rows = []
    for j in range(2, 5):
        n = sylvester[j]
        p = min(q for q in factor(n) if q%4==3)
        m, g = n//p, (-sylvester[j-1])%n
        chi = [0] + [1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]
        require([sum(chi[x]*chi[(x+d)%p] for x in range(p)) for d in range(p)] == [p-1]+[-1]*(p-1), "Paley Gram")
        symbol = [0]*n
        for x in range(1,n):
            if x%p:
                symbol[x] = chi[x%p]
        for packet in packets(n,g):
            if symbol[packet[0][0]]==0:
                for x,v in packet:
                    symbol[x] = v
        require(all(symbol[x]==-symbol[-x%n] and symbol[x]==symbol[g*x%n] for x in range(n)), "deterministic orientation")
        require(all(symbol[x]==chi[x%p] for x in range(n) if x%p), "fixed outer ports")
        inner = [symbol[p*t] for t in range(m)]
        require(sum(inner)==0 and all(inner[t]==-inner[-t%m] for t in range(m)), "inner regularity")
        deterministic_rows.append([j,n,p,m,m*m*p])
    report["deterministic_family_j_N_p_M_norm_squared"] = deterministic_rows

    period_rows = []
    for k in range(1, 7):
        modulus = 19 ** k
        x, count = iterate(2, 1, modulus), 1
        while x != 2:
            x = (x * x - x + 1) % modulus
            count += 1
            require(count <= modulus, "not a cycle at source")
        predicted = 4 if k <= 2 else 12 * 19 ** (k - 3)
        require(count == predicted, "19-adic period")
        period_rows.append([k, modulus, count])
    report["19_adic_periods"] = period_rows
    require(iterate(2, 12, 19**4) == 2 + 17 * 19**3, "first nonzero return defect")
    require(3 * 5 * 13 * 85 % 361 == 330, "four-step derivative")
    require(pow(7, 3, 19) == 1 and 7 != 1, "derivative order")
    require([x for x in range(19) if (x*x-x+1)%19 == 0] == [8, 12], "roots")
    require([iterate(2, i, 19) for i in range(4)] == [2, 3, 7, 5], "root-avoidance orbit")

    exact_rows = []
    for n, g in [(7, 4), (13, 9), (19, 7), (43, 36), (121, 3), (1807, 1764)]:
        H, ps = subgroup(n, g), packets(n, g)
        require(len(H)%2 == 1 and all(gcd(h-1, n)==1 for h in H if h!=1), "Frobenius guard")
        c = covariance(n, ps)
        require(c == [n-1]+[-1]*(n-1), "exact isotropic covariance")
        if n <= 43:
            require(exhaustive_covariance(n, ps) == [x * 2**len(ps) for x in c], "independent covariance")
        exact_rows.append([n, g, len(H), len(ps)])
    report["isotropic_packets_N_g_order_bits"] = exact_rows

    # Complete small-modulus universe, including trivial groups and composites.
    # Every group satisfying the guard embeds into F_p^* for p|N, hence is cyclic.
    guarded_census = 0
    for n in range(3,102,2):
        seen = set()
        for g in range(1,n):
            if gcd(g,n)!=1:
                continue
            H = tuple(sorted(subgroup(n,g)))
            if H in seen:
                continue
            seen.add(H)
            if len(H)%2==0 or any(gcd(h-1,n)!=1 for h in H if h!=1):
                continue
            require(covariance(n,packets(n,g)) == [n-1]+[-1]*(n-1), "complete guarded census")
            guarded_census += 1
    report["all_guarded_subgroups_odd_N_3_through_101"] = guarded_census

    hostile_rows = []
    for n, g in [(21, 16), (49, 2), (121, 4), (1331, 3)]:
        ps = packets(n, g)
        c = covariance(n, ps)
        failures = [d for d in range(n) if c[d] != (n-1 if d==0 else -1)]
        require(failures, "missing hostile")
        hostile_rows.append([n, g, len(failures), failures[0], c[failures[0]]])
        if n in (49,121,1331):
            p = 7 if n==49 else 11
            require(all(c[d]==-p for d in range(n) if gcd(d,n)==1), "native unit-shift defect")
    report["ramified_hostiles_N_g_bad_shifts_first_value"] = hostile_rows

    tower_rows = []
    for p, base, expected_d, e in [(1093,2,364,2), (3511,2,1755,2), (11,3,5,2), (11,4,5,1)]:
        require(order(base,p,p-1)==expected_d, "prime order")
        require(pow(base, expected_d, p**e)==1 and pow(base, expected_d, p**(e+1))!=1, "exact valuation")
        for k in range(1,5):
            require(order(base,p**k,(p-1)*p**(k-1)) == expected_d*p**max(0,k-e), "order tower")
        odd_d = expected_d
        while odd_d%2 == 0:
            odd_d //= 2
        g = pow(base, expected_d//odd_d, p**e)
        H = subgroup(p**e,g)
        require(len(H)==odd_d and all(gcd(h-1,p**e)==1 for h in H if h!=1), "isotropy window")
        tower_rows.append([p,base,expected_d,e,odd_d,(p**e-1)//(2*odd_d)])
    report["order_towers_p_base_d_e_oddpart_bits_at_e"] = tower_rows

    # Actual source residues, not merely endpoint periods: n=(u*2^a-1)/3.
    source_rows = []
    for u in (1,5,7,11):
        for k in range(1,6):
            modulus = 3**(k+1)
            period = 2*3**k
            eligible = [a for a in range(period) if (u*pow(2,a,modulus)-1)%3==0]
            residues = {((u*pow(2,a,modulus)-1)//3)%3**k for a in eligible}
            require(len(eligible)==3**k and residues==set(range(3**k)), "source residue bijection")
            for phase in eligible[:1]:
                blocks = {((u*pow(2,phase+6*b,modulus)-1)//3)%3**k for b in range(3**(k-1))}
                require(len(blocks)==3**(k-1), "six-step address period")
            source_rows.append([u,k,len(residues)])
    report["Collatz_source_residue_counts"] = source_rows

    exponents = {2:5,3:7,5:3,7:4,11:1,13:2}
    def height(b):
        return min(exponents.get(p,0)//a for p,a in factor(b).items())
    require([height(x) for x in (4,2,5,3,7,13)]==[2,5,3,7,4,2], "inherited divisor-height repair")
    require(3*5*13*9==(3511-1)//2, "representative-dependent coincidence")
    require((3*5*13*85-3*5*13*9)%19==0, "derivative representative guard")
    report["status"] = "ALL EXACT CHECKS PASSED; no universal Collatz inference"
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
