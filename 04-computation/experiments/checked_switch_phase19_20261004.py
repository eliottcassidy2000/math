"""Exact route switches and a missing finite-field phase with a Paley geometry.

Run: python3 -B 04-computation/experiments/checked_switch_phase19_20261004.py
Standard library only. No assert statements; checks survive -O.
The imported arithmetic primitives have independent replay controls below.
"""
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from math import gcd, isqrt, lcm

from creative_decoder_20260925 import Word, Join, base_joins, odd_step
from cyclotomic_depth_towers_20261004 import BinaryField


CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(message)


def independent_step(n, sign=1):
    z, exponent = 3*n+sign, 0
    while z % 2 == 0:
        z //= 2
        exponent += 1
    return z, exponent


def replay(n, word, sign=1):
    for exponent in word:
        n, actual = independent_step(n, sign)
        check(actual == exponent, "independent valuation mismatch")
    return n


def compositions(total, length):
    for cuts in combinations(range(1, total), length-1):
        ends = (0,)+cuts+(total,)
        yield tuple(ends[i+1]-ends[i] for i in range(length))


def catalog(depth=6):
    """Keep every word in the declared finite universe, not only its minimum."""
    banks, rows = {}, []
    for length in range(1, depth+1):
        bank = defaultdict(list)
        for total in range(length, 2*length+5):
            for ks in compositions(total, length):
                w = Word.make(ks)
                address = w.C*pow(2**w.K, -1, 3**length) % 3**length
                bank[address].append(w)
        check(len(bank) == 2*3**(length-1), "inverse address coverage")
        banks[length] = bank
        rows.append({"length": length, "max_total": 2*length+4,
                     "words": sum(map(len, bank.values())),
                     "addresses": len(bank),
                     "largest_min_total": max(min(w.K for w in ws)
                                               for ws in bank.values())})
    return banks, rows


def first_reset_switch(n, sign=1):
    """An unbounded word schema, not a bounded orbit-search shortcut."""
    x, r = n, 0
    seen = set()
    while True:
        # Only the minus root 1 can stay in the initial exponent-one run.
        if x in seen:
            return None
        seen.add(x)
        y, a = odd_step(x, sign)
        if a != 1:
            break
        r += 1
        x = y
    if r == 0 or a < 3:
        return None
    m = (n-sign)//2
    w, v = (1,)*r+(a,), (1,)*(r-1)+(2, a-2)
    check(0 < m < n and m % 2, "reset rank guard")
    check(replay(n, w, sign) == replay(m, v, sign), "reset common future")
    return Join(n, m, r+1, r+1, "reset-switch")


def available(n, mode, banks):
    moves = list(base_joins(n))
    if mode == "baseline":
        return moves
    if mode in {"reset", "switch"}:
        reset = first_reset_switch(n)
        if reset:
            moves.append(reset)
    z = n
    for j in range(1, 7):
        z, _ = odd_step(z)
        if z < n:
            moves.append(Join(n, z, j, 0, "lookahead"))
        if mode == "switch":
            for k, bank in banks.items():
                for w in bank.get(z % 3**k, ()):
                    m = w.inverse(z)
                    if m is not None and m < n:
                        moves.append(Join(n, m, j, k, "address-switch"))
    return moves


def compile_interval(limit, mode, banks):
    """Increasing-source DAG. Missing entries use an explicitly counted oracle.

    The oracle is bounded actual first-descent search, not a theorem that
    this search always halts. Every learned word is closed under suffixes.
    """
    heights, learned, uses, seeds = {1: 0}, {}, defaultdict(int), []
    for n in range(3, limit+1, 2):
        moves = available(n, mode, banks)
        for w in learned.values():
            m = w.forward(n)
            if m is not None and m < n:
                moves.append(Join(n, m, w.L, 0, "learned-port"))
        if not moves:
            z, ks = n, []
            while z >= n:
                z, a = odd_step(z)
                ks.append(a)
                check(len(ks) <= 1000, "finite search cap")
            seeds.append(n)
            for start in range(len(ks)):
                w = Word.make(ks[start:])
                check(2**w.K > 3**w.L, "contracting suffix")
                learned[w.ks] = w
            moves.append(Join(n, z, len(ks), 0, "new-port"))
        choice = min(moves, key=lambda j: (j.a+max(heights[j.y]-j.b, 0),
                                           j.y, j.a, j.b, j.label))
        choice.verify()
        heights[n] = choice.a+max(heights[choice.y]-choice.b, 0)
        z = n
        for _ in range(heights[n]):
            z = independent_step(z)[0]
        check(z == 1, "independent root-certificate replay")
        uses[choice.label] += 1
    return {"mode": mode, "odd_sources": len(heights), "seed_count": len(seeds),
            "seed_sources": seeds, "ports": len(learned), "uses": dict(uses),
            "max_odd_steps": max(heights.values())}


def route_controls(banks):
    tested = 0
    for sign in (1, -1):
        for n in range(3, 20001, 2):
            j = first_reset_switch(n, sign)
            if j:
                j.verify(sign)
                tested += 1
    for a in range(2, 202, 2):
        n = 2**a-1
        j = first_reset_switch(n)
        check(j is not None and (j.y, j.a, j.b) == (2**(a-1)-1, a, a),
              "even Mersenne schema")
    for n in (7, 27):
        check(first_reset_switch(n) is None, "reset-two positive hostile")
    for n in (1, 5, 17):
        check(first_reset_switch(n, -1) is None, "minus basin minimum hostile")
    debt_examples, debt_count = [], 0
    for n in range(3, 20001, 4):
        r = (n+1 & -(n+1)).bit_length()-2
        t = (n+1) >> (r+1)
        a = 1+((3**(r+1)*t-1) & -(3**(r+1)*t-1)).bit_length()-1
        if a != 2:
            continue
        z = 3**r*t-1
        b = 1+(z & -z).bit_length()-1
        M, h = z >> (b-1), b-2
        check(b >= 3 and M % 2 == 1, "reset-two debt coordinates")
        Y = replay(n, (1,)*r+(2,))
        check(Y == 3*2**h*M+1, "reset-two debt transport")
        check(replay((n-1)//2, (1,)*(r-1)+(b,)) == M,
              "smaller source debt route")
        exponent3 = 1
        while h >= 3:
            Y = replay(Y, (2,))
            h -= 2
            exponent3 += 1
            check(Y == 3**exponent3*2**h*M+1, "debt normalization")
        target = independent_step(Y)[0]
        if h == 2:
            z = 3**(exponent3+1)*M+1
            check(target == z//(z & -z), "debt even boundary")
        else:
            check(h == 1 and target == 3**(exponent3+1)*M+2,
                  "debt odd boundary")
        debt_count += 1
        if len(debt_examples) < 5:
            debt_examples.append({"source": n, "smaller": (n-1)//2,
                                  "M": M, "h_at_boundary": h,
                                  "power3": exponent3, "boundary_output": target})
    # Negative control: changing 3 to 5 invalidates the reset switch.
    def u5(n):
        z = 5*n+1
        return z//(z & -z)
    check([u5(n) for n in (13, 33, 83)] == [33, 83, 13], "5n+1 cycle")
    check(u5(u5(19)) != u5(u5(9)), "do not transfer plus-three join to plus-five")
    # Direct enumeration of words tests both ternary guard and exact replay.
    inverse_tests = 0
    for k, bank in banks.items():
        for address, words in bank.items():
            z0 = address if address % 2 else address+3**k
            for w in words:
                for h in (0, 1):
                    z = z0+2*3**k*h
                    m = w.inverse(z)
                    if m is not None:
                        check(replay(m, w.ks) == z, "catalog inverse guard")
                        inverse_tests += 1
    return {"signed_reset_replays": tested, "even_Mersenne_exponents": [2, 200],
            "catalog_inverse_replays": inverse_tests,
            "reset_two_debt_replays": debt_count, "debt_examples": debt_examples,
            "hostiles": ["plus 7,27 reset exponent 2", "minus roots 1,5,17",
                         "5n+1 cycle 13,33,83"]}


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def prime_factors(n):
    result = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return result


def multiplicative_order(x, p):
    order = p-1
    for q in prime_factors(order):
        while order % q == 0 and pow(x, order//q, p) == 1:
            order //= q
    check(pow(x, order, p) == 1, "prime-field order upper")
    check(all(pow(x, order//q, p) != 1 for q in prime_factors(order)),
          "prime-field exact order")
    return order


def deeper_phase_controls():
    rows = []
    for p in (19, 5779, 87211):
        check(prime(p) and p % 24 == 19, "exact phase prime")
        o2 = multiplicative_order(2, p)
        o3 = multiplicative_order(3, p)
        o32 = multiplicative_order(3*pow(2, -1, p) % p, p)
        o4 = multiplicative_order(4, p)
        group_order = lcm(o2, o3)
        # In a cyclic group the join of subgroups has lcm order.
        translation_exponents = [a for a in range(o2) if pow(2, a, p) == 3]
        rows.append({"p": p, "order2": o2, "order3": o3, "order3_over2": o32,
                     "torus_layers": (p-1)//6,
                     "orientation_and_layer_period": o2//gcd(o2, 3),
                     "Frobenius4_layer_period": o4//gcd(o4, 6),
                     "group_2_3_order": group_order,
                     "odd_step_register_for_unit_slope": group_order//o2,
                     "translation_branch_exponents": translation_exponents})
    check(rows[-1]["odd_step_register_for_unit_slope"] == 323,
          "next binary prime requires odd-step count modulo323")
    check(rows[-1]["translation_branch_exponents"] == [],
          "offset translation does not persist at87211")
    return rows


def torus_layers(p):
    squares = {i*i % p for i in range(1, p)}
    omega = next(w for w in range(2, p) if pow(w, 3, p) == 1)
    h3 = {1, omega, omega*omega % p}
    left, all_edges, rows = set(squares), set(), []
    while left:
        d = min(left)
        steps = {d*h % p for h in h3}
        left -= steps
        faces = {tuple(sorted((x, (x+s) % p, (x+s+t) % p)))
                 for x in range(p) for s in steps for t in steps if s != t}
        incidence = defaultdict(list)
        for face in faces:
            for edge in combinations(face, 2):
                incidence[edge].append(face)
        check(all(len(fs) == 2 for fs in incidence.values()), "edge has two faces")
        check(len(faces) == 2*p and len(incidence) == 3*p, "torus counts")
        check(not (all_edges & incidence.keys()), "edge layers disjoint")
        all_edges.update(incidence)
        # Cell links, not induced graph neighborhoods; these differ at p=7.
        for x in range(p):
            link = defaultdict(set)
            for face in faces:
                if x in face:
                    a, b = (v for v in face if v != x)
                    link[a].add(b); link[b].add(a)
            check(len(link) == 6 and all(len(vs) == 2 for vs in link.values()),
                  "cellular link has six degree-two vertices")
            seen, pending = set(), [next(iter(link))]
            while pending:
                v = pending.pop()
                if v not in seen:
                    seen.add(v); pending.extend(link[v]-seen)
            check(len(seen) == 6, "cellular link is one C6")
        # Connectivity, independently of the quotient-lattice proof.
        seen, pending = set(), [next(iter(faces))]
        while pending:
            face = pending.pop()
            if face not in seen:
                seen.add(face)
                for edge in combinations(face, 2):
                    pending.extend(incidence[edge])
        check(len(seen) == 2*p, "each layer is connected")
        rows.append({"steps": sorted(steps), "V": p, "E": 3*p, "F": 2*p})
    check(len(all_edges) == p*(p-1)//2, "layers cover complete graph")
    return {"p": p, "omega": omega, "layers": rows}


def field_and_geometry():
    f = BinaryField(3)
    f.irreducibility_check()
    M, order_b = 2**18-1, 13797
    check(M == 27*7*19*73, "full order factorization")
    check(f.power(3, order_b) == 1, "b order upper")
    check(all(f.power(3, order_b//p) != 1 for p in (3, 7, 73)), "b exact order")
    seed = next(x for x in range(2, 50) if f.power(x, order_b) != 1)
    eta = f.power(seed, order_b)
    beta = f.mul(3, eta)
    check(seed == 11 and f.power(eta, 19) == 1 and eta != 1, "phase restorer")
    exclusions = {str(p): f.power(beta, M//p) for p in (3, 7, 19, 73)}
    check(f.power(beta, M) == 1 and all(x != 1 for x in exclusions.values()),
          "completed generator exact order")
    # Independent polynomial long multiplication, no imported mul or power.
    modulus = (1 << 18) | (1 << 9) | 1
    def mul(x, y):
        z = 0
        for i in range(y.bit_length()):
            if (y >> i) & 1:
                z ^= x << i
        while z.bit_length() >= modulus.bit_length():
            z ^= modulus << (z.bit_length()-modulus.bit_length())
        return z
    x, phase_values = 1, []
    for i in range(19):
        phase_values.append(x)
        check(x == f.power(eta, i), "independent phase arithmetic")
        x = mul(x, eta)
    check(x == 1 and len(set(phase_values)) == 19, "phase census")
    x = 1
    for _ in range(M):
        x = mul(x, beta)
        if x == 1:
            break
    check(_+1 == M, "independent exhaustive primitive period")
    geometry = [torus_layers(p) for p in range(7, 101) if p % 12 == 7 and prime(p)]
    row19 = next(row for row in geometry if row["p"] == 19)
    steps = [set(row["steps"]) for row in row19["layers"]]
    perm = [steps.index({4*x % 19 for x in s}) for s in steps]
    check(perm == [1, 2, 0], "fourth-power three-layer rotation")
    qr = {x*x % 19 for x in range(1, 19)}
    for a in range(19):
        for b in range(19):
            if a != b:
                check(((b-a) % 19 in qr) != ((2*b-2*a) % 19 in qr),
                      "Frobenius reverses Paley orientation")
    # Eisenstein multiplication (a+b omega)(c+d omega).
    def emul(x, y):
        a, b = x; c, d = y
        return a*c-b*d, a*d+b*c-b*d
    check(emul((1, -1), (5, 2)) == (7, -1), "Eisenstein factorization")
    check((5+2*7) % 19 == 0 and 25-10+4 == 19, "torus ideal kernel")
    check(3**3-2**3 == 19 and 3*pow(2, -1, 19) % 19 == 11,
          "Collatz clock equals missing phase modulus")
    branch_rows = []
    for a in range(1, 55):
        inv = pow(2, -a, 19)
        multiplier = 3*inv % 19
        preserved = multiplier in qr
        check(preserved == bool(a % 2), "exponent parity orientation")
        positive_multiplier = multiplier if preserved else -multiplier
        layer_perm = [steps.index({positive_multiplier*x % 19 for x in s})
                      for s in steps]
        if a <= 6:
            branch_rows.append({"a_mod6": a % 6, "orientation_preserved": preserved,
                                "layer_permutation": layer_perm})
        else:
            reference = branch_rows[(a-1) % 6]
            check(layer_perm == reference["layer_permutation"], "six-phase law")
        for x in range(19):
            new_index = inv*(3*x+1) % 19
            phase_update = f.power(f.mul(eta, f.power(phase_values[x], 3)), inv)
            check(phase_update == phase_values[new_index], "Collatz phase semiconjugacy")
    check(pow(2, 13, 19) == 3, "exponent-thirteen translation")
    for x in range(19):
        check((pow(2, -13, 19)*(3*x+1)-x) % 19 == 13,
              "translation branch exact increment")
    offset_triple = []
    for a in (1, 7, 13):
        inv = pow(2, -a, 19)
        values = list(range(19))
        for period in range(1, 20):
            values = [inv*(3*x+1) % 19 for x in values]
            if values == list(range(19)):
                break
        check(period == (19 if a == 13 else 3), "two rotations and one translation")
        offset_triple.append({"exponent": a, "permutation_order": period})
    check(2**12-2**6+1 != 57 and 2**6-2**3+1 == 57,
          "cyclotomic exponents: Phi18(2)=Phi6(8)")
    return {"field_dimension": 18, "modulus_bits": modulus, "order_b": order_b,
            "seed_bits": seed, "eta_bits": eta, "beta_bits": beta,
            "full_order": M, "prime_order_exclusions": exclusions,
            "Frobenius4_layer_permutation": perm, "Collatz_branches": branch_rows,
            "same_mod6_different_orders": offset_triple,
            "geometry": geometry}


def main():
    banks, rows = catalog()
    result = {"status": "FINITE-EXACT; universal Collatz coverage OPEN",
              "catalog": rows, "route_controls": route_controls(banks),
              "field_geometry": field_and_geometry(),
              "deeper_phase_controls": deeper_phase_controls()}
    result["compiler"] = []
    for mode in ("baseline", "lookahead", "reset", "switch"):
        row = compile_interval(10000, mode, banks)
        result["compiler"].append(row)
        print(mode, {k: v for k, v in row.items() if k != "seed_sources"}, flush=True)
    residual = {"odd_debt_height": 0, "even_debt_height": 0}
    phases = set()
    for n in result["compiler"][-1]["seed_sources"]:
        r = (n+1 & -(n+1)).bit_length()-2
        t = (n+1) >> (r+1)
        z = 3**(r+1)*t-1
        check((z & -z) == 2, "every residual seed has reset exponent2")
        z = 3**r*t-1
        h = (z & -z).bit_length()-2
        residual["odd_debt_height" if h % 2 else "even_debt_height"] += 1
        phases.add(n % 19)
    residual["source_phases_mod19"] = sorted(phases)
    result["residual_seed_types"] = residual
    result["checks"] = CHECKS
    out = Path(__file__).resolve().parents[2]/"05-knowledge/results"
    (out/"checked_switch_phase19_20261004.json").write_text(
        json.dumps(result, indent=2)+"\n")
    print("catalog", rows)
    print("controls", result["route_controls"])
    print("field", {k: v for k, v in result["field_geometry"].items() if k != "geometry"})
    print("torus primes", [r["p"] for r in result["field_geometry"]["geometry"]])
    print("deeper phases", result["deeper_phase_controls"])
    print("residual seed types", residual)
    print("explicit checks", CHECKS)


if __name__ == "__main__":
    main()
