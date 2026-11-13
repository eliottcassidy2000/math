"""Exact finite checks for the five-atom Hamiltonian-count fibre library.

Run with Python 3; no third-party dependencies are required.  Mathematical
all-order connectivity is proved separately by a normal form in PROOF.md.
"""

from __future__ import annotations

from collections import deque
from itertools import combinations, permutations
from math import comb, factorial
import json
from pathlib import Path


ATOMS = {
    "S": (0,),
    "T": (2, 4, 1),
    "Q": (6, 12, 8, 1),
    "U": (14, 28, 24, 16, 1),
    "V": (14, 20, 24, 18, 1),
}
MOVES = ((1, -2, 0, 1, 0), (2, -1, -1, 0, 1), (1, 1, -1, -1, 1))
BRIDGE_BEFORE = (44, 17, 2, 6, 13, 30, 63)
BRIDGE_AFTER = (110, 16, 2, 6, 13, 30, 62)


def need(value, message):
    if not value:
        raise RuntimeError(message)


def validate(rows):
    n = len(rows)
    for i in range(n):
        need(not (rows[i] >> i) & 1, "loop")
        need(rows[i] < (1 << n), "outside vertex")
        for j in range(i):
            need(((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) == 1, "not tournament")


def strong(rows):
    n = len(rows)
    for start in range(n):
        reached = 1 << start
        previous = -1
        while previous != reached:
            previous = reached
            for v in range(n):
                if (reached >> v) & 1:
                    reached |= rows[v]
        if reached != (1 << n)-1:
            return False
    return True


def counts(rows):
    """Hamiltonian path count of every induced subset, by terminal-vertex DP."""
    validate(rows)
    n = len(rows)
    ending = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            smaller = mask ^ (1 << v)
            if smaller == 0:
                ending[mask][v] = 1
            else:
                ending[mask][v] = sum(ending[smaller][u] for u in range(n)
                                      if (smaller >> u) & 1 and (rows[u] >> v) & 1)
        h[mask] = sum(ending[mask])
    return h


def permutation_count(rows):
    return sum(all((rows[p[i]] >> p[i+1]) & 1 for i in range(len(rows)-1))
               for p in permutations(range(len(rows))))


def ordinal(words):
    result = []
    offset = 0
    total = sum(len(ATOMS[word]) for word in words)
    for word in words:
        block = ATOMS[word]
        later = ((1 << total) - 1) ^ ((1 << (offset + len(block))) - 1)
        result.extend((r << offset) | later for r in block)
        offset += len(block)
    return tuple(result)


def c1_context(rows):
    """Add b,a with old->b, b->a, a->old: repository strong C_1 context."""
    n = len(rows)
    return tuple(r | (1 << n) for r in rows) + (1 << (n+1), (1 << n) - 1)


def profile(rows):
    """Unordered path-cover profile d_k; ordered c_k = k! d_k.

    Choose the block containing the least vertex to count each set partition
    once, and weight each block by its induced Hamiltonian count.
    """
    h = counts(rows)
    n = len(rows)
    d = [[0] * (n + 1) for _ in h]
    d[0][0] = 1
    for mask in range(1, 1 << n):
        first = mask & -mask
        rest = mask ^ first
        subset = rest
        while True:
            block = subset | first
            other = mask ^ block
            for k in range(1, n+1):
                d[mask][k] += h[block] * d[other][k-1]
            if subset == 0:
                break
            subset = (subset - 1) & rest
    return d[-1]


def polynomial(rows):
    # pi(z) = sum_k d_k (z)_k, where bins may be empty.
    d = profile(rows)
    n = len(rows)
    coefficients = [0] * (n+1)
    falling = [1]
    for k in range(n+1):
        for i, a in enumerate(falling):
            coefficients[i] += d[k] * a
        new = [0] * (len(falling)+1)
        for i, a in enumerate(falling):
            new[i] -= k * a
            new[i+1] += a
        falling = new
    return coefficients


def p2(rows):
    h = counts(rows)
    full = len(h)-1
    return sum(h[mask] * h[full ^ mask] for mask in range(len(h)))


def marked_parity(rows, root):
    """Independent enumeration of all simple directed paths starting at root."""
    h = counts(rows)
    full = len(h)-1
    eo = [0, 0]
    def visit(mask, last):
        eo[mask.bit_count() % 2] += h[full ^ mask]
        for nxt in range(len(rows)):
            if not (mask >> nxt) & 1 and (rows[last] >> nxt) & 1:
                visit(mask | (1 << nxt), nxt)
    visit(1 << root, root)
    return [h[-1], *eo]


def module_coefficients(rows, port):
    """Exact coefficients of ordered path covers of a substituted module.

    Enumerate the external vertex order. Each gap is forced occupied, forced
    empty, or optional. Selecting k gaps contributes one ordered k-path cover.
    """
    outside = [v for v in range(len(rows)) if v != port]
    coefficients = [0] * (len(rows)+1)
    for path in permutations(outside):
        forced = 0
        optional = int((rows[port] >> path[0]) & 1) + int((rows[path[-1]] >> port) & 1)
        possible = True
        for left, right in zip(path, path[1:]):
            direct = (rows[left] >> right) & 1
            through = ((rows[left] >> port) & 1) and ((rows[port] >> right) & 1)
            if not direct:
                if not through:
                    possible = False
                    break
                forced += 1
            elif through:
                optional += 1
        if possible:
            for k in range(optional+1):
                coefficients[forced+k] += comb(optional, k)
    return coefficients


def substitute(rows, ports, module):
    atoms = [module if v in ports else (0,) for v in range(len(rows))]
    offsets = []
    total = 0
    for atom in atoms:
        offsets.append(total)
        total += len(atom)
    output = [0]*total
    for v, atom in enumerate(atoms):
        for i, row in enumerate(atom):
            output[offsets[v]+i] = row << offsets[v]
            for w, other in enumerate(atoms):
                if (rows[v] >> w) & 1:
                    output[offsets[v]+i] |= ((1 << len(other))-1) << offsets[w]
    return tuple(output)


def bridge_audit():
    before, after = BRIDGE_BEFORE, BRIDGE_AFTER
    for rows in (before, after):
        need(counts(rows)[-1] == 25 and permutation_count(rows) == 25, "literal bridge count")
    need(before[6] == 63 and strong(tuple(r & 63 for r in before[:6])), "source plus strong six")
    need(strong(after), "strong seven endpoint")
    changed = [(i,j) for j in range(1,7) for i in range(j)
               if ((before[i] >> j) & 1) != ((after[i] >> j) & 1)]
    need(changed == [(0,1),(0,6)], "literal two-arc bridge")
    m6square = (124,97,114,118,2,16,48)
    permutation = (0,4,5,6,1,2,3)
    need(all(((m6square[i] >> j) & 1) == ((after[permutation[i]] >> permutation[j]) & 1)
             for i in range(7) for j in range(7)), "M6-square isomorphism")
    need(marked_parity(before,1) == [25,12,15], "before marked state")
    need(marked_parity(after,1) == [25,10,15], "after marked state")
    tt2, tt3 = (2,0), (6,4,0)
    transitive = tuple(tuple(((1 << m)-1) ^ ((1 << (i+1))-1) for i in range(m))
                       for m in range(8))
    modules = transitive + (ATOMS["T"], ATOMS["Q"])
    module_profiles = [(module,[factorial(k)*d for k,d in enumerate(profile(module))])
                       for module in modules]
    coefficients = []
    direct_checks = 0
    for port in range(7):
        a = module_coefficients(before,port)
        b = module_coefficients(after,port)
        delta = [y-x for x,y in zip(a,b)]
        need(delta[1] == 0, "singleton neutrality")
        need(delta[2] != 0 and all(x*delta[2] >= 0 for x in delta[2:]), "sign-definite inflation obstruction")
        for module, ordered in module_profiles:
            for rows, coeff in ((before,a),(after,b)):
                expected = sum(c*d for c,d in zip(coeff,ordered))
                actual = counts(substitute(rows,{port},module))[-1]
                need(actual == expected, "external-gap compiler versus independent subset DP")
                direct_checks += 1
        coefficients.append({"port": port, "before": a, "after": b, "delta": delta})
    two_port = []
    for pair in combinations(range(7),2):
        deltas = []
        for module in (tt2,tt3):
            deltas.append(counts(substitute(after,set(pair),module))[-1]
                          - counts(substitute(before,set(pair),module))[-1])
        need(deltas[0] != 0 and deltas[1] != 0, "unexpected two-port diagonal survivor")
        two_port.append({"ports": pair, "TT2_TT3_deltas": deltas})
    return {
        "provenance": "absorption/corridor.json codes1162->33931; orbits M6-square comparison",
        "before": {"rows": before, "H":25, "P2":p2(before), "polynomial":polynomial(before),
                   "root1_h_E_O":marked_parity(before,1)},
        "after": {"rows": after, "H":25, "P2":p2(after), "polynomial":polynomial(after),
                  "root1_h_E_O":marked_parity(after,1)},
        "M6_square_rows":m6square, "M6_square_to_after_permutation":permutation,
        "one_port_response_coefficients":coefficients,
        "independent_module_DP_checks":direct_checks,
        "modules": "transitive orders0..7, cycle3, strong4Q",
        "two_port_identical_transitive_probes":two_port,
    }


def signature(x):
    s, t, q, u, v = x
    return (s + 3*t + 4*q + 5*u + 5*v, t + 2*u + v, q+v)


def states(n, alpha, beta):
    result = []
    for v in range(min(alpha, beta)+1):
        for u in range((alpha-v)//2 + 1):
            t = alpha-2*u-v
            q = beta-v
            s = n-3*alpha-4*beta+u+2*v
            if s >= 0:
                result.append((s, t, q, u, v))
    return result


def normal_path(x):
    n, alpha, beta = signature(x)
    vmax = min(alpha, beta)
    route = [x]
    while x[4] < vmax:
        move = MOVES[1] if x[1] >= 1 else MOVES[2]
        x = tuple(a+b for a, b in zip(x, move))
        need(min(x) >= 0, "normal form left nonnegative fibre")
        need(signature(x) == (n, alpha, beta), "normal form left A fibre")
        route.append(x)
    while x[1] >= 2:
        x = tuple(a+b for a, b in zip(x, MOVES[0]))
        need(min(x) >= 0, "normal form left nonnegative fibre")
        need(signature(x) == (n, alpha, beta), "normal form left A fibre")
        route.append(x)
    return route


def response(x):
    return 2**(x[0] + 2*x[1] + 5*x[2] + 2*x[3] + 2*x[4]) * 3**(x[1] + x[3] + 3*x[4]) * 7**x[3]


def inverse_response(n, alpha, beta, value):
    """Recover a multiplicity vector from its desired two-bin response."""
    if value <= 0:
        return None
    remainder = value
    valuations = []
    for prime in (2, 3, 7):
        valuation = 0
        while remainder % prime == 0:
            remainder //= prime
            valuation += 1
        valuations.append(valuation)
    if remainder != 1:
        return None
    e2, e3, e7 = valuations
    u = e7
    twice_v = e3-alpha+u
    if twice_v < 0 or twice_v % 2:
        return None
    v = twice_v // 2
    x = (n-3*alpha-4*beta+u+2*v, alpha-2*u-v, beta-v, u, v)
    if min(x) < 0 or e2 != n-alpha+beta-u-3*v:
        return None
    return x


def components(fibre, moves):
    remaining = set(fibre)
    result = []
    while remaining:
        start = min(remaining)
        queue = deque([start])
        remaining.remove(start)
        component = [start]
        while queue:
            x = queue.popleft()
            for move in moves:
                for sign in (-1, 1):
                    y = tuple(a + sign*b for a, b in zip(x, move))
                    if y in remaining:
                        remaining.remove(y)
                        queue.append(y)
                        component.append(y)
        result.append(sorted(component))
    return sorted(result)


def word_stats(word, cover=False):
    rows = ordinal(word)
    out = {"word": word, "rows": list(rows), "n": len(rows), "H": counts(rows)[-1],
           "P2": p2(rows)}
    out["C1_H"] = counts(c1_context(rows))[-1]
    need(out["C1_H"] == out["H"] + out["P2"], "strong context formula mismatch")
    product_h = 1
    product_p2 = 1
    for atom in word:
        product_h *= counts(ATOMS[atom])[-1]
        product_p2 *= p2(ATOMS[atom])
    need(out["H"] == product_h, "ordinal H product mismatch")
    need(out["P2"] == product_p2, "ordinal P2 product mismatch")
    if cover:
        out["polynomial_ascending"] = polynomial(rows)
    return out


def main():
    atom_data = {}
    for name, rows in ATOMS.items():
        h = counts(rows)[-1]
        need(h == permutation_count(rows), "independent atom permutation count")
        need(strong(rows), "atom must be strongly connected")
        atom_data[name] = {"rows": list(rows), "n": len(rows), "H": h,
                          "P2": p2(rows), "polynomial_ascending": polynomial(rows)}
    need([atom_data[s]["H"] for s in ATOMS] == [1, 3, 5, 9, 15], "wrong atom counts")

    exchanges = []
    for before, after in (("TT", "SU"), ("TQ", "SSV"), ("QU", "STV")):
        pair = [word_stats(before), word_stats(after)]
        need(pair[0]["n"] == pair[1]["n"] and pair[0]["H"] == pair[1]["H"], "failed exchange")
        exchanges.append(pair)
    helper_trace = [word_stats(w) for w in ("SQU", "SUQ", "TTQ", "TSSV", "SSTV")]
    need({(s["n"], s["H"]) for s in helper_trace} == {(10, 45)}, "helper trace broke exact target")
    hostile = states(9, 2, 1)
    need(len(hostile) == 2, "expected two-state hostile fibre")
    need(len(components(hostile, MOVES[:2])) == 2, "two obvious moves should fail")
    need(len(components(hostile, MOVES)) == 1, "third move should repair")

    checked_fibres = 0
    checked_states = 0
    disconnected_two_move = 0
    maximum_normal_steps = 0
    # Exact universe: every nonempty A fibre with 0<=alpha,beta<=12, 0<=n<=60.
    for alpha in range(13):
        for beta in range(13):
            for n in range(61):
                fibre = states(n, alpha, beta)
                vmax = min(alpha, beta)
                nmin = 3*alpha + 4*beta - (alpha-vmax)//2 - 2*vmax
                need(bool(fibre) == (n >= nmin), "minimum order formula")
                if not fibre:
                    continue
                checked_fibres += 1
                checked_states += len(fibre)
                endpoints = set()
                for x in fibre:
                    route = normal_path(x)
                    endpoints.add(route[-1])
                    maximum_normal_steps = max(maximum_normal_steps, len(route)-1)
                    need(inverse_response(n, alpha, beta, response(x)) == x, "response inversion")
                need(len(endpoints) == 1, "nonunique normal form")
                need(len(components(fibre, MOVES)) == 1, "three-move disconnection")
                if len(components(fibre, MOVES[:2])) > 1:
                    disconnected_two_move += 1
    minimality = []
    for i, (n, alpha, beta) in enumerate(((6, 2, 0), (7, 1, 1), (9, 2, 1))):
        f = states(n, alpha, beta)
        need(len(f) == 2, "indispensability fibre has unexpected states")
        need(len(components(f, tuple(m for j, m in enumerate(MOVES) if j != i))) == 2,
             "move not indispensable")
        minimality.append({"signature": [n, alpha, beta], "states": f})
    output = {
        "status": "FINITE-EXACT",
        "atoms": atom_data,
        "moves": MOVES,
        "exchanges": exchanges,
        "helper_trace": helper_trace,
        "hostile_H45_n9": {"states": hostile, "two_move_components": components(hostile, MOVES[:2]),
                            "three_move_components": components(hostile, MOVES)},
        "indispensability": minimality,
        "response_geometry": {
            "parameterization": "(s,t,q,u,v)=(n-3*alpha-4*beta+u+2*v,alpha-2*u-v,beta-v,u,v)",
            "P2": "2^(n-alpha+beta)*3^alpha*(7/6)^u*(9/8)^v",
            "move_P2_ratios": ["7/6", "9/8", "27/28"],
            "hostile_desired_C1_count": 2735,
            "hostile_inverse": inverse_response(9, 2, 1, 2735-45),
        },
        "cross_lane_H25_bridge": bridge_audit(),
        "finite_normal_form_check": {
            "universe": "0<=alpha,beta<=12; 0<=n<=60; every nonnegative integer A-fibre",
            "fibres": checked_fibres, "states": checked_states,
            "two_move_disconnected_fibres": disconnected_two_move,
            "maximum_normal_steps": maximum_normal_steps,
        },
    }
    text = json.dumps(output, indent=2, sort_keys=True) + "\n"
    Path(__file__).with_name("receipt.json").write_text(text, encoding="utf-8")
    print(json.dumps({"atoms": atom_data, "exchanges": exchanges,
                      "finite_normal_form_check": output["finite_normal_form_check"]}, indent=2))


if __name__ == "__main__":
    main()
