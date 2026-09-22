"""Exhaustive order-four tournament/XOR audit, with two named order-five hostiles.

Standard library only; all checks remain active under python -O.
Bits in lexicographic edge order mean smaller vertex points to larger vertex.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path


CHECKS = []
EDGES = list(combinations(range(4), 2))
PERMS = list(permutations(range(4)))


def require(test, label):
    if not test:
        raise RuntimeError(label)
    CHECKS.append(label)


def bit(mask, i, j):
    value = (mask >> EDGES.index(tuple(sorted((i, j))))) & 1
    return value if i < j else 1-value


def relabel(mask, p):
    result = 0
    for i, j in EDGES:
        u, v = p[i], p[j]
        direction = bit(mask, i, j)
        if u > v:
            u, v, direction = v, u, 1-direction
        result |= direction << EDGES.index((u, v))
    return result


def cutmask(vertex_bits):
    return sum((((vertex_bits >> i) ^ (vertex_bits >> j)) & 1) << k
               for k, (i, j) in enumerate(EDGES))


CUTS = sorted({cutmask(v) for v in range(16)})


def signature(mask):
    return tuple(bit(mask, 0, i) ^ bit(mask, 0, j) ^ bit(mask, i, j)
                 for i, j in combinations((1, 2, 3), 2))


def fixed_path_gauge(mask):
    signs = [0]
    for i in range(3):
        signs.append(signs[-1] ^ bit(mask, i, i+1) ^ 1)
    return mask ^ cutmask(sum(x << i for i, x in enumerate(signs)))


def path_count(adj):
    n = len(adj)
    direct = sum(all(adj[p[i]][p[i+1]] for i in range(n-1))
                 for p in permutations(range(n)))
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for used in range(1, 1 << n):
        for last in range(n):
            ways = dp.get((used, last), 0)
            if not ways:
                continue
            for nxt in range(n):
                if not used >> nxt & 1 and adj[last][nxt]:
                    key = used | 1 << nxt, nxt
                    dp[key] = dp.get(key, 0)+ways
    other = sum(dp.get(((1 << n)-1, v), 0) for v in range(n))
    require(direct == other, f"independent Hamilton-path counts on {adj}")
    return direct


def strongly_connected(adj):
    n = len(adj)
    reach = [[bool(adj[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] |= reach[i][k] and reach[k][j]
    return all(map(all, reach))


def parity(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def determinant(matrix):
    result = 0
    for p in permutations(range(len(matrix))):
        product = parity(p)
        for i, j in enumerate(p):
            product *= matrix[i][j]
        result += product
    return result


def invariants(mask):
    adj = [[0 if i == j else bit(mask, i, j) for j in range(4)] for i in range(4)]
    degrees = sorted(map(sum, adj))
    kind = {(0,1,2,3):"transitive", (1,1,1,3):"source_cycle",
            (0,2,2,2):"sink_cycle", (1,1,2,2):"strong"}[tuple(degrees)]
    triangles = sum(adj[i][j]*adj[j][k]*adj[k][i]+adj[i][k]*adj[k][j]*adj[j][i]
                    for i,j,k in combinations(range(4),3))
    degree_triangles = 4-sum(d*(d-1)//2 for d in degrees)
    require(triangles == degree_triangles, f"mask {mask} triangle degree identity")
    h = path_count(adj)
    require(h == 1+2*triangles, f"mask {mask} order-four odd-cycle identity")
    cycles = sum(all(adj[p[i]][p[(i+1)%4]] for i in range(4))
                 for p in PERMS if p[0] == 0)
    strong = strongly_connected(adj)
    require(strong == (kind == "strong") == (cycles == 1), f"mask {mask} strong/cycle classification")
    s = [[0 if i == j else 2*adj[i][j]-1 for j in range(4)] for i in range(4)]
    pf = s[0][1]*s[2][3]-s[0][2]*s[1][3]+s[0][3]*s[1][2]
    det = determinant(s)
    require(det == pf*pf == 1+8*(triangles % 2), f"mask {mask} Pfaffian determinant/parity")
    autos = sum(relabel(mask,p) == mask for p in PERMS)
    require(autos == (3 if kind in ("source_cycle","sink_cycle") else 1), f"mask {mask} automorphisms")
    return {"mask":mask,"edge_bits":[bit(mask,i,j) for i,j in EDGES],"type":kind,
            "outdegrees_sorted":degrees,"directed_triangles":triangles,"hamilton_paths":h,
            "hamilton_cycles_mod_rotation":cycles,"automorphisms":autos,"pfaffian":pf,
            "determinant":det,"isomorphism_representative":min(relabel(mask,p) for p in PERMS),
            "switch_signature":list(signature(mask)),"path_gauge":fixed_path_gauge(mask)}


def main():
    require(len(CUTS) == 8, "cut-space has eight elements")
    rows = [invariants(mask) for mask in range(64)]
    by_type = dict(sorted(Counter(r["type"] for r in rows).items()))
    require(by_type == {"sink_cycle":8,"source_cycle":8,"strong":24,"transitive":24}, "all64 type census")
    require(len({r["isomorphism_representative"] for r in rows}) == 4, "four isomorphism classes")
    fixed_path = [m for m in range(64) if all(bit(m,i,i+1) for i in range(3))]
    fixed_cycle = [m for m in fixed_path if bit(m,3,0)]
    require(len(fixed_path) == 8 and len(fixed_cycle) == 4, "path and cycle family cardinalities")
    path_table = []
    for m in fixed_path:
        a,b,c = bit(m,0,2),bit(m,1,3),bit(m,0,3)
        require(rows[m]["hamilton_paths"] == 5-2*c*(a+b), f"path mask {m} exact formula")
        require(rows[m]["directed_triangles"] == 2-c*(a+b), f"path mask {m} triangle formula")
        path_table.append({"a":a,"b":b,"c":c,"mask":m,"type":rows[m]["type"],
                           "hamilton_paths":rows[m]["hamilton_paths"]})
    cycle_rep = rows[fixed_cycle[0]]["isomorphism_representative"]
    rotation = (1,2,3,0)
    diagonal_a = 1 << EDGES.index((0,2))
    diagonal_b = 1 << EDGES.index((1,3))
    for m in fixed_cycle:
        a,b = bit(m,0,2),bit(m,1,3)
        rotated = relabel(m,rotation)
        require((bit(rotated,0,2),bit(rotated,1,3)) == (1-b,a), f"cycle mask {m} affine rotation")
        require(rows[m]["isomorphism_representative"] == cycle_rep and rows[m]["hamilton_paths"] == 5,
                f"cycle mask {m} common type/path count")
        require(rows[m]["pfaffian"] == -(2*a-1)*(2*b-1), f"cycle mask {m} diagonal XOR Pfaffian")
        require(signature(m) == (a,1-b,1-a), f"cycle mask {m} affine signature plane")
        require(relabel(relabel(m,rotation),rotation) == m^diagonal_a^diagonal_b,
                f"cycle mask {m} half-turn equals both flips")
    require(all(not all(relabel(m,p) == m^diagonal_a for m in fixed_cycle) for p in PERMS),
            "no uniform vertex relabelling implements one diagonal flip")
    require(diagonal_a not in CUTS and diagonal_b not in CUTS and diagonal_a^diagonal_b not in CUTS,
            "diagonal toggles are not vertex cuts")
    require(len({signature(m) for m in fixed_cycle}) == 4, "isomorphic fixed-cycle states have distinct switching classes")
    orbits = []
    seen = set()
    for m in range(64):
        if m in seen:
            continue
        orbit = sorted(m^cut for cut in CUTS)
        seen.update(orbit)
        kinds = dict(sorted(Counter(rows[x]["type"] for x in orbit).items()))
        require(len({signature(x) for x in orbit}) == 1, f"switching orbit {m} cycle signatures")
        require(len({rows[x]["determinant"] for x in orbit}) == 1, f"switching orbit {m} determinant")
        require(len(set(orbit)&set(fixed_path)) == 1, f"switching orbit {m} unique path gauge")
        require(len({fixed_path_gauge(x) for x in orbit}) == 1, f"switching orbit {m} gauge construction")
        require(sum(rows[x]["hamilton_paths"] for x in orbit) == 24, f"switching orbit {m} inherited mean three")
        require(kinds in ({"strong":4,"transitive":4},{"sink_cycle":4,"source_cycle":4}),
                f"switching orbit {m} exact type mixture")
        orbits.append({"representative":m,"members":orbit,"types":kinds,
                       "signature":list(signature(m)),"determinant":rows[m]["determinant"],
                       "path_representative":fixed_path_gauge(m)})
    require(len(orbits) == 8 and seen == set(range(64)), "complete switching partition")
    for m in range(64):
        for n in range(64):
            require((signature(m) == signature(n)) == ((m^n) in CUTS), f"complete gauge invariant pair {m},{n}")
        for p in PERMS:
            n = relabel(m,p)
            require(rows[n]["pfaffian"] == parity(p)*rows[m]["pfaffian"], f"Pfaffian relabelling {m},{p}")
    hostile = 63^cutmask(1 << 1)
    require(rows[63]["type"] == "transitive" and rows[hostile]["type"] == "strong",
            "switching loses strong connectivity and H")
    hsum = sum(r["hamilton_paths"] for r in rows)
    h2sum = sum(r["hamilton_paths"]**2 for r in rows)
    require(hsum == 192 and h2sum == 768, "Hamiltonian presentation multiplicity totals")
    require(Fraction(sum(rows[m]["hamilton_paths"] for m in fixed_path),8) == Fraction(h2sum,hsum) == 4,
            "fixed-path sampling has H-size bias")
    # Named order-five controls only, not an exhaustive order-five universe.
    near_transitive = [[0 if i == j else int(i < j) for j in range(5)] for i in range(5)]
    near_transitive[0][4],near_transitive[4][0] = 0,1
    regular = [[int((j-i)%5 in (1,2)) for j in range(5)] for i in range(5)]
    for adj in (near_transitive,regular):
        require(all(adj[i][(i+1)%5] for i in range(5)) and strongly_connected(adj), "named order-five same-cycle control")
    h5 = [path_count(adj) for adj in (near_transitive,regular)]
    require(h5 == [9,15], "fixed-cycle isomorphism collapse fails already at five vertices")
    source = Path(__file__)
    result = {"status":"PASS","source_sha256_lf":sha256(source.read_bytes().replace(b"\r\n",b"\n")).hexdigest(),
              "scope":"Complete labelled order-four universe; two explicitly named order-five hostiles only.",
              "universe":{"vertices":[0,1,2,3],"edge_order":EDGES,"masks":[0,63],
                          "label_permutations":24,"vertex_switches_distinct":8,"gauge_pair_checks":4096,
                          "order_five_controls":["transitive except 4->0","cyclic regular increments 1,2"]},
              "isomorphism_type_counts":by_type,"fixed_path_table":path_table,"fixed_cycle_masks":fixed_cycle,
              "switching_orbits":orbits,"all64":rows,
              "switching_hostile":{"transitive_mask":63,"switch_vertex":1,"strong_mask":hostile},
              "path_presentation":{"H_sum":hsum,"H_squared_sum":h2sum,"unconditioned_H_mean":3,
                                   "fixed_path_H_mean":4,"fixed_cycle_H_mean":5},
              "order_five_path_counts":h5,"checks_passed":len(CHECKS),
              "check_families":["independent path DP and permutation counts","direct and degree triangle counts",
                                "Leibniz determinant vs Pfaffian","all64 and all24 relabellings",
                                "all8 cuts and all4096 gauge pairs","named positive and hostile controls"]}
    source.with_suffix('.json').write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":"PASS","checks_passed":len(CHECKS),"all64_types":by_type,
                      "output":source.with_suffix('.json').name}))


if __name__ == '__main__':
    main()
