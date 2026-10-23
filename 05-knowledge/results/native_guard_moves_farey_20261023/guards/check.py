"""Native merged-metagraph harmonic characters and actual path-clock footprints.

All arithmetic is exact; inputs retain their original, different mask gauges.
"""
from array import array
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd, lcm
from pathlib import Path
import csv
import gzip
import json
import sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
CHECKS = 0


def need(ok, detail):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(detail)


def small_atlas(n):
    path = RESULTS / "merged_metagraph_black_quotient_bridges_20261020" / f"census_n{n}.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    data["input_sha256"] = sha256(raw).hexdigest()
    data["chords"] = [(n-x, n-y) for x, y in data["tiles"]]
    return data


def literal_laplacian(data):
    N = len(data["nodes"])
    L = [[0]*N for _ in range(N)]
    for color in ("blue_edges", "black_edges"):
        for i, j, weight in data[color]:
            if i != j:
                L[i][i] += weight
                L[j][j] += weight
                L[i][j] -= weight
                L[j][i] -= weight
    return L


def matrix_from_small_mask(data, mask):
    n = data["n"]
    A = [[int(i < j) for j in range(n)] for i in range(n)]
    for bit, (i, j) in enumerate(data["chords"]):
        A[i][j] = 1-((mask >> bit) & 1)
        A[j][i] = 1-A[i][j]
    return A


def small_path_mask(data, A, path):
    return sum((1-A[path[i]][path[j]]) << bit for bit, (i, j) in enumerate(data["chords"]))


def phase(data, mask, potential):
    """Requires authenticated mask in the small atlas's backward-tile gauge."""
    full = len(data["pi"])-1
    return (potential[data["pi"][mask]]-potential[data["pi"][full ^ mask]]) % 8


def inverse(A):
    n = len(A)
    B = [[Fraction(x) for x in row]+[Fraction(i == j) for j in range(n)]
         for i, row in enumerate(A)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if B[i][k])
        B[k], B[pivot] = B[pivot], B[k]
        d = B[k][k]
        B[k] = [x/d for x in B[k]]
        for i in range(n):
            if i != k and B[i][k]:
                d = B[i][k]
                B[i] = [x-d*y for x, y in zip(B[i], B[k])]
    return [row[n:] for row in B]


def inverse_characters(data):
    L = literal_laplacian(data)
    inv = inverse([row[1:] for row in L[1:]])
    bank, certificates = [], []
    for j in range(len(inv)):
        den = lcm(*(row[j].denominator for row in inv))
        current = den//gcd(den, 8)
        lift = [Fraction(0)]+[8*current*row[j] for row in inv]
        need(all(x.denominator == 1 for x in lift), "integral inverse-column character lift")
        x = [int(v) for v in lift]
        image = [sum(a*b for a, b in zip(row, x)) for row in L]
        expected = [0]*len(L)
        expected[0] = -8*current
        expected[j+1] = 8*current
        need(image == expected, "exact primal current certificate for harmonic character")
        residue = tuple(v % 8 for v in x)
        bank.append(residue)
        certificates.append(dict(current_node=j+1, inverse_column_denominator=den,
                                 current_multiplier=current, lifted_potential=x,
                                 potential_mod8=residue))
    return bank, certificates


def minimal_footprints(paths, values):
    n = len(paths[0])
    conflicts = {sum((a != b) << i for i, (a, b) in enumerate(zip(paths[s], paths[t])))
                 for s in range(len(paths)) for t in range(s) if values[s] != values[t]}
    passing = [mask for mask in range(1 << n) if all(mask & difference for difference in conflicts)]
    cost = min(mask.bit_count() for mask in passing)
    best = [[i for i in range(n) if mask >> i & 1] for mask in passing if mask.bit_count() == cost]
    witnesses = []
    if cost:
        for clocks in combinations(range(n), cost-1):
            found = next((dict(clocks=clocks, path_a=paths[s], path_b=paths[t],
                               value_a=values[s], value_b=values[t])
                          for s in range(len(paths)) for t in range(s)
                          if values[s] != values[t] and all(paths[s][i] == paths[t][i] for i in clocks)), None)
            need(found is not None, "literal witness against every smaller fixed-clock bank")
            witnesses.append(found)
    for clocks in best:
        lookup = {}
        for path, value in zip(paths, values):
            answer = tuple(path[i] for i in clocks)
            need(answer not in lookup or lookup[answer] == value, "direct upper footprint table")
            lookup[answer] = value
    flag_masks = [(mask, sum(1 << i for i in range(n)
                            if (i and mask >> (i-1) & 1) or (i < n-1 and mask >> i & 1)))
                  for mask in range(1 << (n-1))]
    flags = [mask for mask, exposed in flag_masks if all(exposed & difference for difference in conflicts)]
    flag_cost = min(mask.bit_count() for mask in flags)
    endpoints = (1 << 0) | (1 << (n-1))
    hostile = next(((list(paths[s]), list(paths[t]), values[s], values[t])
                    for s in range(len(paths)) for t in range(s)
                    if values[s] != values[t] and paths[s][0] == paths[t][0]
                    and paths[s][-1] == paths[t][-1]), None)
    need((hostile is None) == all(endpoints & difference for difference in conflicts),
         "literal endpoint conflict and footprint test agree")
    return dict(minimum_owner_clocks=cost, all_minimum_owner_clock_sets=best,
                lower_bound_witnesses=witnesses,
                minimum_arc_flags=flag_cost,
                all_minimum_arc_flag_sets=[[i for i in range(n-1) if mask >> i & 1]
                                           for mask in flags if mask.bit_count() == flag_cost],
                endpoints_determine=hostile is None, endpoint_hostile=hostile)


def small_controls(n):
    data = small_atlas(n)
    bank, certificates = inverse_characters(data)
    generators = sorted(set(bank))
    N = len(data["nodes"])
    span = {tuple([0]*N)}
    pending = list(span)
    while pending:
        x = pending.pop()
        for y in generators:
            z = tuple((a+b) % 8 for a,b in zip(x,y))
            if z not in span:
                span.add(z)
                pending.append(z)
    # Independent complete modular kernel, lifting every binary correction.
    L = literal_laplacian(data)
    kernel = {tuple([0]*N)}
    modulus = 1
    level_sizes = []
    for _ in range(3):
        lifted = set()
        for x in kernel:
            for bits in product((0,1), repeat=N-1):
                z = (0,)+tuple(a+modulus*b for a,b in zip(x[1:],bits))
                if all(sum(a*b for a,b in zip(row,z)) % (2*modulus) == 0 for row in L):
                    lifted.add(z)
        kernel = lifted
        modulus *= 2
        level_sizes.append(len(kernel))
    need(span == kernel, "inverse-current bank generates full independently lifted mod8 kernel")
    bank = sorted(kernel)
    full = len(data["pi"])-1
    for x in bank:
        currents = [0]*len(data["nodes"])
        for mask, source in enumerate(data["pi"]):
            currents[source] += phase(data, mask, x)
            need((phase(data, mask, x)+phase(data, mask ^ full, x)) % 8 == 0,
                 "actual complement reverses the port phase")
        need(all(v % 8 == 0 for v in currents), "literal directed port conservation at every merged node")
    reports = []
    for row in data["classes"]:
        A = matrix_from_small_mask(data, row["representative"])
        perms = list(permutations(range(n)))
        paths = [p for p in perms if all(A[i][j] for i, j in zip(p, p[1:]))]
        need(len(paths) == row["H"], "complete actual Hamiltonian path count")
        masks = [small_path_mask(data, A, p) for p in paths]
        need(all(data["class_by_mask"][t] == row["class"] for t in masks),
             "actual source class remains authenticated")
        aut = [g for g in perms if all(A[i][j] == A[g[i]][g[j]] for i in range(n) for j in range(n))]
        need(len(aut) == row["aut"], "literal full automorphism group")
        values = [tuple(phase(data, t, x) for x in bank) for t in masks]
        targets = [data["pi"][full ^ t] for t in masks]
        lookup = dict(zip(paths, zip(values, targets)))
        for p, (value, target) in lookup.items():
            for g in aut:
                gp = tuple(g[v] for v in p)
                need(lookup[gp] == (value, target), "actual Aut source orbit has zero phase footprint")
            opposite = [[A[j][i] for j in range(n)] for i in range(n)]
            reverse = tuple(reversed(p))
            t = small_path_mask(data, opposite, reverse)
            need(tuple(phase(data, t, x) for x in bank) == value, "converse reversal retains merged phase")
        report = dict(ordinary_class=row["class"], representative_mask=row["representative"],
                      merged_source=data["pi"][row["representative"]], H=len(paths), aut=len(aut),
                      phase_values=len(set(values)), target_owners=len(set(targets)),
                      phase_footprint=minimal_footprints(paths, values),
                      target_footprint=minimal_footprints(paths, targets))
        reports.append(report)
    return dict(n=n, atlas_sha256=data["input_sha256"],
                literal_merged_nodes=len(data["nodes"]), character_bank_mod8=bank,
                full_grounded_kernel_sizes_mod2_mod4_mod8=level_sizes,
                primal_character_certificates=certificates, sources=reports)


def order8_controls():
    folder = RESULTS / "metagraph_order8_20261023"
    potential_file = RESULTS / "metagraph_order8_electrical_20261023/receipt.json"
    raw = potential_file.read_bytes()
    inherited = json.loads(raw)["order8"]
    x = inherited["higher_layers"][-1]["lifts"][0]
    need(len(x) == 3528 and inherited["root"] == 3527 and x[3527] == 0,
         "original order8 potential and original node gauge")
    nodes = list(csv.DictReader((folder / "merged_nodes.tsv").open(), delimiter="\t"))
    mid = [None]*6880
    for row in nodes:
        for key in ("class_a", "class_b"):
            mid[int(row[key])] = int(row["merged_id"])
    raw_owner = gzip.decompress((folder / "owner_map.u16.gz").read_bytes())
    owner = array("H")
    owner.frombytes(raw_owner)
    if sys.byteorder != "little":
        owner.byteswap()
    need(len(owner) == 1 << 21 and all(v is not None for v in mid), "native8 owner map")
    edge_raw = gzip.decompress((folder / "merged_edges.tsv.gz").read_bytes())
    image = [0]*len(x)
    frozen_relabel_image = [0]*len(x)
    transported_relabel_image = [0]*len(x)
    perm = list(range(len(x)))
    perm[1238], perm[2242] = 2242, 1238
    transported = [0]*len(x)
    for i in range(len(x)):
        transported[perm[i]] = x[i]
    selected_edges = {}
    for row in csv.DictReader(edge_raw.decode().splitlines(), delimiter="\t"):
        i, j, blue, black = [int(row[k]) for k in ("source", "target", "blue", "black")]
        weight = blue+black
        need(black % 2 == 0, "original literal paired-black multiplicity")
        d = weight*(x[i]-x[j])
        image[i] += d
        image[j] -= d
        a, b = perm[i], perm[j]
        d = weight*(x[a]-x[b])
        frozen_relabel_image[a] += d
        frozen_relabel_image[b] -= d
        d = weight*(transported[a]-transported[b])
        transported_relabel_image[a] += d
        transported_relabel_image[b] -= d
        if i in (1238, 2242) and j in (2246, 2315):
            selected_edges[i, j] = black
    need(all(v % 256 == 0 for v in image), "original full native harmonic equations modulo256")
    need(all(v % 8 == 0 for v in transported_relabel_image), "transported character is relabeling invariant")
    need(any(v % 8 for v in frozen_relabel_image), "frozen character is not invariant under arbitrary node relabeling")
    terminals = (1238, 2242, 2246, 2315)
    need([x[i] for i in terminals] == [97, 18, 65, 212], "original four terminal lift coordinates")
    delta = [2*(x[2246]-x[2315]), -2*(x[2246]-x[2315]),
             2*(x[1238]-x[2242]), -2*(x[1238]-x[2242])]
    need([v % 8 for v in delta] == [2, 6, 6, 2], "original mod8 switch defect")
    for t in (-1, 0, 1, 2):
        need(all(t*v % 8 == 0 for v in delta) == (t == 0), "all declared feasible switches exactly rejected")
    pairs = list(combinations(range(8), 2))
    chords = [(i, j) for i in range(8) for j in range(i+2, 8)]
    full = (1 << 21)-1
    reports = []
    for source, expected_H in zip(terminals, (409,407,523,527)):
        row = nodes[source]
        canon = int(row["canon_a"])
        ci = int(row["class_a"])
        A = [[0]*8 for _ in range(8)]
        for bit, (i, j) in enumerate(pairs):
            A[i][j] = (canon >> (27-bit)) & 1
            A[j][i] = 1-A[i][j]
        paths = [p for p in permutations(range(8)) if all(A[i][j] for i, j in zip(p, p[1:]))]
        need(len(paths) == expected_H, "original selected source full Hamiltonian count")
        masks = [sum(A[p[i]][p[j]] << bit for bit, (i, j) in enumerate(chords)) for p in paths]
        need(all(owner[t] == ci for t in masks), "original8 forward-bit source gauge authenticated")
        targets = [mid[owner[full ^ t]] for t in masks]
        values = [(x[source]-x[target]) % 8 for target in targets]
        aut = [g for g in permutations(range(8)) if all(A[i][j] == A[g[i]][g[j]] for i, j in pairs)]
        lookup = dict(zip(paths, zip(values, targets)))
        for p in paths:
            for g in aut:
                need(lookup[tuple(g[v] for v in p)] == lookup[p], "native8 zero phase footprint on actual aut orbit")
        phase_footprint = minimal_footprints(paths, values)
        target_footprint = minimal_footprints(paths, targets)
        expected_banks = {tuple(k for k in range(8) if k not in (i,j))
                          for i,j in combinations(range(8),2) if j-i <= 2}
        need(set(values) == set(range(8)), "original character attains all eight residues")
        need(phase_footprint["minimum_owner_clocks"] == target_footprint["minimum_owner_clocks"] == 6,
             "original phase and target have exact six-clock minima")
        need({tuple(s) for s in phase_footprint["all_minimum_owner_clock_sets"]} == expected_banks,
             "every minimum original phase bank has close two-hole structure")
        need(phase_footprint["all_minimum_arc_flag_sets"] == [[0,2,4],[0,2,6],[0,4,6],[2,4,6]],
             "exact four minimum original phase flag banks")
        reports.append(dict(merged_source=source, ordinary_source=ci, canon=canon, H=len(paths), aut=len(aut),
                            phase_histogram=dict(sorted(Counter(values).items())),
                            target_owners=len(set(targets)), phase_footprint=phase_footprint,
                            target_footprint=target_footprint))
    padded_A = [[0]*8 for _ in range(8)]
    for i in range(8):
        for j in range(i+1, 8):
            if i == 0 or j == 7:
                value = 1
            else:
                ci, u = divmod(i-1, 3)
                cj, v = divmod(j-1, 3)
                value = int((v-u) % 3 == 1) if ci == cj else int((v-u) % 3 in (1, 2))
            padded_A[i][j] = value
            padded_A[j][i] = 1-value
    padded_rows = []
    for a, b in product(range(3), repeat=2):
        path = (0,)+tuple(1+(a+i) % 3 for i in range(3))+tuple(4+(b+i) % 3 for i in range(3))+(7,)
        valid = all(padded_A[u][v] for u, v in zip(path, path[1:]))
        need(valid == ((b-a) % 3 in (0, 1)), "actual padded-cell junction with fixed endpoint owners")
        mask = sum(padded_A[path[i]][path[j]] << bit for bit, (i, j) in enumerate(chords)) if valid else None
        source = mid[owner[mask]] if valid else None
        target = mid[owner[full ^ mask]] if valid else None
        if valid:
            cell_difference = (b-a) % 3
            need(source == 3517 and target == (532 if cell_difference == 0 else 1290)
                 and (x[source]-x[target]) % 8 == (1 if cell_difference == 0 else 5),
                 "pinned padded source retains original nonconstant native character")
        padded_rows.append(dict(a=a,b=b,path=path,valid=int(valid),normalized_forward_mask=mask,
                                source_merged=source,target_merged=target,
                                phase_mod8=(x[source]-x[target]) % 8 if valid else None))
    return dict(original_potential_sha256=sha256(raw).hexdigest(), owner_bytes_sha256=sha256(raw_owner).hexdigest(),
                edge_bytes_sha256=sha256(edge_raw).hexdigest(), root=3527, original_potential_modulus=256,
                used_phase_modulus=8, terminals=terminals, terminal_residues=[x[i] % 8 for i in terminals],
                switch_defect_mod8=[v % 8 for v in delta], four_black_weights={str(k):v for k,v in selected_edges.items()},
                frozen_relabel_nonzero_defect_count=sum(v % 8 != 0 for v in frozen_relabel_image), sources=reports,
                pinned_padded_two_cells=dict(adjacency=padded_A,actual_pins=[[0,0],[7,7]],rows=padded_rows))


def two_cyclic_cells():
    data = small_atlas(6)
    bank, certificates = inverse_characters(data)
    A = [[0]*6 for _ in range(6)]
    for i in range(6):
        for j in range(i+1, 6):
            ci, u = divmod(i, 3)
            cj, v = divmod(j, 3)
            A[i][j] = int((v-u) % 3 == 1) if ci == cj else int((v-u) % 3 in (1, 2))
            A[j][i] = 1-A[i][j]
    rows = []
    for a, b in product(range(3), repeat=2):
        path = tuple((a+i) % 3 for i in range(3))+tuple(3+(b+i) % 3 for i in range(3))
        valid = all(A[u][v] for u, v in zip(path, path[1:]))
        need(valid == ((b-a) % 3 in (0, 1)), "actual two-cell junction guard")
        mask = small_path_mask(data, A, path) if valid else None
        rows.append(dict(a=a, b=b, path=path, valid=int(valid), normalized_mask=mask,
                         source_ordinary=data["class_by_mask"][mask] if valid else None,
                         source_merged=data["pi"][mask] if valid else None,
                         target_merged=data["pi"][(len(data["pi"])-1) ^ mask] if valid else None))
    best = max(range(len(bank)), key=lambda j:len({phase(data, row["normalized_mask"], bank[j])
                                                  for row in rows if row["valid"]}))
    x = bank[best]
    for row in rows:
        row["phase_mod8"] = phase(data, row["normalized_mask"], x) if row["valid"] else None
        row["valid_times_phase_mod8"] = row["phase_mod8"] if row["valid"] else 0
        if row["valid"]:
            delta = (row["b"]-row["a"]) % 3
            need(row["source_merged"] == 25 and row["target_merged"] == (32 if delta == 0 else 11)
                 and row["phase_mod8"] == 4*delta, "two actual source electrical sectors")
    need(len({row["source_ordinary"] for row in rows if row["valid"]}) == 1,
         "all valid block phases share the exact native source")
    residue_image = [sum(a*b for a, b in zip(row, x)) for row in literal_laplacian(data)]
    need(all(v % 8 == 0 for v in residue_image), "short residue primal certificate")
    def iso(U, V):
        return next((p for p in permutations(range(6)) if all(U[i][j] == V[p[i]][p[j]]
                     for i in range(6) for j in range(6))), None)
    actual_targets = []
    for row in rows:
        if not row["valid"]:
            continue
        path = row["path"]
        arcs = set(zip(path, path[1:]))
        target = [[A[i][j] if (i,j) in arcs or (j,i) in arcs else (i != j)-A[i][j]
                   for j in range(6)] for i in range(6)]
        target_mask = (len(data["pi"])-1) ^ row["normalized_mask"]
        target_class = data["class_by_mask"][target_mask]
        representative = data["classes"][target_class]["representative"]
        witness = iso(target, matrix_from_small_mask(data, representative))
        need(witness is not None, "independent actual target-isomorphism witness")
        row["target_ordinary"] = target_class
        row["target_isomorphism_to_atlas_representative"] = witness
        if row["a"] == 0:
            actual_targets.append(target)
    need(iso(actual_targets[0], actual_targets[1]) is None
         and iso(actual_targets[0], list(map(list, zip(*actual_targets[1])))) is None,
         "two native electrical sectors have different merged target owners")
    automorphisms = [p for p in permutations(range(6)) if all(A[i][j] == A[p[i]][p[j]]
                    for i in range(6) for j in range(6))]
    need(len(automorphisms) == 3, "literal two-cell source symmetry")
    return dict(n=6, cells=2, regular_source=False, adjacency=A, atlas_sha256=data["input_sha256"],
                potential_mod8=x, primal_character_certificate=certificates[best],
                residue_primal_current_div8=[v//8 for v in residue_image], actual_source_automorphisms=automorphisms,
                rows=rows, distinct_valid_targets=sorted({row["target_merged"] for row in rows if row["valid"]}),
                distinct_valid_phases=sorted({row["phase_mod8"] for row in rows if row["valid"]}))


def paired_slot_controls(rows, shared_clocks):
    lookup = {(row["a"], row["b"]):row for row in rows}
    first, second = lookup[0,0], lookup[0,1]
    target_bank = sorted((first["target_merged"], second["target_merged"]))
    states = []
    for b, c in product(range(3), repeat=2):
        left, right = lookup[0,b], lookup[0,c]
        valid = left["valid"] and right["valid"]
        accepted = valid and sorted((left["target_merged"], right["target_merged"])) == target_bank
        if accepted:
            need(all(left["path"][j] == right["path"][j] == first["path"][j] for j in shared_clocks),
                 "literal shared owner-clock pins survive coupled assignment")
            need((left["phase_mod8"]+right["phase_mod8"]) % 8
                 == (first["phase_mod8"]+second["phase_mod8"]) % 8,
                 "actual aggregate electrical contribution preserved")
            states.append(dict(b=b,c=c,left_path=left["path"],right_path=right["path"],
                               ordered_targets=[left["target_merged"],right["target_merged"]],
                               ordered_phases=[left["phase_mod8"],right["phase_mod8"]]))
    need([(s["b"],s["c"]) for s in states] == [(0,1),(1,0)], "two exact aggregate-guard assignments")
    need(sum(states[0][key] != states[1][key] for key in ("b","c")) == 2,
         "sharp budget two in actual phase-coordinate menu")
    need(states[0]["ordered_phases"] != states[1]["ordered_phases"],
         "per-port phase pins correctly reject the aggregate swap")
    return dict(shared_owner_clock_pins=[[j,first["path"][j]] for j in shared_clocks],
                fixed_first_cell_phase=0,target_multiset=target_bank,
                phase_multiset=sorted([first["phase_mod8"],second["phase_mod8"]]),
                minimum_phase_coordinates_changed=2,states=states,
                scope="labelled observer-slot reassignment; same unlabelled path and line multiset")


def main():
    reports = [small_controls(n) for n in range(3, 6)]
    eight = order8_controls()
    two = two_cyclic_cells()
    paired = dict(native6=paired_slot_controls(two["rows"], (0,1,2)),
                  original8_pinned=paired_slot_controls(eight["pinned_padded_two_cells"]["rows"], (0,1,2,3,7)))
    result = dict(status="FINITE-EXACT original native electrical characters", small_native=reports,
                  original_eight=eight, two_cyclic_cells=two, paired_observer_slots=paired, checks=CHECKS)
    text = json.dumps(result, indent=2, sort_keys=True)+"\n"
    target = HERE / "receipt.json"
    if "--write" in sys.argv:
        target.write_text(text, encoding="utf-8")
    else:
        need(json.loads(target.read_text(encoding="utf-8")) == json.loads(text), "frozen native guard receipt")
    print(json.dumps(dict(status="PASS", checks=CHECKS,
                          clock_costs={r["n"]:Counter(s["phase_footprint"]["minimum_owner_clocks"] for s in r["sources"])
                                       for r in reports},
                          original8_clock_costs=[s["phase_footprint"]["minimum_owner_clocks"] for s in eight["sources"]])))


if __name__ == "__main__":
    main()
