#!/usr/bin/env python3
"""
LRC(14) dihedral endpoint atlas for OPEN-Q-108.

This script refines the uniform-fattening gauntlet by keeping the endpoint
geometry that the speed-load tournament deliberately forgot.

Local dihedral principle:
  For a speed v, the danger set D_v={t: ||v t|| <= 1/14} is the pullback of
  the base moat under t -> v*t.  Its endpoints are

      (14*k - 1)/(14*v), (14*k + 1)/(14*v)      (mod 1),

  and they are permuted by the dihedral clock action on the v-gon:
  rotations k -> k+r and reflection k -> -k, left <-> right.  For a whole
  core C, the common rotation usually dies, but the global reflection
  t -> 1-t survives and pairs safe components.

Endpoint-mouth formula:
  If a safe component begins at the right endpoint of speed a, center k/a,
  and ends at the left endpoint of speed b, center l/b, then

      length = (a*(14*l - 1) - b*(14*k + 1)) / (14*a*b).

  The numerator is the dihedral mouth determinant.  OPEN-Q-108 can be read as
  a lower bound on the total surviving determinant mass after these local
  clock mouths are paired by reflection.

Tournament Analysis contract:
  * Primary vertices: reflection orbits of safe components ("dihedral mouths").
  * Pairwise observable: total safe measure contributed by a mouth orbit.
  * Switch/gauge: orient orbit A -> B when A contributes more measure; ties use
    the left-endpoint order as the Hamiltonian path.
  * Secondary comparison: edge flips against the raw determinant-sum gauge,
    which forgets denominator normalization.
  * Fingerprints: score histograms, directed cycles, SCC sizes, Hamiltonian
    path counts, and measure-vs-determinant edge flips.

Assumption challenge:
  This session explicitly rejects the default "vertices are runners" quotient
  for this proof step.  Alternatives considered were runners, omitted AP
  speeds, safe components, reflected component orbits, endpoint events,
  omitted-clock teeth, residues mod 14, wall-crossing events, and proof
  obligations.  The chosen vertex set is reflected safe-component orbits: it
  preserves exact measure and endpoint adjacency, but destroys most runner
  identity except the two boundary owners.  The retained LRC predicate is the
  exact safe measure meas(G_C); the destroyed data is the full cover history
  inside danger arcs away from the component boundaries.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations

import lrc14_uniform_fattening_gauntlet_codex as base


N = 14
AP13 = tuple(range(1, N))
DROP6 = tuple(v for v in AP13 if v != 6)


@dataclass(frozen=True, order=True)
class Endpoint:
    x: F
    v: int
    side: str
    k: int


def mod1(x):
    return x % 1


def endpoint_events(v, q=N):
    events = []
    for k in range(v):
        left = mod1(F(q * k - 1, q * v))
        right = mod1(F(q * k + 1, q * v))
        events.append(Endpoint(left, v, "L", k))
        events.append(Endpoint(right, v, "R", k))
    return events


def endpoint_map(core, q=N):
    by_x = defaultdict(list)
    for v in core:
        for event in endpoint_events(v, q):
            by_x[event.x].append(event)
    return {x: tuple(sorted(events)) for x, events in by_x.items()}


def reflect_endpoint(event):
    return Endpoint(
        mod1(1 - event.x),
        event.v,
        "L" if event.side == "R" else "R",
        (-event.k) % event.v,
    )


def reflect_interval(interval):
    a, b = interval
    return (mod1(1 - b), mod1(1 - a))


def event_label(event):
    return f"{event.v}{event.side}{event.k}@{event.x}"


def component_records(core, q=N):
    core = tuple(sorted(core))
    _, comps = base.safe_components(core, q)
    events = endpoint_map(core, q)
    records = []
    for idx, (a, b) in enumerate(comps):
        left_events = events.get(a, ())
        right_events = events.get(b, ())
        left_candidates = tuple(e for e in left_events if e.side == "R") or left_events
        right_candidates = tuple(e for e in right_events if e.side == "L") or right_events
        mouths = []
        for left in left_candidates:
            for right in right_candidates:
                det = (b - a) * q * left.v * right.v
                mouths.append((left, right, det))
        records.append(
            {
                "idx": idx,
                "interval": (a, b),
                "length": b - a,
                "left_events": left_events,
                "right_events": right_events,
                "mouths": tuple(mouths),
            }
        )
    return records


def reflection_orbits(records):
    index = {rec["interval"]: i for i, rec in enumerate(records)}
    seen = set()
    orbits = []
    for i, rec in enumerate(records):
        if i in seen:
            continue
        refl = reflect_interval(rec["interval"])
        j = index.get(refl)
        members = [i]
        if j is not None and j != i:
            members.append(j)
        for member in members:
            seen.add(member)
        members = tuple(sorted(members))
        measure = sum(records[m]["length"] for m in members)
        det_sum = sum(records[m]["mouths"][0][2] for m in members if records[m]["mouths"])
        left_order = min(records[m]["interval"][0] for m in members)
        orbits.append(
            {
                "members": members,
                "measure": measure,
                "det_sum": det_sum,
                "left_order": left_order,
            }
        )
    return sorted(orbits, key=lambda item: item["left_order"])


def orbit_tournament(orbits, key):
    n = len(orbits)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            ai = orbits[i][key]
            aj = orbits[j][key]
            if ai > aj or (ai == aj and orbits[i]["left_order"] < orbits[j]["left_order"]):
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


def tournament_fingerprint(adj):
    scores = base.tournament_scores(adj)
    cycles, transitive = base.directed_triangles(adj)
    return {
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_cycles": cycles,
        "transitive_triples": transitive,
        "scc_sizes": sorted((len(c) for c in base.sccs(adj)), reverse=True),
        "hamiltonian_paths": base.hamiltonian_path_count(adj),
    }


def missing_arc_hits(t, missing, q=N):
    hits = []
    for m in missing:
        ks = []
        for k in range(m):
            center = F(k, m)
            dist = abs(t - center)
            dist = min(dist, 1 - dist)
            if dist <= F(1, q * m):
                ks.append(k)
        if ks:
            hits.append((m, tuple(ks)))
    return tuple(hits)


def ap_delete_rows():
    rows = []
    for missing in AP13:
        core = tuple(v for v in AP13 if v != missing)
        measure, comps = base.safe_components(core, N)
        records = component_records(core, N)
        orbits = reflection_orbits(records)
        shadow = defaultdict(F)
        for rec in records:
            a, b = rec["interval"]
            mid = (a + b) / 2
            for hit in missing_arc_hits(mid, (missing,), N):
                shadow[hit] += rec["length"]
        owner_path = []
        for rec in records:
            if rec["mouths"]:
                left, right, _ = rec["mouths"][0]
                owner_path.append(f"{left.v}-{right.v}")
            else:
                owner_path.append("?")
        rows.append(
            {
                "missing": missing,
                "core": core,
                "measure": measure,
                "components": len(comps),
                "orbits": len(orbits),
                "shadow": dict(sorted(shadow.items())),
                "owner_path": owner_path,
                "orbit_records": orbits,
                "component_records": records,
            }
        )
    return sorted(rows, key=lambda row: row["measure"])


def fmt(fr):
    return f"{fr} = {float(fr):.9f}"


def print_drop6_decomposition():
    records = component_records(DROP6, N)
    orbits = reflection_orbits(records)
    print("Drop-6 detailed dihedral mouth decomposition")
    print(f"  core={DROP6}")
    print(f"  total meas(G_C)={fmt(sum(rec['length'] for rec in records))}")
    for oi, orbit in enumerate(orbits):
        print(
            f"  orbit {oi}: members={orbit['members']}; "
            f"orbit_measure={fmt(orbit['measure'])}; det_sum={orbit['det_sum']}"
        )
        for member in orbit["members"]:
            rec = records[member]
            left, right, det = rec["mouths"][0]
            print(
                "    component {}: {} -> {}; interval=({}, {}); "
                "length={}; determinant={}".format(
                    member,
                    event_label(left),
                    event_label(right),
                    rec["interval"][0],
                    rec["interval"][1],
                    rec["length"],
                    det,
                )
            )
    adj_measure = orbit_tournament(orbits, "measure")
    adj_det = orbit_tournament(orbits, "det_sum")
    print("  mouth-orbit tournament by measure:")
    print(f"    {tournament_fingerprint(adj_measure)}")
    print(
        "  measure-vs-raw-determinant edge flips="
        f"{base.edge_flips(adj_measure, adj_det)}"
    )


def perturbation_trade_scan(limit=180):
    drop_measure, drop_comps = base.safe_components(DROP6, N)
    best = None
    best_with_loss = None
    max_loss = None
    min_surplus = None
    count_missing6 = 0
    count_loss = 0
    for missing in combinations(AP13, 2):
        stem = [v for v in AP13 if v not in missing]
        for w in range(N, limit + 1):
            if w in stem:
                continue
            core = base.primitive(stem + [w])
            if len(core) != 12:
                continue
            measure, comps = base.safe_components(core, N)
            retained = base.intersect_measure(comps, drop_comps)
            lost = drop_measure - retained
            outside = measure - retained
            surplus = outside - lost
            rec = {
                "measure": measure,
                "missing": missing,
                "replacement": w,
                "core": core,
                "retained": retained,
                "lost": lost,
                "outside": outside,
                "surplus": surplus,
            }
            if best is None or measure < best["measure"]:
                best = rec
            if 6 in missing:
                count_missing6 += 1
                if lost > 0:
                    count_loss += 1
                    if best_with_loss is None or measure < best_with_loss["measure"]:
                        best_with_loss = rec
                    if max_loss is None or lost > max_loss["lost"]:
                        max_loss = rec
                if min_surplus is None or surplus < min_surplus["surplus"]:
                    min_surplus = rec
    return {
        "drop6_measure": drop_measure,
        "limit": limit,
        "best": best,
        "count_missing6": count_missing6,
        "count_loss": count_loss,
        "best_with_loss": best_with_loss,
        "max_loss": max_loss,
        "min_surplus": min_surplus,
    }


def print_trade_record(label, rec):
    if rec is None:
        print(f"  {label}: none")
        return
    print(
        f"  {label}: missing={rec['missing']}, replacement={rec['replacement']}, "
        f"meas={fmt(rec['measure'])}"
    )
    print(
        "    retained_drop6={}, lost_drop6={}, new_outside={}, "
        "surplus_over_drop6={}".format(
            rec["retained"], rec["lost"], rec["outside"], rec["surplus"]
        )
    )
    print(f"    core={rec['core']}")


def print_ap_delete_atlas():
    rows = ap_delete_rows()
    print("AP13 deletion atlas through dihedral mouths")
    for row in rows:
        shadow_items = ", ".join(
            f"{hit}:{length}" for hit, length in row["shadow"].items()
        )
        print(
            f"  delete {row['missing']:2d}: meas={fmt(row['measure'])}; "
            f"components={row['components']}; reflection_orbits={row['orbits']}; "
            f"owner_path={'/'.join(row['owner_path'])}"
        )
        print(f"    omitted-clock shadow: {shadow_items}")
    measures = [{"measure": row["measure"], "left_order": F(row["missing"], 1)} for row in rows]
    adj = orbit_tournament(measures, "measure")
    print("  deletion-compression tournament (smaller safe measure ranks higher after reversal):")
    # orbit_tournament points larger -> smaller; reverse all edges for compression ranking.
    rev_adj = [[False] * len(adj) for _ in range(len(adj))]
    for i in range(len(adj)):
        for j in range(len(adj)):
            if i != j:
                rev_adj[i][j] = adj[j][i]
    print(f"    order_by_measure={[row['missing'] for row in rows]}")
    print(f"    {tournament_fingerprint(rev_adj)}")


def sample_core_summary(name, core):
    measure, _ = base.safe_components(tuple(sorted(core)), N)
    records = component_records(core, N)
    orbits = reflection_orbits(records)
    print(f"{name}: core={tuple(sorted(core))}")
    print(f"  measure={fmt(measure)}; components={len(records)}; reflection_orbits={len(orbits)}")
    for orbit in orbits:
        labels = []
        for member in orbit["members"]:
            rec = records[member]
            if rec["mouths"]:
                left, right, det = rec["mouths"][0]
                labels.append(f"{left.v}{left.side}->{right.v}{right.side}:det{det}")
        print(
            f"  orbit members={orbit['members']}; measure={orbit['measure']}; "
            f"mouths={', '.join(labels)}"
        )


def main():
    print("=" * 78)
    print("LRC(14) dihedral endpoint atlas for OPEN-Q-108")
    print("=" * 78)
    print("Endpoint formula: R(a,k) -> L(b,l) has length")
    print("  (a*(14*l - 1) - b*(14*k + 1))/(14*a*b).")
    print("The global reflection t -> 1-t pairs every safe component orbit.")
    print()

    print_ap_delete_atlas()
    print()
    print_drop6_decomposition()
    print()

    print("Two-delete/one-replacement trade around the drop-6 mouths")
    trade = perturbation_trade_scan(limit=180)
    print(f"  drop6 base measure={fmt(trade['drop6_measure'])}")
    print(
        f"  rows with missing speed 6 tested={trade['count_missing6']}; "
        f"rows that damage old drop6 mouths={trade['count_loss']}"
    )
    print_trade_record("best total row", trade["best"])
    print_trade_record("minimum surplus over drop6 among missing-6 rows", trade["min_surplus"])
    print_trade_record("best row that actually damages drop6 mouths", trade["best_with_loss"])
    print_trade_record("largest damage to old drop6 mouths", trade["max_loss"])
    print()

    print("Selected non-AP / perturbed decompositions")
    sample_core_summary(
        "best two-delete/one-replacement",
        (1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 20),
    )
    sample_core_summary(
        "damage-heavy replacement w=69",
        (1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 69),
    )
    sample_core_summary(
        "sporadic AP-neighbor",
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13),
    )
    print()

    print("Conclusion")
    print("  The current worst core is not just AP13 minus 6; it is a D6")
    print("  missing-clock shadow.  Its measure is the sum of two reflected")
    print("  endpoint-mouth orbits, both inside the omitted speed-6 moat:")
    print("    2*(1/728) + 2*(5/1848) = 7/858.")
    print("  In the tested two-delete/one-replacement family, any attempt to")
    print("  damage those old hexagon mouths forces enough new reflected mouth")
    print("  mass elsewhere; the smallest surplus over drop6 is 1/980 and occurs")
    print("  without damaging the old mouths.  This suggests the next proof target:")
    print("  a dihedral mouth-exchange inequality for coordinated growth.")


if __name__ == "__main__":
    main()
