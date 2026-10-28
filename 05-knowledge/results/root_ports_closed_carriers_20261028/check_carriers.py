#!/usr/bin/env python3
"""Exact individualization/refinement with ALL attaining transports retained.

This is an implementation/measurement of a standard canonization mechanism,
not a new graph-isomorphism theorem or a polynomial all-order complexity claim.
"""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import importlib.util
import json
import argparse
import gzip
import hashlib

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('s87', HERE.parent / 'session87_four_vertices_enumeration_20261028/check_generation.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def stable(rows, cells=None):
    if cells is None:
        degree = [r.bit_count() for r in rows]
        cells = tuple(tuple(v for v, d in enumerate(degree) if d == k) for k in sorted(set(degree)))
    while True:
        masks = [sum(1 << v for v in cell) for cell in cells]
        result = []
        for cell in cells:
            groups = {}
            for v in cell:
                sig = tuple((rows[v] & mask).bit_count() for mask in masks)
                groups.setdefault(sig, []).append(v)
            result.extend(tuple(groups[s]) for s in sorted(groups))
        result = tuple(result)
        if len(result) == len(cells):
            return result
        cells = result


def refinement_orders(rows, cells=None):
    cells = stable(rows, cells)
    split = next((i for i, c in enumerate(cells) if len(c) > 1), None)
    if split is None:
        yield tuple(v for cell in cells for v in cell)
        return
    cell = cells[split]
    for v in cell:
        new = cells[:split] + ((v,), tuple(w for w in cell if w != v)) + cells[split+1:]
        yield from refinement_orders(rows, new)


class Canon(old.ParentCanonicalizer):
    def maps(self, rows):
        need(self.max_order is None or len(rows) <= self.max_order, 'parent-only size guard')
        self.requests += 1
        if rows in self.cache:
            return self.cache[rows]
        self.searches += 1
        best = None; transports = set()
        for order in refinement_orders(rows):
            self.orders_tested += 1
            code = old.encode(rows, order)
            if best is None or code < best:
                best = code; transports = {order}
            elif code == best:
                transports.add(order)
        answer = best, tuple(sorted(transports))
        self.cache[rows] = answer
        return answer


def produce(parents, n):
    engine = Canon(max_order=n)
    normalized = []
    for raw in parents:
        code, _ = engine.maps(raw)
        rows = old.matrix(n, code)
        repeated, group = engine.maps(rows)
        need(code == repeated, 'canonical form idempotence')
        normalized.append((code, rows, group))
    normalized.sort()
    prep_leaves = engine.orders_tested
    count = Counter(); children = []
    for code, rows, group in normalized:
        for mask in old.orbit_masks(n, group):
            count['pointed_mask_orbits'] += 1
            child = old.extend(rows, mask)
            # Stable refinement never reorders the initial degree cells.
            # Reject here before paying for refinement of an impossible root.
            if mask.bit_count() < max(r.bit_count() for r in child):
                count['stable_selector_rejected'] += 1
                continue
            cells = stable(child)
            selected = cells[-1]
            if n not in selected:
                count['stable_selector_rejected'] += 1
                continue
            key = code, mask
            for v in selected:
                if v == n:
                    continue
                count['competitor_key_queries'] += 1
                if engine.rooted(*old.delete_root(child, v)) < key:
                    count['complete_key_rejected'] += 1
                    break
            else:
                children.append(child)
    return children, dict(order=n+1, accepted_classes=len(children),
                          preparation_leaves=prep_leaves,
                          competitor_leaves=engine.orders_tested-prep_leaves,
                          canon_searches=engine.searches,
                          production_max_canonicalized_order=n, **dict(count))


def build_hostile(last_converse):
    rows = [0]*12
    for u in range(12):
        for v in range(u+1, 12):
            b, q = divmod(u, 4); c, r = divmod(v, 4)
            arc = (c-b) % 3 == 1 if b != c else q == 0 or q > 0 and (r-q) % 3 == 1
            if b == c == 2 and last_converse:
                arc = not arc
            rows[u if arc else v] |= 1 << (v if arc else u)
    return tuple(rows)


def actual_automorphisms(rows):
    """Independent constraint search: assign images, preserving every arc."""
    degree = [r.bit_count() for r in rows]; n = len(rows)
    result = []; p = []
    def rec():
        u = len(p)
        if u == n:
            result.append(tuple(p)); return
        for v in range(n):
            if v in p or degree[v] != degree[u]:
                continue
            if all(((rows[u] >> k) & 1) == ((rows[v] >> p[k]) & 1) for k in range(u)):
                p.append(v); rec(); p.pop()
    rec()
    return set(result)


def strong_controls():
    # Odd quadratic residues orient Paley7; T5 uses steps 1 and 2.
    atoms = [(n, tuple(sum(1 << ((v+d) % n) for d in ds) for v in range(n)))
             for n, ds in ((3, (1,)), (5, (1, 2)), (7, (1, 2, 4)))]
    objects = [('C3', atoms[0][1]), ('T5', atoms[1][1]), ('Paley7', atoms[2][1]),
               ('81', build_hostile(False)), ('27', build_hostile(True))]
    out = []
    for name, rows in objects:
        canon = Canon(); code, transports = canon.maps(rows)
        base = transports[0]
        # An attaining order is a map from canonical slots to named source
        # vertices. Compose with the inverse base to obtain actual source autos.
        inverse = {v: i for i, v in enumerate(base)}
        group = {tuple(p[inverse[v]] for v in range(len(rows))) for p in transports}
        independent = actual_automorphisms(rows)
        need(group == independent, ('all transporters retained', name))
        masks = old.orbit_masks(len(rows), group)
        need(len(masks) == sum(2**len(old.permutation_cycles(p)) for p in group)//len(group),
             'independent Burnside mask count')
        out.append(dict(object=name, automorphisms=len(group), canon_leaves=canon.orders_tested,
                        binary_attachment_orbits=len(masks),
                        transporter_0_to_4=sum(p[0] == 4 for p in group) if len(rows) > 4 else None))
    need([r['automorphisms'] for r in out] == [3, 5, 21, 81, 27], 'strong atom groups')
    need(out[-2]['binary_attachment_orbits'] == 176 and out[-1]['binary_attachment_orbits'] == 512,
         '81/27 attachment control')
    return out


def k4_guard_closure():
    edges = list(__import__('itertools').combinations(range(4), 2))
    # Successful wedge pairing from Session87, in edge order 01,02,03,12,13,23.
    pairs = ((2, 0), (4, 3), (5, 1))
    ports = {m for pair in pairs for m in (1 << pair[0], 1 << pair[1], (1 << pair[0]) | (1 << pair[1]))}
    group = list(permutations(range(4)))
    def paths(rows):
        return tuple(p for p in group if all(rows[p[i]] >> p[i+1] & 1 for i in range(3)))
    def line(rows, p):
        # Complement all non-spine arcs in this named Hamilton presentation.
        spine = {frozenset((p[i], p[i+1])) for i in range(3)}
        rs = list(rows)
        for u, v in edges:
            if frozenset((u, v)) not in spine:
                rs[u] ^= 1 << v; rs[v] ^= 1 << u
        return tuple(rs)
    seen = set(ports); sizes = [len(seen)]
    while True:
        nxt = set(seen)
        for c in seen:
            rows = old.matrix(4, c)
            nxt.update(old.encode(rows, p) for p in group)
            nxt.add(c ^ 63)
            nxt.update(old.encode(line(rows, p), tuple(range(4))) for p in paths(rows))
        if nxt == seen:
            break
        seen = nxt; sizes.append(len(seen))
    need(seen == set(range(64)), 'actual native closure of the nine ports')
    guards = {c: paths(old.matrix(4, c)) for c in seen}
    need(len(set(guards.values())) == 64, 'named legality guards separate every source')
    # Exact four-bit readout and its hostile fiber, with named edges retained.
    def obs(c):
        degrees = [sum((c >> i) & 1 for i, e in enumerate(edges) if v in e) % 2 for v in range(3)]
        return tuple(degrees)+(c.bit_count() % 2,)
    need(len({obs(c) for c in ports}) == 9, 'four probes decode original ports')
    need(obs(0) == obs(45) and guards[0] != guards[45], 'minimal restored-kernel hostile')
    return dict(seed_ports=sorted(ports), closure_sizes=sizes, states=64,
                distinct_named_path_guard_profiles=64, lost_kernel_bits=2,
                same_four_probes_different_guards=[0, 45])


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--max-order', type=int, default=8)
    parser.add_argument('--skip-baseline', action='store_true')
    args = parser.parse_args()
    need(1 <= args.max_order <= 9, 'explicit finite scope through nine')
    parents = [()]; reports = []; banks = {0: parents}
    baseline = [()]; independent = old.ParentCanonicalizer()
    for n in range(args.max_order):
        parents, report = produce(parents, n)
        if not args.skip_baseline and n < 8:
            baseline, baseline_report = old.produce(baseline, n, audit=False, selector='neighbor_degrees')
            need({independent.maps(t)[0] for t in parents} == {independent.maps(t)[0] for t in baseline},
                 ('exact class-bank equality with Session87', n+1))
            report['baseline_competitor_queries'] = baseline_report.get('competitor_key_queries', 0)
            report['baseline_preparation_leaves'] = baseline_report['parent_preparation_orders_tested']
            report['baseline_competitor_leaves'] = baseline_report['competitor_orders_tested_after_cache']
        banks[n+1] = parents; reports.append(report)
    need([len(banks[n]) for n in range(1, args.max_order+1)] == [1, 1, 2, 4, 12, 56, 456, 6880, 191536][:args.max_order],
         'known counts are audit-only')
    bank_audit = None
    if args.max_order == 9:
        path = HERE.parent/'relation_tower_frobenius_ports_20261028/tournament_bank.json.gz'
        raw = path.read_bytes()
        reference = json.loads(gzip.decompress(raw))
        need(reference['order'] == 9 and len(reference['rows']) == 191536, 'reference scope')
        # The reference was generated by a different parent-isomorphism engine
        # and full root-switch blocks. Our complete invariant matches both banks.
        canon = Canon()
        actual = {canon.maps(t)[0] for t in parents}
        need(len(actual) == len(parents), 'no isomorphic duplicate at nine')
        expected = {canon.maps(tuple(t))[0] for t in reference['rows']}
        need(len(expected) == len(reference['rows']) and actual == expected, 'full order-nine class bijection')
        bank_audit = dict(classes=len(actual), reference_sha256=hashlib.sha256(raw).hexdigest(),
                          exact_class_bijection=True, reference_uses_different_generation_engine=True)
    # Independent full labelled S_n orbit census, not a known-count comparison.
    for n in range(1, min(args.max_order, 6)+1):
        brute, _ = old.labelled_orbit_census(n)
        canon = Canon()
        need({canon.maps(t)[0] for t in brute} == {canon.maps(t)[0] for t in banks[n]},
             ('literal full orbit census', n))
    return dict(status='FINITE-EXACT', generation=reports, strong_attachments=strong_controls(),
                nine_port_native_closure=k4_guard_closure(), checks=CHECKS,
                order_nine_bank_audit=bank_audit,
                full_labelled_independent_order=min(args.max_order, 6),
                no_polynomial_or_scalar_count_speed_claim=True)


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
