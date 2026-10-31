#!/usr/bin/env python3
"""Independent weak-B/Q25 audit. Owned JSON and compressed integer data only.

No imports from producer or previous audit modules; no external C needed.
Boundary edges are extracted from actual expanded paths by port decoding,
independently of the producer's symbolic boundary builder.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import struct

MASKS = [1585793692427, 1599483269979, 1602705117979,
         3372531545867, 3389445068555, 3526867243787]
PORTS = {1, 3, 8, 10}

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def edges(p):
    return {tuple(sorted((a, b))) for a, b in zip(p, p[1:])}

def delta(a, b):
    return len(edges(a) - edges(b))

def validate(p):
    n = len(p)
    check(sorted(p) == list(range(1, n + 1)), 'not a permutation')
    check(list(p[:2]) == [1, 3], 'wrong protected prefix')
    check(p[-1] == (10 if n % 2 else 8), 'wrong terminal port')
    check(all(math.isqrt(a + b) ** 2 == a + b for a, b in zip(p, p[1:])), 'nonsquare edge')
    check({e for e in edges(p) if set(e) <= PORTS} == {(1, 3)}, 'unexpected all-port edge')
    return {v: j % 2 for j, v in enumerate(p, 1)}

def validate_pair(pair):
    a, b = pair
    check(len(b) == len(a) + 1, 'pair lengths')
    c, d = validate(a), validate(b)
    check(all(c[v] == d[v] for v in a), 'pair phase mismatch')

def enumerate_spectrum(n):
    adjacency = [0] * (n + 1)
    for square in (s * s for s in range(2, math.isqrt(2 * n) + 1)):
        for a in range(max(1, square - n), min(n, (square - 1) // 2) + 1):
            b = square - a
            adjacency[a] |= 1 << (b - 1)
            adjacency[b] |= 1 << (a - 1)
    goal = 10 if n % 2 else 8
    spectrum = Counter()
    nodes = 0
    def dfs(last, unused, odd_mask, length):
        nonlocal nodes
        nodes += 1
        if not unused:
            if last == goal:
                spectrum[odd_mask] += 1
            return
        if last == goal:
            return
        available = unused | (1 << (last - 1))
        scan = unused
        while scan:
            bit = scan & -scan
            scan -= bit
            v = bit.bit_length()
            if (adjacency[v] & available).bit_count() < (1 if v == goal else 2):
                return
        reached = todo = 1 << (last - 1)
        while todo:
            bit = todo & -todo
            todo -= bit
            new = adjacency[bit.bit_length()] & available & ~reached
            reached |= new
            todo |= new
        if reached != available:
            return
        options = adjacency[last] & unused
        parity = (length + 1) % 2
        while options:
            bit = options & -options
            options -= bit
            v = bit.bit_length()
            if (v == 8 and parity) or (v == 10 and not parity):
                continue
            dfs(v, unused ^ bit, odd_mask | (bit if parity else 0), length + 1)
    dfs(3, ((1 << n) - 1) ^ 5, 1, 2)
    return spectrum, nodes

def phase_audit(catalogue):
    s42, nodes42 = enumerate_spectrum(42)
    s43, nodes43 = enumerate_spectrum(43)
    restricted = Counter()
    for mask, count in s43.items():
        restricted[mask & ((1 << 42) - 1)] += count
    common = sorted(s42.keys() & restricted.keys())
    check(common == MASKS, 'complete six-mask census mismatch')
    multiplicities = [(s42[m], restricted[m]) for m in common]
    check(sum(a * b for a, b in multiplicities) == 20, 'first pair count')
    source = json.loads(Path(catalogue).read_text())
    witnesses = source['pairs'][-1]['witnesses']
    check(sorted(w['common_odd_mask'] for w in witnesses) == common, 'catalogue disagrees')
    for w in witnesses:
        validate_pair((w['shorter'], w['longer']))
        check(sum(1 << (v - 1) for v in w['shorter'][::2]) == w['common_odd_mask'], 'catalogue witness mask')
    hostile = {}
    p = list(witnesses[0]['shorter'])
    bad_duplicate = p[:]
    bad_duplicate[2] = bad_duplicate[3]
    bad_nonsquare = [1, 3, 2] + [v for v in range(4, 43) if v != 8] + [8]
    for name, bad in [('duplicate', bad_duplicate), ('wrong_prefix', list(reversed(p))), ('nonsquare', bad_nonsquare)]:
        try:
            validate(bad)
        except RuntimeError:
            hostile[name] = True
        else:
            hostile[name] = False
    try:
        validate_pair((witnesses[0]['shorter'], witnesses[1]['longer']))
    except RuntimeError:
        hostile['separately_valid_phase_mismatch'] = True
    else:
        hostile['separately_valid_phase_mismatch'] = False
    check(all(hostile.values()), 'hostile control accepted')
    obstructions = []
    for m in common:
        n = 46 if m in (MASKS[0], MASKS[1], MASKS[3]) else 47
        def color(v):
            return (m >> (v - 1)) & 1 if v <= 42 else v % 2
        ns = [s * s - n for s in range(2, math.isqrt(2 * n) + 1)
              if 1 <= s * s - n < n and color(s * s - n) != color(n)]
        expect = [35] if m == MASKS[1] else ([18] if n == 46 else [34])
        check(ns == expect, 'degree-one fixed-phase obstruction')
        obstructions.append(dict(mask=m, failing_order=n, internal_vertex=n, neighbors=ns))
    port_edges = [(a, b) for a in sorted(PORTS) for b in sorted(PORTS)
                  if a < b and math.isqrt(a + b) ** 2 == a + b]
    check(port_edges == [(1, 3), (1, 8)], 'port-square algebra')
    check(math.isqrt(1 + 10) ** 2 != 11, 'odd endpoint cannot close')
    return dict(onset_counts=[dict(n=42, paths=sum(s42.values()), masks=len(s42), search_nodes=nodes42),
                              dict(n=43, paths=sum(s43.values()), masks=len(s43), search_nodes=nodes43)],
                common_masks=common, shorter_longer_multiplicities=multiplicities,
                pair_count=20, obstructions=obstructions, square_edges_inside_ports=port_edges,
                hostile_controls=hostile)

def owned_audit(operator_dir, base_dir, small_path):
    operator_dir, base_dir = Path(operator_dir), Path(base_dir)
    input_hashes = {}
    def read(path, expected=None):
        raw = Path(path).read_bytes()
        sha = digest(raw)
        if expected is not None:
            check(sha == expected, 'source pin mismatch: ' + str(path))
        resolved = Path(path).resolve()
        if resolved.is_relative_to(operator_dir.resolve()):
            key = 'session124_mixed_radix_square_rays_20261031/' + resolved.relative_to(operator_dir.resolve()).as_posix()
        elif resolved.is_relative_to(base_dir.resolve()):
            key = 'square_sum_B_base_interval_20261019/' + resolved.relative_to(base_dir.resolve()).as_posix()
        elif resolved == Path(small_path).resolve():
            key = 'session124_square_prefix_threshold42_20261031/SMALL_BASES.json'
        else:
            raise RuntimeError('unrecognized owned input kind')
        input_hashes[key] = sha
        return raw
    deps = json.loads(read(operator_dir / 'DEPENDENCIES.json'))
    base_deps = json.loads(read(operator_dir / 'BASE_DEPENDENCIES.json'))
    cert_pins = json.loads(read(operator_dir / 'CERTIFICATES.json'))
    small = json.loads(read(small_path, deps['small_base_sha256']))
    manifest = json.loads(read(base_dir / base_deps['manifest_file'], base_deps['manifest_sha256']))
    check(manifest == base_deps['chunks'], 'base manifest content mismatch')
    # Sequential byte decoding, not the producer's closed-form random-access index.
    bases = {n: tuple(tuple(p) for p in small[str(n)]) for n in range(42, 53)}
    decoded_chunks = []
    for chunk in manifest:
        if chunk['first'] > 1061:
            continue
        compressed = read(base_dir / chunk['file'], chunk['sha256'])
        check(len(compressed) == chunk['bytes'], 'compressed byte count')
        raw = gzip.decompress(compressed)
        check(digest(raw) == chunk['decoded_sha256'] and len(raw) == chunk['decoded_bytes'], 'decoded chunk pin')
        check(raw[:8] == b'SSBAS1\0\0', 'base header magic')
        first, last = struct.unpack('<II', raw[8:16])
        check((first, last) == (chunk['first'], chunk['last']), 'base header interval')
        cursor = 16
        for n in range(first, last + 1):
            pair = []
            for size in (n, n + 1):
                p = [1]
                for root in raw[cursor:cursor + size - 1]:
                    p.append(root * root - p[-1])
                cursor += size - 1
                pair.append(tuple(p))
            if n <= 1061:
                bases[n] = tuple(pair)
        check(cursor == len(raw), 'trailing or missing chunk bytes')
        decoded_chunks.append(chunk['file'])
    check(sorted(bases) == list(range(42, 1062)), 'base interval coverage')
    base_stats = []
    for n, pair in bases.items():
        validate_pair(pair)
        base_stats.append((delta(*pair), n))
    certificates = {}
    for family in ('A', 'B'):
        for parity in (0, 1):
            for digit in range(12, 37):
                name = f'{family}/paired/B_p{parity}_d{digit}.json'
                certificates[family, parity, digit] = json.loads(read(operator_dir / name, cert_pins[name]))
    def lift(family, n, parents):
        parent = (n - 12) // 25
        digit = n - 25 * parent
        output = []
        for extension, word in enumerate(certificates[family, parent % 2, digit]['words']):
            result = []
            members = {}
            heads = []
            tails = []
            singles = []
            for atom in word:
                kind, c = atom[:2]
                if kind == 'single':
                    singles.append(c)
                    result.append(c)
                elif kind == 'head':
                    heads.append(c)
                    result.append(25 - c)
                else:
                    member, orientation = atom[2:]
                    check(kind in ('block', 'tail') and orientation in (-1, 1), 'unknown atom')
                    members[c] = member
                    block = [25 * v + (c if j % 2 else -c) for j, v in enumerate(parents[member])]
                    if kind == 'tail':
                        tails.append(c)
                        block = block[1:]
                    result.extend(block[::orientation])
            check(sorted(singles) == list(range(1, 13)), 'singleton partition')
            check(sorted(members) == list(range(-12, 13)), 'complete signed offset partition')
            check(sorted(heads) == sorted(tails), 'head/tail partition')
            check(all(e == int(25 + (-1) ** (parent + 1) * c <= digit + extension) for c, e in members.items()), 'forced member threshold')
            output.append(tuple(result))
        return tuple(output)
    @lru_cache(maxsize=32)
    def pair(family, n):
        if n <= 1061:
            return bases[n]
        p = lift(family, n, pair(family, (n - 12) // 25))
        validate_pair(p)
        return p
    def boundary(p):
        port_centers = PORTS | {0}
        return {e for e in edges(p) if all((v + 12) // 25 in port_centers for v in e)}
    boundaries = {}
    family_results = {}
    for family in ('A', 'B'):
        states = []
        for parity in (0, 1):
            parent = 42 + parity
            for digit in range(12, 37):
                n = 25 * parent + digit
                p = pair(family, n)
                B0, B1 = map(boundary, p)
                check(len(B0) == len(B1) == 61, 'boundary count not 61')
                boundaries[family, parity, digit, 0] = B0
                boundaries[family, parity, digit, 1] = B1
                k = len(B0 - B1)
                check(delta(*p) == delta(*bases[parent]) + k, 'within-pair decomposition')
                states.append(dict(parity=parity, digit=digit, boundary_deletions=k))
        noncarry = [len(boundaries[family, p, d, 1] - boundaries[family, p, d + 1, 0])
                    for p in (0, 1) for d in range(12, 36)]
        carry = {p: len(boundaries[family, 1 - p, 36, 1] - boundaries[family, p, 12, 0]) for p in (0, 1)}
        check(carry == {0: 0, 1: 0}, 'nonzero pinned carry boundary deletion')
        for parity in (0, 1):
            check(boundaries[family, 1 - parity, 36, 1] == boundaries[family, parity, 12, 0], 'carry boundary sets unequal')
        rays = []
        for seed in (44, 45):
            r = seed
            expected = delta(bases[seed - 1][1], bases[seed][0])
            for depth in range(3):
                old, new = pair(family, r - 1), pair(family, r)
                h = delta(old[1], new[0])
                check(h == expected * 25 ** depth, 'exact carry ray')
                canonical = delta(old[0], new[0])
                within = delta(*old)
                check(h - within - 1 <= canonical <= h + within, 'canonical path triangle')
                local_costs = []
                if depth:
                    preceding = (pair(family, r - 2)[1], old[0])
                    following = (new[1], pair(family, r + 1)[0])
                    check(delta(*preceding) <= 61 and delta(*following) <= 61, 'noncarry neighbors')
                    for x in preceding:
                        for y in following:
                            check(delta(x, new[0]) >= h - within - 62, 'short-middle incoming lower bound')
                            check(delta(old[1], y) >= h - delta(*new) - 62, 'long-middle outgoing lower bound')
                            for z in (old[1], new[0]):
                                local_costs.append(max(delta(x, z), delta(z, y)))
                    check(min(local_costs) >= h - max(within, delta(*new)) - 62, 'bank local maximum lower bound')
                rays.append(dict(seed=seed, depth=depth, order=r, gap=h, canonical_deletions=canonical,
                                 within_pair_deletions=within, smallest_bank_two_step_max=min(local_costs) if local_costs else None))
                r = 25 * r + 12
        family_results[family] = dict(states=states, carry_K=carry, noncarry_gap_max=max(noncarry),
                                     within_pair_K_range=[min(x['boundary_deletions'] for x in states), max(x['boundary_deletions'] for x in states)],
                                     ray_samples=rays)
    check(family_results['A']['within_pair_K_range'] == [4, 22], 'A K range')
    check(family_results['B']['within_pair_K_range'] == [3, 22], 'B K range')
    check([family_results[f]['noncarry_gap_max'] for f in ('A', 'B')] == [37, 36], 'noncarry maxima')
    ratio = Fraction(14, 1) / Fraction(91, 2)
    check(ratio == Fraction(28, 91), 'pinned seed45 ratio')
    return dict(base_pair_count=len(bases), base_range=[42, 1061], base_path_max=max(base_stats),
                base_junction_gaps={r: delta(bases[r - 1][1], bases[r][0]) for r in range(43, 47)},
                states_expanded=len(certificates), decoded_chunks=decoded_chunks,
                boundary_extraction='actual expanded edges with both decoded endpoints in {0,1,3,8,10}',
                families=family_results, universal_bank_path_limsup_lower_bound='4/93',
                pinned_seed45_two_step_max_ratio=str(ratio), source_sha256=input_hashes)

RESULTS=Path(__file__).resolve().parents[3]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--operator-dir', default=str(RESULTS/'session124_mixed_radix_square_rays_20261031'))
    ap.add_argument('--base-dir', default=str(RESULTS/'square_sum_B_base_interval_20261019'))
    ap.add_argument('--small-bases', default=str(RESULTS/'session124_square_prefix_threshold42_20261031'/'SMALL_BASES.json'))
    ap.add_argument('--weak-b-catalogue', default=str(RESULTS/'session125_origin32_numeration_20261031'/'pair_taxonomy'/'weak_b_certificate.json'))
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    report = dict(status='PASS', scope='weak-B Hamiltonian path banks, not cycle closures',
                  producer_modules_imported=False, external_C_required=False,
                  phase=phase_audit(a.weak_b_catalogue),
                  lift=owned_audit(a.operator_dir, a.base_dir, a.small_bases),
                  catalogue_sha256=digest(Path(a.weak_b_catalogue).read_bytes()))
    Path(a.output).write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print('PASS independent weak-B/Q25 audit:', a.output)

if __name__ == '__main__':
    main()
