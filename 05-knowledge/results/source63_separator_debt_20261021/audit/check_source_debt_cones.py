#!/usr/bin/env python3
"""Independent exact audit. Does not import or execute the producer checker.

Published layout: this file is in the packet's audit/ directory. For staging,
pass --packet ../session18-debt-cones. Optional --repository authenticates every
upstream input against its pinned Git object, without modifying the repository.
"""
import argparse
import base64
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import re
import subprocess
import sys
import zlib

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

CHECKS = Counter()


def require(ok, label):
    if not ok:
        raise RuntimeError(label)
    CHECKS[label] += 1


def sha(data):
    return hashlib.sha256(data).hexdigest()


def closed_coefficients(word):
    """Unroll all affine and signed emissions over ONE common denominator."""
    cost = sum(word)
    numerator = positive = debt = 0
    prior = 0
    for i, a in enumerate(word):
        require(type(a) is int and a > 0, 'positive integer valuation')
        base = 3 ** (len(word) - i - 1) * 2 ** prior
        numerator += base
        if a == 1:
            positive += 2 * base
        elif a >= 3:
            debt += (2 ** a - 4) * base
        prior += a
    A, Q = 3 ** len(word), 2 ** cost
    require(A + numerator == Q + positive - debt, 'unrolled signed identity')
    return A, Q, numerator, F(positive, Q), F(debt, Q)


def state(word):
    A, Q, B, P, D = closed_coefficients(word)
    return F(A, Q), P, D


def discover(n, *, odd_limit=None, target_index=None, target_value=None):
    """Run native integer arithmetic, never using a supplied word as driver."""
    x = n
    word = []
    shortcut = ordinary = 0
    first_literal = first_half_literal = first_odd = first_half_odd = None
    target = None
    while True:
        require(x > 0 and x % 2 == 1, 'native positive odd state')
        if x == 1:
            break
        if odd_limit is not None and len(word) == odd_limit:
            break
        require(len(word) < 100000, 'native discovery cap not exhausted')
        y = 3 * x + 1
        ordinary += 1
        a = 0
        while y % 2 == 0:
            y //= 2
            a += 1
            shortcut += 1
            ordinary += 1
            if first_literal is None and y < n:
                first_literal = (shortcut, ordinary, y)
            if first_half_literal is None and 2 * y <= n:
                first_half_literal = (shortcut, ordinary, y)
        word.append(a)
        x = y
        if len(word) == target_index:
            target = x
        if first_odd is None and x < n:
            first_odd = (len(word), shortcut, ordinary, x)
        if first_half_odd is None and 2 * x <= n:
            first_half_odd = (len(word), shortcut, ordinary, x)
        if target_value is not None and x == target_value:
            break
    return dict(word=tuple(word), endpoint=x, target=target,
                shortcut=shortcut, ordinary=ordinary,
                first_literal=first_literal, first_odd=first_odd,
                first_half_literal=first_half_literal,
                first_half_odd=first_half_odd)


def baseline_accepts(j):
    """Reconstruct the inherited sufficient depth-18 criterion for this j."""
    x = (3 * pow(125, j, 2 ** 20) + 1) // 4
    a, b, q = 3, 1, 4
    for _ in range(18):
        if x % 2:
            x = (3 * x + 1) // 2
            a, b = 3 * a, 3 * b + q
        else:
            x //= 2
        q *= 2
        if 2 * (125 * a + b) <= 125 * q:
            return True
    return False


def exponent_with_exact_valuation(P, B, valuation):
    # Invert the isometry one bit at a time, retaining the terminal odd bit.
    residue = 0
    for precision in range(1, valuation):
        modulus = 2 ** (precision + 2)
        candidates = (residue, residue + 2 ** (precision - 1))
        hits = [j for j in candidates
                if (P * pow(125, j, modulus) + B - 2 ** valuation) % modulus == 0]
        require(len(hits) == 1, 'unique terminal-bit exponent inverse')
        residue = hits[0]
    return residue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--packet', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--repository', type=Path)
    args = parser.parse_args()
    packet = args.packet.resolve()
    note = (packet / 'SOURCE63_SEPARATOR_DEBT.md').read_text()
    receipt = json.loads((packet / 'receipt.normal.json').read_text())
    require((packet / 'receipt.normal.json').read_bytes() ==
            (packet / 'receipt.optimized.json').read_bytes(), 'producer normal and optimized receipts equal')
    packet_hashes = {}
    for line in (packet / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ', 1)
        packet_hashes[name] = sha((packet / name).read_bytes())
        require(packet_hashes[name] == digest, 'producer manifest hash')
    require(packet_hashes['SOURCE63_SEPARATOR_DEBT.md'] ==
            '8077d5e3552f93bca4dda19f5314bdadc26556396ea1df8658c4e9ca9466a22a',
            'frozen mathematical note identity')
    raw = (packet / 'inputs/guarded_children_105_20261018g.json').read_bytes()
    require(sha(raw) == 'c08611140beb4db0545e67de0d9ff8f066d6bfbd909b29ca6e9b716432800822',
            'frozen original 63-source input identity')
    inventory = json.loads((packet / 'SOURCE_INVENTORY.json').read_text())
    require(inventory['repository_revision'] == 'c4a7875652a5ea516e93a31c3c6410d2204773b2',
            'pinned inherited repository revision')
    provenance = []
    if args.repository:
        for item in inventory['source_files']:
            data = subprocess.check_output(['git', '-C', str(args.repository), 'show',
                                            inventory['repository_revision'] + ':' + item['path']])
            require(sha(data) == item['sha256'], 'pinned Git object authentication')
            provenance.append(dict(path=item['path'], sha256=sha(data)))
        original = subprocess.check_output(['git', '-C', str(args.repository), 'show',
            inventory['repository_revision'] + ':05-knowledge/results/guarded_children_105_20261018g.json'])
        require(original == raw, 'frozen input byte-identical to original Git source')

    # Parse the displayed source word and each displayed exact fraction.
    match = re.search(r'\bu=\((.*?)\)\.', note, re.S)
    require(match is not None, 'displayed source63 word located')
    U = tuple(int(x) for x in re.findall(r'\d+', match.group(1)))
    fractions = [F(int(n), int(d)) for n, d in re.findall(r'(?<![\w^])([0-9]{3,})\s*/\s*([0-9]{3,})(?!\w)', note)]
    l, p, d = state(U)
    c = d / (l - F(2, 3))
    W = U + (4,)
    A, Q, B, pw, dw = closed_coefficients(W)
    expected_fractions = {
        'lambda': l, 'positive': p, 'debt': d,
        'ordinary cone': d / (l - 1), 'separator cone': c,
        'uniform margin': 1 - c, 'exit debt': dw, 'exit margin': 2 - dw,
        'physical barrier': F(B, Q - A), 'positive barrier': 1 + pw / (1 - F(A, Q))}
    for label, value in expected_fractions.items():
        require(value in fractions, 'displayed exact fraction: ' + label)
    require(len(U) == 33 and sum(U) == 50, 'source63 prefix clocks')
    require(d / (l - 1) > 1 and c < 1 and d < l, 'minimal finite separator repair')
    require(d / (l - 1) <= F(20, 13), 'inherited terminal-packet criterion already sufficient')
    require((A, Q, B) == (3 ** 34, 2 ** 54, 48215861999407177), 'displayed grounded-family carry')
    exact_residue = ((Q - B) * pow(A, -1, 2 * Q)) % (2 * Q)
    require(exact_residue == 63 and 2 * Q == 2 ** 55, 'source63 exact endpoint-odd cylinder')
    original = discover(63)
    require(original['word'] == W + (3, 1, 1, 5, 4), 'independently discovered source63 ROOT word')
    require(original['first_odd'] == (34, 54, 88, 61), 'source63 first odd descent')
    require(original['first_literal'] == (54, 88, 61), 'source63 first literal descent')
    require((len(original['word']), original['shortcut'], original['ordinary']) == (39, 68, 107),
            'source63 first ROOT clocks')
    require((63 % 15, 63 % 35, 61 % 15, 61 % 35) == (3, 28, 1, 26), 'source63 actual prime registers')
    for t in (0, 1, 2, 7, 29):
        n = 63 + 2 ** 55 * t
        run = discover(n, odd_limit=34)
        require(run['word'] == W and run['endpoint'] == 61 + 2 * 3 ** 34 * t, 'exact source63 first-exit lift')
        require(run['first_odd'][:3] == (34, 54, 88) and run['first_literal'][:2] == (54, 88),
                'source63 cylinder preserves first clocks')
        require(n - run['endpoint'] > dw, 'source63 lifted first exit paid')
    lifted_root = discover(63 + 2 ** 69, odd_limit=39)
    require(lifted_root['word'] == original['word'] and lifted_root['endpoint'] == 1 + 2 * 3 ** 39,
            'ROOT-word guard is not ROOT for its positive lift')
    require(3 ** 34 > 2 ** 53, 'no early final-run literal descent coefficient')
    for i in range(1, len(U) + 1):
        require(3 ** i > 2 ** sum(U[:i]), 'every source63 proper odd prefix expands')

    # This family is symbolic: no attempt to allocate 2**H or N_1.
    H = 2 * 3 ** 33
    require(H == receipt['grounded_source63_family']['valuation_period'], 'huge period exact integer')
    for s in range(34):
        h = 2 * 3 ** s
        require(pow(2, h, 3 ** (s + 1)) == 1 and pow(2, h, 3 ** (s + 2)) != 1,
                'lifting controls exact 3-adic valuation')
    require(pow(2, H, 3 ** 34) == 1 and pow(2, H // 2, 3 ** 34) != 1 and
            pow(2, H // 3, 3 ** 34) != 1, 'huge period order exact')
    require(61 * 2 ** 54 - B == 63 * 3 ** 34, 'N zero is exactly original63')
    require(61 * (2 ** 54 - 3 ** 34) - B == 2 * 3 ** 34, 'N family first literal strict gap base')
    for control in receipt['grounded_source63_family']['controls']:
        t = control['t']
        e = 54 + H * t
        require((61 * pow(2, e, 3 ** 34) - B) % 3 ** 34 == 0, 'symbolic N integrality')
        residue = (61 * pow(2, e, 210 * 3 ** 34) - B) % (210 * 3 ** 34)
        require(residue % (3 ** 34) == 0, 'N modular division exact')
        n210 = residue // (3 ** 34)
        require(n210 % 2 == 1, 'symbolic N oddness')
        require([n210 % 15, n210 % 35] == [control['source_mod15'], control['source_mod35']],
                'symbolic N actual prime-register formula')
        require((control['final_valuation'], control['first_literal_shortcut_clock'],
                 control['first_literal_ordinary_clock'], control['first_odd_shortcut_clock'],
                 control['ROOT_shortcut_clock'], control['ROOT_ordinary_clock']) ==
                (4 + H * t, 54, 88, 54 + H * t, 68 + H * t, 107 + H * t), 'symbolic N separated clocks')

    V = (1, 1, 1, 3)
    require(state(V)[0] == F(81, 64) and state(V)[2] == F(1, 2), 'second generator exact data')
    require(state(V)[2] / (state(V)[0] - F(2, 3)) == F(96, 115) < c, 'second generator same cone')
    for k in (1, 2, 3):
        for generators in product((U, V), repeat=k):
            for gaps in product((1, 2), repeat=k - 1):
                word = generators[0]
                for gap, gen in zip(gaps, generators[1:]):
                    word += (1,) * gap + gen
                lm, pp, dd = state(word)
                require(dd <= c * (lm - F(2, 3)), 'mixed-generator separator control')
    require(state(U + U)[2] > c * (state(U + U)[0] - F(2, 3)), 'no-separator cone counterexample')
    hostile = state((3, 1, 1, 1, 1))
    require(hostile[2] >= hostile[0], 'no finite separator repair hostile')
    for row in receipt['macro_controls']:
        prefix = (U + (1,)) * (row['copies'] - 1) + U
        lm, pp, dd = state(prefix)
        require(dd == c * (lm - F(2, 3)), 'all-depth equality finite controls')
        word = prefix + (row['last_valuation'],)
        n = int(row['source'])
        run = discover(n, odd_limit=len(word))
        require(run['word'] == word and run['endpoint'] == int(row['endpoint']), 'macro exact actual source and child')
        require(run['first_odd'][0] == len(word), 'macro actual first odd exit')
        require(n - run['endpoint'] - state(word)[2] > 1 - c, 'macro strict payment margin')
        require((len(word), sum(word), len(word) + sum(word), sum(a >= 3 for a in word)) ==
                (row['odd_steps'], row['shortcut_clock'], row['ordinary_clock'], row['debt_resets']),
                'macro actual clocks and reset count')
        require(run['first_literal'][:2] == (row['first_literal_shortcut'], row['first_literal_ordinary']),
                'macro literal clock')
    require(len(receipt['macro_controls'][-1]['source']) == 1038 and
            receipt['macro_controls'][-1]['odd_steps'] == 2176 and
            receipt['macro_controls'][-1]['debt_resets'] == 257 and
            receipt['macro_controls'][-1]['shortcut_clock'] == 3449, 'displayed largest finite macro control')
    later = discover(55, target_value=53)
    require(len(later['word']) == 39 and state(later['word'])[2] > 2 and later['first_odd'][0] < 39,
            'source55 later-descent hostile is not first exit')
    native4575 = discover(4575)
    w4575 = native4575['word'][:native4575['first_odd'][0]]
    require(state(w4575)[2] > 2, 'source4575 first-exit debt exceeds two')

    rows = json.loads(raw)['rows']
    require(len(rows) == 63 and len({r['j'] for r in rows}) == 63, '63 distinct original exponents')
    residual = [j for j in range(1, 1025) if not baseline_accepts(j)]
    require(residual == sorted(r['j'] for r in rows), 'independent depth18 residual universe exactly retained')
    peaks = []
    source_clocks = []
    total_odd = total_shortcut = total_ordinary = 0
    for row in rows:
        j = row['j']
        n = 125 ** j
        v = tuple(row['prefix'])
        a0 = sum(v)
        c0 = row['c_min']
        ap, qp, bp, pv, dv = closed_coefficients(v)
        P, B0 = 3 * ap, 3 * bp + qp
        require(P == int(row['P']) and B0 == int(row['B']) and a0 == row['prefix_cost'],
                'inherited template exact carrier')
        require(c0 >= 2 and row['guard_bits'] == a0 + c0 - 2 and 0 < j < 2 ** row['guard_bits'],
                'all-height guard quantified least positive anchor')
        require((P * n + B0) % 2 ** (a0 + c0) == 0, 'anchor minimum terminal divisibility')
        require(2 * (P * n + B0) <= 2 ** (a0 + c0) * n, 'anchor proves all-height half inequality')
        for i in range(1, len(v) + 1):
            require(2 * 3 ** i >= 2 ** sum(v[:i]), 'all proper odd checkpoints above half')
        lm, pm, dm = state(v + (c0,))
        R = 1 + pm / (1 - lm)
        require(lm < F(1, 2) and dm < 2 and R < 9 < n, 'all63 exact charge and barrier bounds')
        peaks.append((R, dm, j))
        exported = next(z for z in receipt['templates']['rows'] if z['j'] == j)
        require(F(exported['barrier']) == R and F(exported['debt_at_c_min']) == dm, 'producer row agrees with independent coefficients')
        for cc in (c0, c0 + 1, c0 + 12, c0 + 64):
            valuation = a0 + cc
            J = exponent_with_exact_valuation(P, B0, valuation)
            require((J - j) % 2 ** row['guard_bits'] == 0, 'exact terminal extension preserves original source guard')
            modulus = 2 ** (valuation + 1)
            num = (P * pow(125, J, modulus) + B0) % modulus
            require(num == 2 ** valuation, 'compiled endpoint terminal odd bit')
            # Source modulo 2Q fixes every exact valuation, independently of any real representative size.
            x = pow(125, J, modulus)
            for a in v + (cc,):
                y = 3 * x + 1
                require(y % (2 ** (a + 1)) == 2 ** a, 'compiled word exact native valuation')
                x = y // 2 ** a
            native = []
            for m in (15, 35, 105):
                numerator = (P * pow(125, J, m * 2 ** valuation) + B0) % (m * 2 ** valuation)
                require(numerator % 2 ** valuation == 0, 'prime register full denominator retained')
                observed = numerator // 2 ** valuation
                require(observed == ((P * pow(125, J, m) + B0) * pow(2, -valuation, m)) % m,
                        'prime register exact modular division')
                native.append(observed)
            require(native[0] % 5 == native[1] % 5 and
                    (native[0] + 15 * ((native[1] - native[0]) % 7)) % 105 == native[2], 'native 15/35 compatible CRT')
            lc, pc, dc = state(v + (cc,))
            require(1 + pc / (1 - lc) <= R and dc <= max(F(1), dm) < 2, 'variable terminal valuation charge control')
        run = discover(n, target_index=len(v) + 1)
        word = bytes(run['word'])
        retained = zlib.decompress(base64.b85decode(row['word_zlib_base85']))
        require(word == retained and sha(word) == row['word_sha256'], 'independently discovered complete original ROOT word')
        require(run['endpoint'] == 1, 'actual original 125-power ROOT endpoint')
        require(run['word'][:len(v)] == v and run['word'][len(v)] == row['actual_c'], 'actual original template word prefix')
        child = run['target']
        require(child == int(row['child']), 'actual original child unchanged')
        suffix = word[len(v) + 1:]
        require(sha(suffix) == row['child_word_sha256'] and
                (len(suffix), sum(suffix)) == (row['child_rank'], row['child_cost']), 'actual child suffix identity and clocks')
        require([child % 15, child % 35] == row['child_residues'] and child % 105 == row['child_mod105'],
                'actual original child native registers')
        require([n % 15, n % 35] == row['source_residues'], 'original source native registers')
        require(run['first_half_odd'][0] == len(v) + 1 == row['first_odd_half_index'], 'actual first odd half descent')
        require(run['first_half_literal'][0] == exported['first_literal_half_shortcut_seed'] and
                run['first_half_odd'][1] == exported['odd_half_shortcut_seed'], 'distinct actual literal and odd half clocks')
        require((len(word), sum(word), sum(word) + len(word)) ==
                (row['word_rank'], row['word_cost'], exported['ROOT_ordinary_clock']), 'actual original ROOT clocks')
        total_odd += len(word)
        total_shortcut += run['shortcut']
        total_ordinary += run['ordinary']
        source_clocks.append(dict(j=j, odd=len(word), shortcut=run['shortcut'], ordinary=run['ordinary'],
                                 first_literal_half=run['first_half_literal'][:2],
                                 first_odd_half=run['first_half_odd'][:3], word_sha256=sha(word)))
    require((total_odd, total_shortcut, total_ordinary) == (575652, 1152705, 1728357), 'all63 authentic ROOT total clocks')
    peak = max(peaks)
    peak_debt = max(peaks, key=lambda x: x[1])
    require(peak[2] == peak_debt[2] == 638 and next(r for r in rows if r['j'] == 638)['actual_c'] == 8,
            'common barrier actually attained at638')
    require(peak[0] in fractions and peak_debt[1] in fractions, 'displayed exact common barrier and debt')
    require(str(peak[0]) == receipt['templates']['sharp_common_barrier'] and
            str(peak_debt[1]) == receipt['templates']['largest_debt_at_c_min'], 'exported common extrema')
    require(pow(2, 12, 105) == 1 and all(pow(2, d, 105) != 1 for d in (1, 2, 3, 4, 6)),
            'observer period12 is exact')

    # Parse the displayed residue, rather than copying it from the producer receipt.
    match = re.search(r'j = ([0-9]{20,})\s+mod ([0-9]{20,})', note)
    require(match is not None, 'displayed complement residue parsed')
    J, M = map(int, match.groups())
    L = max(r['prefix_cost'] + r['c_min'] for r in rows)
    require(L == 149 and M == 2 ** 149 and 0 < J < M, 'complement exact modulus and positive class')
    require((3 * pow(125, J, 2 ** 151) + 5) % 2 ** 151 == 0, 'displayed complement actual exponent congruence')
    x = (3 * pow(125, J, 2 ** 151) + 1) // 4
    for _ in range(149):
        require(x % 2 == 1, 'complement 149 actual odd shortcut branches')
        x = (3 * x + 1) // 2
    require(not baseline_accepts(J), 'complement misses baseline depth18 atlas')
    for row in rows:
        require((J - row['j']) % 2 ** row['guard_bits'] != 0, 'complement misses every one of63 guards')
    require(receipt['uncovered_exponent_progression']['residue'] == str(J) and
            receipt['uncovered_exponent_progression']['modulus'] == str(M), 'exported infinite complement identity')

    output = dict(status='PASS: independently audited scoped proof and finite-exact sources; no universal coverage claim',
                  packet_hashes=packet_hashes, upstream_revision=inventory['repository_revision'],
                  upstream_git_authentication=provenance,
                  count=sum(CHECKS.values()), checks=dict(sorted(CHECKS.items())),
                  original_source63=dict(odd=39, shortcut=68, ordinary=107,
                      separator_constant=str(c), first_exit_margin=str(2 - dw)),
                  symbolic_grounded_family_period=H,
                  independent_residual_exponents=residual,
                  original_125_source_clocks=source_clocks,
                  totals=dict(odd=total_odd, shortcut=total_shortcut, ordinary=total_ordinary),
                  sharp_common_barrier=str(peak[0]), maximum_template_debt=str(peak_debt[1]),
                  uncovered_progression=dict(residue=str(J), modulus=str(M), no_half_through=151))
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
