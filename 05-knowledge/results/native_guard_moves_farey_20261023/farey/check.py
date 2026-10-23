"""Exact native incidence, Farey/Berggren, and rational clock controls."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
import json
import sys

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def coupling(w):
    r = [sum(row) for row in w]
    c = [sum(row[j] for row in w) for j in range(len(w[0]))]
    return [[sum(F(w[i][t]*w[j][t], r[i]*c[t])
                 for t in range(len(c)) if c[t])
             for j in range(len(w))] for i in range(len(w))]


def conductance(w):
    return sum(F(w[0][t]*w[1][t], w[0][t]+w[1][t])
               for t in range(len(w[0])) if w[0][t]+w[1][t])


def native_controls():
    chords = ((0, 2), (0, 3), (1, 3))
    target_scores = ((0, 1, 2, 3), (0, 2, 2, 2), (1, 1, 1, 3), (1, 1, 2, 2))
    w = [[0]*4 for _ in range(4)]
    records = []
    for code in range(8):
        arcs = {(i, j) for i in range(4) for j in range(i+1, 4)}
        for k, (i, j) in enumerate(chords):
            if code >> k & 1:
                arcs.remove((i, j)); arcs.add((j, i))
        score = tuple(sorted(sum((i, j) in arcs for j in range(4)) for i in range(4)))
        h = sum(all((a, b) in arcs for a, b in zip(p, p[1:]))
                for p in permutations(range(4)))
        target = target_scores.index(score)
        w[code.bit_count()][target] += 1
        records.append({'word': code, 'target': target, 'Hamiltonian_orders': h})
        need(h == (1, 3, 3, 5)[target], 'Actual target classification')
    expected = [[1, 0, 0, 0], [0, 1, 1, 1], [0, 0, 0, 3], [0, 0, 0, 1]]
    need(w == expected, 'Eight-word incidence')
    p = coupling(w[1:])
    need(p == [[F(11, 15), F(1, 5), F(1, 15)],
               [F(1, 5), F(3, 5), F(1, 5)],
               [F(1, 5), F(3, 5), F(1, 5)]], 'Native three-source matrix')
    need([sum(p[i][j]*[4, -3, -3][j] for j in range(3)) for i in range(3)]
         == [F(8, 15)*v for v in (4, -3, -3)], 'Exact slow native mode')
    w0 = [w[1], [w[2][j]+w[3][j] for j in range(4)]]
    w1 = [[w[0][j]+w[3][j] for j in range(4)],
          [w[1][j]+w[2][j] for j in range(4)]]
    p0, p1 = coupling(w0), coupling(w1)
    need(p0 == [[F(11, 15), F(4, 15)], [F(1, 5), F(4, 5)]], 'Native identical-row lump')
    need(p1 == [[F(3, 5), F(2, 5)], [F(2, 15), F(13, 15)]], 'Complement-paired full incidence')
    need(conductance(w0) == conductance(w1) == F(4, 5), 'Shared Schur conductance')
    need(p0[0][1]+p0[1][0] == F(7, 15), 'Native gap')
    need(p1[0][1]+p1[1][0] == F(8, 15), 'Paired gap')
    need([sum(row) for row in w0] == [3, 4], 'First masses')
    need([sum(row) for row in w1] == [2, 6], 'Paired masses')
    need([r['word'] for r in records if r['word'].bit_count() == 1 and r['target'] == 3]
         == [2], 'Unique one-edge strong-target bridge')
    guarded = [row[:] for row in w0]
    guarded[0][3] -= 1
    need(conductance(guarded) == 0 and coupling(guarded) == [[1, 0], [0, 1]],
         'Single-presentation deletion disconnects')
    # A perfect two-block quotient need not control the full chain.
    double = [[p0[i % 2][j % 2] if i//2 == j//2 else F(0)
               for j in range(4)] for i in range(4)]
    signs = (1, 1, -1, -1)
    need([sum(double[i][j]*signs[j] for j in range(4)) for i in range(4)] == list(signs),
         'Hidden mean-zero invariant copy mode')
    for i in range(4):
        need([sum(double[i][j] for j in block) for block in ((0, 2), (1, 3))]
             == p0[i % 2], 'Same strongly lumped quotient')
    return {'word_records': records, 'incidence': w,
            'nontransitive_lumped_incidence': w0, 'complement_paired_incidence': w1,
            'shared_raw_conductance': '4/5', 'gaps': ['7/15', '8/15'],
            'deleted_bridge_word': 2}


def child(pair, letter):
    m, n = pair
    if letter == 'A':
        return 2*m-n, m
    if letter == 'B':
        return 2*m+n, m
    if letter == 'C':
        return m+2*n, n
    raise ValueError(letter)


def spectral_child(g, letter):
    return {'A': (1-g)/(3-4*g), 'B': (1-g)/(3-2*g), 'C': g/(1+2*g)}[letter]


def decode(pair):
    reverse = []
    m, n = pair
    while (m, n) != (2, 1):
        old_sum = m+n
        if n < m < 2*n:
            m, n = n, 2*n-m
            reverse.append('A')
        elif 2*n < m < 3*n:
            m, n = n, m-2*n
            reverse.append('B')
        elif m > 3*n:
            m, n = m-2*n, n
            reverse.append('C')
        else:
            raise RuntimeError('Unexpected parity/root boundary')
        need(0 < n < m and m+n < old_sum and gcd(m, n) == 1 and (m+n) % 2 == 1,
             'Strict valid ternary parent')
    return ''.join(reversed(reverse))


def cw_word(value):
    a, b = value.numerator, value.denominator
    reverse = []
    while a != b:
        if a < b:
            b -= a; reverse.append('L')
        else:
            a -= b; reverse.append('R')
    return ''.join(reversed(reverse))


def sb_word(value):
    left, right = (0, 1), (1, 0)
    answer = ''
    while True:
        a, b = left[0]+right[0], left[1]+right[1]
        if value == F(a, b):
            return answer
        if value < F(a, b):
            right = (a, b); answer += 'L'
        else:
            left = (a, b); answer += 'R'


def tree_controls():
    coordinates = 0
    longest = 0
    for denominator in range(3, 152, 2):
        for numerator in range(1, (denominator+1)//2):
            if gcd(numerator, denominator) != 1:
                continue
            pair = denominator-numerator, numerator
            word = decode(pair)
            rebuilt = (2, 1)
            for letter in word:
                old_g = F(rebuilt[1], sum(rebuilt))
                rebuilt = child(rebuilt, letter)
                need(F(rebuilt[1], sum(rebuilt)) == spectral_child(old_g, letter), 'Spectral branch formula')
            need(rebuilt == pair, 'Complete rational decoder replay')
            coordinates += 1
            longest = max(longest, len(word))
    level = {(2, 1)}
    seen = set(level)
    counts = [1]
    for depth in range(1, 8):
        new_level = {child(pair, letter) for pair in level for letter in 'ABC'}
        need(len(new_level) == 3**depth and not (new_level & seen), 'Ternary tree has no collisions')
        for pair in new_level:
            need(len(decode(pair)) == depth, 'Tree-depth clock')
        level = new_level; seen |= level; counts.append(len(level))
    need(decode((8, 7)) == 'AAAAAA', 'Actual spectral pair address')
    for j in range(30):
        pair = (j+2, j+1)
        need(decode(pair) == 'A'*j, 'Full consecutive Euclid spine')
        m, n = pair
        triple = (m*m-n*n, 2*m*n, m*m+n*n)
        need(triple[0] == m+n and triple[0]**2+triple[1]**2 == triple[2]**2,
             'Odd-leg denominator identity')
        need(abs(2*n-(m+n)) == abs(2*m-(m+n)) == 1, 'Two Farey adjacencies to one-half')
        need(abs(n*(m+n)-m*(m+n)) == m+n, 'The complementary fractions are not neighbors')
    return {'odd_denominator_bound': 151, 'reduced_coordinates': coordinates,
            'maximum_decoded_depth_in_box': longest, 'tree_level_counts': counts,
            'gap_7_over_15_address': 'AAAAAA', 'euclid_pair': [8, 7],
            'primitive_triple': [15, 112, 113]}


def network_controls():
    rows = {}
    for value in (F(4, 5), F(7, 15), F(8, 15), F(8, 7)):
        word = cw_word(value)
        index = 1
        c = F(1)
        for letter in word:
            index = 2*index+(letter == 'R')
            c = c+1 if letter == 'R' else c/(1+c)
        need(c == value, 'Exact unit series/parallel realization')
        rows[str(value)] = {'CW_word': word, 'heap_index': index,
                            'unit_edges_in_realization': 1+len(word), 'SB_word': sb_word(value)}
    need(rows['4/5']['heap_index'] == 30 and rows['7/15']['heap_index'] == 508
         and rows['8/15']['heap_index'] == 258 and rows['8/7']['heap_index'] == 129,
         'Recovered exact heap addresses')
    need(rows['7/15']['SB_word'] == 'LLRRRRRR' and rows['8/15']['SB_word'] == 'LRLLLLLL',
         'Stern-Brocot paths differ from Calkin-Wilf')
    for a in range(1, 41):
        for b in range(1, 41):
            value = F(a, b)
            c = F(1)
            for letter in cw_word(value):
                c = c+1 if letter == 'R' else c/(1+c)
            need(c == value, 'All bounded rational conductances decode')
    g = F(7, 15)
    clocks = []
    for pi in ((F(3, 7), F(4, 7)), (F(1, 2), F(1, 2)), (F(99, 100), F(1, 100))):
        p = [[(1-g)*(i == j)+g*pi[j] for j in range(2)] for i in range(2)]
        need(all(sum(row) == 1 for row in p), 'Markov reconstruction')
        need(pi[0]*p[0][1] == pi[1]*p[1][0], 'Retained weighted reversibility')
        need(p[0][1]+p[1][0] == g, 'Same gap does not fix masses')
        clocks.append(str(1/p[0][1]))
    need(clocks == ['15/4', '30/7', '1500/7'], 'Actual source hitting clocks differ')
    first_zero = Counter()
    for block in range(15):  # 1111 is the rejected block.
        bits = [(block >> k) & 1 for k in (3, 2, 1, 0)]
        first_zero[bits.index(0)+1] += 1
    need([first_zero[r] for r in range(1, 5)] == [8, 4, 2, 1], 'Four-bit valuation residue counts')
    for r in range(1, 5):
        residue_mass = F(2**(4-r), 15)
        need(F(1, 2**r)/(1-F(1, 16)) == residue_mass, 'Exact periodic valuation mass')
        for k in range(13):
            need(F(1, 2**(4*k+r))/residue_mass == F(15, 16)*F(1, 16)**k,
                 'Retained rejection clock conditional law')
    need(F(1)+4*F(1, 16)/(1-F(1, 16)) == F(19, 15), 'Conditional valuation mean')
    return {'address_table': rows, 'bounded_fraction_pairs': 1600, 'same_gap_hitting_times': clocks,
            'four_bit_first_zero_counts': [8, 4, 2, 1], 'conditional_residue_one_mean': '19/15'}


def main():
    out = {'status': 'PASS', 'native': native_controls(), 'tree': tree_controls(),
           'networks': network_controls(), 'checks': CHECKS}
    text = json.dumps(out, indent=2, sort_keys=True)+'\n'
    target = Path(__file__).with_name('receipt.json')
    if '--write' in sys.argv:
        target.write_text(text, encoding='utf-8')
    else:
        need(json.loads(target.read_text(encoding='utf-8')) == out, 'Frozen receipt')
    print(json.dumps({'status': 'PASS', 'checks': CHECKS, 'reduced_coordinates': out['tree']['reduced_coordinates'],
                      'native_conductance': out['native']['shared_raw_conductance']}))


if __name__ == '__main__':
    main()
