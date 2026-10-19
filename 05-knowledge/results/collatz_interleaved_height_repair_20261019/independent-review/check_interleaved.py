#!/usr/bin/env python3
"""Independent exact controls for the 24-episode first-crossing theorem.

Uses no inputs, network, third-party packages, or assertions. Coefficients use
an odd-position closed sum, and short lawful residues use dyadic lifting with
literal Terras replay. These differ from the author's forward coefficient
recurrence and modular-inverse residue construction. Prints deterministic JSON.
Finite controls complement the symbolic proof audit; they are not a coverage
theorem for arbitrary positive Collatz sources.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
import hashlib
import json

CHECKS = 0
NOTE_SHA256 = '744401a10e3ce93bf3a62c5f5efdd6212f349670cc41675ee2ccf74a64a74d5c'


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def qtext(q):
    return str(q.numerator) + '/' + str(q.denominator)


def step(x):
    return (3*x + 1)//2 if x & 1 else x//2


def replay(x, n):
    bits = []
    for _ in range(n):
        bits.append(str(x & 1))
        x = step(x)
    return ''.join(bits), x


def closed_branch(word):
    # b = sum_{i odd} 2^i * 3^{number of odd positions strictly after i}.
    # Evaluate from right to left, rather than propagating forward carries.
    tail = 1
    b = 0
    for i in range(len(word) - 1, -1, -1):
        if word[i] == '1':
            b += (1 << i)*tail
            tail *= 3
        elif word[i] != '0':
            raise ValueError('Invalid parity bit')
    return tail, b, 1 << len(word)


@lru_cache(maxsize=None)
def dyadic_residue(word):
    # At stage i the two lifts differ by 2^i. Their ith iterates differ by an
    # odd integer, so exactly one matches the next prescribed parity.
    r = 0
    for i, expected in enumerate(word):
        _, y = replay(r, i)
        if str(y & 1) != expected:
            r += 1 << i
        observed, _ = replay(r, i+1)
        require(observed == word[:i+1], 'dyadic lift lost a parity guard')
    return r


K = '110110110110110101011011101010'
Ku = '001111111111111101010110000100'
F, Fu = '10', '01'
H, Hu = '110101000', '001101001'
V, Vu = K*4 + F, Ku*4 + Fu
lower = {'K': K, 'V': V, 'F': F, 'H': H}
upper = {'K': Ku, 'V': Vu, 'F': Fu, 'H': Hu}
P, C0, Q = closed_branch(K)
require((P, C0, Q) == (1162261467, 4912087879, 1073741824), 'K constants')
require(closed_branch(F) == (3, 1, 4), 'F constants')
require(closed_branch(H) == (81, 101, 512), 'H constants')
av, bv, mv = closed_branch(V)
C = Fraction(bv, av-mv)
c = Fraction(C0, P-Q)
E = C+1
require(C == Fraction(83829391903560839836419813590128166431,
                      157489106280555890461927692448373459), 'V center')
require(0 < c < C, 'ordered expansion centers')
require(12798 < 24*E < 12799, '24-episode displacement range')
require(13332 < 25*E < 13333, '25-episode displacement range')

bank_receipt = {}
for name in ('K', 'V', 'F', 'H'):
    word = lower[name]
    a, b, M = closed_branch(word)
    au, bu, Mu = closed_branch(upper[name])
    r, ru = dyadic_residue(word), dyadic_residue(upper[name])
    require((au, bu, Mu) == (a, b+M-a, M), name+' adjacent affine identity')
    require(ru == r+1, name+' independently lawful adjacent guard')
    if name in ('K', 'V'):
        for j in range(1, len(word)+1):
            require(3**word[:j].count('1') > 2**j, name+' prefix multiplier')
        require(Fraction(b, a-M) <= C, name+' expansion center')
    else:
        require(0 < a < M and 0 <= b <= M-a, name+' reset center')
        for j in range(len(word)):
            ap, bp, Mp = closed_branch(word[:j])
            require(ap*M > a*Mp, name+' proper-prefix slope difference')
            require((ap+bp)*M >= (a+b)*Mp,
                    name+' proper-prefix difference at input one')
    bank_receipt[name] = {'length': len(word), 'odd_bits': word.count('1'),
                          'source_residue': r, 'companion_residue': ru}

r = dyadic_residue(K)
r4 = dyadic_residue(K*4)
M4 = 1 << 120
require(r == 815464443, 'minimum lawful K source')
require(r4 == 588943162011782273059357097228697595, 'forced K^4 source')
a4, b4, modulus4 = closed_branch(K*4)
require((P-Q)*b4 == C0*(P**4-Q**4), 'independent K^4 geometric carry')
require(a4 == P**4 and modulus4 == Q**4 == M4, 'K^4 scale')
require(((P-Q)*r4+C0) % M4 == 0, 'alternative K^4 cylinder identity')
require(466400 < r4 < M4 and r4 % 4 == 3, 'source threshold and companion divide by four')
for ell in (1, 2, 3):
    a, b, M = closed_branch(K*ell + F)
    require(a < M, 'short expansion reset contracts')
    require(a*r+b < M*r, 'short expansion reset descends at minimum source')
ah, bh, mh = closed_branch(H)
af, bf, mf = closed_branch(F)
require(ah*mf < af*mh and bh*mf < bf*mh, 'H pointwise below F on positive inputs')

rows = []
classes = {}
thresholds = {}
for d in range(1, 12799):
    value = r4+d
    bits = []
    odd_count = 0
    for j in range(1, 121):
        bit = value & 1
        bits.append(str(bit))
        odd_count += bit
        value = step(value)
        if 3**odd_count < 2**j:
            break
    else:
        raise RuntimeError('No inherited contracting prefix for offset '+str(d))
    word = ''.join(bits)
    a, b, M = closed_branch(word)
    residue = dyadic_residue(word)
    require(j <= 92 and M4 % M == 0, 'repair depth within inherited guard')
    require((r4+d) % M == residue, 'offset lawful guard')
    require(a*(r4+d)+b == M*value, 'literal replay agrees with closed sum')
    require(a*d+b < 466400*(M-a), 'uniform original-source descent threshold')
    for lift in (1, 17):
        x = r4+lift*M4
        observed, actual = replay(x+d, j)
        require(observed == word and M*actual == a*(x+d)+b and actual < x,
                'independent lifted orbit control')
    threshold = Fraction(a*d+b, M-a)
    thresholds[d] = threshold
    # Same row encoding as the published receipt permits a cross-check while
    # coefficient and residue construction remain independently implemented.
    rows.append((d, j, a, b, M, qtext(threshold)))
    if word not in classes:
        classes[word] = {'a': a, 'b': b, 'M': M,
                         'offset_residue': (residue-r4) % M, 'offsets': []}
    classes[word]['offsets'].append(d)

words = sorted(classes)
for left, right in zip(words, words[1:]):
    require(not right.startswith(left), 'first-contracting cylinders are prefix-free')
covered = 0
for word, row in classes.items():
    M, residue = row['M'], row['offset_residue']
    first = residue if residue >= 1 else M
    expected = list(range(first, 12799, M))
    require(expected == row['offsets'], 'complete offset progression in cylinder')
    covered += len(expected)
require(covered == 12798 and len(classes) == 718, 'disjoint complete class cover')
require(max(row[1] for row in rows) == 92, 'maximum continuation depth')
require([row[0] for row in rows if row[1] == 92] == [4], 'unique deepest continuation')
max_threshold = max(thresholds.values())
require(max_threshold == Fraction(810692053755261919, 1738366812781), 'maximum threshold')
require([d for d, q in thresholds.items() if q == max_threshold] == [11808], 'threshold maximizer')
old_max = max(thresholds[d] for d in range(1, 534))
require(old_max == Fraction(116927, 13), 'predecessor threshold control')
require([d for d in range(1, 534) if thresholds[d] == old_max] == [480], 'predecessor maximizer')

hostile = 12852
hostile_word, _ = replay(r4+hostile, 120)
for j in range(1, 121):
    require(3**hostile_word[:j].count('1') >= 2**j, 'hostile inherited prefix cannot contract')
require(12798 < hostile < 25*E, 'hostile lies inside next budget only')

# Controls for parse-dependent first crossing and exact normalization. The
# written induction, not this finite enumeration, covers arbitrary lengths.
def crossing(tokens):
    length = odds = 0
    for i, token in enumerate(tokens):
        word = lower[token]
        length += len(word)
        odds += word.count('1')
        if token in ('F', 'H') and 3**odds < 2**length:
            return i+1, length
    return None


def normalize(tokens):
    result = []
    i = 0
    while i < len(tokens):
        if tuple(tokens[i:i+5]) == ('K', 'K', 'K', 'K', 'F'):
            result.append('V')
            i += 5
        else:
            result.append(tokens[i])
            i += 1
    return tuple(result)


def episodes(tokens):
    return sum(t in ('K', 'V') and (i == 0 or tokens[i-1] in ('F', 'H'))
               for i, t in enumerate(tokens))


require(crossing(('K', 'F', 'F')) == (2, 32), 'cross at individual reset, not run endpoint')
require(crossing(('V', 'F')) == (2, 124), 'internal V reset not explicit boundary')
normalization_cases = changed_cases = merged_cases = 0
for length in range(1, 9):
    for first in ('K', 'V'):
        for suffix in product(('K', 'V', 'F', 'H'), repeat=length-1):
            tokens = (first,)+suffix
            hit = crossing(tokens)
            if hit is None:
                continue
            tokens = tokens[:hit[0]]
            rewritten = normalize(tokens)
            after = crossing(rewritten)
            require(after is not None and after[1] == hit[1], 'normalization preserves first-crossing time')
            require(''.join(lower[t] for t in tokens) == ''.join(lower[t] for t in rewritten),
                    'normalization preserves lower parity word')
            require(''.join(upper[t] for t in tokens) == ''.join(upper[t] for t in rewritten),
                    'normalization preserves companion parity word')
            require(episodes(rewritten) <= episodes(tokens), 'normalization cannot increase episodes')
            normalization_cases += 1
            changed_cases += (rewritten != tokens)
            merged_cases += (episodes(rewritten) < episodes(tokens))
require(changed_cases > 0 and merged_cases > 0, 'normalization controls include actual rewrites and merges')

lam, mu = Fraction(P, Q), Fraction(3, 4)
s = lam**11*mu**3
e = c*s+1-(c+1)*mu**3
w = (K*11+F*3)*24+F
aw, bw, Mw = closed_branch(w)
require((len(w), w.count('1')) == (8066, 5089), 'long carry control length and parity count')
require(1 < s and mu*s**24 < 1, 'long carry first-crossing slope')
carry = mu*e*(s**24-1)/(s-1)+Fraction(1, 4)
require(Fraction(aw, Mw) == mu*s**24 and Fraction(bw, Mw) == carry > 534,
        'closed-sum long carry control')
for copy in range(24):
    for reset in (1, 2, 3):
        require(s**copy*lam**11*mu**reset > 1, 'every earlier explicit long-control reset')
# This extra control authenticates a full 8066-bit cylinder, which the writer's
# numerical long-carry control does not explicitly replay.
rw = (-bw*pow(aw, -1, Mw)) % Mw
observed, yw = replay(rw, len(w))
require(rw > 0 and observed == w and Mw*yw == aw*rw+bw,
        'full long-control cylinder has a positive lawful representative')

family_controls = []
for m, L, j in ((24, 69, 19), (100, 69, 19), (1000, 69, 19)):
    s = lam**L*mu**j
    e = c*s+1-(c+1)*mu**j
    require(j >= 2 and 1 < s < Fraction(4, 3), 'exact floor bracket and remaining reset after rewrite')
    require(e > c and s**m < Fraction(4, 3), 'unbounded-family sufficient inequalities')
    # As the cumulative copy endpoint is >1, and every earlier reset in that
    # copy has larger slope, these inequalities place first crossing at the
    # final extra F. The geometric carry then exceeds (3c/4)m+1/4.
    family_controls.append({'episodes': m, 'K_blocks_per_episode': L,
                            'F_blocks_per_episode': j,
                            'strict_carry_lower_bound': qtext(Fraction(3, 4)*c*m+Fraction(1, 4))})

digest = hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()
require(digest == '10a84fa213b73ab9b182d1a22ecbb99dc74036613c0a17ef80d896eec8fde77b',
        'published independent row digest')
print(json.dumps({
    'status': 'PASS',
    'audited_note_sha256': NOTE_SHA256,
    'explicit_checks': CHECKS,
    'method': 'reverse odd-position closed sums; literal replay; dyadic-lift short guards; prefix-free offset cover',
    'scope': 'bounded first-crossing bank language; no arbitrary-source coverage or general Collatz claim',
    'banks': bank_receipt,
    'bank_center': qtext(C),
    'per_episode_bound': qtext(E),
    'forced_source_residue': r4,
    'forced_source_modulus': M4,
    'positive_offsets': len(rows),
    'first_contracting_classes': len(classes),
    'maximum_repair_steps': max(row[1] for row in rows),
    'maximum_repair_offsets': [row[0] for row in rows if row[1] == 92],
    'maximum_threshold': qtext(max_threshold),
    'maximum_threshold_offsets': [d for d, q in thresholds.items() if q == max_threshold],
    'offset_rows_sha256': digest,
    'repair_time_counts': dict(sorted(Counter(row[1] for row in rows).items())),
    'hostile_25_episode_offset': hostile,
    'normalization_controls': {'maximum_original_macroblocks': 8,
                               'crossing_cases': normalization_cases,
                               'changed_cases': changed_cases,
                               'episode_merging_cases': merged_cases},
    'long_carry_control': {'steps': len(w), 'odd_bits': w.count('1'),
                           'all_72_prior_reset_endpoints_checked': True,
                           'full_guard_replayed': True, 'normalized_carry_strictly_above': 534},
    'additional_family_controls': family_controls,
}, indent=2, sort_keys=True))
