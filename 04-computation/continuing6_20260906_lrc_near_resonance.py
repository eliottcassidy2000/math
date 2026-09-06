"""Exact near-resonance pair geometry and three complete connected-complement clocks.
Only 14904, 14898, 14892 are newly consumed. No producer module is imported.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import gcd, comb
from collections import Counter
import hashlib
import json
import sys

sys.stdout.reconfigure(newline='\n')
G = 0
CLOCKS = (14904, 14898, 14892)
PROFILE_NAME = 'overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
PROFILE_PIN = '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
OLD_NAME = 'continuing5_20260906_lrc_clock16704_complete_certificate.json'
OLD_PIN = 'ddbe0e091d36d54c8f6a7c8ea631bbf363d799b9adaa2ec4fe4cf56250d11a76'


def need(ok, label):
    global G
    G += 1
    if not ok:
        raise ArithmeticError(label)


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def atlas_sum(n):
    p = 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e and (p % 3 != 2 or e > 2):
            return False
        p += 1
    return n == 1 or n % 3 == 2


def geometry(p, q):
    """Literal open-arc intersection at common denominator14pq."""
    L = 14*p*q
    a = [(max(0, 14*k*q-q), min(L, 14*k*q+q)) for k in range(p+1)]
    b = [(max(0, 14*k*p-p), min(L, 14*k*p+p)) for k in range(q+1)]
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        elif b[j][1] < a[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return L, [(out[-1][0]-L, out[0][1])] + out[1:-1]


def arc_count(n, L, intervals, alpha):
    x, d = alpha.numerator, alpha.denominator
    return sum(-((x*L-n*b*d)//(L*d))-(n*a*d-x*L)//(L*d)-1 for a, b in intervals)


def located(n, p, q, literal=False):
    L, intervals = geometry(p, q)
    walls = sorted({n*x % L for interval in intervals for x in interval})
    phases = [F(x, L) for x in walls]
    phases += [F(x+walls[(j+1) % len(walls)]+(L if j == len(walls)-1 else 0), 2*L) % 1
               for j, x in enumerate(walls)]
    rows = []
    for alpha in phases:
        count = arc_count(n, L, intervals, alpha)
        if literal:
            den = n*alpha.denominator
            raw = sum(all(14*min((v*(alpha.numerator+j*alpha.denominator)) % den,
                                 (-v*(alpha.numerator+j*alpha.denominator)) % den) < den
                          for v in (p, q)) for j in range(n))
            need(raw == count, 'native modular predicate at every wall and cell')
        rows.append([str(alpha), count])
    values = [(count, F(alpha)) for alpha, count in rows]
    low, high = min(values), max(values)
    return {'n': n, 'p': p, 'q': q, 'components': len(intervals),
            'wall_count': len(walls), 'minimum': low[0], 'min_phase': str(low[1]),
            'maximum': high[0], 'max_phase': str(high[1]), 'walls_and_cells': rows}


def near_resonance(q, s, endpoint=False):
    """Bounded controls for the separately proved all-parameter closed-gap lemma."""
    R, n = q//14, 7*q-s
    need(q >= 2 and q % 14 and s > 0 and n > 0, 'near-resonance types and full-tooth scope')
    L, intervals = geometry(1, q)
    expected = sorted(((14*k-1) % L, 2) for k in range(-R, R+1))
    need(sorted((a % L, b-a) for a, b in intervals) == expected,
         'full literal pair geometry equals the centered tooth list')
    centers = [(F(1, 2)-F(s*k, q)) % 1 for k in range(-R, R+1)]
    length = F(s, 7*q)
    distances = [min((a-b) % 1, (b-a) % 1) for a, b in combinations(centers, 2)]
    if endpoint:
        need(R >= 1 and q == 14*R+1 and s == 7, 'specified equality family')
        need(s*(14*R+1) == 7*q and min(distances) == length,
             'closed-gap contact at the equality boundary')
        need(sum(d == length for d in distances) == 1, 'only the extreme two gaps touch')
        wanted = 2*R-1
        need(arc_count(n, L, intervals, F(0)) == wanted, 'two open teeth absent at alpha0')
    else:
        need(s*(14*R+1) < 7*q, 'strict near-resonance inequality')
        need(all(d > length for d in distances), 'all closed absence gaps strictly separated')
        wanted = 2*R
    profile = located(n, 1, q, literal=True)
    need(profile['minimum'] == wanted and profile['maximum'] == 2*R+1,
         'exact all-translate extrema match the general gap proof')
    profile.update(s=s, R=R, gap_length=str(length), endpoint=endpoint)
    return profile


def main():
    here = Path(__file__).resolve().parent
    path = here/PROFILE_NAME
    if not path.is_file():
        path = Path('C:/w/s0905/04-computation')/PROFILE_NAME
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == PROFILE_PIN, 'complete inherited profile input pin')
    data = json.loads(raw)
    profiles = {int(k): {(c, tuple(w)) for c, w in row['profiles']} for k, row in data['levels'].items()}
    allowed = data['levels']['6']['gcds']
    need(len(allowed) == 42 and allowed == sorted({c for c, w in profiles[6]}), 'entire allowed state set')
    tests = [(k, I, tuple(j for j in range(7) if j not in I))
             for k in range(1, 7) for I in combinations(range(7), k)]
    need(len(tests) == 126, 'every proper nonempty positional subset, including repeated states')

    def valid(S):
        if gcd(*S) != 1:
            return False
        for k, I, J in tests:
            c = gcd(*(S[i] for i in I))
            word = tuple(sorted(gcd(c, S[j]) for j in J))
            if (c, word) not in profiles[7-k]:
                return False
        return True

    pairs = [(p, q) for p in range(1, 178) for q in range(p+1, 357-p)
             if gcd(p, q) == 1 and atlas_sum(p+q)]
    need(len(pairs) == 5855 and len({p+q for p, q in pairs}) == 94, 'complete strict pair atlas')
    need((1, 355) in pairs and (1, 43) in pairs and (1, 15) not in pairs,
         'consumer and true-atlas equality hostile; generic minimal hostile is separate')
    pair_data = []
    for p, q in pairs:
        L, I = geometry(p, q)
        nums = [2*p] + [min(2*p, p+q-14*k)
                        for k in range(1, (p+q-1)//14+1) for _ in range(2)]
        need(sorted(b-a for a, b in I) == sorted(nums), 'all atlas clipped lengths by literal intersections')
        pair_data.append((p, q, L, nums))

    # No new physical clocks are searched by these local geometric controls.
    lemma_controls = [near_resonance(2, 1), near_resonance(13, 1),
                      near_resonance(43, 6)]
    endpoint_controls = [near_resonance(15, 7, True), near_resonance(43, 7, True)]
    consumers = {7*355-s: near_resonance(355, s) for s in (1, 2, 3)}
    expected = {
        14904: ([1,2,3,4,6,8,9,12,18,23,24,27,36,46,54,72], 170544, 9593, 127, 113, 74, (1,3,6,6,8,46,54)),
        14898: ([1,2,3,6], 120, 61, 10, 12, 12, (2,2,2,3,3,6,6)),
        14892: ([1,2,3,4,6,12,17,34,51], 6435, 681, 39, 38, 32, (1,1,3,6,6,12,17)),
    }
    clock_records = []
    for T in CLOCKS:
        values = [d for d in allowed if T % d == 0]
        want_values, want_total, want_valid, want_types, want_max, want_cond, want_owner = expected[T]
        need(values == want_values, 'exact divisor domain at authorized clock')
        total = 0
        words = []
        maxima, owners = {}, {}
        for S in combinations_with_replacement(values, 7):
            total += 1
            if not valid(S):
                continue
            E = sum(d*((T+7*d-1)//(7*d)) for d in S)-T
            need(7*E == sum(d*((-(T//d)) % 7) for d in S), 'literal ceiling versus residue cost for every valid word')
            words.append([list(S), E])
            for pair in set(combinations(S, 2)):
                if E > maxima.get(pair, -1):
                    maxima[pair], owners[pair] = E, S
        need(total == comb(len(values)+6, 7) == want_total and len(words) == want_valid,
             'complete unpruned multiset universe and full-word filter')
        M = max(E for S, E in words)
        need(M == want_max and len(maxima) == want_types, 'global and all conditional costs')
        need(maxima[6, 6] == want_cond and owners[6, 6] == want_owner,
             'exact compatible owner for the near-resonance pair')
        for pair, S in owners.items():
            need(valid(S) and all(S.count(d) >= pair.count(d) for d in pair),
                 'attaining owner obeys every full profile and pair multiplicity')
            need(sum(d*((T+7*d-1)//(7*d)) for d in S)-T == maxima[pair],
                 'attaining owner literal ceiling cost')
        margin_table = [[a, b, maxima.get((a, b)), owners.get((a, b))]
                        for a, b in combinations_with_replacement(values, 2)]
        es = [e for e in range(1, 7) if T % e == 0]
        rows = []
        histogram = Counter()
        global_survivors = []
        for e in es:
            n = T//e
            for p, q, L, nums in pair_data:
                dp, dq = e*gcd(n, p), e*gcd(n, q)
                pair = tuple(sorted((dp, dq)))
                Csep = e*sum((n*a+L-1)//L-1 for a in nums)
                bound = maxima.get(pair)
                need(gcd(dp, dq) == e, 'native forced margins preserve their sheet gcd')
                if bound is None:
                    mode = 'infeasible'
                elif Csep > bound:
                    mode = 'separate'
                else:
                    mode = 'located'
                    need((e,p,q) == (6,1,355), 'only near-resonance pair needs common location')
                    need(e*consumers[n]['minimum'] > bound, 'strict located credit beats compatible cost')
                histogram[mode] += 1
                if Csep <= M:
                    global_survivors.append([e, p, q, Csep])
                rows.append([e, p, q, dp, dq, Csep, bound, mode])
        need(len(rows) == 5855*len(es) and sum(histogram.values()) == len(rows), 'complete candidate accounting')
        need(global_survivors == [[6,1,355,0]], 'global full-word cost excludes every other atlas edge')
        need(histogram['located'] == 1, 'one located consumer per authorized clock')
        need(T//6 in consumers and 300 > M >= want_cond, 'location alone also beats the global maximum')
        record = {'clock': T, 'allowed_divisors': values, 'word_universe': total,
                  'valid_words_and_costs': words, 'valid_word_count': len(words),
                  'valid_word_sha256': digest(words), 'global_maximum': M,
                  'conditional_maxima': margin_table, 'realized_margin_count': len(maxima),
                  'edge_gcds': es, 'base_count': len(rows), 'all_candidates': rows,
                  'candidate_sha256': digest(rows), 'conditional_first_counts': dict(histogram),
                  'global_separate_count': len(rows)-1, 'global_survivors': global_survivors,
                  'located_profile': consumers[T//6], 'located_credit': 300,
                  'global_surplus': 300-M, 'conditional_surplus': 300-want_cond}
        clock_records.append(record)
        print('CLOCK', T, 'WORDS', total, len(words), 'M', M, 'CONDITIONAL_TYPES', len(maxima),
              'BASE', len(rows), 'GLOBAL_CSEP', len(rows)-1, 'LOCATED1_CREDIT300',
              'M66', want_cond, 'COUNTS', json.dumps(dict(histogram), sort_keys=True))

    # The exact profile and gaps are also invariant under the physical unit multiplier.
    multiplicity_controls = []
    for T in CLOCKS:
        e, n = 6, T//6
        for h in (1, 5, 7):
            if gcd(h, n) != 1:
                continue
            for beta in (F(0), F(1, 17), F(2, 7)):
                den = T*beta.denominator
                raw = sum(all(14*min((e*h*v*(beta.numerator+j*beta.denominator)) % den,
                                    (-e*h*v*(beta.numerator+j*beta.denominator)) % den) < den
                              for v in (1, 355)) for j in range(T))
                L, I = geometry(1, 355)
                normalized = arc_count(n, L, I, (h*beta) % 1)
                need(raw == e*normalized and raw >= 300,
                     'literal common-factor physical grid has multiplicity e, including nontrivial unit')
                multiplicity_controls.append([T, h, str(beta), raw])

    prior = here.parent/'05-knowledge/results'/OLD_NAME if here.name == '04-computation' else here/OLD_NAME
    if not prior.is_file():
        prior = Path('C:/w/s0905/05-knowledge/results')/OLD_NAME
    old = json.loads(prior.read_bytes())['new_scales']
    need(digest(old) == OLD_PIN and old == sorted(set(old)) and len(old) == 8201 and max(old) == 14904,
         'complete current necessary set semantic pin')
    need(set(CLOCKS).issubset(old), 'delete only currently retained clocks')
    new = [t for t in old if t not in CLOCKS]
    need(len(new) == 8198 and max(new) == 14886 and set(old)-set(new) == set(CLOCKS),
         'actual set difference, not assumed clock adjacency')
    certificate = {'scope': 'exactly three new clocks; weak safety in the inherited general connected-complement setting',
                   'status': 'FINITE-EXACT complete certificate; proof and independent audit recorded in companion reports',
                   'profile_sha256': PROFILE_PIN, 'atlas': pairs, 'atlas_sha256': digest(pairs),
                   'clocks': clock_records, 'strict_lemma_controls': lemma_controls,
                   'equality_hostiles': endpoint_controls, 'multiplicity_controls': multiplicity_controls,
                   'old_scales_sha256': OLD_PIN, 'deleted_clocks': list(CLOCKS),
                   'new_scales': new, 'new_scales_sha256': digest(new)}
    dest = here.parent/'05-knowledge/results' if here.name == '04-computation' else here
    output = dest/(Path(__file__).stem+'_certificate.json')
    output.write_bytes((json.dumps(certificate, separators=(',', ':'), sort_keys=True)+'\n').encode())
    print('LEMMA strict s(14R+1)<7q: minimum2R maximum2R+1; equality q14R+1,s7: minimum2R-1')
    print('EQUALITY_ATLAS_HOSTILE q43 n294 minimum5 atalpha0 maximum7')
    print('TOTAL_WORDS', sum(r['word_universe'] for r in clock_records),
          'VALID', sum(r['valid_word_count'] for r in clock_records),
          'BASE', sum(r['base_count'] for r in clock_records))
    print('SCALE_SET count8198 maximum14886 sha256', digest(new))
    print('CERTIFICATE_SHA256', hashlib.sha256(output.read_bytes()).hexdigest())
    print('SCOPE no further clock scan; proper six-component phase remains CITED; LRC14 OPEN')
    print('PASS', G, 'always-active exact gates; actual LF')


if __name__ == '__main__':
    main()
