"""Exact fixed-domain residue tables and three complete connected-complement clocks.
Only 14886, 14880, 14874 are newly consumed. No producer module is imported.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import gcd, comb, lcm
from collections import Counter
import hashlib
import json
import sys

sys.stdout.reconfigure(newline='\n')
G = 0
CLOCKS = (14886, 14880, 14874)
PROFILE_NAME = 'overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
PROFILE_PIN = '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
OLD_NAME = 'continuing6_20260906_lrc_near_resonance_certificate.json'
OLD_PIN = 'c146d73c149c448e744e138aa3bb6f8c286748be81353b8a6605adc47a129117'


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
    need(all(gcd(d,7)==1 for d in allowed), 'all inherited states are septimal units; residue transport domain')
    state_lcm = lcm(*allowed)
    need(state_lcm == 5388292800, 'same inherited arithmetic state modulus')
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

    # These upper bounds are an audited input, not re-proved by the owners.
    global_path = here.parent/'05-knowledge/results/third_20260906_grid_full_words.out' if here.name == '04-computation' else Path('C:/w/s0905/05-knowledge/results/third_20260906_grid_full_words.out')
    global_rows = [json.loads(line[7:]) for line in global_path.read_text().splitlines() if line.startswith('GLOBAL ')]
    global_caps = [[row[0],row[2]] for row in global_rows]
    need(global_caps == [[1,210],[2,270],[3,192],[4,239],[5,197],[6,224]], 'exact proved universal full-word cap input')
    for tau,clock,cost,owner,_ in global_rows:
        need(valid(owner) and clock % 7 == tau and all(clock % d == 0 for d in owner), 'inherited cap owner retains literal profiles and clock')
        need(sum(d*((clock+7*d-1)//(7*d)) for d in owner)-clock == cost, 'inherited cap attainment, separate from its upper-bound proof')

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
    consumers = {7*355-s: near_resonance(355, s) for s in (4, 5, 6)}
    expected = {
        14886: ([1,2,3,6,9,18], 792, 221, 21, 44, 32, (1,2,2,6,6,9,18)),
        14880: ([1,2,3,4,5,6,8,10,12,15,16,20,24,30,32,40,48,60], 346104, 17826, 158, 95, 77, (1,3,6,20,30,32,48)),
        14874: ([1,2,3,6], 120, 61, 10, 20, 20, (1,2,2,3,6,6,6)),
    }
    clock_records = []
    for T in CLOCKS:
        values = [d for d in allowed if T % d == 0]
        state_gcd = gcd(T,state_lcm)
        need(state_gcd == {14886:18,14880:480,14874:6}[T] and values == [d for d in allowed if state_gcd % d == 0],
             'actual clock domain is exactly its gcd with the inherited state modulus')
        want_values, want_total, want_valid, want_types, want_max, want_cond, want_owner = expected[T]
        need(values == want_values, 'exact divisor domain at authorized clock')
        critical_pair = (6, 6*gcd(T//6,355))
        total = 0
        words = []
        residue_maxima = [{} for _ in range(7)]
        residue_owners = [{} for _ in range(7)]
        maxima, owners = {}, {}
        for S in combinations_with_replacement(values, 7):
            total += 1
            if not valid(S):
                continue
            E = sum(d*((T+7*d-1)//(7*d)) for d in S)-T
            need(7*E == sum(d*((-(T//d)) % 7) for d in S), 'literal ceiling versus residue cost for every valid word')
            words.append([list(S), E])
            present_pairs = set(combinations(S, 2))
            for pair in present_pairs:
                if E > maxima.get(pair, -1):
                    maxima[pair], owners[pair] = E, S
            for tau in range(7):
                sevenE = sum(d*((-tau*pow(d,-1,7)) % 7) for d in S)
                need(sevenE % 7 == 0, 'every fixed-domain residue cost is integral')
                Etau = sevenE//7
                if tau == T % 7:
                    need(Etau == E, 'full native ceiling equals the residue table at actual clock')
                for pair in present_pairs:
                    if Etau > residue_maxima[tau].get(pair,-1):
                        residue_maxima[tau][pair], residue_owners[tau][pair] = Etau,S
        need(total == comb(len(values)+6, 7) == want_total and len(words) == want_valid,
             'complete unpruned multiset universe and full-word filter')
        M = max(E for S, E in words)
        need(M == want_max and len(maxima) == want_types, 'global and all conditional costs')
        need(maxima[critical_pair] == want_cond and owners[critical_pair] == want_owner,
             'exact compatible owner for the near-resonance pair')
        for pair, S in owners.items():
            need(valid(S) and all(S.count(d) >= pair.count(d) for d in pair),
                 'attaining owner obeys every full profile and pair multiplicity')
            need(sum(d*((T+7*d-1)//(7*d)) for d in S)-T == maxima[pair],
                 'attaining owner literal ceiling cost')
        margin_table = [[a, b, maxima.get((a, b)), owners.get((a, b))]
                        for a, b in combinations_with_replacement(values, 2)]
        residue_tables = []
        for tau in range(7):
            table = [[a,b,residue_maxima[tau].get((a,b)),residue_owners[tau].get((a,b))]
                     for a,b in combinations_with_replacement(values,2)]
            need(set(residue_maxima[tau]) == set(maxima), 'all residues retain exactly the same feasible forced pairs')
            for pair,word in residue_owners[tau].items():
                need(valid(word) and all(word.count(d)>=pair.count(d) for d in pair), 'residue owner keeps full profiles and forced multiplicity')
            residue_tables.append({'residue':tau,'maximum':max(residue_maxima[tau].values()),'conditional_maxima':table})
        need(residue_tables[T % 7]['conditional_maxima'] == margin_table, 'entire actual conditional table is the correct residue specialization')
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
        noncritical = min((row[5],row[0],row[1],row[2]) for row in rows if row[:3] != [6,1,355])
        need(noncritical[0] == {14886:150,14880:148,14874:150}[T], 'complete least noncritical separate credit')
        need(noncritical[0] > max(row['maximum'] for row in residue_tables), 'fixed domain cap works in every residue on these same complete candidate geometries')
        need(histogram['located'] == 1, 'one located consumer per authorized clock')
        need(T//6 in consumers and 300 > M >= want_cond, 'location alone also beats the global maximum')
        record = {'clock': T, 'allowed_divisors': values, 'state_gcd':state_gcd, 'word_universe': total,
                  'valid_words_and_costs': words, 'valid_word_count': len(words),
                  'valid_word_sha256': digest(words), 'global_maximum': M,
                  'conditional_maxima': margin_table, 'realized_margin_count': len(maxima),
                  'residue_tables': residue_tables, 'critical_pair': critical_pair,
                  'edge_gcds': es, 'base_count': len(rows), 'all_candidates': rows,
                  'candidate_sha256': digest(rows), 'conditional_first_counts': dict(histogram),
                  'global_separate_count': len(rows)-1, 'global_survivors': global_survivors,
                  'located_profile': consumers[T//6], 'located_credit': 300,
                  'global_surplus': 300-M, 'conditional_surplus': 300-want_cond}
        record['least_noncritical_credit'] = noncritical
        clock_records.append(record)
        print('CLOCK', T, 'WORDS', total, len(words), 'M', M, 'CONDITIONAL_TYPES', len(maxima),
              'BASE', len(rows), 'GLOBAL_CSEP', len(rows)-1, 'LOCATED1_CREDIT300',
              'CRITICAL_PAIR', critical_pair, 'M_PAIR', want_cond, 'RESIDUE_M', [r['maximum'] for r in residue_tables], 'COUNTS', json.dumps(dict(histogram), sort_keys=True))

    # New all-height selected-edge consumer, not a scan of another clock.
    L355,I355 = geometry(1,355)
    need(len(I355)==51 and 7*(14*25+1)<7*355 and 6*50>270,
         'all seven strict near-resonance deficits beat the inherited universal cost')
    need(arc_count(2485,L355,I355,F(1,2))==0, 'exact resonance is a genuine local pair-count hostile')
    need(6*gcd(2485,355)==2130 and 2130 not in allowed, 'the exact resonance is impossible for a supplied full-profile edge')
    need(6*sum((2486*(b-a)+L355-1)//L355-1 for a,b in I355)==306>270,
         'monotone separate credit for every n>=2486')
    critical_edge = {'edge':[6,1,355],'minimum_t':14868,'minimum_n':2478,'near_deficits':[1,7],
                     'near_credit':300,'universal_cost':270,'resonance_n':2485,'resonance_count_at_half':0,
                     'forbidden_resonance_margin':2130,'separate_n_min':2486,'separate_credit_min':306}
    need(clock_records[1]['global_maximum']==95 and clock_records[2]['residue_tables'][5]['maximum']==15,
         'same residue alone cannot transport a cost cap across different domains')

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
    need(digest(old) == OLD_PIN and old == sorted(set(old)) and len(old) == 8198 and max(old) == 14886,
         'complete current necessary set semantic pin')
    need(set(CLOCKS).issubset(old), 'delete only currently retained clocks')
    new = [t for t in old if t not in CLOCKS]
    need(len(new) == 8195 and max(new) == 14868 and set(old)-set(new) == set(CLOCKS),
         'actual set difference, not assumed clock adjacency')
    certificate = {'scope': 'exactly three new clocks; weak safety in the inherited general connected-complement setting',
                   'status': 'FINITE-EXACT complete certificate; proof status recorded in companion report',
                   'profile_sha256': PROFILE_PIN, 'atlas': pairs, 'atlas_sha256': digest(pairs),
                   'state_modulus':state_lcm,'inherited_global_caps':global_caps,'global_cap_semantic_sha256':digest(global_caps),
                   'uniform_critical_edge':critical_edge,
                   'clocks': clock_records, 'strict_lemma_controls': lemma_controls,
                   'equality_hostiles': endpoint_controls, 'multiplicity_controls': multiplicity_controls,
                   'old_scales_sha256': OLD_PIN, 'deleted_clocks': list(CLOCKS),
                   'new_scales': new, 'new_scales_sha256': digest(new)}
    dest = here.parent/'05-knowledge/results' if here.name == '04-computation' else here
    output = dest/(Path(__file__).stem+'_certificate.json')
    output.write_bytes((json.dumps(certificate, separators=(',', ':'), sort_keys=True)+'\n').encode())
    print('LEMMA strict s(14R+1)<7q: minimum2R maximum2R+1; equality q14R+1,s7: minimum2R-1')
    print('EQUALITY_ATLAS_HOSTILE q43 n294 minimum5 atalpha0 maximum7')
    print('UNIFORM_SELECTED_EDGE (6,1,355) closes for every t>=14868: near300, resonance forbiddenmargin2130, far306; other edges untouched')
    print('TOTAL_WORDS', sum(r['word_universe'] for r in clock_records),
          'VALID', sum(r['valid_word_count'] for r in clock_records),
          'BASE', sum(r['base_count'] for r in clock_records))
    print('SCALE_SET count8195 maximum14868 sha256', digest(new))
    print('CERTIFICATE_SHA256', hashlib.sha256(output.read_bytes()).hexdigest())
    print('SCOPE no further clock scan; proper six-component phase remains CITED; LRC14 OPEN')
    print('PASS', G, 'always-active exact gates; actual LF')


if __name__ == '__main__':
    main()
