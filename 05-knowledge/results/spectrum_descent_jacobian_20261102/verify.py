"""Independent consequence checks: actual graphs, exact fibres and optimal costs."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import heapq
import json
from descent import (BLOCKS, SEED, ROOT, actual, admissible, apply, clock_lower_bound,
                     greedy, hamilton, height_exponent, jump_optimal, need,
                     optimal, predecessors, rank, single, step, g_values)

HERE = Path(__file__).resolve().parent
LETTERS = [('S', 0), ('S', 1)]+[('D', i) for i in range(16)]
B19 = [38, 60, 56, 49, 1, 16]
SEED_PARENT = [12, 121, 2, 100, 45, 69, 21]


def verify_sources():
    manifest = json.loads((HERE/'sources.json').read_text(encoding='utf8'))
    for entry in manifest['files']:
        digest = sha256((ROOT/entry['path']).read_bytes()).hexdigest()
        need(digest == entry['sha256'], 'pinned source bytes: '+entry['path'])
    return {'status': 'passed', 'pinned_files': len(manifest['files'])}


def forward(Fcap, qcap, hcap):
    distances = {SEED: 0}
    queue = [(0, SEED)]
    while queue:
        cost, state = heapq.heappop(queue)
        if cost != distances[state]:
            continue
        for letter in LETTERS:
            target = step(state, letter)
            if target[0] > Fcap or target[1] > qcap or target[3] > hcap:
                continue
            value = cost+(1 if letter[0] == 'S' else 4)
            if value < distances.get(target, 10**100):
                distances[target] = value
                heapq.heappush(queue, (value, target))
    return distances


def replay(word):
    state = SEED
    for letter in word:
        state = step(state, letter)
    return state


def inverse_controls():
    large = forward(512, 1024, 125)
    for state, distance in large.items():
        a, b = optimal(state), jump_optimal(state)
        need(a is not None and b is not None and a[0] == b[0] == distance,
             'forward Dijkstra / one-step inverse / jump inverse minimum')
        need(replay(b[1]) == state, 'jump word is an actual witness')
        need(distance <= 4*height_exponent(state[3])+state[1].bit_length(),
             'canonical finite clock bound')
    small = forward(63, 128, 25)
    universe = 0
    for h in (1, 5, 25):
        for F in range(1, 64, 2):
            for q in range(129):
                for g in g_values(h):
                    state = (F, q, g, h)
                    if not admissible(state):
                        continue
                    universe += 1
                    a, b = optimal(state), jump_optimal(state)
                    need((a is not None) == (b is not None) == (state in small),
                         'complete bounded membership including negative states')
                    if b is not None:
                        need(a[0] == b[0] == small[state], 'small minimum clock')
    hostiles = []
    for state in ((17,12,0,5), (49,48,0,5), (19,14,0,5), (15,10,0,5)):
        answer = jump_optimal(state)
        options = predecessors(state)
        first = min(options, key=lambda item: (rank(item[0]), item))
        rows = actual(answer[1])
        need(hamilton(rows) == state[0], 'literal hostile output H')
        parent_rows = [row >> 1 for row in rows[1:]]
        parent_paths, parent_covers = literal_parent(parent_rows)
        cut = rows[0] >> 1
        independent = (state[0],
                       sum(bool(Q) and (cut >> Q[0] & 1) for P,Q in parent_covers),
                       sum(cut >> P[-1] & 1 for P in parent_paths), len(parent_paths))
        need(independent == state, 'literal hostile full marked state')
        hostiles.append({'state': state, 'minimum_vertices': answer[0]+1,
                         'word': answer[1], 'rows': rows, 'predecessors': options,
                         'minimum_rank_predecessor': first, 'greedy': greedy(state)})
    need(all(h['greedy'] is None for h in hostiles[:3]), 'three greedy failures')
    need(hostiles[3]['greedy'][0] == hostiles[3]['minimum_vertices'],
         'greedy uses one unnecessary vertex')
    # Large numeric states are audited by the actual serialized generator words;
    # no exponential-size Hamiltonian DP is claimed for these long graphs.
    long_words = []
    for r in (1, 2):
        word = tuple([('S',1)]*13+[('D',2)]*r+[('S',b) for b in [0,1,1,0,0,1]*4])
        state = replay(word)
        answer = jump_optimal(state)
        need(answer is not None and replay(answer[1]) == state, 'large-state inverse')
        need(answer[0] <= sum(1 if k=='S' else 4 for k,a in word), 'large-state cost')
        long_words.append({'source_word': word, 'state': state,
                           'optimal_vertices': answer[0]+1, 'optimal_word': answer[1]})
    return {'positive_window': [512,1024,125], 'positive_states': len(large),
            'complete_window': [63,128,25], 'admissible_states_tested': universe,
            'actual_small_states': len(small), 'hostiles': hostiles,
            'large_state_controls': long_words,
            'one_step_cache': optimal.cache_info()._asdict(),
            'jump_cache': jump_optimal.cache_info()._asdict()}


def determinant(matrix):
    n = len(matrix)
    result = 0
    for p in permutations(range(n)):
        sign = (-1)**sum(p[i] > p[j] for i,j in combinations(range(n),2))
        value = sign
        for i in range(n):
            value *= matrix[i][p[i]]
        result += value
    return result


def smith(matrix):
    divisors = [1]
    for k in range(1,5):
        common = 0
        for ii in combinations(range(4),k):
            for jj in combinations(range(4),k):
                common = gcd(common, determinant([[matrix[i][j] for j in jj] for i in ii]))
        divisors.append(abs(common))
    return [divisors[i+1]//divisors[i] for i in range(4)]


def fibres():
    records = []
    for b in BLOCKS:
        H,E,p,q,s,t,F = (b[key] for key in ('H','E','p','q','s','t','F'))
        M = [[H,p,H-s,F-H], [0,E,0,q], [0,0,1,0], [0,0,0,H]]
        factors = smith(M)
        d = gcd(5,p*q-32*(F-5))
        need(factors == [1,1,d,800//d], 'all-minor Smith certificate')
        need(determinant(M) == H*H*E, 'lifted Jacobian determinant')
        divisor = gcd(H,H-s)
        kernel = (-(H-s)//divisor, 0, H//divisor, 0)
        need(sum(a*x for a,x in zip((H,p,H-s,F-H),kernel)) == 0,
             'primitive response-shear kernel')
        records.append({'cut': b['mask'], 'Smith': factors, 'kernel': kernel})
    target = (15,10,0,5)
    need([s for s,l in predecessors(target)] == [(5,5,5,5),(7,5,3,5),(9,5,1,5)],
         'three-point admissible integer fibre')
    need(jump_optimal((7,5,3,5)) is None, 'middle fibre point is impossible')
    need(jump_optimal((5,5,5,5)) and jump_optimal((9,5,1,5)), 'fibre endpoints realized')
    return {'lifted_D_profiles': records,
            'Smith_histogram': dict(Counter(str(r['Smith']) for r in records)),
            'fibre_target': target, 'fibre_points': [s for s,l in predecessors(target)]}


def literal_parent(rows):
    n = len(rows)
    paths = []
    covers = []
    for order in permutations(range(n)):
        bad = [i+1 for i,(u,v) in enumerate(zip(order,order[1:])) if not(rows[u] >> v & 1)]
        if not bad:
            paths.append(order)
        if len(bad) <= 1:
            for k in bad if bad else range(n+1):
                covers.append((order[:k], order[k:]))
    return paths, covers


def append_zero_block(rows, block):
    # Apex is vertex0; the new parent block follows all old parent vertices.
    n, m = len(rows), len(block)
    out = list(rows)+[0]*m
    for v in range(1,n):
        out[v] |= ((1 << m)-1) << n
    for j in range(m):
        out[n+j] = (block[j] << n) | 1
    return out


def amplified():
    seed = [67 << 1]+[(r << 1) | (not(67 >> v & 1)) for v,r in enumerate(SEED_PARENT)]
    paths,covers = literal_parent(SEED_PARENT)
    h = len(paths)
    g = sum(67 >> p[-1] & 1 for p in paths)
    q = sum(bool(Q) and (67 >> Q[0] & 1) for P,Q in covers)
    F = hamilton(seed)
    need((F,q,g,h) == (449,512,54,117), 'independent deficit-nine seed')
    atoms = []
    for name,rows in (('D',(8,1,3,6)), ('B19',B19)):
        paths,covers = literal_parent(rows)
        H,E = len(paths),len(covers)
        need((H,E) == ((5,32) if name=='D' else (19,256)), 'literal amplifier counts')
        for root in range(len(rows)):
            reached = 1 << root
            while True:
                expanded = reached
                for vertex in range(len(rows)):
                    if reached >> vertex & 1:
                        expanded |= rows[vertex]
                if expanded == reached:
                    break
                reached = expanded
            need(reached == (1 << len(rows))-1, 'amplifier is strong')
        atoms.append({'name':name,'rows':rows,'H':H,'E':E})
    old = ROOT/'05-knowledge/results/square_sum_spectrum_overnight_20261031h'
    banks=[]
    for name in ('strong_profiles256.json','strong_profiles256.ears.json'):
        data=(old/name).read_bytes(); bank=json.loads(data)
        types=sorted({tuple(r['profile'][:2]) for r in bank['profiles']
                      if r['profile'][1] & (r['profile'][1]-1) == 0})
        need(types == [(1,2),(5,32),(19,256)], 'bounded strong dyadic-cover types')
        orders=sorted({(r['profile'][0], r['profile'][1], len(r['rows']))
                       for r in bank['profiles']
                       if r['profile'][1] & (r['profile'][1]-1) == 0})
        need(orders == [(1,2,1),(5,32,4),(19,256,6)], 'dyadic representatives sizes')
        banks.append({'path':str((old/name).relative_to(ROOT)), 'sha256':sha256(data).hexdigest(),
                      'types':types, 'orders':orders, 'accepted_inherited_census':True})
    controls=[]
    for a,b,r in product(range(3),range(2),range(3)):
        if a+b+r == 0:
            continue
        state=(F,q,g,h); rows=list(seed)
        sequence=[((8,1,3,6),5,32)]*a+[(B19,19,256)]*b+[((0,),1,2)]*r
        for block,H,E in sequence:
            oldF,oldq,oldg,oldh=state
            state=(H*oldF+(E-H)*oldq+H*oldg,E*oldq,0,H*oldh)
            rows=append_zero_block(rows,block)
        k=9+5*a+8*b+r
        count=(1 << k)-9*5**a*19**b
        need(state==(count,1 << k,0,117*5**a*19**b), 'all-parameter amplifier formula')
        need(len(rows)==k-a-2*b-1, 'physical vertex clock')
        direct = len(rows) <= 16
        if direct:
            need(hamilton(rows)==count, 'literal amplified tournament Hamilton count')
        controls.append({'a':a,'b':b,'r':r,'k':k,'H':count,'vertices':len(rows),
                         'direct_adjacency_DP':direct})
    # Every cut of each fixed atom: only the empty cut can retain dyadic slope.
    cut_controls=0
    for atom in atoms:
        paths,covers=literal_parent(atom['rows'])
        E=atom['E']
        for cut in range(1 << len(atom['rows'])):
            digit=sum(bool(Q) and (cut >> Q[0] & 1) for P,Q in covers)
            out=E*512+digit*117
            need((out & (out-1)==0) == (cut==0), 'dyadic-slope cut rigidity')
            cut_controls+=1
    return {'seed_state':[F,q,g,h], 'seed_rows':seed, 'atoms':atoms,
            'inherited_census_comparison':banks, 'family_controls':controls,
            'all_cut_controls':cut_controls}


if __name__ == '__main__':
    integrity = verify_sources()
    report={'source_integrity':integrity, 'inverse':inverse_controls(),
            'fibres':fibres(),'amplification':amplified()}
    (HERE/'verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf8')
    print(json.dumps({'status':'passed','positive_optimal_states':report['inverse']['positive_states'],
                      'complete_membership_states':report['inverse']['admissible_states_tested'],
                      'Smith_types':report['fibres']['Smith_histogram'],
                      'actual_family_graphs':sum(r['direct_adjacency_DP'] for r in report['amplification']['family_controls'])}))
