"""Exact paired-black charge compiler for the ORIGINAL native mod-eight bank.

The old dyadic producer is replayed, not inferred from vertex separation.
The resulting fast gate is additive and retains labelled source coordinates.
Literal port allocation and chosen per-slot phases are separate obligations.
"""
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p/'01-canon').is_dir())
RESULTS = REPO/'05-knowledge/results'
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def add_code(terms, codes, cyclic_modulus=4):
    """terms is an iterable of (actual node, signed integer coefficient)."""
    out = [0]*5
    for node, coefficient in terms:
        for j, value in enumerate(codes[node]):
            out[j] += coefficient*value
    return tuple(value % (cyclic_modulus if j == 0 else 2) for j, value in enumerate(out))


def gate_trade(sources, targets, matrix, codes, cyclic_modulus=4):
    """Exact universal mod8 gate for paired-black DISJOINT bipartite trades.

    M is new minus old paired-port incidence. Each entry changes literal
    conductance by 2*M; all integer row/column sums must be zero.
    Returning PASS says nothing about path eligibility or distinct capacity.
    """
    if len(set(sources)) != len(sources) or len(set(targets)) != len(targets):
        raise ValueError('Repeated owner label in a coordinate block')
    if set(sources) & set(targets):
        raise ValueError('Source/target blocks must be disjoint')
    if len(matrix) != len(sources) or any(len(row) != len(targets) for row in matrix):
        raise ValueError('Trade shape')
    if any(sum(row) for row in matrix) or any(sum(row[j] for row in matrix) for j in range(len(targets))):
        raise ValueError('Integer margins must be zero')
    rows = [add_code(zip(targets, row), codes, cyclic_modulus) for row in matrix]
    columns = [add_code(((s, matrix[i][j]) for i, s in enumerate(sources)), codes, cyclic_modulus)
               for j in range(len(targets))]
    return {'passes': all(not any(v) for v in rows+columns),
            'source_row_syndromes': rows, 'target_column_syndromes': columns}


def small_controls():
    """Literal all-character and image checks on composite-ring six-node graphs."""
    reports = []
    for modulus in (2, 4):
        graphs = [([(i, (i+1) % 6, 1+(i % 2)) for i in range(6)]),
                  ([(0, i, 1) for i in range(1, 6)]),
                  ([(i, j, 2) for i, j in combinations(range(6), 2)])]
        for case, edges in enumerate(graphs):
            lap = [[0]*6 for _ in range(6)]
            for a, b, w in edges:
                lap[a][a] += w; lap[b][b] += w
                lap[a][b] -= w; lap[b][a] -= w
            def image(v):
                return tuple(sum(a*b for a, b in zip(row, v)) % modulus for row in lap)
            vectors = [v+(0,) for v in product(range(modulus), repeat=5)]
            harmonics = [v for v in vectors if not any(image(v))]
            images = {image(v) for v in vectors}
            tested = 0
            for a, b, c, d in product((-1, 0, 1), repeat=4):
                matrix = ((a, b, -a-b), (c, d, -c-d), (-a-c, -b-d, a+b+c+d))
                delta = [[0]*6 for _ in range(6)]
                for i in range(3):
                    for j in range(3):
                        delta[i][3+j] = delta[3+j][i] = -2*matrix[i][j]
                by_image = all(tuple(x % modulus for x in row) in images for row in delta)
                literal = all(not any(sum(a*b for a, b in zip(row, x)) % modulus
                                     for row in delta) for x in harmonics)
                need(by_image == literal, 'General three-by-three trade image iff all-bank gate')
                tested += 1
            reports.append({'modulus': modulus, 'graph': case, 'harmonic_vectors': len(harmonics),
                            'balanced_integer_trades': tested})
    return reports


def main():
    files = {'edges': RESULTS/'metagraph_order8_20261023/merged_edges.tsv.gz',
             'electrical_receipt': RESULTS/'metagraph_order8_electrical_20261023/receipt.json',
             'dyadic_producer': RESULTS/'metagraph_order8_electrical_20261023/check.py'}
    hashes = {key: {'path': str(path.relative_to(REPO)).replace('\\', '/'),
                    'hash_mode':'raw' if path.suffix == '.gz' else 'LF-normalized',
                    'sha256': sha256(path.read_bytes() if path.suffix == '.gz' else
                                     path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
              for key, path in files.items()}
    pinfile = HERE/'SOURCE_PINS.json'
    need('--write' in sys.argv or hashes == json.loads(pinfile.read_text()), 'Frozen original source bytes')
    original = json.loads(files['electrical_receipt'].read_text())['order8']
    blue, black, edges = [], [], []
    for row in csv.DictReader(gzip.decompress(files['edges'].read_bytes()).decode().splitlines(), delimiter='\t'):
        i, j, b, k = (int(row[s]) for s in ('source', 'target', 'blue', 'black'))
        if b: blue.append((i, j, b))
        if k: black.append((i, j, k//2))
        edges.append((i, j, b+k))
    spec = importlib.util.spec_from_file_location('original_native_dyadic', files['dyadic_producer'])
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    rebuilt = producer.calculate(3528, blue, black)
    need(rebuilt == original, 'Replay full actual dyadic filtration and owner-marked lifts')
    dims = [original['critical_two_rank'], original['four_divisible_factors'], original['eight_divisible_factors']]
    need(dims == [3363, 5, 1], 'Exact three-level harmonic dimensions')
    need(sum(dims)-dims[0] == 6, 'Paired-black character quotient has cardinality 64')
    ys = original['radical_lifts_mod4']; z = original['eighth_radical_lifts_mod8'][0]
    x = original['higher_layers'][-1]['lifts'][0]
    codes = [(z[i] % 4,)+tuple(y[i] % 2 for y in ys[1:]) for i in range(3528)]
    need(all(z[i] % 2 == (ys[0][i]+ys[4][i]) % 2 for i in range(3528)), 'Eliminated binary coordinate relation')
    need(codes[3527] == (0,)*5 and len(set(codes)) == 64, 'Rooted charge code is onto all 64 elements')
    need(set(codes) == set(product(range(4), range(2), range(2), range(2), range(2))), 'Exact product carrier')
    full_codes = [(x[i] % 128,)+tuple(y[i] % 2 for y in ys[1:]) for i in range(3528)]
    need(all((x[i]-z[i]-2*(ys[0][i]+ys[3][i])) % 4 == 0 for i in range(3528)),
         'Exact conversion between primitive and modulus-eight coordinates')
    valuations = {1:original['exact_valuation_one'],2:original['exact_valuation_two'],
                  3:original['exact_valuation_three']}
    valuations.update({row['level']:row['rank'] for row in original['higher_layers']})
    need(original['complete_two_primary'] and valuations == {1:3358,2:4,3:0,4:0,5:0,6:0,7:0,8:1},
         'Complete native dyadic exponent and invariant-factor multiplicities')
    full_log_h = sum(level*count for level,count in valuations.items())
    need(full_log_h-dims[0] == 11, 'Paired full dyadic quotient has cardinality 2048')
    subgroup = {(0,)*5}; generating_nodes = []
    def plus_code(a,b):
        return tuple((u+v) % (128 if j == 0 else 2) for j,(u,v) in enumerate(zip(a,b)))
    for node,code in enumerate(full_codes):
        if code in subgroup: continue
        generating_nodes.append(node)
        cycle = [(0,)*5]; value = code
        while any(value):
            cycle.append(value); value = plus_code(value,code)
        subgroup = {plus_code(a,b) for a in subgroup for b in cycle}
        if len(subgroup) == 2048: break
    need(len(subgroup) == 2048 and generating_nodes == [1,2,4,5,8], 'Native charges generate full dyadic product code')
    need(len(set(full_codes)) == 1695, 'Vertex-code partition differs from full additive group')
    for potential, modulus in [(y, 4) for y in ys]+[(z, 8)]:
        current = [0]*3528
        for i, j, w in edges:
            value = w*(potential[i]-potential[j]); current[i] += value; current[j] -= value
        need(not any(v % modulus for v in current), 'Literal full Laplacian authenticity')
    sources = (1238, 2242, 2246, 2315)
    expected = [(3,0,0,0,0), (0,0,0,1,0), (3,0,0,1,1), (2,0,1,0,1)]
    need([codes[s] for s in sources] == expected, 'Actual four-source coordinates')
    relations = []
    for a, b, c in product(range(-4, 5), repeat=3):
        d = -a-b-c
        coefficients = (a,b,c,d)
        actual = not any(add_code(zip(sources, coefficients), codes))
        closed = all(v % 2 == 0 for v in coefficients) and (a+c) % 4 == 0
        need(actual == closed, 'Exact four-source balanced relation congruences')
        if actual and any(coefficients) and max(map(abs, coefficients)) <= 2:
            if next(v for v in coefficients if v) > 0:
                relations.append(coefficients)
    unit = [v for v in product((-1,0,1), repeat=4)
            if any(v) and sum(v) == 0 and not any(add_code(zip(sources,v),codes))]
    need(not unit, 'No nonzero original-four unit-capacity column charge')
    minimum_l1 = min(sum(map(abs,v)) for v in relations)
    minimum_relations = [v for v in relations if sum(map(abs,v)) == minimum_l1]
    need(minimum_l1 == 4 and minimum_relations == [(0,2,0,-2),(2,0,-2,0)], 'Sharp smallest original-four integer relations')
    full_relations = []
    for a,b,c in product(range(-4,5),repeat=3):
        d = -a-b-c; coefficients = (a,b,c,d)
        actual = not any(add_code(zip(sources,coefficients),full_codes,128))
        closed = all(v % 2 == 0 for v in coefficients) and (13*a-66*b-19*c) % 128 == 0
        need(actual == closed, 'Exact full-dyadic original-four kernel congruence')
        if actual and any(coefficients) and sum(map(abs,coefficients)) <= 8:
            if next(v for v in coefficients if v) > 0: full_relations.append(coefficients)
    need(full_relations == [(4,0,-4,0)], 'Smallest full-dyadic original-four relation has l1 eight')
    bycode = defaultdict(list)
    for node, code in enumerate(codes): bycode[code].append(node)
    representatives = {str(s): [node for node in bycode[codes[s]] if node != s] for s in sources}
    # The original unsafe rectangle, checked by the additive gate.
    original_gate = gate_trade((1238,2246),(881,2762),((1,-1),(-1,1)),codes)
    need(not original_gate['passes'], 'Old chosen-character-safe rectangle fails universal bank')
    need(gate_trade((1238,2246),(881,2762),((2,-2),(-2,2)),codes)['passes'], 'Old rectangle amplitude two pays electrical debt')
    # This algebraic pattern has minimum nonzero balanced trade support four.
    # Native path eligibility and physical capacity are deliberately not inferred.
    expanded = (1238,32)
    targets = next((u,v) for u,v in combinations(range(3528),2)
                   if u not in expanded and v not in expanded and codes[u] == codes[v] and x[u] % 8 == x[v] % 8)
    need(gate_trade(expanded,targets,((1,-1),(-1,1)),codes)['passes'], 'Expanded same-code rectangle algebra passes')
    image_data = json.loads((HERE/'image_witnesses.json').read_text())
    need(image_data['modulus'] == 256 and image_data['input_edges_sha256'] == hashes['edges']['sha256'],
         'Image witness modulus and original edge hash')
    need(len(image_data['solutions']) == 3, 'Both source circuits and shared target difference certified')
    witness_reports = []
    for record in image_data['solutions']:
        charge = dict(record['charge']);potential = record['potential_mod256']
        need(record['coefficient'] == 2 and len(potential) == 3528 and potential[3527] == 0,
             'Grounded original-coordinate image witness')
        current = [0]*3528
        for i,j,w in edges:
            value = w*(potential[i]-potential[j]);current[i] += value;current[j] -= value
        expected = [2*charge.get(i,0) for i in range(3528)]
        need(all((a-b) % 256 == 0 for a,b in zip(current,expected)), 'Full literal mod256 image certificate')
        correction = [[i,(a-b)//256] for i,(a,b) in enumerate(zip(current,expected)) if a != b]
        need(correction == record['integer_correction_div256_nonzero'], 'Exact saved integral image correction')
        need(not any(add_code(charge.items(),full_codes,128)), 'Image witness agrees with full additive code')
        primitive_pairing = sum(c*x[i] for i,c in charge.items())
        need(primitive_pairing == 0, 'Primitive selected character cancels as an exact integer')
        witness_reports.append({'charge':record['charge'],'coefficient':2,'modulus':256,
                                'primitive_integer_pairing':primitive_pairing})
    first = dict(image_data['solutions'][0]['charge']);second = dict(image_data['solutions'][2]['charge'])
    six_charge = {i:first.get(i,0)-second.get(i,0) for i in first.keys() | second.keys()}
    six_charge = {i:c for i,c in six_charge.items() if c}
    six_potential = [(a-b) % 256 for a,b in zip(image_data['solutions'][0]['potential_mod256'],
                                                image_data['solutions'][2]['potential_mod256'])]
    current = [0]*3528
    for i,j,w in edges:
        value = w*(six_potential[i]-six_potential[j]);current[i] += value;current[j] -= value
    need(all((current[i]-2*six_charge.get(i,0)) % 256 == 0 for i in range(3528)), 'Six-source image certificate by linear subtraction')
    need(len(six_charge) == 6 and sum(six_charge.values()) == 0, 'Six-source balanced unit relation')
    primitive_reports = []
    for charge in [first,second,six_charge]:
        terms = sorted(charge.items());proper_balanced = 0
        for mask in range(1,(1 << len(terms))-1):
            selected = [pair for i,pair in enumerate(terms) if mask >> i & 1]
            if sum(coefficient for node,coefficient in selected): continue
            proper_balanced += 1
            need(any(add_code(selected,codes)), 'Every nonempty proper balanced subblock fails even the mod8 bank')
        primitive_reports.append({'charge':terms,'support':len(terms),'proper_balanced_subblocks_rejected':proper_balanced})
    separator = [0]*3528;separator[1965] = 4
    current = [0]*3528
    incident = []
    for i,j,w in edges:
        value = w*(separator[i]-separator[j]);current[i] += value;current[j] -= value
        if i != j and 1965 in (i,j): incident.append(w)
    need(all(w % 2 == 0 for w in incident) and sum(incident) == 926, 'Actual target has only even incident conductances')
    need(not any(value % 8 for value in current), 'Sparse authentic per-slot separating character')
    need(separator[1965]-separator[2234] == 4, 'Full per-slot bank still separates two targets')
    result = {'status': 'FINITE-EXACT additive gate; native eligibility remains a separate certificate',
              'input_sha256': hashes, 'dyadic_dimensions': dims, 'effective_charge_group': 'Z/4 x (Z/2)^4',
              'effective_group_size': 64, 'coordinate_order': ['Z_mod4','Y1_mod2','Y2_mod2','Y3_mod2','Y4_mod2'],
              'codes': codes, 'original_potential_mod256': x,
              'all_dyadic': {'group':'Z/128 x (Z/2)^4','size':2048,'codes':full_codes,
                              'image_witnesses_sha256_lf':sha256((HERE/'image_witnesses.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
                              'coordinate_order':['X_mod128','Y1_mod2','Y2_mod2','Y3_mod2','Y4_mod2'],
                              'generating_nodes':generating_nodes,'distinct_vertex_codes':1695,
                              'two_primary_valuation_multiplicities':valuations,
                              'original_four_kernel':'all coefficients even;13*a-66*b-19*c=0 mod128;d=-a-b-c',
                              'minimum_original_four_relation_l1':8,'minimum_original_four_relations':full_relations,
                              'explicit_image_certificates':witness_reports,
                              'six_source_charge':sorted(six_charge.items()),
                              'primitive_selected_source_blocks':primitive_reports},
              'per_slot_universal_hostile':{'potential_nonzero':[[1965,4]],'targets':[1965,2234],
                                             'target_value_difference_mod8':4,'weighted_degree':926,
                                             'integer_current_div8_nonzero':[[i,value//8] for i,value in enumerate(current) if value]},
              'code_class_size_histogram': dict(sorted(Counter(map(len,bycode.values())).items())),
              'original_four_sources': sources, 'original_four_codes': expected,
              'balanced_original_four_kernel': 'all coefficients even and coefficient_1238+coefficient_2246 divisible by4',
              'minimum_original_four_relations_l1': minimum_l1, 'minimum_original_four_relations': minimum_relations,
              'same_code_source_candidates': representatives, 'original_rectangle_gate': original_gate,
              'algebraic_expanded_rectangle': {'sources': expanded,'targets': targets,'chosen_target_phases_equal':True,
                                               'native_eligibility_asserted':False},
              'small_controls': small_controls(), 'checks_before_replay': CHECKS}
    result = json.loads(json.dumps(result))
    if '--write' in sys.argv:
        (HERE/'SOURCE_PINS.json').write_text(json.dumps(hashes,indent=2,sort_keys=True)+'\n')
        (HERE/'receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        need(result == json.loads((HERE/'receipt.json').read_text()), 'Frozen code/quotient receipt')
    print(json.dumps({'status':'PASS','checks':CHECKS,'effective_group_size':64,
                      'original_four_unit_relation_count':len(unit),'small_controls':len(result['small_controls'])}))


if __name__ == '__main__':
    main()
