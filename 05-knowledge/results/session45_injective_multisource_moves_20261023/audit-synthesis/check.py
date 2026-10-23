#!/usr/bin/env python3
"""Independent synthesis cross-check, not a second native-atlas authentication."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--repo', type=Path, required=True)
p.add_argument('--bundle', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
checks = 0
inputs = {}

def need(ok, why):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(why)

def path(packet, name):
    for base in (a.bundle / packet, a.bundle / ('session45-' + packet)):
        if (base / name).is_file():
            return base / name
    raise FileNotFoundError((packet, name))

def read(packet, name):
    f = path(packet, name)
    raw = f.read_bytes()
    inputs[packet + '/' + name] = hashlib.sha256(raw).hexdigest()
    return json.loads(raw)

def pair(packet, normal, optimized):
    d = read(packet, normal)
    e = read(packet, optimized)
    need(path(packet, normal).read_bytes() == path(packet, optimized).read_bytes(), packet + ' modes differ')
    need(d['status'] in ('PASS', 'FINITE-EXACT PASS', 'PASS: independently certified original filtration, complete marked dyadic tower and shared-row cut/carry laws', 'PASS AFTER EDGEWISE-HYPOTHESIS REPAIR'), packet + ' reported failure: ' + str(d['status']))
    return d


manifest_counts = {}
for packet in ('dyadic-constraints','multisource-trades','source-circuits','injective-completion','balanced-circuits','audit-dyadic','audit-completion','audit-large-circuits','incoming-integration'):
    manifest = read(packet, 'MANIFEST.json')
    entries = manifest.get('files', manifest.get('publish_only'))
    if isinstance(entries, dict):
        entries = [{'path':name, 'sha256':digest} for name,digest in entries.items()]
    manifest_dir = path(packet,'MANIFEST.json').parent
    for entry in entries:
        if 'base' in entry:
            options = [manifest_dir / manifest[key][entry['base']] / entry['path'] for key in ('base_paths_relative_to_manifest','development_base_paths_relative_to_manifest') if key in manifest]
            target = next((f for f in options if f.is_file()), None)
            need(target is not None, 'manifest dependency exists: ' + entry['path'])
        else:
            target = manifest_dir / entry['path']
        payload = target.read_bytes()
        need(hashlib.sha256(payload).hexdigest() == entry['sha256'], 'frozen manifest hash: ' + packet + '/' + entry['path'])
        if 'bytes' in entry:
            need(len(payload) == entry['bytes'], 'frozen manifest size: ' + packet + '/' + entry['path'])
    manifest_counts[packet] = len(entries)

raw = (a.repo / '05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json').read_bytes()
E = json.loads(raw)['order8']
inputs['repo:original_electrical_receipt'] = hashlib.sha256(raw).hexdigest()
X = E['higher_layers'][-1]['lifts'][0]
Y = E['radical_lifts_mod4'][:4]
for pinfile in ('SOURCE_PINS.json', 'WORKSPACE_SOURCE_PINS.json'):
    pins = read('incoming-integration', pinfile)
    for name, digest in pins['files'].items():
        need(hashlib.sha256((a.repo / name).read_bytes()).hexdigest() == digest, 'incoming source pin: ' + name)

incoming = pair('incoming-integration', 'normal.json', 'optimized.json')
need(incoming['checks'] == 7063, 'bounded integration check count')
dyadic = pair('audit-dyadic', 'normal.json', 'optimized.json')
f = dyadic['filtration']
need(f['paired_charge_order'] == 2048 and f['node_label_fibres'] == 1695, 'group and occupied-label sizes')
need(sum(int(n)*c for n,c in f['fibre_size_distribution'].items()) == 3528, 'marked fibre population')
need(sum(f['fibre_size_distribution'].values()) == 1695, 'marked fibre count')

# An independent meet-in-the-middle state reconstruction from original labels.
# This does not import any producer or prior auditor implementation.
cases = read('balanced-circuits', 'candidates.json')['found']
old = pair('audit-large-circuits', 'normal.json', 'optimized.json')
old_cases = {c['slots']: c for c in old['circuits']}
reports = []
expected = {6:(2,6,2), 8:(2,8,2), 10:(8,4,64), 12:(16,4,2048), 14:(68,4,16384), 16:(196,4,65536)}

def sums(sources, cyclic):
    groups = {}
    for mask in range(1 << len(sources)):
        selected = [s for j,s in enumerate(sources) if mask >> j & 1]
        key = (len(selected), tuple(sum(y[s] for s in selected) % 2 for y in Y), sum(X[s] for s in selected) % cyclic)
        groups.setdefault(key, []).append(mask)
    return groups

for k, case in sorted(cases.items(), key=lambda item:int(item[0])):
    plus, minus = case['plus'], case['minus']
    half, n = len(plus), case['slots']
    need(len(minus) == half and 2*half == n and len(set(plus+minus)) == n, 'disjoint balanced source support')
    full = (1 << n)-1
    accepted = {}
    for cyclic in (4,128):
        left, right = sums(plus, cyclic), sums(minus, cyclic)
        accepted[cyclic] = {x | (y << half) for key,L in left.items() for x in L for y in right.get(key,[])}
    H8, H256 = accepted[4], accepted[128]
    need(H256 == {0,full}, 'complete higher-bank relation')
    count, budget, closure = expected[n]
    need(len(H8) == count == old_cases[n]['H8_accepted'], 'independent H8 count')
    need(old_cases[n]['balanced_states'] == math.comb(n,half), 'complete target-margin population')
    eq = [(i,j) for i,j in itertools.combinations(range(n),2) if all((m>>i&1)==(m>>j&1) for m in H8)]
    components = [{i} for i in range(n)]
    for i,j in eq:
        ci = next(c for c in components if i in c)
        cj = next(c for c in components if j in c)
        if ci is not cj:
            ci.update(cj)
            components.remove(cj)
    need(2**len(components) == closure, 'binary equality-projection closure')
    for i,j in itertools.combinations(range(n),2):
        values = {(m>>i&1,m>>j&1) for m in H8}
        need(values in ({(0,0),(1,1)}, {(0,0),(0,1),(1,0),(1,1)}), 'literal pair projection type')
    attained = None
    for bound in range(2,n+1,2):
        reached, todo = {0}, [0]
        while todo:
            m = todo.pop()
            found = {other for other in H8-reached if (other^m).bit_count() <= bound}
            reached.update(found)
            todo.extend(found)
        if full in reached:
            attained = bound
            break
    need(attained == budget == old_cases[n]['H8_minimax_slot_budget'], 'exact independently reconstructed H8 move budget')
    debt = 2*(sum(X[s] for s in plus)-sum(X[s] for s in minus))
    need(debt == old_cases[n]['integer_X_target_debt'], 'integer X debt convention')
    need(old_cases[n]['ordinary_owner_clock_cost'] == 2*n and old_cases[n]['closed_ports_union'] == 8*n, 'clock and resource count conventions')
    if n == 16:
        hostile = (0,1806,6500)
        majority = (hostile[0]&hostile[1]) | (hostile[1]&hostile[2]) | (hostile[2]&hostile[0])
        need(all(m in H8 for m in hostile) and majority == 260 and majority not in H8, 'majority obstruction')
        need((majority & ((1<<half)-1)).bit_count() == (majority>>half).bit_count(), 'majority retains target margin')
        selected = [s for j,s in enumerate(plus+minus) if majority>>j&1]
        signs = [1 if j<half else -1 for j in range(n) if majority>>j&1]
        need(4*sum(sign*y for sign,y in zip(signs,[Y[0][s] for s in selected])) % 8 == 4, 'genuine 2Y0 current separator')
    reports.append({'slots':n, 'H8_states':len(H8), 'H256_states':len(H256), 'pairwise_closure':closure, 'H8_minimax_budget':attained, 'integer_X_debt':debt})

maximum = pair('audit-large-circuits', 'maximum.normal.json', 'maximum.optimized.json')
need(maximum['nine_subsets'] == math.comb(20,9) == 167960, 'nine-subset universe')
need(maximum['eighteen_relations'] == sum(x['relations'] for x in maximum['omitted_pair_counts']) == 1129, 'boundary relation bookkeeping')
need(sum(maximum['proper_subrelation_side_size_counts'].values()) == 1129, 'every boundary relation has a proper witness')
need(maximum['sharp_maximum_support'] == 16 and maximum['twenty_source_fivebit_parity'] != 0, 'scope-specific upper boundary')
need(old['all_pool_union_ports'] == 160 and old['verified_full_image_equations'] == 7, 'full pool and direct witnesses')

census = pair('source-circuits', 'receipt.normal.json', 'receipt.optimized.json')
u = census['universe']
need((u['canonical_sources'],u['actual_source_paths'],u['black_endpoint_arc_menus'],len(u['capacity2_menus'])) == (3528,1102392,700,8), 'finite source census counts')
need(u['genuine_three_source_four_slot_repairs'] == [[2625,2628,1981,2458,2659]] and u['genuine_two_source_four_slot_repairs'] == [], 'scoped genuine repair census')
need(census['joint_integer_X_debt'] == {'2625':256,'2628':-256}, 'three-source integer debt')
need(not census['raw_X_mod512_harmonic'] and census['scaled2X_mod512_harmonic_and_preserved'], 'raw versus scaled harmonic guard')
completion = pair('audit-completion', 'independent.normal.json', 'independent.optimized.json')
q = completion['resource_carrier_quotient_hostile']
need(q['source'] == 43 and q['normalized_word'] == 344062 and q['raw_union_ports'] == 8 and q['normalized_tile_capacity'] == 1, 'raw lift quotient hostile')
need(completion['native_hostile']['pairwise_feasible'] and not completion['native_hostile']['triple_feasible'], 'native Hall hostile')
need(math.comb(13,6) == 1716 < 2048 < math.comb(14,7) == 3432, 'middle-binomial criterion threshold')
need(math.comb(64,2) == 2016 < 2048 < math.comb(65,2) == 2080, 'distinct-label pair-sum threshold')

# Finite controls for the elementary projection obstruction and its capacity escape.
for n in range(1,9):
    for layers in range(1,5):
        W = [(bit,aux) for bit in (0,(1<<n)-1) for aux in range(layers)]
        for bound in range(n):
            reached, todo = {(0,0)}, [(0,0)]
            while todo:
                w = todo.pop()
                for v in W:
                    if (w[0]^v[0]).bit_count() <= bound and v not in reached:
                        reached.add(v); todo.append(v)
            need(all(bit == 0 for bit,_ in reached), 'neutral decoration cannot bridge projected support')
need(len({('p',0),('p',1)}) == 2 and len({p for p,_ in {('p',0),('p',1)}}) == 1, 'forgetting layers can lose original injection')
result = {'status':'PASS', 'checks':checks, 'independent_circuit_reconstructions':reports, 'source_hashes':inputs, 'verified_manifest_entries':manifest_counts, 'scope':'Independent bounded integration replay is retained separately. This checker independently reconstructs the six charge relations and H8 move budgets from original marked labels, cross-checks quoted receipts and pins, and gives finite projection controls. Native path authentication, full filtration, and exhaustive 18-source upper-bound search remain the separate pinned audits, not new claims of replay here.'}
a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','checks':checks,'circuits':reports},indent=2))
