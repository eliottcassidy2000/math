"""Full native n8 Hamiltonian fibres with fixed endpoints and frozen mod8 phase.

Cost is the number of changed ACTUAL OWNER CLOCKS, not cell phases.
All possible simultaneous reorderings in each explicitly finite fibre count.
"""
from array import array
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


class DSU:
    def __init__(self,n):
        self.p = list(range(n))
    def root(self,i):
        while self.p[i] != i:
            self.p[i] = self.p[self.p[i]]
            i = self.p[i]
        return i
    def join(self,i,j):
        i,j = self.root(i),self.root(j)
        if i == j:
            return False
        self.p[j] = i
        return True
    def sizes(self):
        return sorted(Counter(self.root(i) for i in range(len(self.p))).values())


def fiber(ids, path_rows, key):
    paths = [path_rows[i]['path'] for i in ids]
    edges = sorted((sum(a != b for a,b in zip(paths[i],paths[j])),i,j)
                   for i,j in combinations(range(len(ids)),2))
    dsu = DSU(len(ids))
    tree = []
    for cost,i,j in edges:
        if dsu.join(i,j):
            tree.append([ids[i],ids[j],cost])
    need(len(tree) == len(ids)-1, 'spanning tree on every finite accepted fibre')
    budget = max((e[2] for e in tree),default=0)
    below = DSU(len(ids))
    for cost,i,j in edges:
        if cost < budget:
            below.join(i,j)
    need(len(ids) == 1 or len(below.sizes()) > 1, 'budget-minus-one is disconnected')
    upper = DSU(len(ids))
    for cost,i,j in edges:
        if cost <= budget:
            upper.join(i,j)
    need(len(upper.sizes()) == 1, 'sharp budget sufficient')
    for i,j,cost in tree:
        p,q = path_rows[i]['path'],path_rows[j]['path']
        need(p[0] == q[0] and p[-1] == q[-1], 'actual endpoint pins retained')
        need(path_rows[i]['phase_mod8'] == path_rows[j]['phase_mod8'], 'original electrical phase retained')
        need(cost == sum(a != b for a,b in zip(p,q)), 'actual owner-clock cost')
        if len(key) > 3:
            jclock = key[3]
            need(p[jclock:jclock+2] == q[jclock:jclock+2], 'actual indexed arc pin retained')
    cross = [[cost,ids[i],ids[j]] for cost,i,j in edges
             if path_rows[ids[i]]['target_merged'] != path_rows[ids[j]]['target_merged']]
    below_groups = defaultdict(list)
    for i,global_id in enumerate(ids):
        below_groups[below.root(i)].append(global_id)
    return {'pins_and_phase':list(key),'size':len(ids),'path_indices':ids,
            'distinct_targets':len({path_rows[i]['target_merged'] for i in ids}),
            'sharp_clock_budget':budget,'component_sizes_at_budget_minus_one':below.sizes(),
            'components_at_budget_minus_one':sorted(below_groups.values()),
            'spanning_tree':tree,'minimum_cross_target_edge':min(cross) if cross else None}


def paired_rectangles(sources, x, folder, edge_hash):
    by_source = {s['merged_source']:s for s in sources}
    s,r = by_source[1238],by_source[2246]
    def candidates(source, family):
        out = {}
        for f in source[family]:
            for i,j in combinations(f['path_indices'],2):
                a,b = source['path_rows'][i],source['path_rows'][j]
                key = tuple(sorted((a['target_merged'],b['target_merged'])))
                if key[0] == key[1]:
                    continue
                cost = sum(u != v for u,v in zip(a['path'],b['path']))
                witness = (cost,tuple(f['pins_and_phase']),i,j)
                if key not in out or witness < out[key]:
                    out[key] = witness
        return out
    aa,bb = [candidates(z,'nontrivial_indexed_arc_fibres') for z in (s,r)]
    common = sorted(set(aa)&set(bb))
    need(common == [(881,2762),(2650,2843)], 'complete shared-arc rectangle search')
    pa,pb = [candidates(z,'endpoint_phase_fibres') for z in (s,r)]
    need(len(set(pa)&set(pb)) == 14, 'complete endpoint-only rectangle search')
    wanted = {tuple(sorted((source,target))) for source in (1238,2246) for pair in common for target in pair}
    edge_raw = gzip.decompress((folder/'merged_edges.tsv.gz').read_bytes())
    need(sha256(edge_raw).hexdigest() == edge_hash, 'original native coloured-edge bytes')
    weights = {}
    for row in csv.DictReader(edge_raw.decode().splitlines(),delimiter='\t'):
        key = (int(row['source']),int(row['target']))
        if key in wanted:
            weights[key] = (int(row['blue']),int(row['black']))
    need(len(weights) == 8, 'all rectangle edges present in original atlas')
    def dyadic(vs):
        return 1 << min((abs(v)&-abs(v)).bit_length()-1 for v in vs if v)
    reports = []
    for pair in common:
        u,v = pair
        slots = []
        capacities = []
        color_rows = []
        for source,table in ((s,aa),(r,bb)):
            cost,pins,i,j = table[pair]
            target_to_index = {source['path_rows'][k]['target_merged']:k for k in (i,j)}
            before,after = (u,v) if source is s else (v,u)
            slots.append({'source':source['merged_source'],'pins_and_phase':pins,'cost':cost,
                          'before_index':target_to_index[before],'after_index':target_to_index[after],
                          'before':source['path_rows'][target_to_index[before]],
                          'after':source['path_rows'][target_to_index[after]]})
            for target in pair:
                actual = [k for k,row in enumerate(source['path_rows']) if row['target_merged'] == target]
                color = weights[tuple(sorted((source['merged_source'],target)))]
                need(color[0] == 0 and color[1] == 2*len(actual), 'original paired-black actual port count')
                capacities.append(len(actual))
                color_rows.append({'source':source['merged_source'],'target':target,
                                   'blue':color[0],'black':color[1],'ordinary_path_indices':actual})
        need(Counter(slot['before']['target_merged'] for slot in slots) ==
             Counter(slot['after']['target_merged'] for slot in slots), 'selected actual target capacities retained')
        need(all(slot['before']['phase_mod8'] == slot['after']['phase_mod8'] for slot in slots),
             'individual original port phases retained')
        delta = [x[u]-x[v],x[v]-x[u],x[1238]-x[2246],x[2246]-x[1238]]
        paired = [2*z for z in delta]
        need(all(z % 8 == 0 for z in paired), 'paired-black original mod8 Kirchhoff guard')
        need(any(z % 256 for z in paired), 'higher original 256 debt remains')
        need(min(capacities) == 1, 'scarce actual port forbids amplitude two')
        # The aggregate shadow alone would wrongly allow amplitude two in this orientation.
        removals = [weights[tuple(sorted((1238,u)))][1],weights[tuple(sorted((2246,v)))][1]]
        need(min(removals) >= 4, 'aggregate black shadow permits removing four')
        reports.append({'targets':pair,'frozen_values':[x[1238],x[2246],x[u],x[v]],
                        'slots':slots,'original_colored_port_rows':color_rows,
                        'unit_laplacian_defect':delta,'paired_black_laplacian_defect':paired,
                        'highest_unit_dyadic_modulus':dyadic(delta),'highest_paired_dyadic_modulus':dyadic(paired),
                        'hypothetical_amplitude_two_highest_dyadic_modulus':dyadic([2*z for z in paired]),
                        'maximum_distinct_actual_observer_amplitude':1,
                        'amplitude_two_status':'aggregate-capacity-feasible but impossible with distinct actual ports on the scarce target side'})
    # Disjoint target coordinates retain a 64*t defect for each paired rectangle.
    for a in (-1,0,1):
        for b in (-1,0,1):
            div = all((64*t) % 256 == 0 for t in (a,b))
            need(div == (a == b == 0), 'no nonzero literal-capacity combination repairs modulus256')
    return {'endpoint_only_common_pairs':len(set(pa)&set(pb)),
            'shared_arc_common_pairs':common,'rectangles':reports,
            'scope':'balanced selection of distinct actual ports plus converse partners; not a new realization of the full all-path metagraph'}


def main():
    results = HERE.parents[1]
    guard = json.loads((HERE.parent/'guards'/'receipt.json').read_text(encoding='utf-8'))['original_eight']
    folder = results/'metagraph_order8_20261023'
    potential_bytes = (results/'metagraph_order8_electrical_20261023'/'receipt.json').read_bytes()
    need(sha256(potential_bytes).hexdigest() == guard['original_potential_sha256'], 'original potential hash')
    x = json.loads(potential_bytes)['order8']['higher_layers'][-1]['lifts'][0]
    nodes = list(csv.DictReader((folder/'merged_nodes.tsv').open(),delimiter='\t'))
    merged = [None]*6880
    for row in nodes:
        for key in ('class_a','class_b'):
            merged[int(row[key])] = int(row['merged_id'])
    raw = gzip.decompress((folder/'owner_map.u16.gz').read_bytes())
    need(sha256(raw).hexdigest() == guard['owner_bytes_sha256'], 'original owner map hash')
    owner = array('H');owner.frombytes(raw)
    if sys.byteorder != 'little':
        owner.byteswap()
    need(len(owner) == 1 << 21, 'original forward word universe')
    pairs = list(combinations(range(8),2))
    chords = [(i,j) for i in range(8) for j in range(i+2,8)]
    full = (1 << 21)-1
    source_reports = []
    for source in guard['sources']:
        canon = source['canon']
        A = [[0]*8 for _ in range(8)]
        for bit,(i,j) in enumerate(pairs):
            A[i][j] = (canon >> (27-bit)) & 1
            A[j][i] = 1-A[i][j]
        paths = [p for p in permutations(range(8)) if all(A[i][j] for i,j in zip(p,p[1:]))]
        need(len(paths) == source['H'], 'all actual native Hamiltonian paths recovered')
        autos = [g for g in permutations(range(8)) if all(A[i][j] == A[g[i]][g[j]] for i,j in pairs)]
        need(len(autos) == source['aut'] == 1, 'actual native source has trivial automorphism group')
        swap_count = Counter()
        for p in paths:
            for i,j in pairs:
                changed = list(p)
                changed[i],changed[j] = changed[j],changed[i]
                valid = all(A[u][v] for u,v in zip(changed,changed[1:]))
                predecessors = ([p[i-1]] if i else [])+[p[j-1]]
                successors = [p[i+1]]+([p[j+1]] if j < 7 else [])
                local = j-i >= 3 and all(A[z][p[i]] and A[z][p[j]] for z in predecessors) \
                        and all(A[p[i]][z] and A[p[j]][z] for z in successors)
                need(valid == local, 'exact local transposition legality')
                if valid:
                    swap_count[j-i] += 1
                    if j-i == 3:
                        block = p[i:j+1]
                        need(sorted(sum(A[u][v] for v in block if u != v) for u in block) == [1,1,2,2],
                             'three-clock-span switch has genuine strong T4 core')
        path_rows = []
        for path in paths:
            mask = sum(A[path[i]][path[j]] << bit for bit,(i,j) in enumerate(chords))
            need(owner[mask] == source['ordinary_source'], 'exact original source class')
            target = merged[owner[full ^ mask]]
            value = (x[source['merged_source']]-x[target]) % 8
            reverse = tuple(reversed(path))
            revmask = sum(A[reverse[j]][reverse[i]] << bit for bit,(i,j) in enumerate(chords))
            need(merged[owner[revmask]] == source['merged_source'], 'converse actual source')
            need(merged[owner[full ^ revmask]] == target, 'converse actual target')
            need(all(A[j][i] for i,j in zip(reverse,reverse[1:])), 'converse literal Hamiltonian path')
            path_rows.append({'path':path,'forward_mask':mask,'target_merged':target,'phase_mod8':value})
        groups = defaultdict(list)
        arcgroups = defaultdict(list)
        for i,row in enumerate(path_rows):
            p = row['path']
            key = (p[0],p[-1],row['phase_mod8'])
            groups[key].append(i)
            for clock in range(1,6):
                arcgroups[key+(clock,p[clock],p[clock+1])].append(i)
        fibres = [fiber(ids,path_rows,key) for key,ids in sorted(groups.items())]
        # Nontrivial refinements only; all five internal indexed arc clocks tested.
        arcfibres = [fiber(ids,path_rows,key) for key,ids in sorted(arcgroups.items()) if len(ids)>1]
        def strongest(rows):
            return max(rows,key=lambda r:(r['sharp_clock_budget'],r['distinct_targets'],r['size'])) if rows else None
        candidates = [(r['minimum_cross_target_edge'],r) for r in arcfibres if r['minimum_cross_target_edge']]
        best_arc = min(candidates,key=lambda pair:(pair[0][0],-pair[1]['distinct_targets'],-pair[1]['size']))[1] if candidates else None
        need(best_arc is not None, 'actual same-phase cross-target move with shared indexed arc exists')
        report = {'merged_source':source['merged_source'],'ordinary_source':source['ordinary_source'],
                  'canon':canon,'H':len(paths),'aut':len(autos),'path_rows':path_rows,
                  'legal_transpositions_by_clock_span':dict(sorted(swap_count.items())),
                  'endpoint_phase_fibres':fibres,'nontrivial_indexed_arc_fibres':arcfibres,
                  'endpoint_phase_budget_histogram':dict(sorted(Counter(r['sharp_clock_budget'] for r in fibres).items())),
                  'indexed_arc_budget_histogram':dict(sorted(Counter(r['sharp_clock_budget'] for r in arcfibres).items())),
                  'hardest_endpoint_fibre':strongest(fibres),'hardest_arc_fibre':strongest(arcfibres),
                  'best_shared_arc_cross_target_fibre':best_arc}
        if source['merged_source'] == 1238:
            p = (7,1,0,4,2,6,5,3)
            z = (7,1,0,5,2,6,4,3)
            indices = [paths.index(v) for v in (p,z)]
            chosen = next(f for f in arcfibres if f['pins_and_phase'] == [7,3,7,4,2,6])
            need(chosen['path_indices'] == indices and chosen['size'] == 2 and chosen['sharp_clock_budget'] == 2,
                 'complete actual shared-middle-arc fibre')
            need([path_rows[i]['target_merged'] for i in indices] == [881,2762],
                 'original electrical target change')
            need(x[881] % 8 == x[2762] % 8 == 2 and x[1238] % 8 == 1,
                 'equal target voltage residues authenticate retained phase seven')
            report['strong_four_shared_middle_arc_witness'] = {
                'fibre':chosen,'paths':[path_rows[i] for i in indices],
                'changed_owner_clocks':[3,6],'fixed_middle_arc_clock':4,'fixed_middle_arc':[2,6],
                'outside_boundary_owners':[0,3],'source_residue':1,'target_residues':[2,2]}
        source_reports.append(report)
    rectangles = paired_rectangles(source_reports,x,folder,guard['edge_bytes_sha256'])
    return {'status':'FINITE-EXACT','checks':CHECKS,'original_potential_sha256':guard['original_potential_sha256'],
            'owner_bytes_sha256':guard['owner_bytes_sha256'],'sources':source_reports,
            'paired_selected_port_switches':rectangles,
            'scope':'all Hamiltonian paths of four original n8 sources; endpoints and optional one actual indexed internal arc pinned; original frozen port phase retained; cost is actual changed owner clocks, not cell parameters or full-wiring capacity'}


if __name__ == '__main__':
    result = main()
    receipt = HERE/'receipt.json'
    if '--write' in sys.argv:
        receipt.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    else:
        need(json.loads(receipt.read_text(encoding='utf-8')) == json.loads(json.dumps(result)), 'saved receipt')
    print(json.dumps({'status':'PASS','checks':CHECKS,'sources':[
        {k:r[k] for k in ('merged_source','H','endpoint_phase_budget_histogram','indexed_arc_budget_histogram')}
        for r in result['sources']]},indent=2))
