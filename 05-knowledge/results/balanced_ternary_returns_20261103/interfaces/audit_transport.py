#!/usr/bin/env python3
"""Independent direct-edge replay and subset-DP audit of the native repair."""
from pathlib import Path
from itertools import combinations, permutations
import hashlib
import json

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def validate(rows):
    n = len(rows)
    for i, row in enumerate(rows):
        need(0 <= row < (1 << n) and not (row >> i & 1), 'valid row')
    for u, v in combinations(range(n), 2):
        need(((rows[u] >> v) & 1) + ((rows[v] >> u) & 1) == 1, 'tournament pair')


def profiles(rows):
    validate(rows)
    n = len(rows)
    end = [[0]*n for _ in range(1 << n)]
    start = [[0]*n for _ in range(1 << n)]
    counts = [1] + [0]*((1 << n)-1)
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                smaller = mask ^ (1 << v)
                if not smaller:
                    end[mask][v] = start[mask][v] = 1
                else:
                    end[mask][v] = sum(end[smaller][u] for u in range(n)
                                       if rows[u] >> v & 1)
                    start[mask][v] = sum(start[smaller][u] for u in range(n)
                                         if rows[v] >> u & 1)
        counts[mask] = sum(end[mask])
    full = (1 << n)-1
    return [(counts[full], start[full][v], counts[full ^ (1 << v)],
             sum(end[mask][v]*counts[full ^ mask] for mask in range(1, full+1)))
            for v in range(n)]


def native(rows, path):
    n = len(rows)
    need(sorted(path) == list(range(n)), 'spanning native word')
    need(all(rows[u] >> v & 1 for u, v in zip(path, path[1:])), 'directed native word')
    kept = {tuple(sorted((u,v))) for u,v in zip(path,path[1:])}
    answer = list(rows)
    for u, v in combinations(range(n), 2):
        if (u,v) not in kept:
            answer[u] ^= 1 << v
            answer[v] ^= 1 << u
    return tuple(answer)


def restriction(rows, n):
    return tuple(row & ((1 << n)-1) for row in rows[:n])


def relabel(rows, image):
    answer = [0]*len(rows)
    for u, row in enumerate(rows):
        for v in range(len(rows)):
            if row >> v & 1:
                answer[image[u]] |= 1 << image[v]
    return tuple(answer)


def birth(rows, root):
    n = len(rows)
    answer = list(rows) + [((1 << n)-1) ^ (1 << root)]
    answer[root] |= 1 << n
    return tuple(answer)


def arbitrary_birth(rows, mask):
    n = len(rows)
    return tuple(row | (0 if mask >> v & 1 else 1 << n)
                 for v,row in enumerate(rows)) + (mask,)


def all_extension_counts(rows):
    """Direct Hamilton-path deletion polynomial, independent of cut weights."""
    n = len(rows)
    end = [[0]*n for _ in range(1 << n)]
    start = [[0]*n for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1 << v)
                if not rest:
                    end[mask][v] = start[mask][v] = 1
                else:
                    end[mask][v] = sum(end[rest][u] for u in range(n) if rows[u] >> v & 1)
                    start[mask][v] = sum(start[rest][u] for u in range(n) if rows[v] >> u & 1)
    full = (1 << n)-1
    w = [[0]*n for _ in range(n)]
    for mask in range(1, full):
        rest = full ^ mask
        for u in range(n):
            if end[mask][u]:
                for v in range(n):
                    w[u][v] += end[mask][u]*start[rest][v]
    linear = [start[full][v]-end[full][v]+sum(w[u][v] for u in range(n))
              for v in range(n)]
    values = [sum(end[full])]*(1 << n)
    for mask in range(1, 1 << n):
        bit = mask & -mask
        v = bit.bit_length()-1
        rest = mask ^ bit
        values[mask] = values[rest]+linear[v]-sum(w[u][v]+w[v][u] for u in range(n)
                                                 if rest >> u & 1)
    return values


def strong(rows):
    reach = [row | 1 << v for v,row in enumerate(rows)]
    for k in range(len(rows)):
        for v in range(len(rows)):
            if reach[v] >> k & 1:
                reach[v] |= reach[k]
    return all(row == (1 << len(rows))-1 for row in reach)


def main():
    receipt_path = PACKET/'transport/check.json'
    source = json.loads(receipt_path.read_text())
    c = tuple(source['source_rows'])
    need(profiles(c)[2] == (59,17,17,59), 'balanced input')
    # The preceding session's C2 is graph2032 followed by births at x=1,p=6.
    old_b = (18,52,57,51,32,1)
    need(profiles(old_b)[1] == (17,3,3,21), 'literal prior B2 source')
    old_c = birth(birth(old_b,1),6)
    old_to_word = [0,1,4,5,6,7,2,3]
    need(relabel(old_c,old_to_word) == c, 'physical prior names transported')
    exchange = source['source_exchange']
    template = tuple(exchange['template_rows'])
    need(profiles(template)[0] == (59,5,5,61), 'padded ready profile')
    need(template[7] == 0, 'declared global sink')
    target = tuple(exchange['target_rows'])
    need(relabel(template, exchange['template_to_physical']) == target,
         'full template name map')
    need(profiles(target)[2] == (59,5,5,61), 'named ready endpoint')
    # Exhaust every image preserving root0->physical2 and a root out-neighbor
    # as physical3, INCLUDING the padded sink as a permitted child image.
    minimum = 99
    candidates = 0
    by_child = {}
    for child in range(1,8):
        if not (template[0] >> child & 1):
            continue
        child_min = 99
        for perm in permutations(v for v in range(8) if v not in (2,3)):
            image = [None]*8
            image[0],image[child] = 2,3
            for v,w in zip((v for v in range(8) if v not in (0,child)),perm):
                image[v] = w
            candidate = relabel(template,image)
            distance = sum(((c[u] >> v) ^ (candidate[u] >> v)) & 1
                           for u,v in combinations(range(8),2))
            child_min = min(child_min,distance)
            candidates += 1
        by_child[child] = child_min
        minimum = min(minimum,child_min)
    need(candidates == 3600 and minimum == 3, 'complete padded-template minimum')
    need(exchange['min_flips_within_declared_template'] == minimum, 'reported minimum')
    current = c
    checked_rows = 0
    intermediate_port_edge = []
    for step in exchange['program']:
        u,v = step['reversed_pair']
        extended = tuple(step['birth_rows'])
        first = tuple(step['first_native_rows'])
        second = tuple(step['second_native_rows'])
        after = tuple(step['after_deletion_rows'])
        need(restriction(extended,8) == current, 'temporary birth preserves all names and arcs')
        need(native(extended,step['hamilton_words'][0]) == first, 'first native replay')
        need(native(first,step['hamilton_words'][1]) == second, 'second native replay')
        need(restriction(second,8) == after, 'delete only temporary name')
        differences = [(i,j) for i,j in combinations(range(8),2)
                       if ((current[i] >> j) ^ (after[i] >> j)) & 1]
        need(differences == [(u,v)], 'exact intended old arc reversal')
        states = [current,extended,first,second,after]
        counts = [profiles(rows)[0][0] for rows in states]
        need(counts == step['hamilton_counts'], 'independent subset-DP counts')
        intermediate_port_edge.extend(bool(rows[2] >> 3 & 1) for rows in states)
        checked_rows += len(states)
        current = after
    need(current == target, 'native excursion endpoint')
    need(c[2] >> 3 & 1 and target[2] >> 3 & 1, 'port edge restored at endpoints')
    need(not all(intermediate_port_edge), 'receipt must not claim stepwise port edge preservation')
    first = birth(target,2)
    final = birth(first,8)
    need(list(first) == source['branch181']['first_rows'], 'first NR birth')
    need(list(final) == source['branch181']['final_rows'], 'second NR birth')
    need(profiles(final)[8] == (181,59,59,181), 'balanced181 consequence')
    need(not strong(target) and not strong(final) and target[6] == final[6] == 0,
         'the preserved old sink remains a weak-source sidecar')
    joint = source['joint_attachment_data']
    query_checks = 0
    for key,initial in [('source',c),('after_source_exchange',target),
                        ('balanced_child_old_ports',final),('balanced_child_new_ports',final)]:
        ports = joint[key]['physical_ports']
        values = []
        for word in joint['query_words']:
            rows = initial
            for letter in word:
                incoming = sum(1 << ports[i] for i in (0,1) if letter >> i & 1)
                rows = arbitrary_birth(rows,((1 << len(rows))-1) ^ incoming)
            values.append(profiles(rows)[0][0])
            query_checks += 1
        need(values == joint[key]['values'], 'independent literal two-port query births')
    interior_path = PACKET/'interior/check.json'
    interior = json.loads(interior_path.read_text())
    edges = list(combinations(range(8),2))
    neutral = []
    for length in (1,2):
        for changed in combinations(edges,length):
            rows = list(c)
            for u,v in changed:
                rows[u] ^= 1 << v
                rows[v] ^= 1 << u
            if profiles(rows)[0][0] == 59:
                neutral.append((changed,tuple(rows)))
    need(len(neutral) == 5 and all(len(x[0]) == 2 for x in neutral), 'complete neutral flip audit')
    recorded = interior['neutral_arc_reversal_census']['two_arc_hits']
    need(neutral == [(tuple(map(tuple,x['reversed_arcs'])),tuple(x['rows'])) for x in recorded],
         'exact five mutant graphs independently recovered')
    remaining = []
    for (_,rows), recorded_row in zip(neutral,recorded):
        need(strong(rows) == recorded_row['strong'], 'strong mutant type')
        if not (rows[2] >> 3 & 1):
            need(not strong(rows), 'only excluded port-reversing mutant is weak')
            continue
        image = recorded_row['isomorphism_from_source']
        if image is not None:
            need(relabel(c,image) == rows, 'literal source isomorphism')
            continue
        image = recorded_row['isomorphism_from_converse']
        if image is not None:
            reverse = tuple(sum(1 << u for u in range(8) if c[u] >> v & 1) for v in range(8))
            need(relabel(reverse,image) == rows, 'literal converse isomorphism')
            continue
        remaining.append(rows)
    need(len(remaining) == 1, 'one genuinely different port-preserving mutant')
    remaining_rows = remaining[0]
    need(181 not in all_extension_counts(remaining_rows), 'remaining mutant one-birth exclusion')
    cube_hash = hashlib.sha256()
    for cut in range(256):
        child = arbitrary_birth(remaining_rows,cut)
        values = all_extension_counts(child)
        need(181 not in values, 'remaining mutant full two-birth cube exclusion')
        cube_hash.update(json.dumps(values,separators=(',',':')).encode())
    print(json.dumps({'status':'PASS', 'scope':'Independent direct native replay and subset-DP counts; full3600 padded-template port maps, not all possible H59 source classes or a globally shortest native route.',
        'prior_C2_to_word202032':old_to_word,
        'physical_names':['c_2','x','reset_root','reset_child','y','c_1','birth_a','birth_b'],
        'candidate_maps':candidates,'min_flips':minimum,'min_by_template_child':by_child,
        'native_moves_checked':2*len(exchange['program']), 'stage_count_checks':checked_rows,
        'general_two_port_joint_queries_checked':query_checks,
        'endpoint_profiles':{'source':profiles(c)[2],'repair':profiles(target)[2],
                             'final':profiles(final)[8]},
        'port_edge_preserved_at_endpoints_only':True,
        'endpoint_and_balanced_child_are_strong':False,
        'retained_global_sink_name':6,
        'neutral_endpoint_audit':{'one_arc_sets':28,'two_arc_sets':378,'neutral_mutants':5,
            'port_preserving_strong_mutants':4,'source_isomorphic':2,'converse_isomorphic':1,
            'remaining_mutant':remaining_rows,'remaining_full_two_birth_cube':131072,
            'remaining_cube_sha256':cube_hash.hexdigest(),
            'scope':'Together with the separately audited original-source181 exclusion and strong-component cutoff, three is the minimum old-arc distance to an eight-name H59 endpoint retaining2->3 and admitting any induced181 extension.'},
        'input_sha256':hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
        'interior_receipt_sha256':hashlib.sha256(interior_path.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
