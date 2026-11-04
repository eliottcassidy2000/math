#!/usr/bin/env python3
"""Literal path audit of all two-birth extensions and named native repair.

No producer imports are used for the extension census.  Its endpoint tables
enumerate paths directly; the inherited native compiler is used and replayed
only in the separate source-transport certificate.
"""
from pathlib import Path
from itertools import permutations, combinations
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def graph(word):
    n = len(word)+2
    rows = [0]*n
    rows[0] |= 1 << 1
    for j, letter in enumerate(word, 2):
        for apex in (0, 1):
            u, v = (apex, j) if int(letter) >> apex & 1 else (j, apex)
            rows[u] |= 1 << v
        for k in range(j+1, n):
            rows[j] |= 1 << k
    return tuple(rows)

def literal_tables(rows):
    n = len(rows)
    starts = [[0]*n for _ in range(1 << n)]
    ends = [[0]*n for _ in range(1 << n)]
    visits = 0
    def visit(mask, first, last):
        nonlocal visits
        visits += 1
        starts[mask][first] += 1
        ends[mask][last] += 1
        todo = rows[last] & ~mask
        while todo:
            bit = todo & -todo
            todo -= bit
            visit(mask | bit, first, bit.bit_length()-1)
    for v in range(n):
        visit(1 << v, v, v)
    return starts, ends, visits

def hamilton(rows):
    return sum(literal_tables(rows)[1][-1])

def count_dp(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v]=1
    for mask in range(1,1 << n):
        for v,c in enumerate(dp[mask]):
            if c:
                todo=rows[v] & ~mask
                while todo:
                    bit=todo & -todo
                    todo-=bit
                    dp[mask | bit][bit.bit_length()-1]+=c
    return sum(dp[-1])

def joint_queries(rows,ports,suffixes):
    values=[]
    for word in suffixes:
        T=rows
        for mask in word:
            incoming=sum(1 << ports[i] for i in (0,1) if mask >> i & 1)
            T=birth(T,((1 << len(T))-1) ^ incoming)
        values.append(count_dp(T))
    return {'physical_ports':ports,'values':values}

def cuts(rows):
    """All birth counts; mask marks old vertices beaten by the newborn."""
    n = len(rows)
    starts, ends, visits = literal_tables(rows)
    full = (1 << n)-1
    joint = [[0]*n for _ in range(n)]
    for mask in range(1, full):
        right = full ^ mask
        for u, left_count in enumerate(ends[mask]):
            if left_count:
                for v, right_count in enumerate(starts[right]):
                    if right_count:
                        joint[u][v] += left_count*right_count
    values = []
    for out in range(1 << n):
        incoming = [u for u in range(n) if not (out >> u & 1)]
        outgoing = [v for v in range(n) if out >> v & 1]
        value = sum(ends[full][u] for u in incoming)
        value += sum(starts[full][v] for v in outgoing)
        value += sum(joint[u][v] for u in incoming for v in outgoing)
        values.append(value)
    return values, visits

def birth(rows, out):
    n = len(rows)
    return tuple(r | (0 if out >> i & 1 else 1 << n)
                 for i, r in enumerate(rows)) + (out,)

def native_birth(rows, root):
    return birth(rows, ((1 << len(rows))-1) ^ (1 << root))

def profile(rows, root):
    starts, ends, _ = literal_tables(rows)
    full = (1 << len(rows))-1
    counts = [sum(x) for x in ends]
    counts[0] = 1
    F = sum(ends[mask][root]*counts[full ^ mask]
            for mask in range(1, full+1))
    return (counts[full], starts[full][root], counts[full ^ (1 << root)], F)

def ready59():
    # READY §3, fifth exceptional tail: 11^2,01,10.
    letters = ((1,1), (1,1), (0,1), (1,0))
    rows = [0]*7
    rows[1] |= 1 << 0
    rows[2] |= 1 << 0
    rows[1] |= 1 << 2
    for i, (a,b) in enumerate(letters, 3):
        rows[0] |= 1 << i
        for port, bit in ((1,a), (2,b)):
            u,v = (i,port) if bit else (port,i)
            rows[u] |= 1 << v
        for j in range(i+1,7):
            rows[i] |= 1 << j
    return tuple(rows)

def relabel(rows, old_to_new):
    result = [0]*len(rows)
    for u in range(len(rows)):
        for v in range(len(rows)):
            if rows[u] >> v & 1:
                result[old_to_new[u]] |= 1 << old_to_new[v]
    return tuple(result)

def delta_pairs(a,b):
    return tuple((u,v) for u,v in combinations(range(len(a)),2)
                 if ((a[u] >> v) ^ (b[u] >> v)) & 1)

def main():
    source = graph('202032')
    need(profile(source,2) == (59,17,17,59), 'literal balanced input')
    first, visits = cuts(source)
    targets = (179,181)
    hits = {h:[] for h in targets}
    # Deliberately no first-count pruning: independent complete 2^17 cube.
    digest = hashlib.sha256()
    for a in range(256):
        child = birth(source,a)
        second, path_visits = cuts(child)
        visits += path_visits
        digest.update(json.dumps(second,separators=(',',':')).encode())
        for b,h in enumerate(second):
            if h in hits:
                hits[h].append((a,b))
    need(hits[179] == [(128,128),(128,384),(128,448),(192,128),(192,192),(192,448)],
         'all six two-birth positives')
    need(not hits[181] and 181 not in first, 'no 181 extension')
    for a,b in hits[179]:
        need(hamilton(birth(birth(source,a),b)) == 179, 'literal positive replay')

    ready = ready59()
    need(profile(ready,0) == (59,5,5,61), 'inherited exceptional ready source')
    padded = tuple(r | (1 << 7) for r in ready) + (0,)
    need(profile(padded,0) == (59,5,5,61), 'sink padding preserves rooted profile')
    # Retain the original physical root2 and child3, including their edge.
    # Optimize only this declared padded seven-vertex template and port guard.
    choices=[]
    for child in range(3,8):
        domain=[v for v in range(8) if v not in (0,child)]
        image=[v for v in range(8) if v not in (2,3)]
        for perm in permutations(image):
            mapping=[None]*8
            mapping[0], mapping[child] = 2,3
            for v,w in zip(domain,perm):
                mapping[v]=w
            candidate=relabel(padded,mapping)
            differences=delta_pairs(source,candidate)
            choices.append((len(differences),candidate,tuple(mapping),differences))
    distance,target,mapping,differences=min(choices)
    need(profile(target,2) == (59,5,5,61), 'named repair endpoint profile')
    need(source[2] >> 3 & 1 and target[2] >> 3 & 1, 'protected port edge')

    compiler_path=ROOT/'05-knowledge/results/adaptive_metagraph_reformulation_20261103/metagraph_probe.py'
    compiler=load('inherited_native_compiler',compiler_path)
    current=source
    program=[]
    for u,v in differences:
        U,V,W,p,q,case,backward=compiler.compile_arc(current,u,v)
        need(compiler.restrict(U,8) == current, 'temporary birth retains source')
        need(compiler.native(U,p) == V and compiler.native(V,q) == W,
             'two actual legal native moves')
        result=compiler.restrict(W,8)
        need(delta_pairs(current,result) == ((u,v),), 'one named arc changes')
        program.append({'reversed_pair':[u,v], 'birth_rows':U,
                        'first_native_rows':V,'second_native_rows':W,
                        'hamilton_words':[p,q],'delete_name':8,
                        'after_deletion_rows':result,
                        'hamilton_counts':[hamilton(current),hamilton(U),hamilton(V),
                                           hamilton(W),hamilton(result)],
                        'case':case,'backward_gadget':backward})
        current=result
    need(current == target, 'exact endpoint after native excursion')
    first_reset=native_birth(target,2)
    final=native_birth(first_reset,8)
    need(profile(final,8) == (181,59,59,181), 'actual balanced 181 output')
    query_path=ROOT/'05-knowledge/results/adaptive_state_induction_20261103/two_port_envelope.json'
    suffixes=json.loads(query_path.read_text())['rank']['suffixes']
    need(len(suffixes)==17,'seventeen inherited fixed-port birth queries')
    for T in (source,target,final):
        need(count_dp(T)==hamilton(T),'subset DP versus literal path count')
    # The two original ports remain physically named2 and3 at all steps.
    out={'status':'PASS','source_word':'202032','source_rows':source,
         'source_physical_names':['c_2','x','reset_root','reset_child','y','c_1','birth_a','birth_b'],
         'prior_balanced_graph_to_current_labels':[0,1,4,5,6,7,2,3],
         'source_named_ports':{'root':2,'child':3},
         'independent_extension_audit':{
             'mask_convention':'out-mask: the newborn beats exactly these old vertices; producer in-masks are bitwise complements at each source order.',
             'path_enumeration':'Every directed simple path, all endpoint tables; no subset DP and no producer imports.',
             'all_first_masks':256,'all_second_mask_pairs':131072,
             'first_counts_at_most_181':sum(h<=181 for h in first),
             'first_hits':{h:[i for i,v in enumerate(first) if v==h] for h in targets},
             'two_birth_hits':hits,'second_cube_sha256':digest.hexdigest(),
             'simple_path_visits':visits},
         'source_exchange':{'template':'READY exceptional fifth tail, k=2, plus a global sink',
             'template_rows':padded,'template_to_physical':mapping,
             'target_rows':target,'new_root_profile':profile(target,2),
             'min_flips_within_declared_template':distance,
             'candidate_ported_template_maps':len(choices),'changed_pairs':differences,
             'program':program,'temporary_births':distance,'native_moves':2*distance,
             'temporary_deletions':distance,'peak_order_during_exchange':9,
             'preserved':'All eight old physical names; initial/final port edge2->3; endpoint H=59.',
             'not_asserted':'Stepwise H preservation, strong intermediates, or a globally shortest native route.'},
         'branch181':{'first_rows':first_reset,'final_rows':final,
                      'new_named_root':8,'new_named_child':9,
                      'balanced_profile':profile(final,8)},
         'joint_attachment_data':{
             'query_words':suffixes,
             'scope':'Static two-live-port births: bit i means old port i beats the newborn; all other old vertices lose to the newborn. These are not two-apex word suffixes.',
             'source':joint_queries(source,(2,3),suffixes),
             'after_source_exchange':joint_queries(target,(2,3),suffixes),
             'balanced_child_old_ports':joint_queries(final,(2,3),suffixes),
             'balanced_child_new_ports':joint_queries(final,(8,9),suffixes)},
         'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in (compiler_path,query_path)},
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
