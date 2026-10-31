#!/usr/bin/env python3
"""Read-only check of the inherited, NOT single-arc-flip, n=5 metagraph.

Input is git-show copy from:
05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/census_n5.json
Loops and literal black multiplicities are retained. No layout coordinates.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json


p = Path(__file__).with_name('inherited_census_n5.json')
data = json.loads(p.read_text())
n = len(data['nodes'])

# Independent reconstruction from the fixed path 1->2->...->5. A set bit
# reverses one of the six nonconsecutive arcs, in the stated explorer gauge.
vertex_pairs = list(combinations(range(1,6),2))
pair_index = {edge:i for i,edge in enumerate(vertex_pairs)}
relabel_actions = []
for perm0 in permutations(range(1,6)):
    perm = dict(zip(range(1,6), perm0))
    relabel_actions.append([(pair_index[tuple(sorted((perm[u],perm[v])))],
                              perm[u] > perm[v]) for u,v in vertex_pairs])

def relabel(word, action):
    return sum(((((word >> i) & 1) ^ rev) << j)
               for i,(j,rev) in enumerate(action))

def canonical(word):
    return min(relabel(word,a) for a in relabel_actions)

full_word = (1 << len(vertex_pairs))-1
canon_by_mask = []
for mask in range(64):
    word = full_word
    for k,(a,b) in enumerate(data['tiles']):
        if mask >> k & 1:
            word ^= 1 << pair_index[tuple(sorted((a,b)))]
    canon_by_mask.append(canonical(word))
classes = {}
class_by_mask = []
for c in canon_by_mask:
    classes.setdefault(c,len(classes))
    class_by_mask.append(classes[c])
if not (class_by_mask == data['class_by_mask']):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 49')
merged = {}
class_to_merged = {}
for c,idx in classes.items():
    key = min(c,canonical(c ^ full_word))
    merged.setdefault(key,len(merged))
    class_to_merged[idx] = merged[key]
pi = [class_to_merged[c] for c in class_by_mask]
if not (pi == data['pi']):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 57')
blue_reconstructed, black_reconstructed = Counter(),Counter()
for mask in range(32):
    reflected = sum(((mask >> k) & 1) << j
                    for k,j in enumerate(data['reflection_permutation']))
    counter = blue_reconstructed if reflected == mask else black_reconstructed
    counter[tuple(sorted((pi[mask],pi[63-mask])))] += 1
if not (blue_reconstructed == Counter({(u,v):m for u,v,m in data['blue_edges']})):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 64')
if not (black_reconstructed == Counter({(u,v):m for u,v,m in data['black_edges']})):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 65')

blue = [[0]*n for _ in range(n)]
black = [[0]*n for _ in range(n)]
for key, mat in [('blue_edges', blue), ('black_edges', black)]:
    for u, v, m in data[key]:
        mat[u][v] += m
        if u != v:
            mat[v][u] += m

def degree(mat, v):
    return sum(mat[v])+mat[v][v]

signatures = [(degree(blue,v), degree(black,v), blue[v][v], black[v][v])
              for v in range(n)]
blocks = defaultdict(list)
for v, sig in enumerate(signatures):
    blocks[sig].append(v)
blocks = list(blocks.values())
autos = []
tested = 0
for block_images in product(*(list(permutations(block)) for block in blocks)):
    perm = list(range(n))
    for block, images in zip(blocks, block_images):
        for v, w in zip(block, images):
            perm[v] = w
    tested += 1
    if all(mat[u][v] == mat[perm[u]][perm[v]]
           for mat in (blue, black) for u in range(n) for v in range(n)):
        autos.append(perm)

if not (n == 10):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 96')
if not (sum(m for _,_,m in data['blue_edges']) == 8):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 97')
if not (sum(m for _,_,m in data['black_edges']) == 24):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 98')
if not (sum(black[v][v] for v in range(n)) == 4):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 99')
if not ([v for v in range(n) if degree(black,v) == 0] == [0,7,8]):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 100')
if not (autos == [list(range(n))]):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 101')
pb_fibers = [[mask for mask,node in enumerate(pi) if node == v] for v in [0,7,8]]
def distance_profile(triple):
    return sorted((a^b).bit_count() for a,b in combinations(triple,2))
pb_profiles = sorted({tuple(distance_profile(t)) for t in product(*pb_fibers)})
user_profile = distance_profile((7,17,55))
if not (pb_fibers == [[0],[11,46,53],[16]]):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 107')
if not (user_profile == [2,3,3]):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 108')
if not (tuple(user_profile) not in pb_profiles):
    raise RuntimeError('check_inherited_metagraph.py: gate originally at line 109')
out = {'inherited_input_sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
       'independently_reconstructed_from_64_path_normal_forms': True,
       'nodes': n, 'literal_blue_lines': 8, 'literal_black_lines':24,
       'literal_black_loops': 4,
       'pure_blue_nodes': [0,7,8],
       'pure_blue_mask_fibers': pb_fibers,
       'pure_blue_representative_distance_profiles': pb_profiles,
       'user_7_17_55_Hamming_distance_profile': user_profile,
       'signed_coordinate_gauge_can_identify_user_triple_with_PB': False,
       'blue_degrees': [degree(blue,v) for v in range(n)],
       'black_degrees': [degree(black,v) for v in range(n)],
       'color_and_multiplicity_preserving_vertex_automorphisms': autos,
       'degree_signature_compatible_permutations_tested': tested,
       'not_Johnson_or_Petersen': 'Literal loops already obstruct either simple-graph identification.',
       'not_a_faithful_A5_action': 'Colored weighted graph has trivial vertex automorphism group.'}
print(json.dumps(out, indent=2))
