#!/usr/bin/env python3
"""Independent audit: arc-set isomorphism, cyclic integer shifts, pair separation.

Does not import verify_ring_ca.py. All gates survive -O. Reads only generated
JSON/TSV and the bundled inherited atlas beside this file. No sibling-path reads.
"""
from pathlib import Path
from itertools import permutations, combinations
from collections import Counter
import csv
import hashlib
import json

ROOT = Path(__file__).parent
checks = 0

def need(value, message):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(message)

P = list(permutations(range(1,6)))
CELLS = [(5,1),(4,1),(3,1),(5,2),(4,2),(5,3)]
ORDER = [0,1,2,4,5,3]
ALL_ARCS = frozenset((i,j) for i in range(1,6) for j in range(1,6) if i != j)

def arcs(mask):
    e = {(i,j) for i,j in combinations(range(1,6),2)}
    for k,(j,i) in enumerate(CELLS):
        if (mask//(2**k))%2:
            e.remove((i,j))
            e.add((j,i))
    return frozenset(e)


def iso(e,f):
    de = {i:sum((i,j) in e for j in range(1,6)) for i in range(1,6)}
    df = {i:sum((i,j) in f for j in range(1,6)) for i in range(1,6)}
    if sorted(de.values()) != sorted(df.values()):
        return False
    for p in P:
        if all(de[i] == df[p[i-1]] for i in range(1,6)):
            if all((p[i-1],p[j-1]) in f for i,j in e):
                return True
    return False

es = [arcs(x) for x in range(64)]
class_reps, node_reps, cid, pi = [], [], [], []
for e in es:
    matches = [i for i,r in enumerate(class_reps) if iso(e,r)]
    if not matches:
        matches = [len(class_reps)]
        class_reps.append(e)
    need(len(matches)==1,'unique directed isomorphism class')
    cid.append(matches[0])
    matches = [i for i,r in enumerate(node_reps) if iso(e,r) or iso(ALL_ARCS-e,r)]
    if not matches:
        matches = [len(node_reps)]
        node_reps.append(e)
    need(len(matches)==1,'unique converse-merged class')
    pi.append(matches[0])
need(len(class_reps)==12 and len(node_reps)==10,'independent census')

def to_ring(x):
    return sum(((x//(2**k))%2)*(2**i) for i,k in enumerate(ORDER))

def from_ring(y):
    return sum(((y//(2**i))%2)*(2**k) for i,k in enumerate(ORDER))

def shift(y):
    return ((y*2)%64)+(y//32)

def unshift(y):
    return y//2+32*(y%2)

def a(x):
    y=to_ring(x)
    return from_ring(shift(y)^unshift(y))

def b(x):
    return x^a(x)

def r(x):
    y=to_ring(x)
    z=sum(((y//(2**i))%2)*(2**((-i)%6)) for i in range(6))
    return from_ring(z)

c=lambda x:x^63
blues=[r(x)==x for x in range(64)]
results=json.loads((ROOT/'ring_ca_results.json').read_text())
quotients=json.loads((ROOT/'observer_quotients.json').read_text())
need(results['census']['pi']==pi,'all 64 independently reconstructed merged nodes')
need(results['census']['class_by_mask']==cid,'all 64 directed classes')
rows=list(csv.DictReader((ROOT/'transitions.tsv').open(),delimiter='\t'))
need(len(rows)==64,'table contains exactly 64 rows')
for x,row in enumerate(rows):
    expected={'mask':x,'tournament_class':cid[x],'merged_node':pi[x],
              'blue':int(blues[x]),'C':c(x),'R':r(x),
              'A':a(x),'A2':a(a(x)),'B':b(x),'B2':b(b(x)),
              'pi_A':pi[a(x)],'pi_B':pi[b(x)],'pi_C':pi[c(x)]}
    for k,v in expected.items():
        need(int(row[k])==v,'transition '+str(x)+'/'+k)
    need(row['explorer_bits_lsb_first']==''.join(str(x//(2**i)%2) for i in range(6)),
         'explorer bit rendering')
    need(row['perimeter_bits']==''.join(str(to_ring(x)//(2**i)%2) for i in range(6)),
         'perimeter bit rendering')
    need(a(a(a(x)))==a(x),'A cubic identity')
    need(b(b(b(x)))==b(b(x)),'B stabilization identity')
    need(a(c(x))==a(x) and b(c(x))==c(b(x)),'complement response')
    need(a(r(x))==r(a(x)) and b(r(x))==r(b(x)),'reflection response')

controls={'CR':[c,r],'A':[a],'B':[b],'AB':[a,b],'ACR':[a,c,r],
          'BCR':[b,c,r],'ABCR':[a,b,c,r],
          'coordinate_flips':[lambda x,k=k:x^(2**k) for k in range(6)]}
observations={'pi':pi,'pi_blue':list(zip(pi,blues))}


def indistinguishability(obs,fs):
    """Pair-removal recursion, independent of producer's partition refinement."""
    equivalent=[[obs[x]==obs[y] for y in range(64)] for x in range(64)]
    depths=[[None if equivalent[x][y] else 0 for y in range(64)] for x in range(64)]
    depth=0
    while True:
        removed=[]
        for x in range(64):
            for y in range(x):
                if equivalent[x][y] and any(not equivalent[f(x)][f(y)] for f in fs):
                    removed.append((x,y))
        if not removed:
            break
        depth+=1
        for x,y in removed:
            equivalent[x][y]=equivalent[y][x]=False
            depths[x][y]=depths[y][x]=depth
    blocks=[]
    seen=set()
    for x in range(64):
        if x not in seen:
            block=[y for y in range(64) if equivalent[x][y]]
            need(not(seen&set(block)),'equivalence classes partition states')
            blocks.append(block)
            seen.update(block)
    return equivalent,blocks,depth,depths

quotient_audits={}
for on,obs in observations.items():
    quotient_audits[on]={}
    for cn,fs in controls.items():
        eq,blocks,depth,depths=indistinguishability(obs,fs)
        q=quotients[on][cn]
        need(blocks==q['blocks'],'minimal blocks '+on+'/'+cn)
        need(len(blocks)==q['states'],'state count '+on+'/'+cn)
        for x in range(64):
            for y in range(64):
                need(eq[x][y] == (q['class_by_mask'][x]==q['class_by_mask'][y]),
                     'all-pair equality '+on+'/'+cn)
        quotient_audits[on][cn]={'states':len(blocks),
          'maximum_shortest_distinguishing_word_length':depth,
          'separated_pair_depth_histogram':dict(sorted(Counter(depths[x][y]
            for x in range(64) for y in range(x) if depths[x][y] is not None).items()))}

# Audit the complete weighted, colored complement-line graph.
edges={'blue':Counter(),'black':Counter()}
for x in range(64):
    if x<c(x):
        edges['blue' if blues[x] else 'black'][tuple(sorted([pi[x],pi[c(x)]]))]+=1
for color in edges:
    need([[u,v,n] for (u,v),n in sorted(edges[color].items())]==results['census']['edges'][color],
         'independent '+color+' complement edges')

# The literal requested masks are not pure-blue nodes; 17 is a blue tiling in
# a mixed node, which must not be mistaken for a pure-blue node.
need([pi[x] for x in [7,17,55]]==[1,9,6],'literal triplet nodes')
need([blues[x] for x in [7,17,55]]==[False,True,False],'literal triplet blue flags')
need([b(x) for x in [7,17,55]]==[26,63,22],'literal one-step triplet')
need([b(b(x)) for x in [7,17,55]]==[26,63,37],'literal attractor triplet')
need([pi[b(x)] for x in [11,46,53]]==[6,0,9],'B does not descend')
need([pi[a(x)] for x in [11,46,53]]==[7,7,5],'A does not descend')
jb=results['linear_decomposition']['basis_explorer_masks']
need(jb==[37,26,46,10,29,5],'explicit Jordan basis')
need([b(v) for v in jb]==[37,26,0,46,0,29],'Jordan basis B images')
vectors=[0]
for v in jb:
    vectors += [x^v for x in vectors[:]]
need(len(set(vectors))==64,'Jordan basis independence')

# Required portable provenance pin: the two research-lane copies were identical.
path=ROOT/'reference_n5_atlas.json'
data=json.loads(path.read_text())
reference_hash=hashlib.sha256(path.read_bytes()).hexdigest()
need(reference_hash=='c5e68df205d8b9ce5de786e46352db83aebe93d954c9a0d5057bf862afdd143c',
     'bundled inherited atlas exact source hash')
need(data['tiles']==[list(cell) for cell in CELLS],'inherited explorer coordinate gauge')
need(data['reflection_permutation']==[0,3,5,1,4,2],'inherited reflection coordinate gauge')
need(data['pi']==pi and data['class_by_mask']==cid,'inherited atlas classes')
for color in edges:
    need(data[color+'_edges']==results['census']['edges'][color], 'inherited colored graph')
provenance=[{'relative_path':'reference_n5_atlas.json','sha256':reference_hash}]

out={'status':'PASS','checks':checks,
     'independent_methods':['pairwise directed arc-set isomorphism with degree filtering',
       'cyclic integer rotation for CA','pair-removal continuation indistinguishability'],
     'normal_and_optimized_safe':True,
     'inherited_pin_comparisons':provenance,'quotient_audits':quotient_audits,
     'artifact_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
       for name in ['ring_ca_results.json','observer_quotients.json','transitions.tsv']}}
print(json.dumps(out,indent=2,sort_keys=True))
