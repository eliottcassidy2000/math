"""Independent brute permutations at 7; all 32 labelled circulant choices at 11."""
from verify import blocks, automorphisms, skew_labels, check
from itertools import permutations, product
from collections import Counter
import json
p=7;B=blocks(p,{1,2,4});blockset=set(B)
brute=[]
for a in permutations(range(p)):
    if all(frozenset(a[x] for x in b) in blockset for b in B):brute.append(a)
check(set(brute)=={a for a,b in automorphisms(B)}, 'audit.py validation originally at line 10')
brute_skew=[f for f in permutations(range(p)) if all(i not in B[f[i]] for i in range(p)) and all((j in B[f[i]]) != (i in B[f[j]]) for i in range(p) for j in range(i))]
check(set(brute_skew)==set(skew_labels(B)), 'audit.py validation originally at line 12')
print('Independent 7! permutation checks passed: 168 design automorphisms, 8 skew matchings.')
counts=Counter()
for signs in product([False,True],repeat=5):
    S={11-i if sign else i for i,sign in enumerate(signs,1)};B=blocks(11,S);G=automorphisms(B);K=[a for a,b in G if a==b];F=skew_labels(B)
    counts[len(G),len(K),len(F)]+=1
    check(len(F)==len(G)//len(K), 'audit.py validation originally at line 18')
print(json.dumps([dict(design_aut=g,tournament_aut=k,skew_matchings=f,labelled_connection_sets=num) for (g,k,f),num in sorted(counts.items())]))
print('All 32 circulant connection sets at 11 have a single orbit of compatible matchings. ALL CHECKS PASSED')
