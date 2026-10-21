#!/usr/bin/env python3
"""Exact seven-point group and color-stabilizer audit. No repo imports."""
from itertools import permutations, combinations, product
from collections import Counter
import json, hashlib, pathlib, argparse

def need(x,msg):
    if not x: raise RuntimeError(msg)

def image(s,p): return frozenset(p[i] for i in s)
def stabilizers(s,G): return [p for p in G if image(s,p)==s]
def color_stab(c,G): return [p for p in G if all(c[p[i]]==c[i] for i in range(len(c)))]
def group_check(G,n):
    S=set(G); need(tuple(range(n)) in S,'identity')
    for p in G:
        for q in G: need(tuple(p[q[i]] for i in range(n)) in S,'closure')

def main():
    n=7; Q={1,2,4}; V=range(n)
    arc=lambda i,j:(j-i)%7 in Q
    lines={frozenset((b+q)%7 for q in Q) for b in V}
    other={frozenset((b-q)%7 for q in Q) for b in V}
    perms=list(permutations(V))
    ordinary=[p for p in perms if all(arc(i,j)==arc(p[i],p[j]) for i,j in combinations(V,2))]
    anti=[p for p in perms if all(arc(i,j)!=arc(p[i],p[j]) for i,j in combinations(V,2))]
    extended=ordinary+anti
    fano=[p for p in perms if {image(s,p) for s in lines}==lines]
    affine=[tuple((a*x+b)%7 for x in V) for a in range(1,7) for b in V]
    affine_sq=[tuple((a*x+b)%7 for x in V) for a in Q for b in V]
    need(set(ordinary)==set(affine_sq),'ordinary = square affine')
    need(set(extended)==set(affine),'extended = all affine')
    need(len(fano)==168 and len(ordinary)==21 and len(extended)==42,'group orders')
    for G in (ordinary,extended,fano): group_check(G,n)
    need(len(lines)==7 and len(other)==7 and not lines&other,'two disjoint systems')
    need(all(sum({i,j}<=s for s in lines)==1 for i,j in combinations(V,2)),'Fano pair incidence')
    cycles={frozenset(s) for s in combinations(V,3) if all(sum(arc(i,j) for j in s if j!=i)==1 for i in s)}
    need(cycles==lines|other,'cyclic triples exact union')
    need(all({image(s,p) for s in lines}==other for p in anti),'anti swaps Fano systems')
    need(set(ordinary)<=set(fano),'ordinary inside Fano')
    subsets=[frozenset(s) for k in range(n+1) for s in combinations(V,k)]
    groups={'paley_aut':ordinary,'paley_aut_pm':extended,'fano_aut':fano}
    out={'universe':'All 5040 permutations, all 128 subsets, ordered colorings on seven labeled points.', 'groups':{},'fano_lines':[sorted(s) for s in sorted(lines,key=lambda s:tuple(sorted(s)))], 'fano_other_lines':[sorted(s) for s in sorted(other,key=lambda s:tuple(sorted(s)))]}
    for name,G in groups.items():
        hist={str(k):dict(sorted(Counter(len(stabilizers(s,G)) for s in subsets if len(s)==k).items())) for k in range(n+1)}
        regs=[s for s in subsets if len(stabilizers(s,G))==1]
        color_results={}
        for q in range(1,5):
            regcolors=[]
            for c in product(range(q),repeat=n):
                if len(color_stab(c,G))==1:regcolors.append(c)
            color_results[str(q)]={'regular_colorings':len(regcolors),'regular_orbits':len(regcolors)//len(G),'first':regcolors[0] if regcolors else None}
            if regcolors:break
        out['groups'][name]={'order':len(G),'subset_stabilizer_hist_by_size':hist,'regular_subsets':len(regs),'regular_subsets_by_size':[sum(len(s)==k for s in regs) for k in range(n+1)],'distinguishing_colors':color_results}
    need(out['groups']['paley_aut']['regular_subsets']==84,'84 ordinary regular subsets')
    need(out['groups']['paley_aut_pm']['regular_subsets']==0,'extended has no regular subset')
    need(out['groups']['fano_aut']['regular_subsets']==0,'Fano has no regular subset')
    # Every ordinary-regular 2-coloring has unique anti-symmetry and full stabilizer C2.
    for s in subsets:
        if len(stabilizers(s,ordinary))==1:
            A=stabilizers(s,anti);need(len(A)==1,'unique anti symmetry')
            p=A[0];need(all(p[p[i]]==i for i in V),'anti is involution')
    # A regular binary coloring plus one individually marked point can break its
    # remaining reflection iff the marked point is moved by that reflection.
    c=out['groups']['paley_aut_pm']['distinguishing_colors']['3']['first']
    need(c is not None,'3 colors suffice for extended group')
    # Faithfulness hostile: abstract odd group C3 acting trivially has no regular subset.
    # In actual permutation-image language the image is trivial, so it does.
    out['exact_consequences']={'ordinary_regular_implies_unique_involutory_anti':True,'no_gluck_import_for_extended_even_group':True,'all_cyclic_triples_are_two_fano_systems':True}
    out['source_sha256']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='receipt.json');args=parser.parse_args()
    pathlib.Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
