"""Independent exact palette, subtournament, recursion and hostile controls."""
from itertools import combinations, product, permutations
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import argparse
import importlib.util
import json
ROOT = Path(__file__).resolve().parents[3]
CHECKS=0
def need(ok,name):
    global CHECKS
    CHECKS+=1
    if not ok:
        raise RuntimeError(name)

def triples(rows):
    return [sum(1 << i for i in t) for t in combinations(range(len(rows)),3)
            if all(sum(rows[i] >> j & 1 for j in t)==1 for i in t)]
def legal_sets(rows):
    cy=triples(rows)
    return [s for s in range(1 << len(rows)) if not any(s & c == c for c in cy)]

def color_solver(rows):
    good=legal_sets(rows); n=len(rows)
    @lru_cache(None)
    def chi(s):
        if s == 0:
            return 0
        if s in good_set:
            return 1
        pivot=s & -s
        return 1+min(chi(s ^ t) for t in good if t & pivot and t & s == t)
    good_set=set(good)
    return chi((1 << n)-1),max(s.bit_count() for s in good)

P=[sum(1 << j for j in range(7) if (j-i)%7 in (1,2,4)) for i in range(7)]
H=[74,100,25,50,35,5,60]

def weighted():
    cy=triples(P); good=legal_sets(P)[1:]
    @lru_cache(None)
    def dp(d):
        if not any(d):
            return 0
        pivot=next(i for i,x in enumerate(d) if x)
        active=sum(1 << i for i,x in enumerate(d) if x)
        return 1+min(dp(tuple(x-((s >> i)&1) for i,x in enumerate(d))) for s in good
                     if s & active == s and (s >> pivot)&1)
    count=0
    for d in product(range(4),repeat=7):
        expected=max(max(d),(sum(d)+2)//3,
                     max((sum(d[i] for i in range(7) if c >> i & 1)+1)//2 for c in cy))
        need(dp(d)==expected,'all unequal demands in [0,3]^7')
        count+=1
    return {'vectors':count,'range':'every d_i in {0,1,2,3}'}

def hole():
    good=legal_sets(H)
    chi,alpha=color_solver(H)
    need((chi,alpha)==(3,4),'native seven-vertex coloring hole')
    dual=[F(1,2)]+[F(1,4)]*6
    need(sum(dual)==2 and all(sum(dual[i] for i in range(7) if s >> i & 1)<=1 for s in good),
         'exact fractional dual')
    fractional=[37,58,67,92]
    need(all(s in good for s in fractional) and
         all(sum(s >> i & 1 for s in fractional)==2 for i in range(7)), 'fractional cover of cost2')
    two=[11,21,102,120]
    three=[11,19,37,92,102,120]
    one=[1,6,120]
    for k in range(1,65):
        colors=one if k==1 else (two*(k//2) if k%2==0 else three+two*((k-3)//2))
        need(all(s in good for s in colors) and all(sum(s >> i & 1 for s in colors)==k for i in range(7)),
             'all tested demand-ray constructions')
        need(len(colors)==(3 if k==1 else 2*k),'sharp ray cost')
    rank_bound=max(F(s.bit_count(),max(t.bit_count() for t in good if t & s == t))
                   for s in range(1,128))
    need(rank_bound==F(7,4),'sharp unweighted subset-capacity bound is strictly below the fractional optimum')
    # A cover of each degree2 and3 supplies a local feasible point avoiding each chosen prime,
    # but the desired degree-one integral point does not exist.
    return {'rows':H,'dichromatic':chi,'transitive':alpha,'fractional_dichromatic':'2',
            'dual':list(map(str,dual)),'half_weight_supports':fractional,
            'one_palette':one,'two_palette':two,'three_palette':three,
            'ray_formula':'chi_k=0 at k=0; chi_1=3; chi_k=2k for every k>=2',
            'maximum_unweighted_subset_rank_bound':str(rank_bound)}

def substitute(base,children):
    offsets=[0]
    for rows in children:
        offsets.append(offsets[-1]+len(rows))
    out=[0]*offsets[-1]
    for i,child in enumerate(children):
        for v,row in enumerate(child):
            out[offsets[i]+v] |= row << offsets[i]
        for j in range(len(children)):
            if base[i] >> j & 1:
                block=((1 << len(children[j]))-1) << offsets[j]
                for v in range(len(child)):
                    out[offsets[i]+v] |= block
    return out

def recursion_controls():
    cases=0
    # Full choices of which Paley ports get C3, up to three such ports.
    local=[[0],[2,4,1]]
    for mask in range(128):
        if mask.bit_count()>3:
            continue
        children=[local[mask >> i & 1] for i in range(7)]
        rows=substitute(P,children)
        d=[1+(mask >> i & 1) for i in range(7)]
        a=d[:]  # K1 has alpha1 and C3 alpha2.
        predicted=max(max(d),(sum(d)+2)//3,
                      max((sum(d[i] for i in range(7) if c >> i & 1)+1)//2 for c in triples(P)))
        predicted_alpha=max(sum(a[i] for i in range(7) if s >> i & 1) for s in legal_sets(P))
        need(color_solver(rows)==(predicted,predicted_alpha),'literal unequal Paley substitution')
        cases+=1
    h1=color_solver(H)[0]
    need(h1==3,'H-power seed')
    return {'literal_paley_substitutions':cases,'largest_order':13,
            'H_power_dichromatic':'3*2^(r-1), r>=1','H_power_transitive':'4^r',
            'H_power_fractional':'2^r'}

def rank_cliques():
    counts=[]
    for rows in (P,H):
        for k in range(1,5):
            # An order-layer clique has one vertex in each of k layers and all forward arcs.
            lifted=sum(all(rows[word[i]] >> word[j] & 1 for i in range(k) for j in range(i+1,k))
                       for word in permutations(range(7),k))
            direct=sum(s.bit_count()==k for s in legal_sets(rows))
            need(lifted==direct,'rank-layer cliques biject with transitive k-subsets')
            counts.append({'seed':'P7' if rows==P else 'H7','k':k,'count':direct})
    return counts

def hostile():
    spec=importlib.util.spec_from_file_location('regular_control',
        ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py')
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    inherited=mod.hostile_checks()
    all_sets=[]; port_answers=[]
    for flip in (False,True):
        rows=[0]*12
        for u,v in combinations(range(12),2):
            b,s=divmod(u,4); c,t=divmod(v,4)
            arc=(c-b)%3==1 if b!=c else s==0 or (s>0 and (t-s)%3==1)
            if flip and b==c==2: arc=not arc
            rows[u if arc else v] |= 1 << (v if arc else u)
        all_sets.append(legal_sets(rows))
        need(color_solver(rows)==(3,6),'81/27 ordinary coloring and transitivity')
        # Identical added port: 8 -> p, and p -> every other old vertex.
        extended=[row | ((1 << 12) if i==8 else 0) for i,row in enumerate(rows)]
        extended.append(((1 << 12)-1) ^ (1 << 8))
        selected=(1 << 8)|(1 << 9)|(1 << 12)
        answer=1 if selected in legal_sets(extended) else 2
        port_answers.append(answer)
    need(all_sets[0]==all_sets[1] and len(all_sets[0])==547,'entire named palette incidence agrees')
    need(port_answers==[1,2],'one actual port separates the previously identical palette carriers')
    return {'native':inherited,'same_named_transitive_subsets':547,
            'one_port_new_out_mask':3839,'one_port_query_vertices':[8,9,12],
            'one_port_acyclic_color_answers':port_answers,
            'consequence':'All named weighted acyclic-coloring problems coincide; symmetry still differs.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path); args=ap.parse_args()
    result={'weighted':weighted(),'hole':hole(),'substitution':recursion_controls(),
            'rank_layer_cliques':rank_cliques(),'hostile_81_27':hostile(),'checks':CHECKS}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: args.output.write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')
