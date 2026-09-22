"""Exact operation-graph, shifted-arithmetic, and residue controls.

The proofs and quantifiers are in the companion note. No global Collatz
claim, graph-minor theorem, or plane-tiling claim is certified by a census.
"""
from fractions import Fraction
from math import gcd, isqrt
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def plus(a, b):
    return a + b + 1


def times(a, b):
    return a*b + a + b


def additive_arc(a, z):
    return 0 < a < z and z != 2*a


def multiplicative_arc(a, z):
    return 0 < a < z and z % a == 0 and z != a*a


def main():
    for a in range(-1, 13):
        require(plus(a, -1) == a and times(a, 0) == a, "transported identities")
        require(times(a, -1) == -1, "transported absorbing zero")
        for b in range(-1, 13):
            require(plus(a,b)+1 == (a+1)+(b+1), "addition conjugacy")
            require(times(a,b)+1 == (a+1)*(b+1), "multiplication conjugacy")
            for c in range(-1, 13):
                require(plus(plus(a,b),c) == plus(a,plus(b,c)), "addition associativity")
                require(times(times(a,b),c) == times(a,times(b,c)), "multiplication associativity")
                require(times(a,plus(b,c)) == plus(times(a,b),times(a,c)), "distributivity")

    counts=[]
    for n in range(1, 161):
        plus_edges={(a,z) for z in range(1,n+1) for a in range(1,z) if additive_arc(a,z)}
        times_edges={(a,z) for z in range(1,n+1) for a in range(1,z) if multiplicative_arc(a,z)}
        # Independent construction from distinct parent pairs; loops suppressed.
        plus_pairs=set()
        times_pairs=set()
        for a in range(1,n+1):
            for b in range(a+1,n+1):
                if a+b<=n:
                    plus_pairs.update(((a,a+b),(b,a+b)))
                if a*b<=n:
                    if a<a*b:
                        times_pairs.add((a,a*b))
                    if b<a*b:
                        times_pairs.add((b,a*b))
        require(plus_edges==plus_pairs and times_edges==times_pairs, "shadow/parent-pair equality")
        require(len(plus_edges)==n*(n-1)//2-n//2, "additive edge count")
        require(len(times_edges)==sum(n//a for a in range(1,n+1))-n-isqrt(n)+1,
                "multiplicative edge count")
        if n in (1,2,5,6,10,18,40,80,160):
            counts.append({"N":n,"additive_edges":len(plus_edges),"multiplicative_edges":len(times_edges)})

    bipartite=[]
    for name, predicate, left, right in (
        ("additive",additive_arc,(1,2,4),(3,5,6)),
        ("multiplicative",multiplicative_arc,(1,2,3),(6,12,18))):
        for a in left:
            for b in right:
                require(predicate(min(a,b),max(a,b)), "explicit K3,3 edge")
        bipartite.append({"graph":name,"left":left,"right":right})
    for size in range(2,9):
        additive=list(range(size+1,2*size+1))
        exponents=[3**i for i in range(size)]
        multiplicative=[2**e for e in exponents]
        for i in range(size):
            for j in range(i+1,size):
                require(additive_arc(additive[i],additive[j]), "arbitrary clique additive control")
                require(multiplicative_arc(multiplicative[i],multiplicative[j]),
                        "arbitrary clique multiplicative control")

    for multiplier in range(2,31):
        for modulus in range(1,61):
            roots=[q for q in range(modulus) if (multiplier-1)*q % modulus==0]
            count=gcd(multiplier-1,modulus)
            expected=[j*(modulus//count) for j in range(count)]
            require(roots==expected,"scaling fixed-point kernel")
        require([q for q in range(multiplier) if multiplier*q%multiplier==q]==[0],
                "literal modulo-m hostile")
        require(all(multiplier*q%(multiplier-1)==q for q in range(multiplier-1)),
                "repaired modulo-m-minus-one control")

    for a in range(1,21):
        for b in range(1,21):
            modulus=a*b
            for x in range(200):
                quotient,remainder=divmod(x,modulus)
                require(x==modulus*quotient+remainder and 0<=remainder<modulus,
                        "Euclidean address")

    for width in range(1,7):
        modulus=10**width-1
        for x in range(min(modulus,1000)):
            digits=f"{x:0{width}d}"
            rotated=int(digits[1:]+digits[0])
            require((10*x)%modulus==rotated%modulus,"cyclic digit rotation")
    seq=[0]
    for j in range(1,17):
        seq.append(plus(seq[-1],seq[-1]))
        require(seq[-1]+1==2**j,"shifted diagonal/Mersenne identity")
    powers=[1]
    for r in range(1,7):
        powers.append(times(powers[-1],powers[-1]))
        require(powers[-1]+1==2**(2**r),"shifted square/Fermat identity")
    inverse_fibre_checks=0
    for start in range(1,200,2):
        siblings=[start]
        for j in range(1,25):
            siblings.append(4*siblings[-1]+1)
            require(siblings[-1]==4**j*start+(4**j-1)//3,"inverse fibre closed form")
        initial=3*start+1
        image=initial//(initial & -initial)
        for j,node in enumerate(siblings):
            affine=3*node+1
            require(affine==4**j*initial and affine//(affine & -affine)==image,
                    "common accelerated Collatz image")
        for i in range(13):
            for gap in range(1,13):
                require(gcd(siblings[i],siblings[i+gap])==gcd(siblings[i],(4**gap-1)//3),
                        "inverse sibling gcd return law")
                inverse_fibre_checks+=1
    report={
        "status":"PASS", "scope":"Exact finite controls for proved operation identities",
        "universes":{"semiring_variables":[-1,12],"graph_cutoffs":[1,160],
                     "clique_sizes":[2,8],"scaling_multipliers":[2,30],
                     "scaling_moduli":[1,60],"quotient_A_B":[1,20],"quotient_x":[0,199],
                     "inverse_fibre_odd_starts":[1,199],"inverse_fibre_indices":[0,24],
                     "inverse_fibre_gcd_source_indices":[0,12],"inverse_fibre_gaps":[1,12]},
        "graph_edge_counts":counts,"explicit_nonplanarity_witnesses":bipartite,
        "Mersenne_diagonal_from_zero":seq,"shifted_squaring_from_one":powers,
        "modular_fixed_point_count":"gcd(multiplier-1,modulus)",
        "inverse_fibre_gcd_checks":inverse_fibre_checks,
        "inverse_fibre_scope":"Siblings with a common odd Collatz image; not consecutive forward iterates.",
        "tile_interpretation":"Parent-pair grid is complete; multiplicative arc coordinates have density zero.",
        "global_collatz_convergence_proved":False,
    }
    output=Path(__file__).with_suffix('.json')
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print("PASS: transported operations, exact shadows, universal clique controls, modular kernels and address laws")


if __name__=='__main__':
    main()
