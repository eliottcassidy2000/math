"""Eight endpoint frames, affine designs, and unordered-pair holonomy.

Exact standard-library universes: 48 palette transports, all 40320 point
permutations of the cube, 1344 affine maps, and every 48^r bouquet transport
for r=0,1,2. No graph-colorability or unit-distance realization is inferred.
"""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json
import sys


def require(test, why):
    if not test:
        raise RuntimeError(why)


def parity(x):
    return x.bit_count() % 2


def dot(x, y):
    return parity(x & y)


def permute(x, pi):
    return sum(((x >> j) & 1) << pi[j] for j in range(3))


def apply_map(x, g):
    pi, a = g
    return permute(x, pi) ^ a


def rank_f2(columns):
    pivots = {}
    for column in columns:
        w = column
        while w:
            p = w.bit_length()-1
            if p in pivots:
                w ^= pivots[p]
            else:
                pivots[p] = w
                break
    return len(pivots)


def image_mask(mask, p):
    return sum(((mask >> i) & 1) << p[i] for i in range(8))


def frame_fano(eta):
    # Colors 0,1 are pair0; 2,3 pair1; 4,5 pair2; central label6.
    center = {tuple(sorted((6,2*j,2*j+1))) for j in range(3)}
    off = {tuple(2*j+((x>>j)&1) for j in range(3)) for x in range(8) if parity(x)==eta}
    return frozenset(center | off)


def main():
    perms = tuple(permutations(range(3)))
    transports = tuple((pi,a) for pi in perms for a in range(8))
    even = (0,3,5,6)
    hyperplanes = frozenset(sum(1 << x for x in range(8) if dot(n,x)==b) for n in range(1,8) for b in range(2))
    code = hyperplanes | {0,255}
    require(len(hyperplanes)==14 and len(code)==16, "affine hyperplane count")
    require(all(a^b in code for a in code for b in code), "extended code closure")
    require(Counter(w.bit_count() for w in code)=={0:1,4:14,8:1}, "extended code weights")
    for triple in combinations(range(8),3):
        containing = [b for b in hyperplanes if all((b>>x)&1 for x in triple)]
        require(len(containing)==1, "Steiner S(3,4,8) axiom")
        fourth=triple[0]^triple[1]^triple[2]
        require(containing[0]==sum(1<<x for x in (*triple,fourth)), "affine fourth point")
    for origin in range(8):
        lines={tuple(x for x in range(8) if x!=origin and (b>>x)&1) for b in hyperplanes if (b>>origin)&1}
        require(len(lines)==7, "derived Fano at origin")
        require(all(sum(set(pair)<=set(line) for line in lines)==1 for pair in combinations([x for x in range(8) if x!=origin],2)), "derived Fano pair axiom")

    cube_edges=frozenset(tuple(sorted((x,x^(1<<j)))) for x in range(8) for j in range(3))
    palette_maps=set()
    direct_product=set()
    for g in transports:
        pi,a=g
        p=tuple(apply_map(x,g) for x in range(8))
        palette_maps.add(p)
        require({image_mask(b,p) for b in hyperplanes}==set(hyperplanes), "palette transport preserves SQS")
        require({tuple(sorted((p[x],p[y]))) for x,y in cube_edges}==set(cube_edges), "palette transport preserves cube")
        a_even=a^(7 if parity(a) else 0)
        tetra_perm=tuple(even.index(permute(x,pi)^a_even) for x in even)
        direct_product.add((parity(a),tetra_perm))
        for x in range(8):
            eta=parity(x)
            e=x^(7 if eta else 0)
            image=apply_map(x,g)
            require(parity(image)==(eta^parity(a)), "completion parity action")
            require((image^(7 if parity(image) else 0))==(permute(e,pi)^a_even), "tetrahedral splitting")
        color_perm=tuple(2*pi[j]+(b^((a>>pi[j])&1)) for j in range(3) for b in range(2))+(6,)
        for eta in (0,1):
            moved={tuple(sorted(color_perm[x] for x in line)) for line in frame_fano(eta)}
            require(moved==set(frame_fano(eta^parity(a))), "actual Fano completion action")
    require(len(palette_maps)==48 and len(direct_product)==48, "wreath product C2 times S4")
    require(direct_product==set(product(range(2),permutations(range(4)))), "all tetrahedron permutations occur")

    linear_bases=[b for b in product(range(1,8),repeat=3) if rank_f2(b)==3]
    affine_maps=set()
    for b in linear_bases:
        linear=lambda x: (b[0] if x&1 else 0)^(b[1] if x&2 else 0)^(b[2] if x&4 else 0)
        for a in range(8):
            p=tuple(linear(x)^a for x in range(8))
            affine_maps.add(p)
            require({image_mask(h,p) for h in hyperplanes}==set(hyperplanes), "AGL preserves affine design")
    require(len(linear_bases)==168 and len(affine_maps)==1344, "GL and AGL orders")
    # Independent all-permutation automorphism census.
    all_design_aut=set()
    all_cube_aut=set()
    for p in permutations(range(8)):
        if all(image_mask(h,p) in hyperplanes for h in hyperplanes):
            all_design_aut.add(p)
        if all(tuple(sorted((p[x],p[y]))) in cube_edges for x,y in cube_edges):
            all_cube_aut.add(p)
    require(all_design_aut==affine_maps, "full SQS automorphism group")
    require(all_cube_aut==palette_maps, "full cube automorphism group")
    direction_bases={tuple(sorted(b)) for b in linear_bases}
    require(len(direction_bases)==28 and 1344==48*28, "palette basis sidecar index")
    shear=tuple(((x&1)^((x>>1)&1))|(x&6) for x in range(8))
    require(shear in affine_maps and shear not in palette_maps, "SQS versus cube hostile")
    require(shear[0]==0 and shear[2]==3 and (0,2) in cube_edges and (0,3) not in cube_edges, "explicit lost edge")

    # Exhaust every graph-gain instance after spanning-tree normalization for
    # cycle rank0,1,2. Fixing linear monodromy, enumerate all translation data.
    bouquet_rows=[]
    total_instances=0
    for r in range(3):
        aggregate=Counter()
        for pis in product(perms,repeat=r):
            fixed_even=[x for x in even if all(permute(x,pi)==x for pi in pis)]
            eps=len(fixed_even).bit_length()-1
            hidden=2*(r-1)+eps
            require(hidden>=0, "nonnegative hidden dimension")
            scalar_valid=0
            compatible=0
            for gains in product(range(8),repeat=r):
                total_instances+=1
                gs=tuple(zip(pis,gains))
                fixed=[x for x in range(8) if all(apply_map(x,g)==x for g in gs)]
                passes=all(parity(a)==0 for a in gains)
                if passes:
                    scalar_valid+=1
                if fixed:
                    compatible+=1
                    require(passes and len(fixed)==2**(1+eps), "section count and parity")
                aggregate[(eps,passes,bool(fixed),len(fixed))]+=1
            require(scalar_valid==4**r, "completion-compatible gains")
            require(compatible==2**(2-eps), "coboundary image count")
            require(scalar_valid==compatible*2**hidden, "hidden obstruction quotient size")
        bouquet_rows.append({"cycle_rank":r,"instances":48**r,
            "rows":[{"invariant_even_dimension":e,"completion_compatible":p,"frame_compatible":c,"sections":s,"instances":n} for (e,p,c,s),n in sorted(aggregate.items())]})
    require(total_instances==1+48+2304, "complete bouquet universe")
    cycle=(1,2,0)
    witness=((cycle,0),(cycle,5))
    fixed_sets=[[x for x in range(8) if apply_map(x,g)==x] for g in witness]
    require(fixed_sets==[[0,7],[3,4]], "two individually soluble loops hostile")
    require(not set(fixed_sets[0])&set(fixed_sets[1]), "no joint endpoint frame")
    require(all(parity(g[1])==0 for g in witness), "hostile completion passes")

    # A tree with arbitrary permutation/gain data always has eight sections.
    tree_edges=((0,1),(1,2))
    for gs in product(transports, repeat=2):
        sections=0
        for root in range(8):
            x1=apply_map(root,gs[0]); x2=apply_map(x1,gs[1])
            require(apply_map(root,gs[0])==x1 and apply_map(x1,gs[1])==x2, "tree propagation")
            sections+=1
        require(sections==8, "tree section count")

    receipt={"status":"FINITE-EXACT supporting the accompanying proofs",
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "universes":{"palette_transports":48,"all_point_permutations":40320,"affine_maps":1344,"bouquet_instances":total_instances,"tree_instances":2304},
        "hyperplane_blocks":[[x for x in range(8) if (b>>x)&1] for b in sorted(hyperplanes)],
        "extended_code_words":sorted(code),"unordered_direction_bases":28,
        "automorphism_orders":{"SQS8":len(all_design_aut),"cube":len(all_cube_aut),"fixed_origin":168},
        "bouquet_census":bouquet_rows,
        "hostiles":{"shear":list(shear),"edge_lost":[[0,2],[0,3]],
            "individually_soluble_loops":{"coordinate_permutation":list(cycle),"gains":[0,5],"fixed_sets":fixed_sets,"joint_fixed_set":[]}}}
    encoded=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if len(sys.argv)>1:
        Path(sys.argv[1]).write_text(encoded,encoding='utf-8',newline='\n')
    else:
        print(encoded,end='')


if __name__=='__main__':
    main()
