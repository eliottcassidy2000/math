"""Exact direction-coset/inverse-ray clock conjugacy and prime-power tori.

Group labels and cellular faces are retained. This is not a conjugacy of
forward Collatz trajectories. All finite checks survive optimization.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd, isqrt
import json


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


def prime(p):
    return type(p) is int and p >= 2 and all(p % d for d in range(2,isqrt(p)+1))


def prime_factors(n):
    result, d = [], 2
    while d*d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return result


def order(a, modulus, candidate):
    need(gcd(a,modulus) == 1 and pow(a,candidate,modulus) == 1, 'unit return')
    divisors = prime_factors(candidate)
    result = candidate
    for p in divisors:
        while result % p == 0 and pow(a,result//p,modulus) == 1:
            result //= p
    need(all(pow(a,result//p,modulus) != 1 for p in divisors if result % p == 0),
         'all order shortenings excluded')
    return result


def v3(n):
    exponent = 0
    while n % 3 == 0:
        exponent += 1
        n //= 3
    return exponent


@lru_cache(None)
def layers(p, depth=1):
    need(prime(p) and p % 12 == 7 and type(depth) is int and depth >= 1,
         'prime congruent to seven modulo twelve, positive depth')
    modulus = p**depth
    units = {x for x in range(1,modulus) if x % p}
    squares = {x*x % modulus for x in units}
    h = tuple(x for x in sorted(units) if pow(x,3,modulus) == 1)
    need(len(h) == 3 and set(h) <= squares and (-1) % modulus not in squares,
         'retained cubic subgroup and asymmetric square half')
    classes, label, cube = {}, {}, {}
    for x in sorted(squares):
        if x in label:
            continue
        coset = tuple(sorted(x*t % modulus for t in h))
        representative = coset[0]
        classes[representative] = coset
        for y in coset:
            need(y not in label, 'direction cosets disjoint')
            label[y] = representative
        value = pow(x,3,modulus)
        need(value not in cube, 'cube gives an injective layer coordinate')
        cube[value] = representative
        need(all(pow(y,3,modulus) == value for y in coset), 'cube is independent of representative')
    need(len(classes)*6 == len(units) and set(label) == squares, 'complete layer quotient')
    need(set(cube) == {pow(x,3,modulus) for x in squares}, 'exact cubed-square image')
    return modulus, h, classes, label, cube


def layer_orbits(p, depth=1):
    modulus, h, classes, label, cube = layers(p,depth)
    unseen, orbits = set(classes), []
    for start in sorted(classes):
        if start not in unseen:
            continue
        current, cycle = start, []
        while current in unseen:
            unseen.remove(current)
            cycle.append(current)
            successor = label[4*current % modulus]
            need(pow(successor,3,modulus) == 64*pow(current,3,modulus) % modulus,
                 'cubing intertwines fourth-power layer action with inverse block')
            current = successor
        need(current == start, 'layer action is a permutation cycle')
        orbits.append(tuple(cycle))
    return tuple(orbits)


def torsor_source(direction_representative, alpha, modulus):
    need(gcd(alpha,modulus) == 1, 'hub-normalized source torsor must be a unit')
    return (alpha*pow(direction_representative,3,modulus)-1)*pow(3,-1,modulus) % modulus


def normalize_source(source, alpha, modulus):
    need(gcd(alpha,modulus) == 1, 'cannot normalize a collapsed hub')
    return (3*source+1)*pow(alpha,-1,modulus) % modulus


def p_adic_arc(x,y,p,depth):
    """Intrinsic orientation on distinct residues in Z/(p^depth)."""
    difference = (y-x) % p**depth
    if difference == 0:
        return 0
    while difference % p == 0:
        difference //= p
    return 1 if pow(difference,(p-1)//2,p) == 1 else -1


def geometry(p, depth):
    """Enumerate every cellular triangle and edge in every valuation layer."""
    total = p**depth
    edge_owner, link_components, strata = {}, Counter(), []
    for v in range(depth):
        scale, reduced = p**v, p**(depth-v)
        modulus, h, classes, _, _ = layers(p,depth-v)
        omega = h[1]
        need(modulus == reduced and (1+omega+omega*omega) % reduced == 0,
             'Eisenstein direction relation')
        layer_count, edge_count, face_count = 0, 0, 0
        for fibre in range(scale):
            for representative, directions in classes.items():
                edges = {tuple(sorted((x,(x+d) % reduced)))
                         for x in range(reduced) for d in directions}
                need(all(p_adic_arc(fibre+scale*x, fibre+scale*((x+d) % reduced),p,depth) == 1
                         for x in range(reduced) for d in directions),
                     'positive torus directions retain the intrinsic tournament orientation')
                q, wq = representative, representative*omega % reduced
                faces = set()
                for x in range(reduced):
                    faces.add(tuple(sorted((x,(x+q) % reduced,(x+q+wq) % reduced))))
                    faces.add(tuple(sorted((x,(x+wq) % reduced,(x+q+wq) % reduced))))
                need(len(edges) == 3*reduced and len(faces) == 2*reduced,
                     'one torus has exact vertex/edge/face counts')
                boundary_counts, links = Counter(), defaultdict(lambda:defaultdict(set))
                for face in faces:
                    need(len(set(face)) == 3, 'triangle has distinct vertices')
                    a,b,c = face
                    for x,y in ((a,b),(a,c),(b,c)):
                        boundary_counts[tuple(sorted((x,y)))] += 1
                    for x,y,z in ((a,b,c),(b,a,c),(c,a,b)):
                        links[x][y].add(z)
                        links[x][z].add(y)
                need(set(boundary_counts) == edges and set(boundary_counts.values()) == {2},
                     'each layer edge has two incident triangles')
                for x in range(reduced):
                    link = links[x]
                    need(len(link) == 6 and all(len(neighbors) == 2 for neighbors in link.values()),
                         'six degree-two cellular link vertices')
                    seen, stack = set(), [next(iter(link))]
                    while stack:
                        y = stack.pop()
                        if y not in seen:
                            seen.add(y)
                            stack.extend(link[y]-seen)
                    need(len(seen) == 6, 'one connected hexagonal cellular link')
                    link_components[fibre+scale*x] += 1
                for x,y in edges:
                    lifted = tuple(sorted((fibre+scale*x,fibre+scale*y)))
                    need(lifted not in edge_owner, 'all torus layers have disjoint edges')
                    edge_owner[lifted] = (v,fibre,representative)
                layer_count += 1
                edge_count += len(edges)
                face_count += len(faces)
        need(layer_count == (p-1)*p**(depth-1)//6, 'same torus count in every valuation stratum')
        strata.append(dict(valuation=v,tori=layer_count,vertices_per_torus=reduced,
                           edges=edge_count,faces=face_count))
    need(len(edge_owner) == total*(total-1)//2, 'all complete-graph edges partitioned exactly')
    need(set(link_components) == set(range(total)) and set(link_components.values()) == {(total-1)//6},
         'glued cellular link is the stated disjoint union of hexagons')
    outdegrees = Counter()
    for x in range(total):
        for y in range(x+1,total):
            sign = p_adic_arc(x,y,p,depth)
            need(sign in (-1,1) and p_adic_arc(y,x,p,depth) == -sign,
                 'exactly one tournament orientation for every distinct pair')
            first = next(i for i in range(depth) if (x//p**i) % p != (y//p**i) % p)
            digit_difference = (y//p**first-x//p**first) % p
            need(sign == (1 if pow(digit_difference,(p-1)//2,p) == 1 else -1),
                 'least-significant differing digit is the lexicographic Paley factor')
            outdegrees[x if sign == 1 else y] += 1
    need(set(outdegrees.values()) == {(total-1)//2}, 'regular iterated lexicographic tournament')
    return dict(prime=p,depth=depth,vertices=total,strata=strata,edges=len(edge_owner),
                normalized_tori=sum(s['tori'] for s in strata),
                hexagonal_links_per_glued_vertex=(total-1)//6)


def main():
    census, literal_edges = [], 0
    for p in (7,19,31,127,5779,87211):
        modulus,h,classes,label,cubes = layers(p)
        d = order(2,p,p-1)
        block = order(64,p,p-1)
        direction = order(4,p,p-1)
        need(block == d//gcd(d,6), 'inverse-ray block period')
        orbits = layer_orbits(p)
        need(all(len(cycle) == block for cycle in orbits), 'all layer clock orbits have the exact block order')
        need(len(classes) == (p-1)//6 and len(orbits)*block == len(classes), 'visited and unvisited layer count')
        alpha, first = 16 % p, label[1]
        cycle = next(c for c in orbits if first in c)
        current = first
        for b in range(block):
            # Independent actual positive integer, exact odd edge to root1.
            exponent = 4+6*b
            source = (2**exponent-1)//3
            need(source > 1 and source % 2 and 3*source+1 == 2**exponent,
                 'literal first-hit inverse-root edge')
            need(torsor_source(current,alpha,p) == source % p, 'layer torsor matches the actual root ray')
            need(normalize_source(source,alpha,p) == pow(64,b,p), 'normalized inverse source phase')
            need(cubes[pow(64,b,p)] == current, 'cube decoder returns exactly the layer')
            current = label[4*current % p]
            literal_edges += 1
        need(current == first and set(cycle) == {cubes[pow(64,b,p)] for b in range(block)},
             'the actual ray reaches exactly one geometric layer orbit')
        for q in classes:
            source = torsor_source(q,alpha,p)
            successor = torsor_source(label[4*q % p],alpha,p)
            need(successor == (64*source+21) % p, 'full anchored source-torsor conjugacy')
        census.append(dict(prime=p,layers=len(classes),block_order=block,direction_order=direction,
                           clock_orbits=len(orbits),visited_layers=block,unvisited_layers=len(classes)-block,
                           unvisited_orbits=len(orbits)-1,shared_ternary_depth=v3(block),
                           multiplicative_cube_section_exists=(len(classes)%3 != 0)))

    # Geometric direction children and dynamical orbit children are distinct.
    _,_,parent,parent_label,_ = layers(19,1)
    _,_,child,_,_ = layers(19,2)
    fibres = Counter(parent_label[q % 19] for q in child)
    need(set(fibres) == set(parent) and set(fibres.values()) == {19}, 'nineteen child direction layers')
    need(len(layer_orbits(19,1)) == len(layer_orbits(19,2)) == 1,
         'one child clock orbit despite nineteen child layer labels')
    geometric = [geometry(p,k) for p,k in ((7,1),(19,1),(7,2),(19,2))]

    # A chosen direction retains a cubic coordinate that the layer forgot.
    _,h,_,_,_ = layers(19)
    roots = [x for x in range(1,19) if pow(x,9,19) == 1 and pow(x,3,19) == 7]
    need(roots == [4,6,9] and all(order(x,19,18) == 9 for x in roots),
         'no multiplicative section above this order-three cube')
    source0 = (2**18-1)//3
    need(source0 % 19 == 0 and normalize_source(source0,1,19) == 1,
         'zero source residue is harmless when the hub coordinate is a unit')
    collapsed = [(2**(4+6*b)*19-1)//3 for b in range(3)]
    need({n % 19 for n in collapsed} == {6}, 'hub divisible by nineteen collapses the residue clock')
    try:
        normalize_source(collapsed[0],16*19,19)
    except ArithmeticError:
        pass
    else:
        raise ArithmeticError('collapsed zero-over-zero normalization accepted')
    need([((3*n+1)//19)*pow(16,-1,19) % 19 for n in collapsed] == [pow(64,b,19) for b in range(3)],
         'retaining the peeled valuation restores the direction phase')
    # At p13 the square-cube quotient exists, but opposite direction cosets
    # yield the same undirected layer and the square rule is not a tournament.
    q13={x*x % 13 for x in range(1,13)}
    h13={x for x in range(1,13) if pow(x,3,13)==1}
    need(12 in q13 and q13 == h13|{-x % 13 for x in h13}, 'p13 geometry hostile')
    print(json.dumps(dict(status='PROVED torsor and prime-power layer mechanisms; FINITE-EXACT declared controls',
        prime_layer_clocks=census,literal_root_edges=literal_edges,prime_power_geometries=geometric,
        direction_vs_clock_children=dict(prime=19,parent_depth=1,direction_children_per_parent=19,
                                        child_clock_orbits_per_parent=1),
        cube_section_hostile=dict(prime=19,cube=7,cube_roots=roots,root_orders=[9,9,9]),
        normalized_zero_source=source0,collapsed_hub=19,
        boundary='Layer quotient cubing acts on multiplicative directions, not additive vertices or forward Collatz time.'),indent=2))
    print('PASS: exact arithmetic; every assertion remains active under -O')


if __name__ == '__main__':
    main()
