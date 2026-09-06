#!/usr/bin/env python3
"""Exact fixed-pair sampler: uniform Eulerian isomorphism classes, not labels.

Precompute one integer weight per partition; each subsequent sample uses
polynomial-size GF(2) elimination and a random relabelling. No graph
isomorphism or automorphism oracle is used by the sampler. The exhaustive
small-order verifier deliberately uses a separate brute permutation path.
"""
import argparse
from bisect import bisect_right
from collections import Counter
from itertools import permutations
from math import factorial, gcd
import hashlib
import json
import random


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def partitions(n, least=1):
    if n == 0:
        yield ()
    for first in range(least, n+1):
        for tail in partitions(n-first, first):
            yield (first,) + tail


def permutation(parts):
    answer = []
    start = 0
    for length in parts:
        answer.extend(range(start+1, start+length))
        answer.append(start)
        start += length
    return tuple(answer)


def edge_system(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    return edges, {edge:i for i,edge in enumerate(edges)}


def transform(mask, perm, edges, index):
    answer = 0
    for k,(i,j) in enumerate(edges):
        if mask >> k & 1:
            answer |= 1 << index[tuple(sorted((perm[i],perm[j])))]
    return answer


def is_even(mask, n, edges):
    parity = [0]*n
    for k,(i,j) in enumerate(edges):
        if mask >> k & 1:
            parity[i] ^= 1
            parity[j] ^= 1
    return not any(parity)


def fixed_kernel(n, parts):
    """Literal edge orbits and parity rows, then a free-coordinate basis."""
    edges, index = edge_system(n)
    perm = permutation(parts)
    unseen = set(range(len(edges)))
    orbits = []
    while unseen:
        first = min(unseen)
        current = first
        orbit = 0
        while current in unseen:
            unseen.remove(current)
            orbit |= 1 << current
            i,j = edges[current]
            current = index[tuple(sorted((perm[i],perm[j])))]
        require(current == first, "edge orbit must close at its origin")
        orbits.append(orbit)
    rows = [0]*n
    for c,orbit in enumerate(orbits):
        for k,(i,j) in enumerate(edges):
            if orbit >> k & 1:
                rows[i] ^= 1 << c
                rows[j] ^= 1 << c
    pivots = []
    rank = 0
    for c in range(len(orbits)):
        found = next((j for j in range(rank,n) if rows[j] >> c & 1), None)
        if found is None:
            continue
        rows[rank],rows[found] = rows[found],rows[rank]
        for j in range(n):
            if j != rank and rows[j] >> c & 1:
                rows[j] ^= rows[rank]
        pivots.append(c)
        rank += 1
    free = [c for c in range(len(orbits)) if c not in pivots]
    basis = []
    for c in free:
        vector = 1 << c
        for row,pivot in zip(rows,pivots):
            if row >> c & 1:
                vector |= 1 << pivot
        graph = 0
        for j,orbit in enumerate(orbits):
            if vector >> j & 1:
                graph ^= orbit
        basis.append(graph)
    return basis, len(orbits), rank


def closed_dimension(parts):
    b = sum(a//2 for a in parts)
    b += sum(gcd(a,c) for i,a in enumerate(parts) for c in parts[i+1:])
    return b-len(parts)+int(any(a % 2 for a in parts))


def class_size(n, parts):
    z = 1
    for length,count in Counter(parts).items():
        z *= length**count * factorial(count)
    return factorial(n)//z


class OrbitSampler:
    def __init__(self, n):
        require(n >= 1, "positive vertex count")
        self.n = n
        self.types = list(partitions(n))
        self.cumulative = []
        total = 0
        for parts in self.types:
            total += class_size(n,parts) << closed_dimension(parts)
            self.cumulative.append(total)
        self.total = total
        require(total % factorial(n) == 0, "Burnside divisibility")
        self.orbit_count = total//factorial(n)
        self.edges, self.index = edge_system(n)

    def sample(self, rng):
        parts = self.types[bisect_right(self.cumulative,rng.randrange(self.total))]
        basis,b,rank = fixed_kernel(self.n,parts)
        require(len(basis)==closed_dimension(parts)==b-rank, "sample kernel dimension")
        bits = rng.getrandbits(len(basis))
        graph = 0
        for j,vector in enumerate(basis):
            if bits >> j & 1:
                graph ^= vector
        sigma = list(range(self.n))
        rng.shuffle(sigma)
        labelled = transform(graph,sigma,self.edges,self.index)
        canonical = permutation(parts)
        conjugated = [0]*self.n
        for i in range(self.n):
            conjugated[sigma[i]] = sigma[canonical[i]]
        return labelled, tuple(conjugated), parts


def kernel_vectors(basis):
    vectors = [0]
    for vector in basis:
        vectors += [x ^ vector for x in vectors]
    return vectors


def verify_small():
    gates = 0
    report = []
    for n in range(1,6):
        edges,index = edge_system(n)
        states = [x for x in range(1 << len(edges)) if is_even(x,n,edges)]
        perms = list(permutations(range(n)))
        orbit = {x:min(transform(x,p,edges,index) for p in perms) for x in states}
        sizes = Counter(orbit.values())
        totals = Counter()
        sampler = OrbitSampler(n)
        for parts in sampler.types:
            basis,b,rank = fixed_kernel(n,parts)
            literal = {x for x in states if transform(x,permutation(parts),edges,index)==x}
            generated = kernel_vectors(basis)
            require(len(generated)==len(set(generated)), "kernel coordinate bijection")
            require(set(generated)==literal, "literal fixed-space universe")
            require(len(basis)==closed_dimension(parts)==b-rank, "closed rank")
            gates += 3
            for graph in generated:
                totals[orbit[graph]] += class_size(n,parts)
        require(set(totals)==set(sizes), "every graph class is reached")
        require(set(totals.values())=={factorial(n)}, "exact equal class mass")
        require(sampler.orbit_count==len(sizes), "independent class count")
        gates += 3
        # An independent direct enumeration of all fixed pairs verifies
        # the labelled inverse-orbit-size law, not merely class totals.
        for graph in states:
            stabilizer = sum(transform(graph,p,edges,index)==graph for p in perms)
            require(stabilizer*sizes[orbit[graph]]==factorial(n), "labelled fixed-pair law")
            gates += 1
        report.append({"n":n,"labels":len(states),"classes":len(sizes),
                       "class_sizes":sorted(sizes.values()),"fixed_pairs":sampler.total})
    require(report[3]["class_sizes"]==[1,3,4], "hostile: uniform labels not uniform classes")
    gates += 1
    return gates,report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes",type=int,default=30)
    parser.add_argument("--samples",type=int,default=64)
    parser.add_argument("--seed",type=int,default=240300737)
    args = parser.parse_args()
    require(args.nodes >= 1 and args.samples >= 1, "positive nodes and samples")
    gates,small = verify_small()
    for row in small:
        print("EXHAUSTIVE FIXED-PAIR CONTROL", json.dumps(row,sort_keys=True))
    rank_count = 0
    for n in range(1,13):
        for parts in partitions(n):
            basis,b,rank = fixed_kernel(n,parts)
            require(len(basis)==closed_dimension(parts)==b-rank, "independent rank bank")
            gates += 1
            rank_count += 1
    sampler = OrbitSampler(args.nodes)
    rng = random.Random(args.seed)
    samples = []
    for _ in range(args.samples):
        graph,perm,parts = sampler.sample(rng)
        require(is_even(graph,args.nodes,sampler.edges), "sample degree parity")
        require(transform(graph,perm,sampler.edges,sampler.index)==graph, "sample fixed permutation")
        gates += 2
        samples.append([str(graph),list(perm),list(parts)])
    print("INDEPENDENT LITERAL RANK BANK", rank_count)
    print("PREPROCESS",json.dumps({"n":args.nodes,"partitions":len(sampler.types),
          "classes":str(sampler.orbit_count),"fixed_pair_integer_bits":sampler.total.bit_length()},sort_keys=True))
    print("SEEDED SAMPLE CONTROLS",args.samples,"seed",args.seed)
    print("CHECKS",gates)
    payload={"small":small,"rank_bank":rank_count,"n":args.nodes,
             "count":str(sampler.orbit_count),"samples":samples}
    print("SEMANTIC SHA256",hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest())
    print("PROVED exact uniform-orbit law; subexponential preprocessing; expected polynomial time per sample")
    print("NOT uniform labelled graphs; no isomorphism oracle, Markov mixing claim, or Boolean gap bound")


if __name__ == "__main__":
    main()
