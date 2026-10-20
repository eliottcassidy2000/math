"""Exact clone/run growth and alternating-core recurrence.

Run: python 05-knowledge/results/weighted_core_transport_20261020/growth_check.py
Only standard library and the frozen prior packet are used.
All tests remain active with python -O. Default generation ends at order 8.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, factorial
from pathlib import Path
import sys

PRIOR = Path(__file__).resolve().parent.parent / "metagraph_growth_20261020"
sys.path.insert(0, str(PRIOR))
import enumeration as E
import growth_check as G


def need(value, message):
    if not value:
        raise RuntimeError(message)


def compositions(n, b):
    if b == 1:
        yield (n,)
    elif b > 1:
        for head in range(1, n - b + 2):
            for tail in compositions(n - head, b - 1):
                yield (head,) + tail


def expand(q, weights):
    blocks = tuple(i for i, w in enumerate(weights) for _ in range(w))
    return tuple(sum((int(v < u) if blocks[v] == blocks[u]
                      else G.arc(q, blocks[v], blocks[u])) << u
                     for u in range(len(blocks)) if u != v)
                 for v in range(len(blocks)))


def kernel(n, b):
    out = Counter({b: n + b})
    out[b + 1] += (1 << b) - 2 * b
    out[b + 2] += ((1 << b) - 1) * (n - b)
    for s in range(2, n - b + 1):
        out[b + 1 + s] += (1 << b) * comb(n - b, s)
    return +out


def switches(chains, mask):
    return sum(((mask >> u) ^ (mask >> v)) & 1
               for c in chains for u, v in zip(c, c[1:]))


def alternating_masks(chains):
    for phases in product((0, 1), repeat=len(chains)):
        yield sum((phase ^ (j & 1)) << v
                  for c, phase in zip(chains, phases)
                  for j, v in enumerate(c))


def labeled_core(n):
    return sum((-1) ** (n - k) * comb(n - 1, k - 1)
               * (factorial(n) // factorial(k)) * (1 << comb(k, 2))
               for k in range(1, n + 1))


def local_audit():
    count = 0
    for n in range(1, 6):
        for code in range(1 << comb(n, 2)):
            rows = E.matrix(n, code)
            chains = G.clone_chains(rows)
            b = len(chains)
            hist = Counter()
            old_pairs = {(u, v) for c in chains
                         for u, v in zip(c, c[1:])}
            alternate = set(alternating_masks(chains))
            for mask in range(1 << n):
                child = E.extension(rows, mask)
                new_chains = G.clone_chains(child)
                pairs = {(u, v) for c in new_chains
                         for u, v in zip(c, c[1:])}
                d = sum(n in pair for pair in pairs)
                s = switches(chains, mask)
                retained = {(u, v) for u, v in old_pairs
                            if ((mask >> u) & 1) == ((mask >> v) & 1)}
                need({p for p in pairs if n not in p} == retained,
                     "old homogeneous edges split exactly at switches")
                need(len(new_chains) == b + 1 + s - d, "component identity")
                if d:
                    need((s, d) in {(0, 1), (1, 2)}, "clone dichotomy")
                    need(len(new_chains) == b, "clone core unchanged")
                    containing = next(c for c in new_chains if n in c)
                    old = tuple(v for v in containing if v != n)
                    need(old in chains, "clone joins exactly one old chain")
                    need(set(new_chains) - {containing} == set(chains) - {old},
                         "all other chains unchanged")
                else:
                    need(len(new_chains) == b + 1 + s, "strict core growth")
                need((len(new_chains) == n + 1)
                     == (mask in alternate and d == 0), "alternating criterion")
                hist[len(new_chains)] += 1
                count += 1
            need(hist == kernel(n, b), "closed labelled core-size kernel")
    print("EXHAUSTIVE labelled parents orders 1..5:", count,
          "attachments; clone dichotomy, split law, kernel, alternating gate PASS")


def recursive_bank(max_n):
    # No known A000568 values enter generation or its stopping criterion.
    cores = {1: {0: (0,)}}
    classes = {1: {0: (0,)}}
    automorphisms = {1: {0: 1}}
    print("n classes pair_free_cores alternating_trials accepted_trials exact_core_mass")
    print("1 1 1 0 0 1")
    for n in range(2, max_n + 1):
        current = {}
        trials = accepted = 0
        mass = Fraction()
        rooted_multiplicity = Counter()
        for code, parent in classes[n - 1].items():
            chains = G.clone_chains(parent)
            ap = automorphisms[n - 1][code]
            for mask in alternating_masks(chains):
                trials += 1
                child = E.extension(parent, mask)
                if len(G.clone_chains(child)) != n:
                    continue
                accepted += 1
                cc, ac = E.canonical(child, aut=True)
                current[cc] = E.decode(n, cc)
                contribution = Fraction(ac, n * ap)
                mass += contribution
                rooted_multiplicity[cc] += contribution
        need(all(value == 1 for value in rooted_multiplicity.values()),
             "per-child inverse-root mass exactly one")
        need(mass == len(current), "symmetry-corrected core recurrence")
        cores[n] = current
        reps = {}
        auts = {}
        for b, qs in cores.items():
            for q in qs.values():
                for weights in compositions(n, b):
                    child = expand(q, weights)
                    cc, ac = E.canonical(child, aut=True)
                    reps[cc] = E.decode(n, cc)
                    auts[cc] = ac
        classes[n], automorphisms[n] = reps, auts
        labeled = Counter()
        for code, rows in reps.items():
            b = len(G.clone_chains(rows))
            labeled[b] += factorial(n) // auts[code]
        need(sum(labeled.values()) == (1 << comb(n, 2)), "full labelled orbit mass")
        need(labeled[n] == labeled_core(n), "independent signed-binomial core mass")
        for b in range(1, n + 1):
            expected = (labeled_core(b) * (factorial(n) // factorial(b))
                        * comb(n - 1, b - 1))
            need(labeled[b] == expected, "labelled substitution stratum")
        print(n, len(reps), len(current), trials, accepted, mass, flush=True)
    # Comparison values are used only AFTER independent generation.
    known = (0, 1, 1, 2, 4, 12, 56, 456, 6880)
    need(all(len(classes[n]) == known[n] for n in range(1, max_n + 1)),
         "A000568 comparison")
    # A fully separate path-presentation universe checks every class <= 7.
    for n in range(1, min(max_n, 7) + 1):
        fixed = {E.canonical(E.tiling(n, t))
                 for t in range(1 << comb(n - 1, 2))}
        need(fixed == set(classes[n]), "independent fixed-path cube equality")
    print("Complete fixed-path cubes orders 1..7 independently match the generated bank.")
    return cores, classes


def kernel_mass_audit():
    distribution = Counter({1: 1})
    for n in range(1, 20):
        child = Counter()
        for b, count in distribution.items():
            for target, multiplicity in kernel(n, b).items():
                child[target] += count * multiplicity
        need(sum(child.values()) == (1 << comb(n + 1, 2)), "kernel total")
        need(child[n + 1] == labeled_core(n + 1), "kernel / substitution agreement")
        distribution = child
    print("Labelled kernel and formal-series inversion agree through order 20.")
    print("Labelled pair-free core counts n=1..8:",
          [labeled_core(n) for n in range(1, 9)])


def hostile_controls():
    # Same n,b makes the coarse growth kernel agree; target class response need not.
    cycle = E.matrix(3, 5)  # arcs 0->1, 1->2, 2->0
    need(len(G.clone_chains(cycle)) == 3, "C3 control")
    parents = (
        E.extension(cycle, 0),           # new source
        E.extension(cycle, 7),           # new sink
    )
    need(all(len(G.clone_chains(t)) == 4 for t in parents), "pair-free parents")
    profiles = [Counter(E.canonical(E.extension(t, mask))
                        for mask in range(16)) for t in parents]
    need(profiles[0] != profiles[1], "n,b does not determine ordinary target profile")
    print("HOSTILE n=4,b=4 source+C3 versus C3+sink: equal core-size kernels,",
          "different ordinary-class extension profiles.")
    # A single core cannot encode ordinary tau targets with bounded core size.
    for n in range(5, 41):
        transitive = E.matrix(n, (1 << comb(n, 2)) - 1)
        mate = G.tau(transitive)
        need(len(G.clone_chains(mate)) == n, "tau transitive has full core")
    print("HOSTILE tau(transitive chain): 1 core vertex -> n core vertices,",
          "verified n=5..40 (general proof in allocation.md).")
    print("HOSTILE naive rooted count: C3 has one TT2-parent alternating attachment;",
          "1/3 is not 1. Aut(C3)/Aut(TT2)=3 repairs it.")


def cycle_index_hostile(cores):
    # All rigid cores of the same order have exactly the same cycle index.
    # A differing next-core automorphism profile therefore defeats even that state.
    q = 6
    profiles = {}
    for rows in cores[q].values():
        if E.canonical(rows, aut=True)[1] != 1:
            continue
        profile = Counter()
        for mask in range(1 << q):
            child = E.extension(rows, mask)
            if len(G.clone_chains(child)) == q + 1:
                profile[E.canonical(child, aut=True)[1]] += 1
        key = tuple(sorted(profile.items()))
        profiles.setdefault(key, rows)
    need(len(profiles) > 1, "equal cycle indices must fail this hostile")
    first = sorted(profiles.items())[:2]
    print("HOSTILE same complete core cycle index p1^6, different pair-free",
          "child automorphism-order profiles:")
    for profile, rows in first:
        print("  rows=", rows, "child_Aut_histogram=", profile)


def main():
    local_audit()
    kernel_mass_audit()
    hostile_controls()
    cores, _ = recursive_bank(8)
    cycle_index_hostile(cores)
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
