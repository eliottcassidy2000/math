---
id: THM-4622
title: "Choice and the chromatic number of the plane: in ZF the countable real-algebraic plane has an arithmetical, absolute chromatic number chi_alg with 5 <= chi_alg <= chi(R^2) <= 7 (countable de Bruijn-Erdos needs no choice), in ZFC chi(R^2) = chi_alg (Tarski transfer of finite unit-distance graphs plus de Bruijn-Erdos), so choice can only raise the chromatic number of the plane and never lowers it; every countable field plane and every residue colouring is choice-free and absolute; and under the Frucht/Beckman-Quarles lens Aut(UD(R^2)) is the isometry group (KNOWN) while Aut(UD(Q^2)) contains the full symmetric group on its infinitely many isomorphic components (Pythagorean denominators), so the rational plane is as far from rigid as a unit-distance graph can be"
status: >
  KNOWN in substance: the ZFC equality chi(R^2) = chi_alg (de Bruijn-Erdos plus the folklore that
  every finite unit-distance graph has a real-algebraic realization; axiom-sensitivity is discussed
  by Shelah-Soifer 2003/2004), countable de Bruijn-Erdos without choice (the standard Konig argument
  on the definable tree; folklore), the component structure of UD(Q^2) (classical, Woodall 1973
  context) and Beckman-Quarles 1953. PROVED here as written (elementary): the ZF chain, the
  arithmetical/absolute status of chi_alg (the bound 5 rests on the SAT-verified 5-chromatic
  graphs of de Grey / Heule / Parts, finite computations), the absoluteness statement for computably
  presented fields, and the Sym(omega) remark. Not a new bound on chi(R^2). Independently audited
  2026-10-12 (audit B: PASS WITH CORRECTIONS, applied: the 509-vertex graph is Parts's, Heule's is
  510; statement 5 restricted to computably presented L; KNOWN typing for 1 and 6; the density
  assertion replaced by 'not in a line'; Shelah-Soifer's 2004 plane examples noted; MISTAKE-594).
  Inspired by Pinsky (arXiv:2305.11382) and Fernandes-Mezabarba-Rodrigues (arXiv:2609.20140).
session: mac-mini-2026-10-12-logic4
source: 05-knowledge/results/logic_lenses_jc2_collatz_hn_20261012.md
related:
  - THM-4558 (residue colourings of field planes; Lemma R; the table kappa(q))
  - THM-4624 (automorphism census of the finite unit-circle graphs G_q)
  - THM-418 (lattice chromatic dichotomy), THM-4566 (Hilbert class field planes)
  - S. Shelah, A. Soifer, Axiom of choice and chromatic number of the plane, J. Combin. Theory Ser. A 103 (2003); Axiom of choice and chromatic number: examples on the plane, JCTA 105 (2004)
  - J. Parts, Graph minimization, focusing on the example of 5-chromatic unit-distance graphs in the plane, arXiv:2010.12665 (509 vertices); M. Heule, 510 vertices (FoCM 2023)
  - N. G. de Bruijn, P. Erdos, A colour problem for infinite graphs and a problem in the theory of relations (1951)
  - F. S. Beckman, D. A. Quarles, On isometries of Euclidean spaces, Proc. AMS 4 (1953)
  - A. D. N. J. de Grey, The chromatic number of the plane is at least 5, Geombinatorics 28 (2018)
---

# THM-4622 — choice, absoluteness and automorphisms for the chromatic number of the plane

## Setting

`UD(S)` is the unit-distance graph on a subset `S ⊆ C`: `x ~ y` iff `|x - y| = 1`. `R_alg` is the
field of real algebraic numbers; `UD(R_alg^2)` is a countable graph with a computable vertex
enumeration and a decidable edge relation. `χ_alg := χ(UD(R_alg^2))`. For a countable field `L ⊂ C`
closed under complex conjugation, `G(L)` is the field plane of THM-4558.

## Statement

1. **Countable de Bruijn–Erdős in ZF (KNOWN, folklore).** Let `G` be a graph whose vertex set is
   enumerated by `ω`. If every finite subgraph of `G` is `k`-colourable, then `G` is `k`-colourable. No
   form of choice is used.
2. **The real-algebraic plane (ZF).** For each `k`, "`χ_alg ≤ k`" is equivalent to the `Π^0_1`
   arithmetical statement "every finite subgraph of `UD(R_alg^2)` is `k`-colourable"; so `χ_alg` has the
   same value in every transitive model of ZF, and `5 ≤ χ_alg ≤ 7`.
3. **The chain (ZF).** `χ_alg ≤ χ(R^2) ≤ 7`.
4. **The equality (ZFC).** `χ(R^2) = χ_alg`. Hence the Hadwiger–Nelson number of ZFC is the
   arithmetical invariant `χ_alg`; a model of ZF in which `χ(R^2) > χ_alg` would have to use the
   non-algebraic points, and no model has `χ(R^2) < χ_alg`.
5. **Field planes (ZF).** For every conjugation-closed `L ⊆ Q̄` given by a computable enumeration
   (every number field, the Hilbert class fields of THM-4566, the odd-`t` compositum), the statement
   "`χ(G(L)) ≤ k`" is `Π^0_1` and has the same truth value in every transitive model of ZF; for an
   arbitrary countable `L` it is absolute between transitive models `M ⊆ N` with `L ∈ M` countable in
   `M`. The residue colourings of THM-4558 (Lemma R, on number fields) are choice-free: coset
   representatives of the unit-vector group `Γ ⊂ L` are chosen as least elements in a fixed
   enumeration of `L`.
6. **Frucht/Beckman–Quarles.** `Aut(UD(R^2)) = Isom(R^2)` (KNOWN: a unit-distance-preserving map of
   `R^n`, `n ≥ 2`, is an isometry). In contrast, `UD(Q^2)` has infinitely many connected components,
   all isomorphic (translates of the subgroup `Γ_Q` generated by the rational unit vectors; KNOWN,
   classical), and `Aut(UD(Q^2))` contains a copy of the symmetric group on the set of components; it
   has cardinality `2^(ℵ_0)`, and almost all of its elements are not induced by isometries of `R^2`.

## Proofs

**(1).** Let `T` be the tree of proper partial colourings `c : {0, ..., m-1} → {0, ..., k-1}` of the
induced subgraphs on initial segments, ordered by extension. It is finitely branching, and has nodes
at every level by hypothesis. Call a node *good* if it has descendants at every level; the root is
good, and a good node has a good child (otherwise, with finitely many children each dying out at some
finite level, the node would die out at the maximum of those levels). Define `c_0 = ∅` and `c_(m+1)` =
the lexicographically least good child of `c_m`. Every step is a least-element choice from a finite
set, so the sequence is defined by recursion on `ω` without choice, and `∪_m c_m` is a proper
`k`-colouring of `G`.

**(2).** By (1), `UD(R_alg^2)` is `k`-colourable iff all its finite subgraphs are. Finite subsets of
`R_alg^2` are enumerable and the edge relation and `k`-colourability of a finite graph are decidable,
so the statement is `Π^0_1`; arithmetical statements have the same truth value in every transitive
model of ZF (they are decided by the standard `ω`). Upper bound `7`: the hexagonal seven-colouring is
an explicit map. Lower bound `5`: de Grey's 1581-vertex unit-distance graph (or the 510-vertex graph
of Heule, the 509-vertex graph of Parts) is 5-chromatic — a finite, SAT-verified fact — and has a
realization in `R^2`; by (4)'s transfer argument below it has one in `R_alg^2` (in fact the published
coordinates are algebraic).

**(3).** A colouring of `R^2` restricts to `R_alg^2`; the seven-colouring is explicit.

**(4).** In ZFC, de Bruijn–Erdős gives `χ(R^2) = max_G χ(G)` over finite subgraphs `G ⊆ UD(R^2)`. For
such a `G` on `n` points, the existential sentence
`∃ (x_i, y_i)_(i ≤ n) [ ∧_(ij ∈ E(G)) (x_i - x_j)^2 + (y_i - y_j)^2 = 1 ∧ ∧_(i<j) (x_i, y_i) ≠ (x_j, y_j) ]`
holds in the real closed field `R`, hence (Tarski: real closed fields are elementarily equivalent) in
`R_alg`. So `G` is a subgraph of `UD(R_alg^2)` (possibly with extra unit distances, which only raises
the chromatic number of the host) and `χ_alg ≥ χ(G)`. Hence `χ_alg ≥ χ(R^2)`, and (3) gives equality.

**(5).** For computably presented `L ⊆ Q̄` the graph `G(L)` is computable (algebraic arithmetic and
equality are decidable), so (2) applies verbatim. For an arbitrary countable `L ∈ M`, finite subsets of
`L` and `k`-colourability of finite graphs are absolute between `M` and `N`, and (1) holds in both, so
"`χ(G(L)) ≤ k`" is absolute. A countable `L` carries a choice-free enumeration; the components of `G(L)`
are the cosets `x_0 + Γ`, and Lemma R colours `x ∈ x_0 + Γ` by `φ(red_𝔓(x - x_0))` with `x_0` the least
element of its coset and `φ` a proper colouring of a finite Cayley graph (a finite object).

**(6).** The Beckman–Quarles theorem is CITED. For `Q^2`: `UD(Q^2)` is the Cayley graph of `(Q^2, +)` on
the rational unit vectors, so its components are the cosets of `Γ_Q` and translation makes them
pairwise isomorphic. A rational unit vector in lowest terms is `(a/c, b/c)` with `a^2 + b^2 = c^2`,
`gcd(a, b, c) = 1`, and every prime factor of `c` is `≡ 1 (mod 4)` (a primitive hypotenuse is a sum of
two coprime squares of opposite parity). So every element of `Γ_Q` has coordinates whose denominators
are products of primes `≡ 1 (mod 4)`, and the points `(3^(-k), 0)`, `k ≥ 0`, lie in pairwise distinct
cosets: infinitely (countably) many components. Choose (least) base points `p_C` in each component
`C`; for a permutation `π` of the components, `x ↦ x - p_C + p_(π C)` (`x ∈ C`) is an automorphism.
These form a copy of `Sym(ω)`, of cardinality `2^(ℵ_0)`. An automorphism that translates two
components by different vectors `t_1 ≠ t_2` is not an isometry: the condition
`|x + t_1 - y - t_2| = |x - y|` says that `w = x - y` lies on the line `{w : |w + (t_1 - t_2)| = |w|}`, while
the set of differences `x - y` (`x ∈ C_1`, `y ∈ C_2`) is a coset of `Γ_Q ⊇ Z^2`, not contained in a line.  ∎

## Reading

* **Where choice can act.** The only ZF-sensitive object is `UD(R^2)` itself, and only upward:
  `χ(R^2) ≥ χ_alg` is a ZF theorem, while `χ(R^2) = χ_alg` is proved with choice. Whether some model
  of ZF has `χ(R^2) > χ_alg` is OPEN (the Shelah–Soifer question for the unit-distance graph; their
  2003 example of axiom-dependence is a different distance graph on `R`, and their 2004 paper,
  *Axiom of choice and chromatic number: examples on the plane*, JCTA 105, gives distance graphs on
  `R^2` with the same axiom-dependence — so the open case is specifically the unit distance). Every bound the repo proves on
  field planes, every residue colouring, every SAT certificate and every `κ(q)` value lives in the
  absolute, countable world.
* **Pinsky's lesson transported.** Frucht's theorem fails in permutation models because a graph
  realizing a group must be *supported*; here the analogous support is a finite set of algebraic
  numbers: every finite certificate (a 5-chromatic subgraph, a residue colouring) has finite support
  in `R_alg`, and the plane's chromatic number is decided by such certificates.
* **Bases and choice.** Fernandes–Mezabarba–Rodrigues show that "every vector space has a basis"
  implies AC even without Foundation. The Hamel-basis style colourings behind Shelah–Soifer's
  examples are therefore the natural place where choice enters distance-graph colourings; the
  unit-distance graph of the plane sidesteps them because its lower bounds are finite and its upper
  bound is a tiling.
