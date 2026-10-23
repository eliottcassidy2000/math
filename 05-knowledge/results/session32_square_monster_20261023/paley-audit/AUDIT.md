# Session 32 independent Paley–Monster normalizer audit

Date: 2026-10-23 UTC. Mathematical/source and producer-file audit: **PASS**
for the precisely qualified bridge below. Producer review is recorded in
section 8. This audit imports no producer code,
edits no producer source and makes no Git mutation.

## 1. The genuine bridge and its scope

For a Monster element g of class 7A, put K=<g> and N=N_M(K). The sourced
structures are

    C_M(K) = C7 × He,
    N = (7:3 × He):2,
    N ≅ (7:6) ×_{C2} (He:2).

Here the final expression is a fiber product, not a direct product. Its
two parity homomorphisms have kernels 7:3 and He. Projection onto 7:6 is
surjective and its kernel is He. Consequently

    N/He ≅ 7:6 ≅ AGL(1,7).

This is a concrete group quotient, not a numerical coincidence. The
orientation-preserving index-two subgroup maps onto 7:3. The subgroup
He is normal: it is the derived subgroup of the normal centralizer
C7 × He, since He is perfect.

An alternative elementary proof of the quotient avoids interpreting the
dot/colon notation. N/He has order 42 and a normal subgroup K of order 7;
its conjugation quotient on K has order six because the centralizer of K
has order 7|He|. Thus this quotient maps faithfully onto Aut(K)=C6. If x
lifts a generator of that C6 into N/He, then x^6 lies in K and is fixed
by x. A generator multiplier of F7* fixes no nonidentity element of K,
so x^6=1. Therefore x has order six and splits the extension, with its
faithful multiplier action. This gives C7:C6 directly.

## 2. Primary-source verification

- Robert A. Wilson, *The odd-local subgroups of the Monster*, J. Austral.
  Math. Soc. Series A 44 (1988), 1–16, Theorem 7 and §10, printed p.13.
  The 7A cyclic-subgroup normalizer is listed as (7:3 × He):2.
  https://doi.org/10.1017/S1446788700031323
- Jianbei An and Robert A. Wilson, *The Alperin Weight Conjecture and
  Uno's Conjecture for the Baby Monster B, p Odd*, LMS J. Comput. Math.
  7 (2004), 120–166, §3, printed p.124 (PDF p.5). They construct this
  subdirect product using the actions of 7:6 on seven points and He:2
  on 2058 points, matching the two index-two parities.
  https://doi.org/10.1112/S1461157000001078
- Thomas Breuer, Kay Magaard and Robert A. Wilson, *Some steps in the
  verification of the ordinary character table of the Monster group*,
  dated 12 December 2024, §4.13, printed p.20, explicitly gives both
  centralizer and normalizer and the subdirect product of 7:6 and He.2.
  Section 5 records the GAP construction with factors 7:3, 7:6, He, He.2.
  https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf

All three were checked in accessible primary text, rather than relying
on the encyclopedia entry or on a group-shape label alone. The audit
does not claim to reconstruct the Monster or rerun its character-table
verification. These are cited group-theoretic inputs.

## 3. Two seven-element actions that must not be identified

The affine action of N/He on F7 is transitive and faithful. It is the
left-coset action of 7:6 on the seven cosets of a complement C6, or,
equivalently, the N action on cosets of the preimage of that complement.
After choosing an origin and generator of the translation subgroup,
its elements act as x ↦ ax+b. The choice is not canonical.

In contrast, conjugation by N on the seven elements of K fixes the
identity. Its kernel is C_M(K)=C7 × He, and its image is only C6. In
affine coordinates,

    (x ↦ ax+b) (x ↦ x+k) (x ↦ ax+b)^(-1) = (x ↦ x+ak).

The translation b disappears. This action is k ↦ ak, not k ↦ ak+b.
It has the two orbits {0} and F7*. Any claim that conjugation on K
supplies the transitive seven-point Paley action is false.

## 4. Exact tournament symmetry

Define P7 on F7 by x→y iff y−x belongs to Q={1,2,4}. Squares preserve Q
and nonsquares exchange Q with −Q. Thus x↦ax+b preserves every arc if
a∈Q and reverses every arc if a∉Q.

There are no other automorphisms. Translations are transitive. An
automorphism fixing 0 must preserve its outgoing set Q, whose induced
tournament is the directed triangle 1→2→4→1. There are at most three
restrictions to Q. An automorphism fixing Q pointwise also fixes every
other vertex, since their adjacency signatures to Q distinguish them.
The three square multipliers realize the bound. Therefore

    Aut(P7) = F7 : Q ≅ 7:3, of order 21,
    Aut^±(P7) = AGL(1,7), of order 42.

Here Aut^± means all arc-preserving and globally arc-reversing
permutations, not the automorphism group of the underlying graph K7.
The odd coset contains 21 anti-automorphisms, only seven of which are
involutions. These are x↦−x+b; the remaining fourteen have order six.
The full N acts on the tournament only in this signed sense. Its
index-two subgroup (7:3 × He) acts by genuine automorphisms.

## 5. Independent finite-exact counts

`independent_checks.py` searches all 7! permutations, and independently
checks Hamiltonian paths by a subset dynamic program. Explicit exception
checks remain live under Python optimization. The normal and `-O`
receipts are identical.

- Directed Hamiltonian paths: H(P7)=189, with 27 ending at each vertex.
- Genuine automorphisms: 21; anti-automorphisms: 21.
- Fixed-path presentations: 189/21=9, each represented by exactly 21
  labelings. The automorphism action on ordered Hamiltonian paths is
  free, since an automorphism fixing a path fixes its vertex sequence.
- Blue means that the presentation is fixed by converse followed by
  reversal of the fixed path. There are three blue presentations,
  represented by 63 paths. There are six black presentations, giving
  fundamental paired-black degree three under the established convention.
- Each of the seven affine anti-involutions fixes nine paths under
  anti-automorphism followed by path reversal. The fourteen order-six
  anti-automorphisms fix none. The signed action on the 189 paths has
  orbit sizes 21,21,21,42,42,42. This signed action is not free.
- Additional controls: 14 directed 3-cycles, 42 directed 5-cycles and
  24 directed 7-cycles, counted modulo cyclic rotation, and seven
  disjoint odd-cycle pairs. Numerically 1+2·80+4·7=189.

The edge-bitstrings of the three blue presentations are in the receipt,
using pairs (i,j) in lexicographic order, 0≤i<j≤6. The count 24 is a
Hamiltonian-cycle count and must not be confused with H(P7)=189.

## 6. CFT citation boundary

Bae, Harvey, Lee, Lee and Rayhaun, *Conformal Field Theories with Sporadic
Group Symmetry*, arXiv:2002.02970v2, §3.3, proves existence of a P(7)
parafermion subalgebra of V^natural through the lattice embedding chain.
It says the commutant is expected to admit a Held action, and presents
evidence for the 7A lift of the fusion automorphism. Its proposed seven
P(3) subalgebras and the detailed Held-VOA picture remain prospective
in this discussion. These statements are not all proved embedding or
equivariant-construction theorems. The HTML equations are (328)–(330);
the published numbering around (3.173)–(3.175) differs.
https://arxiv.org/html/2002.02970v2
https://doi.org/10.1007/s00220-021-04207-7

No CFT claim is needed for the elementary affine/Paley quotient theorem.
The action constructed here is a finite-group action, not a claimed
construction of a VOA, a fusion category, a twined trace or a modular
transformation law for tournament counts.

## 7. What this does not prove about H=7 and H=21

1. The integers seven and twenty-one here are a vertex count and a group
   order. They are not Hamiltonian-path counts of P7, whose count is189.
2. The universal divisibility |Aut(T)| divides H(T) only implies that
   a tournament with automorphism group of order21 has H divisible by21.
   It does not rule out H=21. P7 supplies one value189, not a universal
   obstruction for all tournaments.
3. The quotient action belongs to a proper Monster normalizer, not to
   the whole Monster. A nontrivial seven-point Monster action would
   contradict simplicity and |M|>7!: it would have to be faithful.
4. He is the kernel of the affine action. Its presence as a kernel does
   not yield new Hamiltonian-path constraints or a faithful Held action
   on the seven vertices.
5. No map from all tournament classes to this local Monster structure,
   preserving Hamiltonian-path number or its absence, is supplied. No
   explanation or proof of the forbidden values {7,21} follows here.

The positive result is exactly the signed affine symmetry quotient and
its compatible Paley orientation. Any broader Hamiltonian-gap or
moonshine interpretation needs a further map and a preserved observable.

## 8. Final producer-file review

The producer's `PALEY_MONSTER_BRIDGE.md`, `check.py`, `exact.json` and
`optimized.json` were read and independently reviewed. No mathematical
correction was required. The auditor reran the producer checker both
normally and with `-O`; each replay byte-matches its saved receipt. The
auditor's separate checker also validates all nine producer presentations
and all three blue presentations by their full adjacency bitstrings,
not only their counts. Its final receipt pins all four producer files by
SHA-256. The producer checker hash is

    5821f709524b50c3acd416c1993e770dbb98ffd8a6ffc72e16b97cd6f981975b

The producer's additional regular-subset count is independently checked:
there are 21 ordinary regular subsets of each size 2,3,4,5, totaling84;
there are no signed regular subsets among all128 subsets. Its sharp
transitivity statements follow from the established group orders and
affine formulas. The distinction between global converse and complement
of nonpath arcs is correct; no metagraph-wiring theorem is being claimed.

The producer cites a further primary source, An–Wilson, *The Alperin
weight conjecture and Uno's conjecture for the Monster M, p odd*, LMS
J. Comput. Math. 13 (2010), 320–356, Proposition4.1, proof Case1 and
Table1 on printed p.324. This exact locator was independently verified:
both C7 × He and (7:3 × He):2 occur there. The year is2010, even though
the DOI contains2009. https://doi.org/10.1112/S1461157009000059

The cited local repository inheritance paths exist. This audit validates
the specific Paley symmetry/count claims needed here, not all results
in those larger inherited theorem files. The producer's CFT wording
correctly distinguishes the established parafermion embedding from the
expected Held action and evidence for a7A lift. The forbidden-H and
full-Monster-action boundaries remain intact.
