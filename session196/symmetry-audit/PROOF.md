# Independent audit: prime Hamiltonian counts and the first symmetric 7p value

**Status:** PROVED consequences of the established odd-cycle collection formula,
with explicitly finite-exact certificates for the small orbit types below.
**Scope:** All finite tournaments; `H(T)` counts directed Hamiltonian vertex
words, and `Aut(T)` means orientation-preserving automorphisms.
**Priority:** No claim of literature novelty is made.

## 1. Minimal imported dependency and basic lemmas

The only imported counting identity is the existing proved OCF in
`01-canon/theorems/THM-002-ocf.md`: `H(T)` is the sum of `2^k` over collections
of `k` vertex-disjoint nontrivial odd directed cycles, including the empty
collection. In particular,

- `H(T)` is positive and odd;
- `H(T) >= 1 + 2 c3(T)`, where `c3` counts cyclic triples;
- for disjoint induced subtournaments `A,B`, `H(T) >= H(A) H(B)`.

The last claim follows by combining any odd-cycle collection in `A` with any
one in `B`. It does not require either subset to be a module or component.
There is also a completely independent path proof: merge one Hamiltonian
word of `A` and one of `B` by repeatedly outputting the head that beats the
other head, then append the remainder. The result is a Hamiltonian word.
Deletion of the other alphabet recovers each input word, so this is an
injection. Every consecutive output pair is either consecutive in one input
word, or was compared as the two heads when the first was output.

Induced monotonicity can be proved without OCF. Given a Hamiltonian word of
`T-v`, insert `v` just before the first vertex that it beats, or at the end
if there is none. The preceding vertex, if any, beats `v`; the new word is
a Hamiltonian word. Deleting `v` recovers the input.

The strong components of a tournament are linearly ordered, all cross-arcs
point forward, and every Hamiltonian word traverses components in that order.
Thus `H(T)` is the product of their Hamiltonian counts. Automorphisms preserve
each strong component, since the component order has no nontrivial symmetry.

## 2. The prime-orbit lower bound

**Theorem.** If a prime `q >= 5` divides `|Aut(T)|`, then

`H(T) >= q(q^2 + 11)/12`.

**Proof.** By Cauchy's theorem choose an automorphism `g` of order `q`. At least
one of its vertex orbits `O` has size `q`. The restriction of `g` is a transitive
cyclic automorphism of the induced tournament `U=T[O]`. Thus every vertex of
`U` has outdegree `(q-1)/2`.

For every tournament on `m` vertices,

`c3 = binom(m,3) - sum_v binom(d_plus(v),2)`,

because a transitive triple has a unique source. Substituting the common
degree gives `c3(U)=q(q^2-1)/24`.

The automorphism group acts freely on Hamiltonian vertex words: fixing one
word fixes every vertex. Therefore `|Aut(U)|` divides `H(U)`, and in particular
`q` divides it. OCF gives

`H(U) >= 1 + q(q^2-1)/12`.

For primes `q>=5`, `a=(q^2-1)/12` is an integer. Rounding the strict lower bound
`qa+1` to the next multiple of `q` gives `H(U)>=q(a+1)`. Induced monotonicity
finishes the proof. Notice that only one orbit of an order-`q` element was
used; no claim about a transitive action on all vertices is needed. QED.

**Corollary (prime-H rigidity).** If `H(T)=p` is prime and `p>=5`, then `T` is
rigid. Indeed `|Aut(T)|` divides `p` by the same free action, so any nontrivial
group has order `p`. The lower bound above is strictly larger than `p`.

**Boundary `p=3`.** Every tournament with `H=3` has exactly one cyclic triple:
OCF bounds the count by one, while a tournament with no cyclic triple is
transitive and has `H=1`. Write its vertices as `C`. Every outside vertex
is uniform toward `C`, since a mixed relation with a directed triangle
creates another cyclic triple. Let `A` dominate `C`, and let `B` be dominated
by `C`. Both `A,B` are transitive; moreover `A -> B`, since a contrary arc
would create a triangle through a vertex of `C`. Consequently

`T = transitive(A) -> C3 -> transitive(B)`.

Conversely every such tournament has `H=3` and automorphism group `C3`.
There are exactly `n-2` unlabelled examples at order `n>=3`, and
`n!(n-2)/3` labelled examples. Thus the exception is structural, not merely
one exceptional graph. There are no tournaments with `H=2`.

**Enumeration consequence.** For each prime `p>=5`, the number of labelled
order-`n` tournaments with `H=p` is exactly `n!` times the number of unlabelled
ones. This is a genuine symmetry-stratum statement about A000568, not an
identification of A000568 with Hamiltonian counts.

## 3. Structural lemmas used for the small boundaries

A module is a vertex set toward which every outside vertex is uniform. If an
order-`q` automorphism has exactly one nontrivial orbit, that orbit is a module:
all other vertices are fixed and must be uniform toward the orbit.

Every vertex of a nontrivial strong tournament lies on a directed triangle.
For a vertex `v`, write `A` for its out-neighbors and `B` for its in-neighbors.
If no arc points from `A` to `B`, then vertices in `A` cannot reach `v` or `B`,
contradicting strong connectivity.

**Rooted strong-four lemma.** In a strong tournament of order at least four,
every prescribed vertex `r` belongs to an induced strong four-vertex tournament.
Take a triangle `C` through `r`. If an outside vertex has mixed relations to
`C`, adjoining it works. Otherwise `C` is a module. In the strong quotient
obtained by contracting `C`, a triangle through the contracted vertex yields
two outside vertices `a,b` with `C -> a -> b -> C`. Any two vertices of `C`,
including `r`, together with `a,b` induce a strong four-vertex tournament.

**Strong-five containment.** Every strong tournament of order at least five
contains an induced strong five-vertex tournament. Start with a strong four
from the preceding lemma. A mixed outside vertex extends it. If none exists,
that four-set is a module; a two-vertex cyclic ear exists as above. A triangle
inside the four-set together with that ear is a strong five-set.

**Consequence `H != 7`.** The finite certificate checks all 1,024 labelled
five-vertex tournaments: the 544 strong ones have counts `9,11,13,15`, with
multiplicities `240,120,120,64`. Every strong tournament of order at least five
therefore has `H>=9`; the smaller nontrivial strong possibilities have counts
`3` and `5`. SCC multiplicativity then excludes `H=7` at every order.
This does not assume any all-order exclusion of 21.

## 4. Precisely bounded finite certificates

`verify_symmetry_threshold.py` uses its own direct adjacency encoding and
subset dynamic program, importing only primitives from the audit's standalone
`verify_prime_rigidity.py`. It checks every strong five-vertex case and every
circulant of orders `3,5,7` against literal permutation counting. All rooted
strong-four substitutions are computed independently by subset DP.

For an order-`q` orbit, label its vertices by powers of the cyclic generator.
The induced tournament has one chosen step from each pair `{s,-s}`. Thus the
following are complete finite universes, not samples:

| Orbit order | Connection sets | Orbit H | C3 quotient inflation H | Rooted strong4 quotient inflation H |
| --- | ---: | --- | --- | --- |
| 3 | 2 | 3 | 15 | 27,33 |
| 5 | 4 | 15 | 135 | 255,345 |
| 7 | 8 | 175 (six),189 (two) | 2163,2205 | 4151,4221,5789,5859 |

The 32 order-eleven connection sets have minimum `H=92411`, checked by subset
DP. For the two possible orientations of a second cyclic three-orbit and all
eight synchronized cross-arc patterns, the induced six-vertex tournament has
`H=9` for either uniform cross-pattern and `H=45` for every mixed pattern.
All 16 of those cases are also checked by literal permutations.

## 5. Stronger consequences and hostile small cases

**All prime-square strata except 9 are rigid.** If `H=p^2` with prime `p>=5`,
nontrivial symmetry would have a prime divisor `p`. For `p>=13`, the analytic
lower bound exceeds `p^2`. For `p=7,11`, the exact orbit minima exceed `p^2`.
For `p=5`, two five-orbits already cost `15^2=225`; otherwise the unique
five-orbit is a module. If it is its own strong component, `15` divides `H`.
If its strong component is larger, the contracted strong quotient contains a
triangle through it, and induced monotonicity gives `H>=135`. Both options
exclude `H=25`.

The identical argument proves that every tournament with `H=35` is rigid:
a prime symmetry divisor would be `5` or `7`; a seven-orbit costs at least
175, and the preceding five-orbit argument excludes 35. Nothing about weak
versus strong ambient tournaments has been suppressed: divisibility applies
when the orbit is its own component; monotonicity applies otherwise.

**Every hypothetical H=21 tournament is rigid.** The prime-orbit bound excludes
symmetry of order seven. Suppose an order-three automorphism exists. OCF
supermultiplicativity implies there are at most two nontrivial three-orbits,
since three disjoint cyclic triples cost `3^3=27`.

With one orbit, it is a `C3` module. With two orbits, the complete six-vertex
certificate gives 45 unless their cross-arcs are uniform, so their cross-arcs
must be uniform. Each orbit is then a module relative to the other orbit and
to every fixed vertex.

Consider a strong component containing such a module, contracting any of the
three-orbits it contains. If the quotient has at least four vertices, the
rooted strong-four lemma and the certified C3 inflation minimum give `H>=27`
(select one representative from any other orbit in that four-set). If the
quotient has three vertices and contains one C3 module, its count is 15,
which cannot divide 21. If it contains both C3 modules, the three blockwise
Hamiltonian words of the quotient already give at least `3*3*3=27`.
The remaining possibility is that each three-orbit is its own strong component.
Two such components force `9|21`, impossible. One forces the product of all
other component counts to equal 7, impossible by the all-order H7 argument.

This establishes only a rigidity condition at 21, not the absence of every
21-count tournament. It therefore does not presume the target spectrum result.

## 6. Exact first nonrigid 7p value

**Theorem.** Among all finite tournaments whose Hamiltonian count is `7p` for
a prime `p`, the smallest count admitting nontrivial automorphisms is

`4151 = 7 * 593`.

**Lower bound.** Let `H=7p<4151`. Parity rules out `p=2`; the preceding section
handles `p=3,5,7`. For `p>=11`, the prime-orbit bound is greater than `7p` for
an orbit of order `p`. Since the automorphism-group order divides `7p`, any
nontrivial group therefore supplies an order-seven element.

Each seven-orbit has `H>=175`. Two such orbits would give `H>=175^2>4151`, so
there is exactly one, called `M`, and it is a module. Its own count is 175 or
189. Let `S` be its ambient strong component and `Q=S/M` its contracted quotient.

- If `Q` has one vertex, then `H(S)=175` or `189`, neither of which divides
  `7p` for prime `p`.
- A strong quotient cannot have two vertices.
- If `Q` has three vertices, `H(S)=2163=7*309` or `2205=7*315`, neither of
  which divides `7p` for prime `p`.
- If `Q` has at least four vertices, it contains a strong four through the
  contracted vertex. Inflating that vertex by `M` gives an induced
  subtournament with at least 4151 Hamiltonian words, a contradiction.

This proves the all-order lower bound.

**Sharpness.** The ten-vertex tournament with outgoing-neighbor bitmasks

`[240,225,195,135,142,156,184,512,255,383]`

is strong, has `H=4151` by exact subset DP, and has automorphism order exactly
seven by complete score-class-preserving permutation enumeration. Its first
seven vertices form a cyclic module. The prime factor 593 is verified by
trial division through 24. Hence the lower bound is attained.

## 7. Reproduction and scope ledger

Run from the audit directory:

```sh
python3 verify_prime_rigidity.py --max-n 6 --output prime_rigidity_labelled_n6.json
python3 verify_symmetry_threshold.py
```

The first command enumerates all 33,867 labelled tournaments of orders 1--6,
checks oddness and the triangle lower bound for all of them, and computes
complete automorphism groups for every prime-H stratum and for H=1 or21.
It verifies prime-H rigidity, the H3 boundary, and the exact finite spectra.

The all-order theorems are proved by the arguments above; the n<=6 audit is
only an independent hostile check. The finite circulant/substitution tables
are load-bearing finite-exact lemmas with fully specified universes. No
extrapolation from the order-six census and no unpublished global H21
nonattainability are used.
