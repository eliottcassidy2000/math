# From runner circles to a four-dimensional symplectic spectrum constructor

**PROVED constructions and reduction; FINITE-EXACT certificates. Independently reviewed.**
Session 192, 2026-11-05. The contribution is a target-preserving geometric lift
and an infinite fixed-dimensional construction class. It does not settle the
full `{7,21}` spectrum, establish two-apex completeness, or prove anything new
about Lonely Runner. No literature-priority claim is made.

## 1. Inheritance and the connection contract

The relevant, initially unrelated source is the actual half-turn comparator
used in the repository's Lonely Runner geometry. Points on a circle orient
`u -> v` when the directed angular difference lies in `(0,1/2)` turns. With
no coincident or antipodal ties, this is the sign of a two-dimensional
alternating determinant. Its output is necessarily locally transitive.

Closest proved sources, cited by ID and filename:

- THM-402, `round-tournaments-are-2-dichromatic`: the circular comparator
  has transitive in- and out-neighborhoods. Its chromatic statistic alone
  does not supply arbitrary Hamiltonian counts.
- THM-468, `tournament-determinant-floor`: local orders are exactly
  cut-switches of transitive tournaments and the determinant-floor class.
  The underlying classification is credited there to earlier literature.
- THM-4602, `tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips`:
  the exact rank-two / local-transitivity correspondence, with its prior-art
  qualification and its warning about transferring only sign information.
- Session 191, [two-apex grammar](../../session191_tournament_growth_20261104/two_apex/TWO_APEX.md), the four-letter physical grammar.
  Its 499 witnesses through 1001 were independently checked there. This
  packet checks them again through the **new geometric relation**, rather
  than importing the old transfer state as a path counter.

The freshest inherited Session 191 route is pinned at
<https://github.com/eliottcassidy2000/math/blob/795788ef710366dcfc0b8cd84fb64644f164c9e5/00-navigation/S191.md>.
The generated witness bank is included as `inherited_prefix_witnesses.json`; the earlier guide links the generating scan but did not include this JSON file. Its entries are independently checked here. No incoming unreviewed theorem is used.

The source-to-target map is explicit:

    circular determinant signs -> sum of two planar determinants
    -> arbitrary two-apex word tournament -> ordered products.

It preserves every labelled arc and therefore every Hamilton path, not just
one count, residue, score, or determinant. The circle's actual time parameter,
runner speeds, distance-to-integers inequalities, and ownership data do **not**
transfer. No rank-four configuration is asserted to be a physical LRC state.
The second symplectic plane is the sidecar that restores arbitrary two-apex
incidences. The cheapest hostile is H=13, a genuinely attainable target which
no rank-two configuration can produce at any order.

Live board: circle signs; forbidden local vortices; two-apex letters;
integer sign margins; affine charts and joins; target factorization.
The productive change is dimension and labelled relation, not replacing the
spectrum question by another scalar resemblance.

## 2. A fully explicit rank-four compiler

Use coordinates **`(q1,p1,q2,p2)`**, and the alternating form

    omega(u,v) = u1*v2-u2*v1 + u3*v4-u4*v3.

Its matrix is `J2 direct-sum J2`, with `J2=[[0,1],[-1,0]]`.
This is not the block-coordinate convention `(q1,q2,p1,p2)`.
Orient `u -> v` iff `omega(u,v)>0`.

Let `w=c_1...c_m` be any word over `{0,1,2,3}`. The two-apex tournament
has names `a,b,x_1,...,x_m`, the arc `a -> b`, the chain `x_i -> x_j`
for `i<j`, and

    bit 0 of c_i = [a -> x_i],
    bit 1 of c_i = [b -> x_i].

Set

    K = 3m+1,
    a = (1,0,K,0),             b = (1,0,0,K),
    x_i = (1,3i,s_i,t_i),
    s_i = -1 if bit 1 is set, and +1 otherwise,
    t_i = +1 if bit 0 is set, and -1 otherwise.

The four sign pairs `(s,t)` are respectively

    0: (+1,-1),  1: (+1,+1),  2: (-1,-1),  3: (-1,+1).

**Theorem 1.** These integer vectors realize exactly the word tournament.
Every pair has absolute symplectic pairing at least one, and the full
alternating pairing matrix has rank at most four.

**Proof.** Every type of pair has an explicit value:

    omega(a,b)   = K^2 > 0,
    omega(a,x_i) = 3i + K*t_i,
    omega(b,x_i) = 3i - K*s_i,
    omega(x_i,x_j) = 3(j-i)+s_i*t_j-t_i*s_j.

Since `K>3i`, the apex signs are `t_i` and `-s_i` respectively.
For `i<j` the last expression is at least `3-2=1`. There are no ties,
and integrality gives the unit margin. If V has the vectors as columns,
the pairing matrix is `V^T J V`, of rank at most four. QED.

This is an all-word theorem with no length, run-count, or clock restriction.
All points lie in the same affine chart `q1=1`; this additional coordinate
property is essential for the closure in Section 4. The construction changes
K when m changes. It does not claim that one fixed numerical configuration
supports every infinite continuation with all old vector coordinates frozen.
The **old tournament relation**, not its particular vector coordinates, is
what an induced word extension retains.

## 3. The first exact rank boundary occurs at H=13

Here H counts directed Hamiltonian vertex words, with H(empty)=1.
Call an alternating realization rank at most r if its real alternating
pairing matrix has rank at most r, has no off-diagonal zero, and has the
tournament's arc signs. Every rank-two alternating matrix factors as planar
determinants, and every such determinant representation has rank at most two.

**Theorem 2.** Among all tournaments with H=13, the least possible rank of
an alternating sign realization is exactly four.

**Lower bound.** For a planar determinant tournament, the vectors beaten by
one fixed vector all lie in one open half-plane. Their angular order is
transitive; the same holds for the vectors beating it. Thus every rank-two
realization is locally transitive, directly, without a numerical rank test.

A nontransitive locally transitive tournament is strong. Otherwise a directed
triangle lies in one strong component and a vertex in any other component
either dominates or is dominated by that triangle, contradicting local
transitivity. A tournament with H=13 is nontransitive, since a transitive
tournament has H=1. Therefore a hypothetical rank-two H13 tournament is strong.

Busch's cited strong-tournament bound `H >= 5^((n-1)/3)` implies `n<=5`,
because `5^5=3125>13^3=2197`. The bound and primary DOI
<https://doi.org/10.37236/1141> are inherited from the audited Session 190/191
packets; no sharper residue-dependent floor is assumed here.

Complete independent enumeration of all labelled tournaments at orders 1..5
in `verify.py` gives the following local-order supports:

    n=1: {1};       n=2: {1};       n=3: {1,3};
    n=4: {1,5};     n=5: {1,9,11,15}.

The local-order counts are `1,2,8,48,384`. Each retained H is independently
checked by literal permutations and subset DP. Hence none has H=13. This is
an all-order conclusion from a stated size cutoff and a complete finite gate,
not an extrapolation from a bounded random search. Alternating ranks are even,
so the lower bound is four.

**Upper bound and named certificate.** Apply Theorem 1 to `w=120`. The five
vectors are

    a  = (1,0,10,0),       b  = (1,0,0,10),
    x1 = (1,3,1,1),        x2 = (1,6,-1,-1),
    x3 = (1,9,1,-1).

Their outgoing bit rows are `(6,8,26,17,3)` and their Hamilton count is 13.
The pairings among `a,b,x1,x2` have Pfaffian

    100*3 - 13*16 + (-4)*(-7) = 120 != 0,

so this actual pairing matrix has rank four as well. This finite relation is
checked by both path-count algorithms. QED.

The obstruction is at the **target value**, not merely at this particular
presentation: no other order, marking, angular configuration, or rank-two
realization can rescue H13. Conversely, a local vortex in one H=9 presentation
does not show that H=9 requires rank four, since a different round H9 witness
exists. The proof deliberately keeps those quantifiers separate.

This does not say that ordinary matrix rank of the +/-1 adjacency matrix is
four, or that the signs are a rank-four oriented-matroid chirotope. The object
is the sign pattern of an alternating form, specifically defined above.

## 4. Ordered products stay in the same four-dimensional affine chart

**Theorem 3.** Suppose each of finitely many tournaments is represented by
integer vectors `(1,p,q,r)` in this chart. Their ordered join has another
integer representation in the same symplectic dimension and with unit margin.

Choose common bounds `|p|<=B`, `|q|,|r|<=C` for all blocks, and set

    L = 2B+2C^2+1.

For block j, replace every point `(1,p,q,r)` by `(1,p+jL,q,r)`.
Internal pairings are unchanged. Between a point in block i and one in block
j>i, the pairing is at least

    (j-i)L - 2B - 2C^2 >= 1.

Thus every earlier block dominates every later block. The Hamilton count of
the ordered join is the product of the block counts, because every Hamilton
word must traverse the blocks in that order, and arbitrary Hamilton words
inside the blocks concatenate uniquely. This proves both claims.

The common affine chart is load-bearing: arbitrary rank-four realizations
without a shared positive linear functional are not asserted to satisfy this
join construction. Theorem 1 supplies exactly the needed chart. Projecting
its configuration onto the first plane destroys it: a and b both project to
`(1,0)`, producing a tie. Discarding the second plane is not a target-preserving
quotient.

## 5. An infinite explicit rank-four target family

**Theorem 4.** Every positive odd integer whose prime factors are at most
139, other than 7 and 21, is the Hamiltonian-path count of a tournament with
an integral rank-at-most-four alternating realization. The construction
returns every vertex vector and has positive integer margins.

The finite input is the following generator bank:

1. Every odd prime p<=139 other than 7.
2. `49,63,343`.
3. Every `7p` with p prime, `5<=p<=139`, p!=7.

Every generator is at most 973 and appears in the inherited word-witness
prefix through 1001. This packet independently builds its geometry and
counts its paths; the actual words are in `receipt.json`.

To construct N, let e be its exponent of 7.

- If e=0, use the ordinary prime factors as generators.
- If e>=2, express e as `2u+3v`, with v=0 or1, and use u copies of49 and v
  copies of343. Use the remaining ordinary prime factors separately.
- If e=1 and another prime p>=5 occurs, use one generator7p and the remaining
  ordinary prime factors.
- If e=1 and only powers of3 remain, write `N=7*3^a`. The excluded cases are
  a=0,1. For a>=2, use63 and a-2 copies of3.
- For N=1 use one vertex.

Apply Theorem 1 to every chosen generator word and Theorem 3 to their ordered
join. This proves the exact target count. It is not merely modular coverage.

**Size and height.** The inherited full bank through1001 has chain length at
most8, so every block has at most10 vertices, B<=24, and C<=25. One may use
the uniform translation step1299. With r generator factors the construction
has at most10r vertices and coordinate height at most
`max(25,1299(r-1)+24)`. Since every generator is at least3,
`r<=log_3 N` for N>1. The dimension is fixed at four while both order and
integer coordinate height grow only logarithmically with N. Factoring in this
fixed smoothness class requires only division by the finite prime bank.

The original integer counts follow from existing two-apex witnesses and the
standard SCC product. The strengthened conclusion is that **all these targets
fit one fixed four-dimensional, unit-margin geometry**, with an explicit
compiler and closure theorem. No previously unknown global integer H is
claimed. The finite prefix itself, even the larger inherited scan through
230519, is not a tail theorem.

## 6. Reproduction, controls, and next precise question

Run:

    python3 verify.py > normal.stdout
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout

The program uses no third-party packages or optimization-removable assertions.
It independently checks all499 geometric witnesses through1001, every word
of length at most5 for its exact arc signs, every labelled tournament of order
at most5 for the rank-two finite gate, the explicit H13 Pfaffian, every
139-smooth odd input through100001 for the factorization algorithm, and named
joined examples. Joined examples of order at most14 receive an independent
full subset-DP count; the much larger example uses the proved join law with
independently verified block counts and checks **every physical cross arc**.
The receipt distinguishes these verification modes.

This opens a specific new construction lane: search the affine symplectic
rank-four chamber space for target-preserving patches **beyond the two-apex
word subfamily**, rather than forcing circular comparator moves to realize
H13 or insisting on fixed-end word growth. The current construction has not
proved a patch calculus, recursive ternary closure, or spectrum completeness.
The rank-two obstruction tells us precisely why the old planar lens is too
small; the explicit lift identifies a larger real carrier that retains every
arc and already supports an infinite constructive target class.
