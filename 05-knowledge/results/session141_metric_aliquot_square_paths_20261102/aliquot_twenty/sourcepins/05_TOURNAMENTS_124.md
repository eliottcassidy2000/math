# The 1–2–4 bridge: weighted divisor downsets and tournament cycle complexes

**Status: PROVED elementary carrier theorem and marked finite-field identities;
CITED inherited Hamiltonian-spectrum results; FINITE-EXACT independent controls.**
No external novelty claim is made for the elementary translations. The useful
connection is a shared weighted counting language with explicit realizability
conditions. It is not an identification of aliquot untouchables with missing
Hamiltonian-path counts.

## 1. Inheritance and the exact question

The closest proved mechanism is the odd-cycle formula
`H(T)=I(Omega(T),2)`, retained in
[THM-002 / ocf](../../../../01-canon/theorems/THM-002-ocf.md).
Here `H` counts Hamiltonian paths as vertex words. `Omega(T)` has one vertex
for each literal directed odd cycle of length at least three, modulo rotation
only; two cycles conflict when their original vertex sets intersect. Distinct
cycles on the same support must remain distinct. The primary formulation is
[Irving–Omar, Corollary 20](https://arxiv.org/html/2412.10572v3), which credits
Grinberg–Stanley. Its complement formulation is equivalent here: reversing
each vertex word identifies the path counts of a tournament and its converse.
The canonical hostile is `I(K3,2)=7`: a perfectly valid graph
count whose graph cannot be the whole `Omega(T)` of a tournament.

The prior [Camion–Busch / polyhedra / Collatz note](../../camion_busch_gaps_polyhedra_collatz_20261002.md)
already explored `7=1+2+4`, `21=1+4+16`, their projective-plane counts, and the
Paley/Collatz occurrence of `{1,2,4}`. The current
[S57 aliquot packet, section 4](../../s57_untouchables_root27_20261026/README.md)
already supplies the weighted divisor and projective-boundary interpretations.
Those observations are inherited here, not rediscovered results.

The corrected near miss is treating three notions as one: ordinary path count
`H`, path-orbit count `tc=H/|Aut|`, and a numerical polynomial evaluation after
forgetting its coefficients. The least-used sidecar is the monomial support:
which exponents or cycle collections are actually permitted. The concept board
is weighted downsets / squarefree faces / flagness / tournament realization /
automorphism quotients / Frobenius markings. The anchor is the exact counting
interface, the niche is its semiprime boundary, and the wildcard is the marked
1–2–4 Frobenius orbit.

### Current 7/21 status

| Object | Current statement | Scope |
|---|---|---|
| `H(T)` | Neither 7 nor 21 occurs at any order | **CITED**, [THM-1370 / h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) |
| Full ordinary H spectrum | Whether every other positive odd integer occurs | **OPEN**, [current spectrum note](../../rotation_owners_spectrum_20261020e/spectrum.md) and [THM-4094 / matching-deficit-and-two-prime-lane](../../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md) |
| `tc(T)` | Neither 7 nor 21 occurs through order 11 | **CITED FINITE-EXACT plus structural reduction**, [THM-4790 / presentation-fibres-7-and-21-absent-through-eleven-vertices](../../../../01-canon/theorems/THM-4790-presentation-fibres-7-and-21-absent-through-eleven-vertices.md) |
| `tc(T)` at arbitrary order | Exclusion of 7 and 21 for non-rigid classes | **OPEN**, [HYP-9352(2)](../../../hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md) |
| Aliquot sum `s(n)` | `s(8)=7` and `s(18)=21` | **PROVED by direct divisor lists** |

The [Session66 tc-scope appendix](../../session66_root27_aliquot_graph_interfaces_20261026/session66-tc-scope/REPORT.md)
correctly separated H from tc at its pinned source. Its through-nine / reported
order-ten verification boundary predates the independently audited THM-4790;
it should not replace that later through-eleven result.

The ordinary H omission argument uses multiplication across strong components
and the strong minima `3,5,9,15,25,45,75,125,...`. The exact residue-class floor
and recurrence `f(n+3)=5f(n)` are classical
[Busch, Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).
Together with the small strong spectra they leave no strong factor 7 or 21.
This packet cites that established global proof; its new finite checks do not
purport to reprove it. The tc argument needs additional automorphism data:
Paley7 has `(H,|Aut|,tc)=(189,21,9)`. Excluding H=7,21 does not exclude the
larger products `7|Aut|` or `21|Aut|`.

## 2. A common carrier, with a sharp flag boundary

Let `n>=2` have distinct prime factors `p_1,...,p_r` and exponent vector
`a=(a_1,...,a_r)`. Its proper-divisor exponent set is

    D(n) = {e: 0<=e_i<=a_i for all i} minus {a}.

This is a multicomplex: a set of nonnegative exponent vectors closed downward.
Its weighted enumerator is

    A_n(X_1,...,X_r) = sum_(e in D(n)) product_i X_i^(e_i)
                     = product_i(1+X_i+...+X_i^(a_i)) - product_i X_i^(a_i).

The divisor bijection `e -> product_i p_i^(e_i)` proves

    A_n(p_1,...,p_r)=s(n),       A_n(1,...,1)=tau(n)-1.       (1)

For the graph comparison, a vertex is a degree-one monomial **actually in
D(n)**. Coordinates whose degree-one monomial is absent are discarded as
ghost coordinates. A graph realization below means literal equality of
monomial supports on these retained coordinates, not merely equality of an
evaluation, and not a new bijection that treats every divisor as a vertex.

**Theorem (prime-factor multiplicity boundary).** `D(n)` is the independence
complex of a simple graph on these active coordinates if and only if

    Omega_arith(n) = sum_i a_i <= 2.                       (2)

**Proof.** For a prime n, D(n) consists only of the zero vector, the
independence complex of the graph with no vertices. For `n=p^2`, it has the
zero vector and one degree-one vector, hence is the complex of one isolated
vertex. For `n=pq` with distinct primes, it has the zero vector and two
degree-one vectors but no pair, hence is the independence complex of K2.

Conversely, suppose total multiplicity is at least three. If some `a_i>=2`,
then `2e_i` belongs to D(n): it cannot be the deleted top vector unless
`n=p_i^2`, which has multiplicity two. Thus D(n) contains a nonsquarefree
monomial and cannot be a simplicial complex of subsets. If all `a_i=1`,
there are at least three active coordinates and D(n) contains every proper
subset. Every pair is allowed, so any graph realizing it would have no
edges and would allow the full set as well. That full set was deleted.
Equivalently, this boundary simplex is not flag. These exhaust the cases. QED.

For `n=1` the proper-divisor set is empty and s(1)=0. This void boundary is
kept separate: every graph independence polynomial contains the empty-set
term one.

The surviving exact weighted graph maps are therefore the familiar easy
aliquot sources:

    n=p       -> empty graph, value 1;
    n=p^2     -> one vertex of weight p, value 1+p;
    n=pq      -> K2 with weights p,q, value 1+p+q.           (3)

These formulas are inherited elementary identities. The theorem identifies
exactly where this particular graph-on-prime-coordinates translation stops.
A different encoding with more vertices or higher-order constraints is not
excluded. The obstruction is the omitted coordinate: multiplicity for a
repeated prime, and a higher-order forbidden face for three or more distinct
primes.

## 3. Why 1+2+4 can be arithmetic without being a tournament count

For 8 the multicomplex is `{1,X,X^2}` and

    s(8)=1+2+4=7.

The coefficients `(1,1,1)` cannot be a simplicial face profile: one face of
size two requires two different singleton faces. Two independent cycle
vertices instead have polynomial `1+2X+X^2`, with value nine at X=2.

Evaluation alone forgets this obstruction. The different graph polynomial
`1+3X` also evaluates to seven, and it is the independence polynomial of K3.
Thus there are two separate gates:

1. Arithmetic multicomplex -> squarefree flag complex: the literal prime
   support may fail, as it does for 8.
2. Arbitrary graph -> tournament cycle-conflict graph: K3 passes the first
   gate but fails the second by the inherited global H theorem.

For twenty-one the arbitrary graph P4 passes the graph gate:

    I(P4,X)=1+4X+3X^2,       I(P4,2)=21.

It likewise cannot be the **whole** cycle-conflict graph of a tournament.
This says nothing about its occurrence as an induced subgraph. On the
arithmetic side, `18` has proper divisors `1,2,3,6,9`, summing to 21.
Its multivariate support is `{1,X,Y,XY,Y^2}`; evaluating both variables at
two would give thirteen instead. The prime weights are essential.

The low coefficient constraints can be checked without a tournament census.
Three independent vertices already contribute `(1+2)^3=27`. Hence an
independence polynomial evaluating to 7 or 21 has no coefficient above degree
two. For seven, `1+2a+4b=7` and `b<=binom(a,2)` force `(a,b)=(3,0)`.
For twenty-one the possible pairs are `(4,3),(6,2),(8,1),(10,0)`.
Tournament realizability, rather than these arithmetic equations, eliminates
the remaining possibilities.

The two image gaps 2 and 5 of the aliquot map do not map automatically to
these H gaps. Conversely, every positive odd number is `I(K_m,2)=1+2m` for
some unrestricted graph. Omissions arise only after a carrier's realizability
conditions are imposed. This is the common research pattern; no reduction
or equivalence between the two image-membership problems is supplied.

## 4. Field-size grouping retains an inner exponent coordinate

The point counts `7=1+2+4` and `21=1+4+16` are those of projective planes
over F2 and F4. This counts nonzero triples modulo scalar multiplication;
it does not itself impose a tournament or an aliquot preimage.
For a prime p and positive integer r, splitting a divisor exponent uniquely
as `br+c`, with `b=0,1,2` and `0<=c<r`, gives the exact identity

    s(p^(3r)) = (1+p+...+p^(r-1)) (1+p^r+p^(2r)).          (4)

Thus the F_(p^r) plane count is the outer factor. The inner residue c is a
real sidecar, with total weight `1+p+...+p^(r-1)`. At p=2,r=2,

    s(64)=63=(1+2)(1+4+16)=3*21.

Discarding the odd-exponent divisor layer would incorrectly replace 63 by
21. Conversely, the ordinary tournament H spectrum is not closed under
division: H=63 occurs, although H=21 and H=7 do not. The checker includes a
literal eight-vertex H=63 tournament whose Omega is K31, independently
confirming the inherited example in
[THM-344 / n8-h63-complete-omega](../../../../01-canon/theorems/THM-344-n8-h63-complete-omega.md).
It reaches 63 as `1+2*31`, not by realizing a forbidden factor seven.

## 5. A marked 1–2–4 orbit links Collatz, Paley7 and F8

Use the **ordinary** Collatz map `C(n)=n/2` for even n and `3n+1` for odd n.
On its three-point cycle only,

    1 -> 4 -> 2 -> 1,          C(x)=4x mod7.

Multiplication by four is an order-three automorphism of Paley7, whose
quadratic-residue difference set is `{1,2,4}`. The induced Paley triangle
on these three vertices has arcs `1->2->4->1`, opposite to the displayed
Collatz traversal. This exact restriction does not define a quotient of
Collatz on all residues: `2==9 mod7`, but `C(2)=1` and `C(9)=28==0 mod7`.
Parity is missing from the residue-only state.

To compare with the seven nonzero vectors of F2^3, choose

    F8 = F2[a]/(a^3+a+1),
    r in Z/7  -> a^r,      r=0,...,6,
    bit-coordinate labels = (1,2,4,3,6,7,5).

This is an explicit multiplicative marking, not an additive identification
of Z/7 with F2^3. Residue doubling becomes Frobenius `x->x^2`; multiplication
by four becomes its square. The marked residue orbit `{1,2,4}` becomes
`{a,a^2,a^4}={2,4,6}`. Their sum is zero, so these are the three nonzero
vectors in a two-dimensional F2 subspace: one Fano line.

The marking and rank distinction both matter. Replacing the generator a by
`b=a^3` sends the same residue labels to `{3,5,7}`. Their sum is one rather
than zero, and they form a basis of F8 over F2. This second Frobenius orbit
is a **normal basis**; the first, lying on a Fano line and having rank two,
is not. Both orbits have length three. Therefore an orbit-length statement
does not permit transferring a normal-basis necklace dictionary to the
Fano-line orbit. The checker verifies the two ranks, the chosen power table,
and the failed global Collatz residue quotient.

## 6. Finite controls and one recovered wording correction

Run from this directory, using only Python's standard library:

```sh
python3 -B check_tournaments_124.py --output normal.replay.json
python3 -B -O check_tournaments_124.py --output optimized.replay.json
cmp normal.replay.json optimized.replay.json
cmp normal.json normal.replay.json
```

Both frozen outputs are byte-identical and pass **8,494 explicit gates**.
The universe is:

- Every integer `2<=n<=1000`: direct proper-divisor sums versus weighted
  exponent sets, the product formula, the F1 count, and full reconstruction
  of the candidate graph independence complex. Exactly 467 satisfy the flag
  criterion and 532 fail. The all-order statement rests on the proof above.
- All 1,099 labelled tournaments on one through five vertices: literal
  Hamiltonian permutations versus subset DP; all literal directed odd cycles
  and the conflict-graph independence polynomial; full permutation
  automorphisms, divisibility, and the finite H/tc gap controls.
- Ordered stacks of one, two and three triangles: OCF profiles `(1+X)^k`
  and counts `3^k`, testing genuine higher-degree cycle collections.
- Literal witnesses H=31 at order six, H=63 at order eight, H=127 at order
  seven, and Paley7. Each has independent word enumeration, DP, OCF and full
  automorphism count. Paley7 additionally rejects collapsing cycles to their
  vertex supports. The Mersenne controls reject a permanent Mersenne-gap law.
- Exact coefficient, graph, exponent-grouping, Collatz-cycle and marked-F8
  positive and hostile controls described above.

A source wording issue was found during the H=127 control. The finite-record
clause of [THM-4748 / two-reversed-arcs-and-the-mersenne-side](../../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md)
said that the classes with H=31,63,127,255 have least feedback-arc number three.
Read as a statement about every class, this is false. The explicit row-mask
tournament

    (106,60,89,48,97,68,10)

has H=127, automorphism order one, OCF profile `(1,45,9)`, and minimum
feedback-arc number **five**, checked over all 7! vertex orders. What survives
is the theorem's minimum over **all realizations of a specified Mersenne
value**: three, attained by its fan construction at a possibly larger order.
The reviewer notified the integrating root; this packet does not silently
modify the inherited proof. The explicit distinction is between a minimum
over a statistic's fibre and a uniform assertion about every member of it.

The next mathematical targets retain their original types: all-order tc
exclusion remains open; ordinary H completeness still needs new strong
realizations; and a broader divisor-to-cycle correspondence would need both
the multiplicity/higher-face sidecar and a tournament-realization theorem.
