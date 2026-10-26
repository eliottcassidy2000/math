# Marked Frobenius orbits, sharp symmetry loss, and the two Collatz sheets

2026-10-26. **PROVED** elementary statements A1-A6 below, with the existence
of normal bases **CITED** through THM-4804. **FINITE-EXACT** controls have
explicit bounds. The proposed global one-cycle/three-cycle classification
and a causal chess-tempo explanation remain **OPEN**, not consequences of
the modular identities. No theorem about convergence is claimed.

## 1. What is inherited, and what marking adds

**THM-4804**,
[module shattering and normal bases](../../../01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md),
already supplies the user's finite-field dictionary. With a chosen normal
basis `beta_i=alpha^(2^i)` of `K=F_(2^m)`, write

    x_w=sum_i w_i beta_i.

Frobenius rotates the binary coordinates. If `d` is the least period of
`w`, it is also `[F_2(x_w):F_2]`; the carousel tournament has

    Aut(T_w)=C_3^m semidirect C_(m/d),
    C_(m/d)=Gal(K/F_2(x_w)).                                      (1)

These statements classify a decoration orbit. They do not identify
Frobenius with the actual complement-line operation.

**A1 (PROVED, marked-port law).** A marked block position `i` corresponds
to the pair `(x_w,beta_i)`. This pair has Frobenius orbit length exactly
`m`, regardless of `d`; its stabilizer is trivial. Indeed the conjugates
of the normal element are distinct, so only the identity power of
Frobenius fixes `beta_i`. Equivalently `F_2(x_w,beta_i)=K`.

On the tournament side a marked block therefore kills the entire residual
cyclic factor in (1). A marked singleton inside that block leaves the
internal order `3^m`; a marked triangle vertex leaves `3^(m-1)`. These are
different rooted objects. Marking the whole ordered Hamilton path kills
**every** automorphism: a permutation preserving its ordered list fixes
all vertices.

There is also a rooted classification count. A fixed decoration orbit of
period `d` has `d` distinct block-rooted types: its product with the `m`
basis-vector positions has `dm` pairs, and every diagonal Frobenius orbit
has size `m`. Thus its ordinary tournament has `d` singleton-rooted and
`d` triangle-vertex-rooted types. The reciprocal-stabilizer check is

    d/3^m + d/3^(m-1) = 4d/3^m = |V(T_w)|/|Aut(T_w)|.

This is the `e=m` case of the incoming [aligned phase theorem, section 6](../session68_tournament_composition_classification_20261026/session68-arithmetic/REPORT.md):
periods `d,e` have `gcd(d,e)` joint rotation orbits, each of length
`lcm(d,e)`. Independently classifying the decoration and its port forgets
these relative phases. These are ordinary rooted classes; converse merging
would require its additional reflection action.

It follows that there is no automorphism-equivariant choice of one
Hamilton path for a tournament with nontrivial automorphism group. If
`s(T)` were such a choice, every automorphism would fix that path, a
contradiction. This concerns a choice on the unmarked object; it does not
forbid canonical labelled encodings or enumerating all path orbits.

## 2. The sharp 4p threshold is a different cost of symmetry breaking

The integrated [fibre theorem F2](../tournament_composition_classification_20261026/fibre/RESULT.md)
proves `v_p(|Aut(T)|)` invariant in a fixed labelled `C3` fibre whenever
`n<4p`, sharply for every odd prime `p`. Its proof retains the child-type
stabilizer. Only children of size at least four can change ordinary type
without changing their cyclic triples, and fewer than `p` such child slots
are fixed pointwise by every `p`-subgroup of the size-preserving quotient
group. Their changed colours cannot change its `p`-part.

**A2 (PROVED, realization in the decoration field).** At `m=p`, a constant
word gives `x_w in F_2`, with residual group `C_p`; a single nonzero bit
gives a normal-basis vector of degree `p`, with trivial residual group.
Changing `A+=C3->1` to `A-=1->C3` in that slot keeps the labelled cyclic
triples and removes one factor `p` from the group. Its vertex cost is `4p`.
Thus the first observer-preserving change in unmarked symmetry is realized
by a degree jump `1 -> p`.

The marked-port law A1 does not itself have a `4p` lower bound: marking a
vertex of `Q_p` already removes its translations at order `p`. A mark
changes the classified object; F2 compares unmarked tournaments inside
one observer fibre. Conflating these two operations would erase the
meaning of the sharp threshold.

For general odd `m`, the cyclic contribution has valuation
`v_p(m)-v_p(d)`. A degree change `d->d'` therefore changes this contribution
by `v_p(d)-v_p(d')`. The full path-response distribution from R1-R4 is
additional information: field degree determines the size normalization,
not which complement targets or attachment ports occur.

**Reflection boundary.** In normal coordinates the reflection
`R(beta_i)=beta_(-i)` satisfies `R F R=F^(-1)`. For odd `m>=3`, it cannot
be a field automorphism: all field automorphisms are powers of `F` and
commute, whereas `F!=F^(-1)`. It is an `F_2`-linear map. Also
`sum_i beta_i=1`, because the trace is a nonzero element of `F_2` (zero
would be a dependence among the basis vectors). Word complementation is
therefore `x->x+1`; the necklace converse is an affine reflected map,
not field multiplication. Normal-basis coordinates require the basis.

## 3. The multiplication-table centre is exact, but sign-symmetric

The precise prior authority is **THM-4792**,
[two sheets modulo odd m](../../../01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md),
which merges the October 5 central-square calculation and the October 25
modular edge-table classification. For `m=2k+1`, the central labels are
`k=-1/2`, `k+1=1/2`, and the four products form

    [ 1/4  -1/4 ]
    [-1/4   1/4 ]  modulo m.                                    (2)

Alternating the **unhalved** affine maps gives
`1/4 --(3x-1)--> -1/4 --(3x+1)--> 1/4`. This is an alternating two-map
cycle. It is not a cycle of either fixed-sign positive-integer Syracuse
map. THM-4792 also proves that negation can be represented by a fixed
halving-clock shift exactly when `-1` lies in `<2> mod m`; this fails,
for example, at `m=7` and `m=15`. The central identity (2) still holds there.

**A3 (PROVED, modular sign obstruction).** For every odd `m` and every fixed
clock `a>=0`, let

    f_(epsilon,a)(x)=(3x+epsilon)2^(-a) mod m,  epsilon=+1 or -1.

Negation `N(x)=-x` satisfies `N f_(+,a) N=f_(-,a)`. Consequently the two
finite functional graphs have identical cycle lengths and multiplicities,
even if `3|m` and the maps are not permutations. This is immediate on every
edge; no probabilistic or asymptotic argument is involved.

At clock zero, for example, **both signs** have cycle lengths `(1,4)` modulo
5, `(1,6)` modulo 7, `(1,5,5)` modulo 11 and `(1,3,3,3,3)` modulo 13.
All four tables have the same central-square mechanism. Thus that mechanism
alone cannot force different positive cycle counts for the two signs.
Sorting or relabelling the vertices also cannot change a directed graph's
cycle inventory. An edge-priority rule would be extra dynamics to specify.

## 4. The missing arithmetic port is the actual valuation clock

Define the positive odd maps

    U_epsilon(n)=(3n+epsilon)/2^v2(3n+epsilon).

On nonzero odd integers negation exchanges the signs, but exchanges the
positive and negative domains as well. This leaves room for different
positive-cycle behavior; it does not explain its complete inventory.
The floor/ceiling translation of THM-4792 concerns different one-step maps
and does not conjugate these positive odd accelerated maps.

**A4 (PROVED).** Fix any odd modulus `m`, any residue `r`, either sign, and
any `a>=1`. Infinitely many positive odd integers satisfy

    n=r mod m,    v2(3n+epsilon)=a.                             (3)

*Proof.* Solve `3n+epsilon=2^a mod 2^(a+1)`; since 3 is invertible modulo
that power of two, this specifies one odd residue. CRT combines it with
`n=r mod m`. All positive representatives in the resulting progression
have the exact same valuation. ∎

So an odd-modulus vertex does not determine its physical halving clock.
Even the next odd residue can differ: modulo 5, `1` and `11` have the same
source residue but `U_+(1)=1`, `U_+(11)=17`; on the minus side use
`1` and `31`, with outputs `1` and `23`. A finite table retaining all
clock-labelled edges is a valid observation, but selecting the actual edge
requires the integer/dyadic sidecar. This is the point where the proposed
modular explanation still owes a global argument.

## 5. What edge order really controls: an exact cycle criterion

This is the inherited affine-cycle mechanism, recovered explicitly from
**THM-4474**, [strategy-cube cycle densities](../../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md),
and the [minus-sheet ladder, S2-S3](../collatz_mod6_20260922_w6_tao_minus_sheet_ladder.md).
The latter already contains the seven-cycle word and numerator below.

Let `a=(a_0,...,a_(r-1))` be positive integers, `A=sum a_i`,
`A_j=sum_(i<j) a_i`, and put

    D=2^A-3^r,
    C(a)=sum_(j=0)^(r-1) 3^(r-1-j) 2^A_j.                       (4)

**A5 (PROVED, exact word gate).** The rational periodic start for that
valuation word is `n=epsilon C(a)/D`. It is a positive odd integer cycle
with precisely those valuations iff

    epsilon D>0  and  D divides C(a).                        (5)

*Proof.* Composing the affine steps gives
`2^A n_r=3^r n_0+epsilon C(a)`, proving the unique periodic candidate.
Both `C` and `D` are odd and `D!=0`. For the cyclically rotated word let
`C_j` denote its numerator. Directly from (4),
`2^(a_j) C_(j+1)=3 C_j+D`. If `D|C_0`, oddness of `D` implies `D|C_j`
for every rotation. The integers `n_j=epsilon C_j/D` are odd, have the
specified common sign, and satisfy `3n_j+epsilon=2^(a_j)n_(j+1)`.
Oddness of the latter quotient makes every valuation exact. This proves
sufficiency as well as necessity. Primitive words modulo rotation classify
primitive cycles; a repeated state would repeat its valuation word. ∎

Every clock word therefore belongs to precisely one **positive rational**
sheet, chosen by the sign of `D`. The plus sheet requires `2^A>3^r`; the
minus sheet requires `2^A<3^r`. The integer gate is additional. Word order
matters through `C`, even when all clock counts and their total agree:

| Valuation word | `3^7-2^11` | `C` | Minus-sheet rational start |
|---|---:|---:|---:|
| `(1,1,1,2,1,1,4)` | 139 | 2363 | 17 |
| `(1,1,1,1,1,2,4)` | 139 | 2123 | `2123/139` |

The first is the known integer cycle
`17->25->37->55->41->61->91->17`; the second has the same valuation multiset
and lies on the same positive rational sheet, but fails integrality.
This supplies a precise role for ordered edges. It does not derive a
one-versus-three total from which edge is displayed first.

The finite census here exhausts **106,761** positive valuation words with
`r<=8` and `A<=18`, or **15,907** primitive rotation classes. It recovers
only `(1)` for `U_+` and `(1)`, `(5,7)`, `(17,25,37,55,41,61,91)` for `U_-`
within that universe. Completeness for arbitrary `r,A` remains **OPEN**.
Even a complete classification of cycles would not alone rule out divergent
nonperiodic positive orbits. The standard Collatz conjecture asks for every
positive orbit to reach 1; see [Lagarias's primary bibliography](https://arxiv.org/abs/math/0309224).

The chess-tempo analogy can suggest studying the positive-domain boundary
and the first-passage order. It is not a proof of (5) for all words or of
cycle completeness. Both accelerated maps already fix the least positive
odd integer 1; no advantage follows just from listing that vertex first.

## 6. Bounded exact controls

**A6 (PROVED, a precise finite-field bridge for clock words).** Fix a prime
power `q`, a bijection from the clock alphabet `{1,...,q}` to `F_q`, and a
normal basis of `F_(q^r)` over `F_q`. Encoding the clock word in that basis
identifies rotation with Frobenius; primitive clock necklaces have field
degree `r`. The integer-cycle predicate (5) is invariant under rotation by
the numerator recurrence, so it selects a well-defined subset of these
Frobenius orbits. It still requires the **integer** decoding of the letters,
their sum `A`, and the carry divisibility; it is not supplied by field degree
or field multiplication. For `q=4`, both seven-letter words in the table
have field degree seven, but only the first passes the integer gate. This
is an actual bridge between the two necklace settings and an explicit
counterexample to discarding the arithmetic sidecar. Unbounded valuations
are not covered by one fixed finite alphabet.

`arithmetic.py` constructs irreducible binary polynomials and explicit normal
bases for `m=3,5,7,9`, testing every field element, its Frobenius degree,
the marked-pair orbit and a literal failure of reflection multiplicativity.
It tests both fixed-clock modular maps for odd `m=3..101`, `a=0..6`
(18,200 edges per sign), and 124,800 explicit CRT lifts for clocks 1..8.
The cycle-word universe is stated above, with exact valuation verification
on every accepted integer cycle. The all-order proofs do not rest on these
finite ranges. No finite cycle census is a proof of the global conjectures.
