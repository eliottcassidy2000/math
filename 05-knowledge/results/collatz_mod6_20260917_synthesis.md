# Collatz from the mod-6 rows: the 3-adic mirror, signed cycles, and what the peripheral scaffolding actually proves

**Status: PROVED scoped elementary and analytic statements; FINITE-EXACT
censuses with independent audits; CITED background. The Collatz conjecture,
completeness of the signed cycle catalogs, the twin-prime and Goldbach
conjectures, and LRC(14) remain OPEN. The user's geometric blueprint and its
"peripheral scaffolding" are REFUTED as a proof (inherited audit, extended
here). No literature-priority claim is made for the classical parts.**

Session `collatz-mod6-20260917` (machine `mac-mini`, worked 2026-09-17 and
2026-09-21). Twelve lane notes carry the proofs and exact tables; this note is
the cross-lane synthesis and the honest connection ledger. Lane notes:

| Lane | Note | Script |
|---|---|---|
| Row braid and the factorization 63 | [row_braid_typing](collatz_mod6_20260917_row_braid_typing.md) | `collatz_mod6_20260917_row_braid_typing.py` |
| Evens also go to `3n+1`: the graph `E` | [extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md) | `collatz_mod6_20260917_extended_collatz_scc.py` |
| The greedy 3-adic map `G` and its duals | [three_adic_g_map](collatz_mod6_20260917_three_adic_g_map.md) | `collatz_mod6_20260917_three_adic_g_map.py` |
| Reverse tree as a three-row automaton | [reverse_tree_pieces](collatz_mod6_20260917_reverse_tree_pieces.md) | `collatz_mod6_20260917_reverse_tree_pieces.py` |
| Bang, `-7/4`, `-29/16` | [zsigmondy_triad](collatz_mod6_20260917_zsigmondy_triad.md) | `collatz_mod6_20260917_zsigmondy_triad.py` |
| Triangles and the semicircle chart | [pythagorean_semicircle](collatz_mod6_20260917_pythagorean_semicircle.md) | `collatz_mod6_20260917_pythagorean_semicircle.py` |
| Sandwich mixed-cell drift | [sandwich_bias](collatz_mod6_20260917_sandwich_bias.md) | `collatz_mod6_20260917_sandwich_bias.py` |
| Sandwich cell ordering and the "4-tournament" | [cell_ordering_scale](collatz_mod6_20260917_cell_ordering_scale.md) | `collatz_mod6_20260917_cell_ordering_scale.py` |
| Divisor balance `F = alpha S + beta U` | [divisor_balance_family](collatz_mod6_20260917_divisor_balance_family.md) | `collatz_mod6_20260917_divisor_balance_family.py` |
| Floor sums and odd functions | [floor_sums_odd_functions](collatz_mod6_20260917_floor_sums_odd_functions.md) | `collatz_mod6_20260917_floor_sums_odd_functions.py` |
| Typing of every "three" and "two" | [wild_typing](collatz_mod6_20260917_wild_typing.md) | `collatz_mod6_20260917_wild_typing.py` |
| Audit of the pasted scaffolding | [scaffolding_audit](collatz_mod6_20260917_scaffolding_audit.md) | `collatz_mod6_20260917_scaffolding_audit.py` |

## Inheritance and concept board

This session continues the two `arithmetic_braids` sessions of 2026-09-17
([first synthesis](arithmetic_braids_20260917_synthesis.md),
[second synthesis](arithmetic_braids2_20260917_synthesis.md)) and the
[blueprint audit of 2026-09-21](collatz_blueprint_20260921_synthesis.md). It
inherits without re-proving: the three rows `6j+1 -> 9j+2`, `6j+3 -> 9j+5`,
`6j+5 -> 9j+8`; the inverse-fibre braid `R(n)=4n+1` with `F(R(n))=4F(n)`;
the cycle gate `n_0 = bB/(2^K-3^L)` and the minimal parameter
`q=|2^K-3^L|/gcd(B,|2^K-3^L|)`; `F=S+U` iff `N in {p,p^3,p^2qr}`; the CRT
sandwich identity; the floor identities detecting squarefreeness and
primality; the Fano/Hamming/E8 sign-gauge map; the exact affine word model
`W(x)=(2^r x - B)/3^m` and the descent certificate
`T^L(n)<n iff 3^L n+B<2^K n`.

Anchor: the user's mod-6 rows and the question "what arrows appear if evens
also go to `3n+1`". Niche: the sign of the sandwich mixed-cell discrepancy and
the scale dependence of the cell order. Wildcard: the triple
`{63, -7/4, -29/16}`. Closest proved mechanism: the inverse-fibre odometer.
Canonical hostile: a long all-one halving word (`n=2^(L+1)-1`), and its new
3-adic mirror `m = 1 mod 3^j`. Corrected near miss: promoting residue-level
coverage to control of one integer's whole orbit. Least-used sidecar: the
ordered exponent word together with its carry, which this session finds is
the *same* object in the 2-adic forward and 3-adic inverse readings.

Research cards used: "Separate unbounded local support from a
height-bounded modular cover", "Search the statement before the method",
"Type every analogy and every implication", "Attack a proposed bound before
extending it", and "Test structured adversaries, not only random samples"
([META-PATTERNS](../../00-navigation/META-PATTERNS.md)).

The six live concepts and how each new object changed them:

| Concept | Retained coordinate | Changed by this session |
|---|---|---|
| Three rows mod 6 | `2^k` exponent class mod 6 | one factorization `63=2^6-1=7*9` produces rows, exponent classes, and the period-3 braid; Wieferich primes bound the reading |
| Inverse fibre / reverse tree | ordered halving word, carry `B` | rows of children are a 3-periodic automaton; the minimal child is a section |
| Nondeterministic graph `E` | arrows `n->n/2`, `n->3n+1` for all `n` | new arrows are one backward `R`-step; multiples of 3 transient; conjecture: one giant SCC |
| Greedy 3-adic map `G` | residue word mod `3^(J+1)` | exact Markov chain, invariant measure, drift `log(2/3)`, tail `(7/9)^(J-1)`; the 3-adic mirror of Terras |
| Signed parameter `b` | `q`, Catalan solutions | the three `3n-1` cycles are two Catalan cycles plus one sporadic cycle |
| Quadratic parameter `c` | Thue unit, parabola `c=-7/4-s^2` | the `n=3` primitive-divisor failures over all of `Q` are exactly `c in {0,-1,-2,-7/4}` |

## 1. One factorization behind the three rows

**PROVED** ([row_braid_typing](collatz_mod6_20260917_row_braid_typing.md),
[zsigmondy_triad](collatz_mod6_20260917_zsigmondy_triad.md)). The user's
observation that the powers of two in the image rows `2,5,8 mod 9` have
exponents `1,5,3 mod 6` is the statement `ord_9(2)=6`, i.e. `9 | 2^6-1 = 63`.
The same number governs the braid: `R^3(n)=64n+21`, `R^3(n)-n=21(3n+1)`,
so `R^3` fixes the source modulo `42` and the target `F(n)` modulo `63`.
The factor `3` sees only the parity gate `k` odd, the factor `7` sees
`k mod 3`, and `9` sees both. Bang's exception `2^6-1` is exactly Catalan's
`3^2-2^3=1`: `2^6-1=(2^3-1)(2^3+1)=7*9`, so no prime has `ord_p(2)=6`, and
that is why the exponent law lives modulo the composite `9`. The reading
"rows exist because the multiplier is 3" is separate: for `pn+1` the powers
of two occupy `ord_(p^2)(2)/ord_p(2)` rows, which is `p` exactly for
non-Wieferich `p`; at `p=1093` and `3511` all powers of two sit in one row.

The trunk of the tree, i.e. the odd predecessors of powers of two, is the
`R`-orbit of `1`: `(4^j-1)/3 = 1,5,21,85,341,1365,5461,...`. Every term with
`j>=2` has a primitive prime (Zsigmondy for `4^j-1` has no exception because
`4-1=3` and `4+1=5` is not a power of two), and at prime index `p>=5` the term
is a base-2 pseudoprime (Cipolla; `341=11*31`, `5461=43*127`, `1398101`).
So the pasted "341 bottleneck" is REFUTED and replaced by a true statement,
which the [wild_typing lane](collatz_mod6_20260917_wild_typing.md) sharpens
into a map: `N_j=(4^j-1)/3=R^(j-1)(1)` satisfies `3N_j+1=4^j`, so the Collatz
halving exponent of the trunk term equals `ord_(N_j)(2)=2j`; `N_j` is a
base-2 pseudoprime iff `6j | 4^j-4`, and this index set is closed under
`j -> N_j`, so every prime `p>=5` starts an infinite tower
`p -> N_p -> N_(N_p) -> ...` of pseudoprimes on the trunk (`|W|=2321` for
`j<=20000`, `61` composite indices). Squarefreeness of `N_j` obeys the
lifting-the-exponent law `p^2 | N_j iff ord_p(4) | j and (p | j or p
Wieferich)`; the guess "`N_j` squarefree iff `gcd(j,N_j)=gcd(j,3)`" is
REFUTED at `j=182=ord_1093(4)`.

The identity `c(c+1)^2=-63/64` at `c=-7/4`, hence `f^3(0)=c/2^6=-7/256`, is
exact, and so is the Thue-unit reading `4rho^2-7 = rho^(-14)` in the plastic
field. These are two different Diophantine equations (`a(a+4)^2=-63` and
`F(a,b)=+-1`) whose small solutions both display `63=7*3^2`. No map between
`63` and `-7/4` beyond that numeral was found (SCOPE; both lanes agree).

## 2. Evens also go to `3n+1`: the graph `E` and its 3-adic mirror

Let `E` be the graph on positive integers with arrows `n -> n/2` for even
`n` and `n -> 3n+1` for **all** `n`.

**PROVED** ([extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md)).
The arrows absent from Collatz are exactly `2j -> 6j+1`. For `v=6j+1>=7`
the least odd `T`-predecessor is `n_0=(4v-1)/3 = 8j+1`, and the new
predecessor is `R^(-1)(n_0)=(n_0-1)/4=2j`: **the new arrows extend each
`1 mod 6` inverse fibre by exactly one backward `R`-step**, while `5 mod 6`
and `3 mod 6` targets receive nothing. Multiples of `3` are transient
(no arrow enters `3Z` from outside; inside `3Z` only halving).

Two questions separate cleanly:

* `Q1`: every `n` reaches `1` in `E` (implied by Collatz).
* `Q2`: `1` reaches every `m` not divisible by `3` in `E`.

`Q2` is equivalent to reducing `m` to `1` by the inverse moves `m -> 2m`,
`m -> (m-1)/3` (any parity), i.e. by compound moves `m -> (2^k m-1)/3` with
`2^k m in {4,7} mod 9`. **Conjecture (E-SCC):** the non-multiples of `3`
form a single strongly connected component of `E`; equivalently `Q1` and `Q2`.

The **greedy 3-adic map** `G(m)=(2^k m-1)/3`, `k` minimal, is a self-map of
`{m : 3 does not divide m}`. **PROVED:** its residues mod 9 form an exact
Markov chain (from `{1,2,4,5}` the next residue is uniform on `{1,4,7}`,
from `{7,8}` uniform on `{2,5,8}`); the stationary law is `2/9` on each of
`1,4,7` and `1/9` on each of `2,5,8`; `E[k]=1`, so the stationary drift is
exactly `log(2/3)`; the invariant measure is
`mu=(2/3)Haar(1+3Z_3)+(1/3)Haar(2+3Z_3)`. The first `J` letters depend only
on `m mod 3^(J+1)` (sharp), and the tilted matrix
`[[5/3,1/3],[10/3,2/3]]` (determinant `0`, trace `7/3`) gives
`E[2^(K_J)] = 3(7/3)^(J-1)` exactly and the tail bound

    #{m<=X : 3 does not divide m, sigma(m)>J} <= (7/9)^(J-1) (2X/3 + 2*3^J)

for the greedy stopping time `sigma`. This is the 3-adic mirror of the
Terras--Everett theorem: same shape (a residue-determined finite prefix, a
negative drift, exponential decay of the non-stopped density), opposite
completion. The [G lane](collatz_mod6_20260917_three_adic_g_map.md)
sharpens it: `(Z_3^x, G)` is topologically conjugate to the six-state
subshift of finite type above (entropy `log 3`, `mu` its Parry measure);
for every tilt `u` the exact exponential moment is
`E[u^(K_J)] = (u+1)(u^2+2)/6 * ((u^2+u+1)/3)^(J-1)`, so the tilted matrix
has rank one for all `u` and `(7/9)^(J-1)` is its `u=2` Chernoff bound;
the exact residue-level stopping densities are monotone with
`d_12 = 530885/531441 = 0.998954` (Collatz's Terras density at `J=12` is
`0.9448`); the residue stopping time equals the true stopping time for
every `m` with `sigma(m)<=12` (proved per level through the exact
threshold `B_i/(2^(K_i)-3^i)`) and for all `m<=10^6`. **FINITE-EXACT:**
every `m<=10^7` coprime to `3` reaches `1` under `G` (maximum `93` steps,
attained at `m=8751065` and `9454814`; worst peak ratio `133.03` at
`m=4847486`, whose word prefix
`3 2^5 0 3 2^9 0` is forced by `v_3(G(m)-1)=6` and `v_3(G^8(m)-1)=10`);
`Q1` holds to `10^6`; `E` restricted to `[1,2000]` has exactly `74` simple
cycles of length `<=40`, every one using an even-to-`3n+1` arrow. The peak
rate `4/3` per step is optimal (`K_n <= 2n+[k_1=3]`), and two consecutive
`k=3` letters are impossible because residue `5 mod 9` is entered only from
`{7,8}`.

The hostile families mirror too: for Collatz the prefix `n=2^(L+1)-1`
grows for `L` steps; for `G` the prefix `m = 1 mod 3^j` forces `j-1` steps
of factor `4/3` and then one of `1/3`, net `4^(j-1)/3^j`, which exceeds `1`
exactly when `j>=5` (`244 -> 325 -> 433 -> 577 -> 769 -> 256`), and the
family `m=3^j+1` lands on `4^(j-1)` with `peak/m ~ (4/3)^(j-1) -> infinity`. In
both readings the obstruction to a proof is identical: a residue-determined
prefix cannot control the whole orbit of one fixed integer. The
generalization `G_b(m)=(2^j m-b)/3` keeps `E[j]=1` for every odd `b`
coprime to `3`, so every member of the family has drift `log(2/3)`, exactly
as every `3n+b` map has forward drift `log(3/4)`. Its cycles obey the gate
`m_0 (2^J-3^L) = b B'(w)` with `B'` the inherited carry of the *reversed*
word, the scaling lemma `G_b(bm)=bG_1(m)` gives every `G_b` the three
universal cycles `b{1}`, `b{-1}`, `b{-4,-11}` (the `3n+b` maps have four),
and a bounded census (`|b|<=49`, `|m|<=10^5`) shows `#cycles(G_b)` for
`b=1,5,7,...,49` equal to `3,5,4,5,7,6,6,8,7,4,4,7,4,6,5,6,8` against the
signed `3n+b` counts `4,9,5,7,13,...`: the two families share the carry and
the gate polynomial but not their branch selection, so there is no
cycle-for-cycle bijection (SCOPE).

**Connection contract (2-adic/3-adic duality).** Source: the affine word
`W(x)=(2^r x-B)/3^m` of the blueprint audit. Target: the same word read as
a `G`-path. Map: reverse the arrows and exchange the roles of `2` and `3`
in the guard. Preserved: the carry `B` and the integrality gate. Lost:
the branch-selection rule (forced valuation forward, minimal admissible
`k` backward). Sidecar: the full residue of `m` modulo `3^(J+1)`. Decisive
test: `m=244` (grows under `G`) versus `n=31` (grows under `T`).

## 3. Signed parameters: why `3n-1` has three cycles

**PROVED** (inherited [signed_cycles](arithmetic_braids2_20260917_signed_cycles.md);
typed here in [wild_typing](collatz_mod6_20260917_wild_typing.md)). Negation
conjugates `3n+b` to `3n-b`; cycles come in negatives. A halving word is a
cycle word at parameter `b` iff `q | b`, `q` the reduced denominator of
`B/(2^K-3^L)`. When `b=2^K-3^L` itself, every composition of `K` into `L`
positive parts is a cycle word, so the number of cycles jumps with the
number of necklaces. Catalan's equation `2^K-3^L=+-1` has only
`(K,L) in {(1,1),(2,1),(3,2)}`, so `3n+1` receives exactly one such cycle
(`{1}`) and `3n-1` exactly two (`{1}` and `{5,7}`); its seven-cycle
`17,25,37,55,41,61,91` is *sporadic* (`2^11-3^7=-139` divides the carry).
So "the ones with three cycles" decompose as two Catalan cycles plus one
sporadic cycle; for `b=+-5` the nine cycles are `5x` the `3n-1` cycles plus
five primitive ones, two of them necklace cycles of `5=2^5-3^3`. The census
of cycle counts `c(b)` for `gcd(b,6)=1`, `|b|<=99`, from starts
`|n|<=2*10^5` (`c(+-1)=4`, `c(+-5)=9`, `c(+-7)=5`, `c(+-11)=7`,
`c(+-13)=13`, maximum `c(65)=19`) attains
`{4,5,6,7,8,9,11,12,13,15,17,19}` and misses `{1,2,3,10,14,16,18}` below
`20`; `c(b)>=c(1)=4` on the known catalog is PROVED (`T_b(bn)=bT_1(n)` scales
the Collatz cycles), the superadditivity `c(b_1b_2)>=c(b_1)+c(b_2)-c(1)` is
PROVED for true counts, the number of fixed points is
`#{K>=1 : (2^K-3) | b}`, and `7=c(11)=c(37)=c(49)` is attained, so the
spectrum has no relation to the tournament holes `{7,21}` (SCOPE). The user's `3n-5` is exactly `T_(-5)(n)=T(n-2)`: it is
conjugate to `m -> T(m)-2`, not to Collatz, its rows are Collatz's rows
precomposed with `r -> r-2 mod 6`, and its nine cycles are the negatives
of the nine `b=5` cycles.

## 4. The reverse tree as a three-row automaton

**PROVED** ([reverse_tree_pieces](collatz_mod6_20260917_reverse_tree_pieces.md),
[row_braid_typing](collatz_mod6_20260917_row_braid_typing.md)). The
`j`-th child of an internal node `u` lies in row `base(u mod 18)+4j mod 6`,
i.e. `rho(2^h u mod 9)` with `4 -> row 1`, `7 -> row 5`, `1 -> row 3`; rows
cycle `1,5,3` with exact period `ord_9(4)=3`; row-`3` children are exactly
the leaves (odd multiples of `3`, which have no `T`-preimage). Because
`ord_27(4)=9`, the first nine children of every internal node hit each
residue class mod `9` exactly once, so the infinite fibre is uniform
(`1/3` leaves, `1/6` per internal class); along a fixed index path the rows
of the first `D` descendants depend on `u mod 3^(D+1)` and not on
`u mod 3^D`, the tree-side twin of the `G` law. All `500,000` odd
`u<=10^6` reach `1` (maximum depth `195` at `837799`), and by depth the
residues mod `9` are uniform to within truncation noise; the fibre split is
`1/2 : 1/2` between the internal classes mod `3`, not `G`'s `2/3 : 1/3`,
because the tree counts every child once while `G` skips leaves and takes
the minimal admissible child. The expected number
of depth-`k` descendants is exactly `(4/3)^k`; a split-independent modulus-9
difference inequality yields the rigorous but weak exponent `0.246`, while
the literature has `0.84` (Krasikov--Lagarias 2003, CITED). The
minimal-child map `m(u)=(2^(h0+1)u-1)/3` is an injective section
(`T o m = id`) with image the odd `n` not congruent to `5 mod 8`, unique
fixed point `1`, and exactly `(2/3)^L` of the classes mod `3^(L+1)` survive
`L` minimal-child steps; `G` and `m` coincide on odd `u` iff
`u = 1,2,8 mod 9`. The tree is the free object on the guarded branches
`D(x)=2x` and `E(x)=(2x-1)/3`; Collatz is its surjectivity onto the odd
integers. The user's "show how the three pieces fit" is
therefore already exact; what is not exact is *surjectivity*, and the
relaxed version `Q2` shows precisely where the difficulty sits: dropping
the parity guard turns the reverse tree into `E`, gives it the certificate
`G`, and still leaves a 3-adic hostile family. Reachability **from** `1`
(`Q2`) and reachability **to** `1` (Collatz) are exchanged by reversing the
arrows; the deterministic Collatz graph is the subgraph of `E` that keeps
only the forced branch.

## 5. Quadratic dynamics and triangles

**PROVED / CITED** ([zsigmondy_triad](collatz_mod6_20260917_zsigmondy_triad.md),
[pythagorean_semicircle](collatz_mod6_20260917_pythagorean_semicircle.md)).

* The user's `x^2+{0,1,2}` graphs are those of `x^2-0`, `x^2-1`, `x^2-2`,
  the three rational post-critically finite quadratics.
* `disc_x Phi_3(x,c) = -(4c+7)^3 (16c^2+4c+7)^2`: `c=-7/4` is the real
  period-3 saddle-node, the collided 3-cycle has multiplier exactly `1`, and
  the unmarked 3-cycle curve is the parabola `c=-7/4-s^2`; `-29/16` is
  `s=+-1/4`, `-2` is `s=+-1/2`.
* The third critical numerator is `a*F(a,b)`, `F=a^3+2a^2b+ab^2+b^3`, every
  prime of `F` primitive, so failure is the unit equation `F(a,b)=+-1` in
  the plastic field. PARI's unconditional Thue solver gives the complete
  solution set `+-{(1,0),(0,1),(-1,1),(-2,1),(-7,4)}`, verified
  independently by the session lead: **over all of `Q` the parameters whose
  third critical term has no new prime are exactly `c in {0,-1,-2,-7/4}`**.
  This closes the "all-height third-numerator unit equation" left open in
  the [first braids synthesis](arithmetic_braids_20260917_synthesis.md).
* Microcosm incidence is finite and exact: `-7/4` has period `1,2,3` for
  `c=-77/16,-37/16,-29/16` and no period `4`; the pattern `f^n(0)=c/2^(2^n-2)`
  at the real period-`n` parabolic parameter holds for `n<=3` and fails at
  `n=4` (REFUTED as a law).
* For a primitive triple, `c+b=(m+n)^2` and `c-b=(m-n)^2` are odd squares
  while `c+a=2m^2`, `c-a=2n^2` never are; the odd pair `(s,t)=(m+n,m-n)` is
  the complete identity and `s^2=c+b` is one of its two coordinates, with
  exactly `phi(s)/2` triples over each `s`. The `(m,n)` triangle is the
  half-angle triangle. The literal reading "hypotenuse `k^2+1`, altitude
  `sqrt k`" has no rational member (`3-4-5` has altitude `12/5`); the reading
  that fits is the Euclid column `(k^2-1,2k,k^2+1)`. The semicircle chart
  `(e/d,l,theta)=(tan^2 theta, sin(2theta)/2, theta)` is monotone on
  `(0,pi/4]`; primitive triples are its rational points, approaching the
  degenerate end along the Euclid column and the isosceles end along the
  Pell orbit `(m,n)->(2m+n,m)` (hypotenuses `5,29,169,985,...`, canon
  THM-3341). Angle doubling is Gaussian squaring (THM-3341), its hypotenuse
  sequence is the multiplicative squaring forest of the summand note, and a
  triple is a double iff its hypotenuse is a perfect square. The appearance
  of `29` both as `-29/16` and as the second Pell hypotenuse is a coincidence
  (REFUTED as a map).

## 6. The sandwich matrix: drift, order, and the false tournament

**PROVED / FINITE-EXACT / HEURISTIC** ([sandwich_bias](collatz_mod6_20260917_sandwich_bias.md),
[cell_ordering_scale](collatz_mod6_20260917_cell_ordering_scale.md)).

* No tournament on four vertices has score profile "max > two equal
  middles > min" (the `n=4` score sequences are `(0,1,2,3),(0,2,2,2),
  (1,1,1,3),(1,1,2,2)`); the user's picture is a weak order with one tie,
  and the only intrinsic relation among cells is the mirror `k -> -k`,
  which orients nothing. The tournament reading is cosmetic.
* `SS > {PS,SP} > PP` holds for every `K in [762,10^7]` (last violations at
  `K=761` and `K=160`); the mirror difference `N_PS-N_SP` never settles
  (`992` sign changes below `10^7`, final `+139` against `sqrt(N)=1125`).
* The class mechanism is exact: a right semiprime is same-class-or-square,
  a left semiprime is mixed-class, and
  `S_1(x)-S_5(x) = (1/2)[sum_(5<=p<=x/5) chi(p) pi_chi(x/p) + pi'(sqrt x)]`.
  The independence prediction for the drift (`36.6, 113.2, 230.0` at
  `K=10^5,10^6,10^7`) splits into a prime-Chebyshev part (`55-59%`) and a
  prime-square/semiprime-class part; the pointwise inequality is REFUTED
  (`K=1453`), the drift is confirmed at `10^6` and only partially at
  `10^7` by a pooled-gap control. The `chi`-sums by `Omega` alternate in
  sign (`-363,+524,-426,+265` at `6*10^7`), the finite-scale face of the
  Ford--Sneed/Meng alternation.
* Landau's leading term is wrong at these scales; a Selberg--Delange product
  model reproduces the tie structure and puts the `CC > SS` crossover far
  beyond `10^7` (HEURISTIC).

## 7. Divisor balance and floor sums

**PROVED / FINITE-EXACT** ([divisor_balance_family](collatz_mod6_20260917_divisor_balance_family.md),
[floor_sums_odd_functions](collatz_mod6_20260917_floor_sums_odd_functions.md)).
`F = alpha S + beta U` has finitely many exponent profiles for every
`(alpha,beta)` except `(1,0)` (`F=S` iff squarefree or `p^2`) and `(2,0)`
(`F=2S` iff `p` or `p^3 m`); every cell with `alpha+beta>=1` contains the
prime power `p^(alpha+beta+1)`; `(1,1)` is the unique finite cell among the
prime-cube cells on which `Omega` is injective. The `p^2qr` solutions
overtake the primes at `n=145119` and lead from `147925` on, with
`#{p^2qr<=x} ~ P(2) x loglog x/log x`, `P(2)=0.45224742`. The `k`-free
analogue `F=S_k+U` has four shapes for every `k>=3`; `k=2` is the
exceptional three-shape case. The squarefree divisors of `p^2qr` are the
seven points of the Fano plane, which is a real map to `F_2^3` and, via the
[Fano code note](arithmetic_braids2_20260917_fano_code.md), to an octonion
sign gauge and `E8`; nothing in it encodes Bott periodicity (SCOPE).

The three floor sums obey one odd-function law: for odd `f in Z[x]`,
`sum floor(f(k)/n) = (sum f)/n - (n-1)/2 + Z_f(n)/2` with
`Z_(x^e)(n)=n/n_e-1`, so every odd-power identity holds iff `n` is
squarefree (truth set of density `6/pi^2`), the cube-root sum is its
lattice reciprocal, the bilinear sum equals `(n-1)^2(n-2)/4+(P(n)-2n+1)/2`
with Pillai's gcd-sum `P`, holding iff `n` is prime; and even powers are
different in kind: for a prime `p = 3 mod 4`, `p>3`, `sum floor(k^2/p)`
deviates from the odd law by exactly the class number `h(-p)` (Dirichlet,
CITED; verified against reduced-form counts for all such `p<500`), the
deviation is `0` for `p = 1 mod 4`, and for a general even exponent `e`
the deviation is `h(-p)` whenever `gcd(e,p-1)=2` and `p = 3 mod 4` (the
converse fails: `(p,e)=(499,6)` also has deviation `h(-499)=3`); the
quartic deviation for `p = 5 mod 8` is `-(2/p)(A_0-A_2)` in terms of
index-class sums and has no class-number closed form here (OPEN). The
squarefree and prime characterizations were proved independently in
[floor_reciprocity](arithmetic_braids2_20260917_floor_reciprocity.md); the
`Z_f` law, the Pillai form and the class-number contrast are this session's.

## 8. The pasted scaffolding, audited

The second half of the paste (peaks `1,26,80`, the "valve" `160` for `23`,
`196`, the "non-local descent inequality") is the attachment already audited
by the concurrent opus session in
[collatz_guards_20260921_synthesis](collatz_guards_20260921_synthesis.md):
the exact `23` reset word `(1,1,5)` is the cylinder `n = 23 mod 256`
(`128 T^3(n)=27n+19`), "`5 mod 9` forces the reset" is REFUTED by
`95->143->215->323`, the family `3*2^(2a+1)-1` extending `23` grows at its
first reset for every `a>=4`, the first-reset descent density is the
irrational Beatty number `P=sum_(L>=1) 2^(-floor(L log_2 3))=0.7137...`,
no positive integer orbit keeps `2^(K_j)/3^j` in a bounded strip, and
`196=14^2` is the identity `1+sum p_i=(k+1)^2+2 sum c_i` (odd composites
skipped), with no trajectory consequence. This session does not repeat
those proofs; it cites them.

[scaffolding_audit](collatz_mod6_20260917_scaffolding_audit.md) tests every
remaining claim of the pasted block.

| Claim | Verdict | Exact replacement or witness |
|---|---|---|
| "Shear" identity `T((A+1)(B+1)-1)-T(AB-1)=T(A)+T(B)+A(B^2-1)+B(A^2-1)` | REFUTED | fails at all 49 pairs in `[1,7]^2`; the true left side is `(A+B+1)(2AB+A+B)/2 = T(A)+T(B)+AB(A+B+2)`; the paste's flows solve to `S(A,B)=AB(A+B)+C(A,2)+C(B,2)`, which is its *right* side |
| Tetrahedral and pentatope laws | TRUE (PROVED) | Chu--Vandermonde `C(A+B+d,d)=sum_i C(A+a,i)C(B+b,d-i)` with splits `(2,1)`, `(2,2)`; holds in every dimension |
| `g`-operator: `AgB=AB(AgB)`, `AgB=(AB)^(AB)`, `1gB=B!` | REFUTED | the recursion forces `(AB-1)AgB=0`; the readings disagree at `B=2`; the unique symmetric reading with `1gB=B!` and `AgA=(A!)^2` is `AgB=A!B!` |
| `C(n,2)-(n-1)=T_(n-2)` | TRUE | arithmetic |
| "Hamiltonian-path edges carry 0 bits" | REFUTED | `sum_T h(T)=n! 2^(T_(n-2))` (double count; census `2,12,192,7680,737280`), so the code costs `C(n,2)+log_2(n!/2^(n-1))` bits |
| `63` "recycles primes" | TRUE-with-repair | Bang's exception via `2^3+1=3^2` (Catalan); "a loop switching multiplication and addition edges" defines nothing |
| `341` "stalls" the `4x+1` sequence | REFUTED | `(4^5-1)/3=341=11*31` has two primitive primes; `4^n-1` never lacks one; true fact: Cipolla pseudoprimes on the trunk |
| `-29/16` "cancels translation after three steps" | SCOPE | the content is THM-4139/4146, cited |
| Peaks `1,26,80` of seeds `5,7,23` | REFUTED | shortcut peaks `8,26,80` (full map `16,52,160`); every shortcut peak of an odd seed `>1` is `2 mod 6` (PROVED) |
| `4k+1` square tracking, next node `169/160` | REFUTED | over the `20,341` distinct peaks below `10^5` the `s^2+-1` cells are at chance (`z=1.91`); bases `1` and `3 mod 4` indistinguishable; `160=2^5*5` is on the trunk of `5` |
| `196` = sum of the first twelve primes | REFUTED | the sum is `197`; the guards note gives the true identity with skipped odd composites |
| `3^L n+B_L<2^(K_L) n` | TRUE (inherited) | the descent certificate; for `27` the first positive margin is at `L=37`, `K=59` |
| Unified cycle, Bott periodicity in tetration, Wythoff diffraction, quasicrystal Fourier matrix | SCOPE | no definitions; the only real object is the mod-30 wheel with eight reduced classes |

## 9. Typing every "three" and "two"

Real maps found this session or inherited: rows `<->` `63` (`ord_9(2)=6`);
`-29/16` 3-cycle `<->` `3:4:5` (THM-4146); doubling forest `<->` halving;
squaring forest `<->` angle doubling; `p^2qr` squarefree divisors `<->`
Fano plane; `3n-1` cycles `<->` Catalan solutions; forward halving word
`<->` backward greedy word (same carry). Only the numeral is shared, no map
found (SCOPE): the three `F=S+U` shapes, the three almost-prime classes,
the three semicircle parameters, the three negative cycles versus the three
rows, `{7,21}` versus anything Collatz, Bott's `8` versus `Omega(C^2B)=8`,
`6/pi^2` versus negative cycles (its numerator is `zeta(2)`'s, and the
rows conditioned on squarefreeness have densities `9/pi^2, 6/pi^2, 9/pi^2`).

## 10. What this means for a Collatz proof

Nothing here proves Collatz, and the session's strongest results say why
the pasted routes cannot: every certificate that reads a finite prefix
from a residue, whether 2-adic (forward `T`) or 3-adic (backward `G`), has
an explicit family of integers that grows for that prefix, and ambient
densities (squarefree `6/pi^2`, row densities, Chebyshev drifts) carry no
every-orbit quantifier. The exact open obligation is unchanged from the
blueprint audit: for each fixed odd `n>1`, exhibit a halving word with
positive margin `2^K n - 3^L n - B > 0` that `n` actually realizes.

What the session adds is a cleaner target. The E-SCC conjecture splits the
problem into two halves that are mirror images under arrow reversal and the
exchange `2 <-> 3`; the backward half has a deterministic certificate `G`
with a proved exponentially decaying non-stopping density and an explicit
invariant measure, and its hostile families are as explicit as Collatz's.
A proof strategy that treats one adic side must survive the other; a
strategy that uses the carry `B` must use it on both sides at once, because
`B` is the only coordinate the two readings share.

Cheapest next probes: (i) an exact census of `G`-cycles on the negative
integers and of `G_b` cycles against the signed `3n+b` catalog, looking for
a bijection through the carry; (ii) whether the joint 2-adic/3-adic word of
one integer (forward word to its maximum, backward greedy word from the
maximum) satisfies a conservation law for `B`; (iii) a Krasikov--Lagarias
difference-inequality system on the `E` reverse tree, where the extra
even-to-odd edges should raise the exponent above `0.84`; (iv) HYPOTHESIS:
the bounded-strip exclusion of the
[guards discrepancy note](collatz_guards_20260921_discrepancy.md) (no
positive orbit keeps `2^(K_j)/3^j` in a bounded interval) transfers to `G`
with `3^J/2^(K_J)` in place of `2^(K_j)/3^j`, using the tilted-matrix tail
in place of the halving-cylinder density; the repeated-state step transfers
unchanged because a `G`-cycle multiplies the ratio by a nonunit
`3^L/2^K`. This is unproved here.

## Reproduction and audit scope

Each lane script is run as `python3 04-computation/experiments/<script> >
05-knowledge/results/<lane>.out`; normal and `-O` streams agree modulo
timing stamps. Independent audit scripts exist for every lane
(`*_audit*.py`) with their own outputs. Universes, filters, positive and
hostile controls are stated inside each note. Session lanes were run by
subagents in three waves; the worktree was pruned between the two working
days and the scripts were recovered from the agent transcripts and rerun,
which is why some outputs are dated 2026-09-21.
