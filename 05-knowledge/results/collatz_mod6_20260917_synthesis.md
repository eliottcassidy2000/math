# Collatz from the mod-6 rows: the 3-adic mirror, signed cycles, and what the peripheral scaffolding actually proves

**Status: PROVED scoped elementary and analytic statements; FINITE-EXACT
censuses with independent audits; CITED background. The Collatz conjecture,
completeness of the signed cycle catalogs, the twin-prime and Goldbach
conjectures, and LRC(14) remain OPEN. The user's geometric blueprint and its
"peripheral scaffolding" are REFUTED as a proof (inherited audit, extended
here). No literature-priority claim is made for the classical parts.**

Session `collatz-mod6-20260917` (machine `mac-mini`, worked 2026-09-17,
2026-09-21 and 2026-09-22). Twenty-seven lane notes across six waves carry
the proofs and exact tables; this note is the cross-lane synthesis and the
honest connection ledger. Sections 1--10 cover the first three waves
(twelve lanes, prefix `collatz_mod6_20260917_`), section 11 the fourth
wave (prefix `collatz_mod6_20260921_`), section 12 the fifth
(prefix `collatz_mod6_20260922_`) and section 13 the sixth
(prefix `collatz_mod6_20260922_w6_`). First-wave lane notes:

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
spectrum has no relation to the tournament holes `{7,21}` of the
Hamiltonian-path count (`THM-1370-h-spectrum-omits-7-21-all-n`, where
completeness of the hole list is conjectural; THM-1745 is its
arborescence shadow) (SCOPE). The user's `3n-5` is exactly `T_(-5)(n)=T(n-2)`: it is
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

## 11. Wave four (2026-09-21/22): the triangle encoding, the fruit curve, Pillai clocks, and the pasted "circuit"

Five further lanes, run after the concurrent opus sessions
[odd_square_20260921](odd_square_20260921_synthesis.md) (the edge
encoding, its angle bound and edge count) and
[prime_shells_20260921](prime_shells_20260921_synthesis.md) (the `3P+-H`
elliptic tree) had appeared, and designed to add only what those left open.

| Lane | Note | Script |
|---|---|---|
| Berggren transport of Collatz-family edges | [berggren_edge_transport](collatz_mod6_20260921_berggren_edge_transport.md) | `collatz_mod6_20260921_berggren_edge_transport.py` |
| Fruit curve rank and positive multiples | [fruit_rank_positive_multiples](collatz_mod6_20260921_fruit_rank_positive_multiples.md) | `collatz_mod6_20260921_fruit_rank_positive_multiples.py` |
| Pillai clocks and the cycle gate | [pillai_convergents_cycle_gates](collatz_mod6_20260921_pillai_convergents_cycle_gates.md) | `collatz_mod6_20260921_pillai_convergents_cycle_gates.py` |
| `G` on the negatives and the joint carry | [g_negatives_joint_carry](collatz_mod6_20260921_g_negatives_joint_carry.md) | `collatz_mod6_20260921_g_negatives_joint_carry.py` |
| Giuga, AAC, Littlewood, Sarnak typing | [grand_circuit_typing](collatz_mod6_20260921_grand_circuit_typing.md) | `collatz_mod6_20260921_grand_circuit_typing.py` |

**The encoding is a tree language on guards (PROVED).** Write an edge of
the family `E_(a,b,k)` as `a x + b = 2^k y` with `x,y` odd, coprime and
distinct, encoded by the primitive triangle with roots `(max,min)`. In root
coordinates the Berggren children are `B1:(s,t)->(s+2t,t)`,
`B2:(2s+t,s)`, `B3:(2s-t,s)`. The complete transport table says how each
child moves the guard `(a,b,k)`: for a `k=1` edge (`y>x`) the children are
the `7x+b` edge at `k=1`, the braid step `(4x+b, y)` at `k=3` in the same
family, and the `3x-b` edge `(2x+b, y)` at `k=2`; for `k>=2` (`y<x`) the
second and third children are edges of the families `(2^(k+1)+3)x+b` and
`(2^(k+1)-3)x-b` at the same `k`, and the first child belongs to no
odd-multiplier family. Hence no child of a `k>=2` Collatz edge keeps the
multiplier `3`, every `k=1` edge has exactly two legal children, and the
only legal children beyond these are five sporadic small readings
(`3->5`, `7->5`, `3->1`), so the legal `3x+-1` sub-forest of the Berggren
tree consists of depth-one stars plus one root cluster of seven pairs. The
whole inverse fibre of a target `y` lies on one `B1`-ray with `t=y`:
`(x_(k+2),y)=B1^(2^(k-1))(x_k,y)` for `k>=2`, preceded by one `B2` step
from `k=1`; and consecutive orbit triangles are never parent and child
(`0` of `1,849,203` pairs below `10^5`). So "triangle descent" is a legal
Collatz move only through the guard table, never along an orbit. This
extends THM-3756 and the odd-square note by the transport of `(a,b,k)`.

**The fruit curve (PROVED by PARI 2-descent).** `rank E(Q)=1` for
`y^2=x^3+109x^2+224x`, torsion `Z/6` generated by `(56,728)`, `G=(-4,28)`
of infinite order with an index in `E(Q)/torsion` free of primes below
`10^4` (index one is CITED from Bremner--Macleod), canonical height
`1.51870...`; the
repaired triple is `9G` with digits `(81,80,79)`; among `mG+kT`,
`1<=m<=80`, the positive fruit triples occur exactly at odd `m` (the
generator lies on the bounded real component, all torsion on the identity
component) with `m in {9,17,43,51,77}` for even `k` and
`{13,21,39,47,55,73}` for odd `k`; fruit digits grow as
`(3/2) m^2 h(G)/log 10 = 0.98935 m^2`. The user's "ternary tree isomorphic
to the primitive triples" is the labelled `3P+-H,3P` tree of the
prime-shells note, which preserves no arithmetic and no positivity: the
children of the positive `9G` are all non-positive.

**Pillai clocks (PROVED/FINITE-EXACT).** `|2^K-3^L| <= 3^L/(4L)` forces
`K/L` to be a convergent of `log_2 3`, `<= 3^L/(2L)` a convergent or
extreme intermediate fraction (the weaker `< 3^L/L` is refuted as a
criterion by `5/3`); with Eliahou's inequality and the `2^68` verification
every nontrivial positive cycle of length up to `1.49*10^10` has `K/L` a
convergent. The gate census on all `40` convergent and intermediate clocks
with `L<=13` finds `q=1` words only at the repeated unit-gap clocks and the
seven rotations of `(1,1,1,2,1,1,4)` on `11/7` (the `3n-1` seven-cycle);
the clock `19/12` has minimum `q=23` (`7153=23*311`, the user's `23` once
more, numerology). The user's `2^10=10^3+24` is complete-listed
(`23` pairs with `|2^a-10^b|<=100`) and typed: the Collatz equation is
`2^K-3^L=Delta` with `Delta | bB`, not base ten.

**`G` on the negatives and the joint carry (PROVED/FINITE-EXACT).** Every
`m` in `[-10^7,-1]` coprime to `3` reaches `{-1}` or `{-4,-11}` under `G`
(no escapes, no new cycles; `56.5%` to the two-cycle); by gate enumeration
of all `1,398,100` greedy words of length `<=10`, the only `G`-cycles of
length `<=10` of either sign are `{1}`, `{-1}`, `{-4,-11}`, so the mirror
table reads `T`: one positive and three known negative cycles (the
negative list is complete for `L<=7`), `G`: one positive and two negative,
the `(K,L)=(3,2)` row splitting the greedy `{-4,-11}` from the all-odd
`{-5,-7}`. The worst negative excursion ratio is again
`133.033` (at `m=-2391484`, the same `(K_17,17)=(34,17)` accident as
`4847486`). The joint carry of a start `n`, its orbit maximum `M` and the
greedy return from `M` to `1` satisfies only the conservation identity
`3^(L_1) n + B_1 = 2^(K_1-K_2)(3^J+B_2)` and word-forced parities: an
exhaustive invariant search on `(B_1,B_2,K_1,K_2,L_1,J)` over all odd
`n<=10^5` finds NONE beyond it. The greedy return retraces `n`'s own orbit
exactly while every pre-peak iterate `p` satisfies `G(u)=p`, which is the
residue criterion "`p = 1 mod 3` and `k<=3`, or `p = 2 mod 3` and `k=1`"
(eleven of sixteen classes mod `48`); at the peak `M = 5 mod 12` and
`B_2` is odd. The `G`-analogue of the guards' bounded-strip theorem is
PROVED on both sheets: on the positive sheet by positivity alone
(`sum 3^s/2^(K_s) <= 3m_0`, so `3^j/2^(K_j) -> 0`), on the negative sheet
by a complete transfer of the guards' source-density argument with the
tilted-matrix tail in place of the halving cylinders. This settles the
hypothesis (iv) of section 10 (stated there before this lane ran).

**The pasted "grand circuit" (verdict table in the typing note).** The
four conjectures are stated exactly with sources. The user's Giuga
criterion `p^2(p-1) | n-p` is PROVED equivalent to the
Borwein--Borwein--Borwein--Girgensohn form, and no composite below `10^5`
satisfies it. For AAC the norm equation forces `v_3(u)=0` and `v_2(u)<=1`
for every `p`, so no `2`- or `3`-adic bridge to this session can exist; no
violation below `3000`. Sarnak's hypothesis is void: Collatz on `Z_2` is
the full `2`-shift (entropy `log 2`) and `G` is a six-state shift of
entropy `log 3`; a `2`-adic point whose parity sequence encodes the Mobius
indicator turns the proposed "orthogonality" into a two-point Chowla sum,
which is itself open. Littlewood involves two numbers and Collatz only
one; if `log_2 3` has unbounded partial quotients (records `1,2,3,5,23,55`
in the first sixty) the conjecture is trivial for every pair containing
it, and `phi` enters the audited notes only inside refuted or SCOPE items
(the blueprint energy's Wythoff word, the Zeckendorf remark of the
summand note). None of the three proposed
formalizations has a map; the honest formalizable bridges remain the
Pillai clocks through the cycle gate and the `E`-graph/`G` mirror.

## 12. Wave five (2026-09-22): the portrait of a counterexample, the minus-sheet control, and the third blueprint

| Lane | Note | Script |
|---|---|---|
| Counterexample portrait | [counterexample_portrait](collatz_mod6_20260922_counterexample_portrait.md) | `collatz_mod6_20260922_counterexample_portrait.py` |
| Minus-sheet positive control | [minus_sheet_positive_control](collatz_mod6_20260922_minus_sheet_positive_control.md) | `collatz_mod6_20260922_minus_sheet_positive_control.py` |
| Paley `T_7`, two Fano planes, octonions | [paley_fano_octonion_design](collatz_mod6_20260922_paley_fano_octonion_design.md) | `collatz_mod6_20260922_paley_fano_octonion_design.py` |
| Block-matrix spectrum audit | [block_spectrum_audit](collatz_mod6_20260922_block_spectrum_audit.md) | `collatz_mod6_20260922_block_spectrum_audit.py` |
| Compression entropy and Lean audit | [compression_and_lean_audit](collatz_mod6_20260922_compression_and_lean_audit.md) | `collatz_mod6_20260922_compression_and_lean_audit.py` |

**Portrait of a counterexample (PROVED from CITED inputs; no contradiction
found).** Every proved necessary condition on a hypothetical nontrivial
positive `3n+1` cycle and on a hypothetical divergent orbit is collected
in two tables and tested on the `3n-1` sheet through the sheet criterion
`T_+(-n)=-T_-(n)`. The gate `q=1`, the halving-word cylinders, the
stopping-time densities and the bounded-strip exclusions are SHEET-BLIND;
only the sign law (`Delta>0` for a positive `3n+1` cycle, `Delta<0` on the
minus sheet), the convergent theorem's magnitude hypothesis, and the `2^68`
verification are SIGN-SPECIFIC. From the `2^68` bound, Legendre and
Hardy--Wright 171 alone, every nontrivial positive `3n+1` cycle has
`L > 14,878,203,146`; the minus-sheet mirror from a new `10^7` census is
`L > 2738` for any fourth positive `3n-1` cycle. No two conditions in the
portrait contradict each other, which is the honest statement of where a
proof stands: the conjunction of everything proved is still satisfiable.

**The minus sheet as a positive control (PROVED).** The parity-word map
modulo `2^J` is a bijection on both sheets and the prefix-descent counts
of the two sheets coincide for every `J`; hence any "density-one descent"
argument proves the same statement for `3n-1`, where three cycles exist,
and cannot be a proof. Coefficient descent implies actual descent on the
minus sheet (and conversely on the plus sheet); the Berggren child `B3`
is a bijection from `k=1` edges of one sheet to `k=2` edges of the other,
the exact coupling of the sheets through the tree; the budget
`sum q_i = 3n_0` is exhausted by the three cycles (inherited from the
glued-lines note). Basins to `10^7`: no escape.

**Paley `T_7`, two Fano planes and the octonions (PROVED).** A regular
tournament on `n` vertices has `(n^3-n)/24` cyclic triples, so the Paley
tournament on `F_7` has `14`, not the pasted `21` (its arc count; `21` is
also a value the Hamiltonian-path count omits, canon
`THM-1370-h-spectrum-omits-7-21-all-n`, an ID collision with
`THM-1370-elliptic-...`). The `14` cyclic triples are exactly the two
disjoint cyclic Steiner triple systems `dev{0,1,3}` and `dev{0,1,5}`; on
the first the Paley orientation is `s->s+1->s+3->s`, the octonion rule
`e_r e_(r+1)=e_(r+3)`, and the algebra it defines is alternative and a
composition algebra (checked by construction); the second plane carries
the opposite table under `x->-x`. For every prime `p = 3 mod 4` the cyclic
triples of the Paley tournament form a `2-(p,3,(p+1)/4)` design by
arc-transitivity; `p=7` is the only case where it splits into two Steiner
systems. `h(T_7)=189` is the Hamiltonian-path count already recorded in
the braids-two Fano notes; the `2640` labelled regular tournaments on
seven vertices split `240/720/1680` by automorphism group order
`21/7/3` with `h=189/175/171`.

**The block matrix (REFUTED).** The pasted block layout is square only
for `m=n=1`; the consistent `(m+n+2)`-vertex reading with a one-directional
cross-link is block triangular, so its characteristic polynomial is
`x^2 chi_A chi_B` and a source or sink only appends zero eigenvalues.
Transitive blocks are nilpotent, but every non-transitive tournament block
has Perron root `rho >= 1` (all `32,994` labelled non-transitive
tournaments on at most six vertices), and more generally a nonnegative
integer matrix has `rho in {0} cup [1,infinity)`; one-directional bipartite
blocks are nilpotent and two-directional ones have real spectrum
`+-sqrt(mn)`, so "purely imaginary pairs" occur only for a signed skew
block; Perron--Frobenius refutes "`Re(lambda)<0` for every trajectory" for
every `0/1` matrix. Brauer--Gentry's bounds `Re(lambda) >= -1/2` and
`|lambda| <= (n-1)/2` are re-proved in one line each and every non-Perron
eigenvalue of a regular tournament has real part exactly `-1/2`.
Kuratowski minors appear trivially (`K_{m+n}`) and carry no map to any
integer orbit (SCOPE).

**Compression and Lean (PROVED/REFUTED).** `8/pi^2` is the density of
coprime pairs among odd pairs (exact count `81,058,757/10^8`), its binary
entropy is `0.7002786534` bits (not `0.704`, and not zero); no lossless
code of arbitrary streams beats one bit per bit (pigeonhole), a suffix
lookahead is a bijection, and the tournament layout costs `C(n,2)` bits
for `T_(n-2)` payload bits; consecutive odd Collatz values are always
coprime, so a "coprime mask" on orbit edges has density `1`. The pasted
Lean fails at its first line in the repo's toolchain (no Mathlib build
exists), and with Mathlib names stubbed every remaining error is the
paste's own (undefined identifiers, an unbound `n`, a one-point
`FanoPlane`, an equivalence between a structure and a set, four
`sorry`s). The audited package `CollatzBlueprintAudit` builds with
`35` theorems and axioms at most `propext`/`Quot.sound`; a correct minimal
descent-certificate statement typechecks in core Lean with `decide`-closed
witnesses `n=3` (`t=6,K=4,L=2,B=5`) and `n=7`.

## 13. Wave six (2026-09-22): the summand holes `{1,4,6}`, the square-filtered summand graph, Tao on the minus sheet, and the root asymmetry

| Lane | Note | Script |
|---|---|---|
| Summand closure and the square filter | [summand_closure_square_filter](collatz_mod6_20260922_w6_summand_closure_square_filter.md) | `collatz_mod6_20260922_w6_summand_closure_square_filter.py` |
| Square-sum Hamiltonicity | [square_sum_hamiltonicity](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md) | `collatz_mod6_20260922_w6_square_sum_hamiltonicity.py` |
| Tao on the minus sheet, the `Delta=4` ladder | [tao_minus_sheet_ladder](collatz_mod6_20260922_w6_tao_minus_sheet_ladder.md) | `collatz_mod6_20260922_w6_tao_minus_sheet_ladder.py` |
| Sign-specific probes | [sign_specific_probes](collatz_mod6_20260922_w6_sign_specific_probes.md) | `collatz_mod6_20260922_w6_sign_specific_probes.py` |
| Lean paste audit | [lean_paste_audit_w6](collatz_mod6_20260922_w6_lean_paste_audit_w6.md) | `collatz_mod6_20260922_w6_lean_paste_audit_w6.py` |

**Two different "threes" (PROVED).** THM-2422's `{1,4,6}` is the closed
hole module of the distinct-summand closure of the seeds `{2,3}`: three
missing elements, not chains. The square-filtered summand graph `Q_n`
(`x~y` iff `x+y` is a square, `x != y`) really does have three chains for
`4<=n<=12`, but they are founded by the vertices born isolated, `1`, `2`
and `4`: `{1,3,6,8,10}`, `{2,7,9}`, `{4,5,11,12}`, the orbits of the
reflections `x->4-x`, `9-x`, `16-x`; `Q_n` is a linear forest exactly for
`n<=12`, vertex `13` merges the chains of `1` and `4` (squares `16`, `25`)
and vertex `14` absorbs `{2,7,9}`, and `Q_n` is connected for every
`n>=14` because for `n>=6` the square `(floor(sqrt n)+1)^2` lies strictly
inside `(n,2n)` (as `(1+sqrt 2)^2 = 5.83`), so the newest vertex always
attaches. The degree law is
`deg_n(m) = floor(sqrt(m+n)) - floor(sqrt m) - [2m square]`; `Q_n` has a
leaf iff `3<=n<=30`, the threshold set solely by vertex `18` (`36=2*18` is
the excluded diagonal, and the square `49` reaches `18` only at `n=31`).
There is no dyadic or fractal law in `Q_n`: the `2j^2` leaf scar exists
only for `j<=3`, the born-isolated and born-leaf sets are finite, degrees
grow like `(sqrt2-1)sqrt n`. And the two threes cannot be identified:
no additive target filter `T` has born-isolated set `{1,4,6}` (founders
`4` and `6` would force `5` to be a founder; exhaustive over all `2^12`
subsets of `{2..13}`), so the closure holes are intrinsically a
labelled-fibre phenomenon. The Collatz target `(3n+1)/2` is never a
square (`2k^2 = 1 mod 3` is impossible), so the square filter and the
Collatz arrow select disjoint targets: the analogy preserves only the
additive hyperedge with its excluded diagonal (SCOPE for any consequence).

**Square-sum Hamiltonicity (FINITE-EXACT and PROVED obstructions).** Path
existence for `1<=n<=40` is exactly `n in {1,15,16,17,23} cup [25,40]`,
path counts up to reversal agree with OEIS A090460 term by term, cycles are
absent at `30,31` and present from `32`; the `Q_15` path is unique, forced
by its eleven degree-two vertices with leaves `8` and `9`, and omits exactly
the edge `{1,3}` of square `4` (the pasted "avoid `4`" is refuted: `4` is
interior with degree two). `Q_18` fails by its three leaves `{16,17,18}`;
`n=19..22` fail by one-round forced-edge certificates and `n=24` by a
two-round one, all written as proofs from the two general obstructions
(three leaves; a cut set `S` with more than `|S|+1` components). OEIS
A090461 records the "all `k>=25`" conjecture as proved with Hamiltonian
cycles for `k>=32` (Gerbicz 2018, CITED via the OEIS comments and the
archived Mersenneforum thread: a 49-fold blow-up of "nice pairs" of
chains); the simpler self-similar step is re-proved here by a finite
junction check: any chain of `Q_n` from `1` to `3` (`n` odd) blows up to a
Hamiltonian cycle of `Q_(25n+12)` from `1` to `3`, which from the base
`n=35` gives the family `(71*25^m-1)/2` (verified at `887`, `22187`,
`554687`). This is the honest form of the user's "microcosm generating the
macrocosm": a finite base plus an explicit self-similar embedding by an
odd square. The pasted horizon table is corrected entry by entry: component
counts `3,2,1` at `12,13,14`; the `18--22` failures are leaf and cut
obstructions, not parity; `23` succeeds because `18` can be an endpoint
next to `7`.

**Tao on the minus sheet (PROVED from CITED inputs; the full theorem stays
OPEN).** Tao's paper (arXiv:1909.03562, fetched and text-converted) makes
no remark about `3n-1`; the pasted "dyadic martingale renewal process" and
"entropy decrement" descriptions are refuted (the paper's renewal process is
two-dimensional and serves only the Fourier decay). What is proved here:
the minus-sheet iteration is `Syr_-^n(N) = 3^n 2^(-|a|) N - F_n(a)` with
Tao's own offset polynomial and the identical valuation-word law, hence
`Syrac_-(Z/3^nZ) = -Syrac_+(Z/3^nZ)` and Propositions 1.14 and 1.17
(fine-scale mixing, characteristic-function decay) hold verbatim on
`3n-1`; whether the whole almost-all theorem transfers is OPEN (UPDATE 2026-09-22: now CITED-PROVED by Gonçalves--Greenfeld--Madrid 2025, arXiv 2111.06170, Thm 1.3, which covers `3N-1`; see [barrier atlas](collatz_procgen_20260922_barrier_atlas.md)) (Sections 3
and 5 of the paper were not re-run with the sign flipped), and "almost-all
bounded values implies a single root" is refuted by the three minus
basins (odd `n<=10^6`: `0.327/0.324/0.349`). The "`Delta=4` prime ladder"
`3,7,11,17` has differences `4,4,6`, and primes in arithmetic progression
with difference `4` have length at most three; "`3` is a first nontrivial
prime target" is false since `3` has no Syracuse preimage on either sheet.

**The root asymmetry (PROVED).** `1` is a fixed point of `(an+b)/2^k` iff
`a+b` is a power of two; in the summand reading the plus arrow
`n -> n+(n+1)/2` has companion `(n+1)/2`, which fixes `1` (the excluded
diagonal `1+1`), while the minus companion `(n-1)/2` sends `1` to `0`,
which is not a positive summand. This asymmetry is confined to `n=1` and
is reading-dependent; the count of `(3,b,k)` edges with hypotenuse `<=X`
has the sheet-blind main term `0.507819 sqrt X` and the sheets differ by at
most `2 log_2(3 sqrt(2X)+1)`. The carry is an odd word function times `b`
(`v_2(B_L)=0` always, `B_L(w,-1)=-B_L(w,+1)`), the `E`-graph picture is
sheet-blind (`Q1_-` and `Q2_-` hold to `10^6`, the minus cycles being
escaped through even arrows), and the general **word-function theorem**
says why: any invariant that is a function of the parameter `b` and the
parity word is sheet-blind, so a sign-specific invariant must be an
*order* statement. The one such statement found is the residue-of-minimum
law on the known cycles (plus-sheet cycle minimum `3 mod 4` and maximum
`1 mod 4`, minus sheet the reverse), recorded as an observation on four
cycles, not a theorem.

**Lean (core Lean, `decide`, no Mathlib, no `sorry`).** The pasted
`square_sum_graph` is not loopless (values `2, 8, 18, 32` are self-adjacent);
`Q_14` connected and `Q_12`, `Q_13` disconnected with component counts
`3,2,1,1,1,1` for `n=12..17` are machine-checked; the pasted
`delta_four` statement is true but `15` and `21` are composite; the
pasted `SignSpecificCertificate` is inhabited for every `n>=1` (`L=0,
K=1, B=0`) and mentions no orbit, so it guards nothing.

## Reproduction and audit scope

Each lane script is run as `python3 04-computation/experiments/<script> >
05-knowledge/results/<lane>.out`; normal and `-O` streams agree modulo
timing stamps. Independent audit scripts exist for every lane
(`*_audit*.py`) with their own outputs. Universes, filters, positive and
hostile controls are stated inside each note. Session lanes were run by
subagents in three waves; the worktree was pruned between the two working
days and the scripts were recovered from the agent transcripts and rerun,
which is why some outputs are dated 2026-09-21.
