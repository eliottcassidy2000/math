# Arithmetic free bits: the dual check, the erased circuit, and the prime guard

**Session85 — 2026-10-28. Status: PROVED elementary statements + FINITE-EXACT realizations.**
The coding-theory mechanisms are classical, not new discoveries. New work here
is the explicit arithmetic enlargement, sharp guarded-continuation comparison,
and repair of the marginal/sign boundary. No unbounded prime-supply theorem is
claimed.

## Results first

1. **A circuit does not give the recovery equation.** In a parity-check
   representation, a primal circuit is a smallest *ambiguity*: toggling all its
   bits preserves the syndrome. Recovery equations live in the **dual code**.
   This explains how to get a bit from other bits without restricting the
   mechanism to ordinary graph cycles.
2. The seven-prime bank extends explicitly to **Hamming [15,11,3]**. Every simple
   binary matroid of rank at most four is a restriction of this one arithmetic
   bank. No conjecture is used for this finite statement.
3. There is also a different, explicit sixteen-prime bank whose square-sigma
   selections are **extended Hamming [16,11,4]**, with **five** independent checks
   and dual **RM(1,4)**. This is a precise 5/15 connection through an affine lift,
   not an inference from matching integers.
4. Perfect one-bit error correction is not perfect one-prime multiplication.
   Hamming rank `r` needs at most `r` legal additions, and this is sharp precisely
   when the unused primes have basis columns.
5. Hamming7 cosets agree on selected-bit marginals through order3, and Hamming15
   through order7. Adding the separate equal sign law to the same *joint*
   observation loses one level: joint sign+bit equality is guaranteed only
   through orders2 and6. The first separating observations are explicit below.
6. Within Hamming7, a different target—the overall selection-parity bit—is
   recoverable from **any Fano line**, never from two coordinates, and from
   **every five-coordinate set**. Four-point caps fail. This is an exact7/5
   access structure, not an ordinary threshold code.

## Inheritance and scope

Read-only source pin:
`83b601158076669edc8729ce461f150da1c3877e`, especially
`05-knowledge/results/incidence_attractors_20261028/THEORY.md`, A3, and its
`AUDIT.md`. The source checkout was preserved; this lane writes only its own
scratch directory.

Closest mechanism: square-class parity of multiplicative sigma on pairwise
coprime atoms. Canonical hostile: selecting both2 and4 invalidates the formal
factorwise sigma calculation. Corrected near miss: Hamming7 has **14 minimal
circuits**, not merely its seven lines; equal small-bit marginals and a separate
sign law do not imply equal joint laws. Least-used relevant sidecar: the dual
check space and a specified set of erased coordinates.

The live concept board is: arithmetic square classes; primal/dual matroids;
erasure fibres; guarded addition versus toggling; affine extension; low-order
joint observations. The map preserves the square/non-square sigma predicate,
not integer sizes, sigma magnitudes, or aliquot trajectories.

## 1. Exact criterion for a genuinely redundant bit

Let `H` have columns `h_i` in `F2^r`, rank `r`, and assume the state is promised
to satisfy `Hx=t`. Write `E` for erased positions and `O` for observed ones.
The residual equation is

    H_E x_E = t + H_O x_O.                                      (1)

If the observed values come from a valid state, its solution set is an affine
space with exactly

    2^(|E| - rank H_E)                                          (2)

points. Consequently:

* **All erased bits are recoverable iff the erased columns are independent.**
* An individual erased bit `j` is recoverable iff
  `h_j not in span{h_i : i in E\{j}}`.
* Equivalently, a row functional `a` exists with `a·h_j=1` and
  `a·h_i=0` for all other erased coordinates. It supplies

      x_j = a·t + sum_{i in O} (a·h_i) x_i.                    (3)

Proof: subtract two solutions. Ambiguities are exactly `ker H_E`. A kernel
vector changes coordinate `j` iff `h_j` is a sum of the other erased columns.
When it is not, extend a basis of their span by `h_j` and define a separating
linear functional. Rank-nullity proves (2).

For a uniform coset state, the conditional Shannon entropy is exactly
`H(X_E | X_O) = |E| - rank H_E` bits. Thus the constraints remove
`rank H_E` bits of freedom; this does **not** mean that every individual erased
bit is determined. For example, a three-column circuit has nullity1 and all
three of its coordinates are ambiguous.

The promise is load-bearing. If membership `Hx=t` is itself unknown, assuming
it to infer a bit would be circular. There are `n-r` independent data bits,
not `n` independent bits plus `r` newly created information bits.

### Primal circuits versus dual recovery checks

A nonempty minimal support `D` with `sum_{i in D} h_i=0` is a **primal circuit**.
Its indicator is a codeword: it produces two indistinguishable completions
whenever all of `D` is erased. A **dual** word `aH`, instead, gives a valid
parity equation. Minimal nonzero dual supports are cocircuits of the column
matroid. They are different objects, even when a graphic picture blurs them.

Literal hostile: columns1,2,3 of the seven-prime bank are primes7,2,5.
Their product70 has sigma144. Their three selected bits are111, whose parity
is1. Hence their primal circuit emphatically does not assert an even parity
check on those three bits.

### A different target: one global bit recoverable from any Fano line

On Hamming7, let `ell(x)=sum_i x_i` be overall selection parity, equivalently
the sign bit of `mu=lambda`. Every line's complement is a dual weight-four
support. Therefore

    ell(x) = sum_{i in L} x_i for every Fano line L.           (4a)

This does not contradict the circuit warning: line parity is the **variable
secret bit**, not a fixed-zero check. In fact the even subcode is precisely
`C^perp`, so `ell` identifies the one-dimensional quotient `C/C^perp`.

A coordinate set `B` recovers `ell` **iff B contains a Fano line**. A linear
functional on a code is determined by its `B`-projection iff it differs from
a functional supported on `B` by a dual word. Here the possible supports of
`1+d`, for `d in C^perp`, are the seven lines and the full seven-set, proving
the criterion. For every unauthorized `B`, each nonempty observation fibre
has equally many secret-zero and secret-one states: its remaining linear
freedom contains a parity-flipping direction. Thus failure means zero
information about this uniform secret, not merely imperfect reconstruction.

There are7 authorized triples,28 authorized four-sets, and all21 five-sets,
all7 six-sets, and the seven-set. The seven failing four-sets are the four-point
caps, the complements of lines. Every five-set contains a line: a line-free
set `B` and its translate `a+B` by any `a in B` are disjoint subsets of `F2^3`,
so `2|B|<=8`. This proves a genuine five-of-seven **guarantee**, while some
triples already succeed; it is not a pure five-threshold access structure.

## 2. The seven-prime bank, correctly decoded

Use square-class coordinates `(2,3,5)` and column order1 through7:

    columns: 1  2  3   4   5   6   7
    primes:  7  2  5  19  89  59  29.

For a squarefree product of the selected primes, sigma is square iff `Hx=0`.
The code is Hamming [7,4,3], with enumerator `1+7z^3+7z^4+z^7`.
All seven weight-three and all seven weight-four supports are minimal circuits:
**14 total**. A weight-four zero support cannot contain a weight-three zero
support, since its one remaining nonzero column would then be zero.

The dual has seven nonzero words, each of weight4. Thus every bit can be
recovered from a suitable three others, and no two other bits suffice. For
example, on the square-sigma code:

    x_7 = x_5 + x_89 + x_29.                                  (4)

Here subscripts denote actual prime labels; addition is XOR. Taking primes
`5,89,59,29` as four free data bits, the other two formulas are

    x_2  = x_5  + x_59 + x_29,
    x_19 = x_89 + x_59 + x_29.

This gives a systematic arithmetic encoder. It does not say every arbitrary
triple of observed coordinates repairs a chosen bit.

All one- and two-erasure sets are recoverable. Of the35 three-erasure sets,
28 are recoverable and seven Fano lines are not. No four-erasure set is uniquely
recoverable. This is independence in the non-graphic Fano matroid, not an
ordinary graph-cycle reconstruction theorem.

## 3. An explicit Hamming [15,11,3] prime realization

Set the square-class basis to `(2,3,5,7)`. For mask `v`, let `d_v` be the product
of its selected basis primes. The following are exact identities
`p_v+1=d_v k_v^2`:

| v | d_v | k_v | p_v |
|---:|---:|---:|---:|
|1|2|2|7|
|2|3|1|2|
|3|6|1|5|
|4|5|2|19|
|5|10|3|89|
|6|15|2|59|
|7|30|1|29|
|8|7|6|251|
|9|14|1|13|
|10|21|2|83|
|11|42|1|41|
|12|35|2|139|
|13|70|9|5669|
|14|105|2|419|
|15|210|2|839|

Every displayed `p_v` is prime, distinct, and the least realizing prime for
that **fixed square class**. Exhausting `1<=k<k_v` proves the last assertion;
no optimum over changing bases or representations is claimed.

The column matrix consists of all nonzero vectors in `F2^4`, so has rank4 and
kernel dimension11. There are exactly2048 square-sigma selections among32768
subsets. The minimum weight is3: one and two distinct nonzero columns cannot
sum to zero, while `{a,b,a+b}` always does.

Every simple binary matroid of rank at most4 can be represented by a subset of
these columns after embedding its representation in `F2^4`. Restricting the
prime bank therefore realizes it exactly. This gives a finite family beyond
graphic matroids, with no infinitude or primality conjecture.

The minimal circuits have sizes3,4,5 and counts35,105,168, respectively:
**308 circuits**. In general, in the full nonzero `F2^r` column bank,

    number of k-circuits = product_{i=0}^{k-2}(2^r-2^i) / k!,
    3 <= k <= r+1.                                           (5)

Choose an ordered independent `(k-1)`-tuple and adjoin its sum. Each circuit
is counted `k!` times by the omitted element and order of the others.

The recoverable erasure counts are

| number erased | Hamming7 | Hamming15 |
|---:|---:|---:|
|0|1|1|
|1|7|15|
|2|21|105|
|3|28 of35|420 of455|
|4|0 of35|840 of1365|
|5 or more|0|0|

The dual simplex code has15 nonzero weight-eight checks. Every individual bit
has minimum repair locality7, versus locality3 in Hamming7. Greater dimension
does not automatically make individual local repair cheaper.

## 4. A sharp arithmetic guard price for “one correction”

In the full nonzero `F2^r` bank, `r>=2`, every syndrome `s!=0` names exactly one
column. Toggling that coordinate gives the unique distance-one square-sigma
selection. If its prime is absent, this is multiplication; if present, it is
division. A code-theoretic bit flip must not silently become an always-legal
coprime multiplication.

Let `U` be the unused columns. The XOR of the full bank is0, hence the current
syndrome is `s=XOR(U)`. A legal addition chooses `T subset U` with `XOR(T)=s`.
Choose a basis of `span(U)`: it represents `s` using at most `rank(U)<=r`
columns. Thus **at most r legal prime additions always suffice**.

Moreover, the minimum equals `r` **iff U is a basis of F2^r**. The forward
direction needs a small argument. If `rank(U)<r`, the preceding bound is
strict. Otherwise choose a basis `B subset U`. A representation requiring
all `r` basis vectors is the only one not already shorter. If `U` has an extra
vector `v`, it is the sum of at least two basis vectors, so replace those by
`v` and shorten the representation. Therefore maximal cost requires `U=B`.
Conversely, XOR of a basis requires every one of its vectors.

The proof uses a simple nonzero column bank; duplicate or zero columns need
separate treatment for the “iff” characterization.

Exact cost distributions, over all source subsets:

| least legal additions | Hamming7 | Hamming15 |
|---:|---:|---:|
|0|16|2048|
|1|56|15360|
|2|28|13140|
|3|28|1380|
|4|0|840|

Concrete Hamming7 hostile: source
`n=5*89*59*29=761395` leaves the basis primes`7,2,19` unused. Its sigma is972000,
and exactly three additions are necessary. Yet removing prime29 produces26255,
whose sigma is32400=`180^2`. One toggle, three legal additions: different
operation languages, different cost.

## 5. Sharp marginal tower, including the sign-law correction

Assume **r>=2**. Let `C_t={x:Hx=t}` carry the uniform law for the full Hamming
bank of length `n=2^r-1`. Every nonzero dual word has weight `w=2^(r-1)`. Consequently all
cosets have the same uniform coordinate marginals through order `w-1`.

Direct proof: after fixing fewer than `w` selected coordinates, the remaining
columns still span `F2^r`. Otherwise a nonzero functional would vanish on all
of them, with its weight-`w` support contained in fewer than `w` positions.
Every selected assignment therefore has the same `2^(n-r-|S|)` completions.
At `w` positions, a dual support gives a parity relation that distinguishes
two cosets with different corresponding syndrome components. The bound is
sharp: orders3/4 for Hamming7, orders7/8 for Hamming15.

For squarefree sources, `mu=lambda=(-1)^(sum x_i)`. Each coset has an equal
number of both signs: translation by any weight-three codeword reverses parity.
This is a statement about the **separate** `(mu,lambda)` law.

The lower boundary is deliberate: at `r=2`, the code is `{000,111}` and the
sign law is balanced in every coset. At `r=1`, the two cosets are singletons
with opposite, deterministic signs; the sign-law assertion is false there.
The verifier explicitly checks both boundary cases.

The joint limit is different. Take a nonzero row functional `a` with `a·t=1`,
let `D=support(aH)`, and observe `S=D^c`, which has `w-1` coordinates. Then

    parity(all selected bits) XOR parity(bits in S) = a·t.     (6)

Thus sign together with these `w-1` selected bits perfectly distinguishes
`C_0` and `C_t`, even though the selected bits alone and sign alone agree.

The joint laws do agree for every `|S|<=w-2`. To prove this, use binary
characters. A character of sign plus the selected bits has exponent
`epsilon*1 + v`, with `support(v) subset S`. Its expectation on a uniform coset
vanishes unless the exponent is a dual word. A nonzero dual word needs `w`
positions; its complement needs `w-1`. Neither is available below `w-1`,
apart from the trivial character. This proves the sharp joint threshold.

For the seven-prime bank and `a=t=1`, `D={1,3,5,7}` in column labels, so the
three observed primes are`2,19,59`. For the fifteen-prime bank the seven observed
primes are the even-column labels`2,19,59,251,83,139,419`. These are explicit
hostiles to upgrading separate marginal agreements into joint agreement.

## 6. Why five checks arise from fifteen nonzero four-bit labels

The fifteen columns above are `F2^4 \ {0}`. Adjoin0 and lift every point `v` to
the affine column `(1,v)` in `F2^5`. This matrix has16 columns and rank5, giving
extended Hamming [16,11,4]. Its dual consists of all evaluations

    v -> a_0 + a_1 v_1 + a_2 v_2 + a_3 v_3 + a_4 v_4,

so is RM(1,4): five coefficients, dimension5, with30 weight-eight words and
one weight-sixteen word. The constant check imposes even cardinality. No
weight-one, -two, or -three word is possible, and an affine plane of four
points supplies weight4.

There are two genuinely different arithmetic realizations.

### 6a. Keep the fifteen primes, but retain the sign guard

Adjoin prime3, whose sigma4 is square and hence has zero column. Square sigma
alone then gives `C_15 x F2`, a [16,12,1] code with4096 words, **not** extended
Hamming. Imposing even selected cardinality, equivalently `mu=lambda=+1`,
cuts it to the2048-word extended code. Prime3 is the overall-parity bit.

### 6b. A new pure-square-sigma prime bank, without an extra sign observer

Use square-class basis `(2,3,5,7,11)`, with
`d_v=2*product of the selected members of (3,5,7,11)` for `v in F2^4`.
Its column is exactly `(1,v)`; the table supplies `p_v+1=d_v k_v^2`:

| v | d_v | k_v | p_v |
|---:|---:|---:|---:|
|0|2|2|7|
|1|6|1|5|
|2|10|3|89|
|3|30|1|29|
|4|14|1|13|
|5|42|1|41|
|6|70|9|5669|
|7|210|2|839|
|8|22|3|197|
|9|66|2|263|
|10|110|1|109|
|11|330|2|1319|
|12|154|9|12473|
|13|462|1|461|
|14|770|1|769|
|15|2310|1|2309|

Again, each is the least prime in its specified square class. Since each
`v_2(p_v+1)` is odd, square sigma automatically forces an even number of
selected primes. The extra parity is built into the arithmetic carrier rather
than supplied by a separate observer. There are2048 square-sigma selections
among65536 subsets, with weight enumerator

    1+140z^4+448z^6+870z^8+448z^10+140z^12+z^16.

All sets of at most three erasures recover. Among four- and five-erasure sets,
1680 and2688 recover, respectively; more than five cannot. Minimum individual
repair locality is7. Its minimal circuits number588:140 of size4 and448 of
size6. The size-four circuits are the affine planes in `F2^4`.

The integer15 also counts the nonzero words of the dimension-four Hamming7
code, but that is a different role from15 columns of its rank-four successor.
The explicit maps above, not the bare integers, license the connection.

## 7. What general prime realization would require

Given distinct coordinate primes `q_1,...,q_r`, prescribe a nonzero vector `v`
by squarefree `d_v=product q_i^(v_i)>1`. A realizing prime must have

    p=d_v k^2-1.                                             (7)

For fixed `d_v`, this polynomial has positive leading coefficient, is
irreducible over the rationals since `d_v` is nonsquare, and has no fixed prime
divisor: at `k=0` its residue is-1 modulo every prime. The Bunyakovsky conjecture
would give infinitely many realizing primes for every such class. Conditional
on that conjecture, every finite loopless binary matroid admits a realization
by distinct primes, including repeated parallel columns.

This is **CONDITIONAL**, not a consequence of the finite tables or of Dirichlet's
theorem. No simultaneous-prime-values conjecture is needed: the finitely many
columns may use different `k` values. The primary discussion in
[Jones–Zvonkin, §3](https://arxiv.org/html/2105.03915) states the relevant
polynomial-primality conditions and distinguishes computations from conjectures.

There is an unconditional obstruction even before such conjectures: a zero
column requires `p+1=k^2`, so `(k-1)(k+1)=p`. The only prime possibility is
`k=2,p=3`. Hence **two loops cannot be realized by two distinct prime atoms**.
Under Bunyakovsky, every finite binary matroid with at most one loop is
realizable; with two loops it is not. This statement concerns distinct prime
atoms, not arbitrary coprime composite or prime-power atoms.

For now the unconditional claim is the explicit finite rank-four simple bank,
the affine sixteen-column bank, and any of their restrictions. Quotienting a
square-class coordinate is not automatically arithmetic matroid contraction:
it forgets a sigma valuation-parity condition and needs a new predicate or
realization.

## 8. Reproduction, independent path, and limits

From the directory containing the published `check_codes.py`:

```sh
python3 check_codes.py > /tmp/codes-normal.json
python3 -O check_codes.py > /tmp/codes-optimized.json
cmp /tmp/codes-normal.json /tmp/codes-optimized.json
```

Both modes pass and yield byte-identical JSON receipts. All98432 subsets of
the three banks are checked arithmetically by enumerating their divisor sums
with a subset-zeta transform and testing integer squares. This path does not
derive the answer from the parity matrix. Every prime is checked by complete
trial division; every claimed least realization exhausts the smaller `k`.

All erasure sets are checked by two implementations: rank, versus zeta counts
of supported codewords. Their support unions independently identify ambiguous
individual bits. The circuit counts, dual distributions, repair locality,
uniform-projection rank criterion, sharp sign hostile, and legal-addition cost
tables are also checked. No removable `assert`, random sample, third-party
package, or inherited producer implementation is used.
Every Fano secret-access set and every associated observation fibre is checked
as well, including perfect privacy of all unauthorized sets.

This is same-agent mathematical work with independently implemented checks,
not external peer review. The elementary code theory belongs to the classical
tradition initiated by [Hamming's 1950 paper](https://zoo.cs.yale.edu/classes/cs323/doc/Hamming.pdf)
([publisher DOI](https://doi.org/10.1002/j.1538-7305.1950.tb00463.x)); the original
paper's indexed primary text was available, but the direct PDF open timed out.
The proofs needed here are written out rather than delegated to that source.
No novelty priority claim is made for the general coding-theory facts.

Nothing here proves unrestricted prime realization, a graph or tournament
composition theorem, aliquot termination, or LRC(14).
