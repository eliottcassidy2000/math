# One fixed integer: exact compatibility, finite support, and endpoint obligations

Status: **PROVED SCOPED RESULTS + EXACT REDUCTIONS + FINITE-EXACT controls;
independently audited.** Geometry and automata lanes independently accepted
Sections 3--6; the automata lane additionally audited Section 4A and
replayed the final script. No Collatz convergence claim. All residue and
landing computations use integers.

## 1. Inheritance and live concept board

The closest successful mechanism is the separation of an exact integer map
from the asymptotic certificate in
[THM-4500, pairing Bellman contraction](../../01-canon/theorems/THM-4500-pairing-bellman-contraction-and-natural-density.md).
Its contracting object is an auxiliary Haar-L1 potential, not the Collatz
orbit. Actual global realizability needs a separate reset/phase argument.

The canonical hostile is
[THM-4471, the fixed-point proof refutation](../../01-canon/theorems/THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md):
the misplaced pairwise alternative also certifies a successor map and the
positive `3n-1` cycle. A correct local inequality need not have the quantifier
needed to iterate it. The corrected near miss is
[THM-4482, ranks are tensions](../../01-canon/theorems/THM-4482-ranks-are-tensions-strategy-cube.md)
and [THM-4483, forced charges](../../01-canon/theorems/THM-4483-forced-charges-two-place-ranks.md):
finite or height-summable counter banks fail; allowing enough private height
can merely encode stopping times. The underused sidecar here is the monotone
least positive-height representative of a compatible 2-adic address.

[THM-4476, thin divergence](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md),
[THM-4495, no-descent Spitzer count](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md),
and [THM-4499, little-o thinness](../../01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md)
bound one hypothetical nonperiodic orbit without excluding it. Their
nonperiodicity, injectivity, endpoint and landing assumptions remain active.

| lane | object / map | preserved predicate | lost information / decisive test |
|---|---|---|---|
| anchor | parity prefixes to least residues | exact finite arithmetic | finite support / test bounded representatives |
| anchor | fixed start to its affine orbit | every exact successor | termination / positive hostile cycles |
| niche | finite bounded-height prefix tree | one common integer witness | a moving height bound / nested compactness |
| niche | dippers to first landing | actual orbit incidence | endpoints / build a critical hover then drop |
| wildcard | uniform poset order to height-selected words | selected permutations | uniform-poset law / two-extreme selection |

The last row is audited independently in the [width lane, Section 4](crossroads_poset_20260926_width.md): at length 6 and
five odd steps, the four expanding words `110111,111011,111101,111110`
have start residues `27,39,47,31`. Cutting starts at 31 selects the extreme
two insertion positions. Their common order relations still permit all four,
so the selected law is not the uniform extension law of another poset on
the same events. Arithmetic selection needs its own sidecar.

## 2. Assumption ledger for one fixed integer

Use the ordinary half-step map `T(n)=n/2` for even n and `(3n+1)/2` for odd n.
Do not make 1 absorbing in parity-word algebra; its ordinary orbit is `1,2,1`.

| statement | status and exact consequence | first invalid upgrade |
|---|---|---|
| one integer n has finite binary input | TRUE | future orbit values have uniformly bounded bit length |
| T is a total computable function | TRUE | the search for a visit to 1 terminates |
| every finite iterate T^k(n) is exactly computable | TRUE | one finite computation decides all future iterates |
| one fixed n determines one infinite parity word | TRUE | that word is statistically typical |
| every finite binary parity word has a unique residue class | TRUE | every infinite word represents a positive integer |
| compatible prefixes determine one 2-adic start | TRUE | its ordinary integer height is finite |
| least residues are nondecreasing | TRUE | they stabilize without a boundedness argument |
| least residues stabilize for a fixed positive n | TRUE, automatically once 2^k>n | its orbit descends; even a nonterminating computable map can have this property |
| a computable infinite address has a short generating program | TRUE | the address has finitely many nonzero binary digits |
| each depth has a bounded representative | TRUE with a depth-dependent bound 2^k | one common bound works at all depths |
| every finite horizon has some successful witness | quantified statement `forall k exists n_k` | one n works at every horizon |
| a positive orbit is bounded by H | implies eventual periodicity | its eventual cycle is necessarily 1,2 |
| a rank is lower bounded and decreases strictly | not sufficient without a well-founded range or a fixed decrement | strict real decreases cannot occur infinitely often |
| an orbit has finite reciprocal sum | compatible with sparse divergence | it is finite or cannot consist of integers |
| an averaged residue law contracts | says something about that law | each deterministic integer trajectory contracts |
| a finite orbit segment has high landing multiplicity | actual finite witness | the same average persists on an entire nonperiodic orbit |

For the rank row, `R(j)=1/(j+1)` along the successor map is a strictly
decreasing nonnegative real rank. A nonnegative integer rank, or a fixed
positive decrease in a lower-bounded real rank, does force termination.
The missing requirement is well-foundedness, not arithmetic computability.

Finite input also does not mean finite-state evolution. Unconditionally,
`T(n)+1<=(3/2)(n+1)`, hence
`T^k(n)+1<=(3/2)^k(n+1)`. This is a computable height bound for each
specified finite horizon; it grows with that horizon. A single finite
program can therefore compute states requiring unboundedly many bits.
Replacing this bound by one independent of k is an additional dynamical
statement, not a consequence of knowing n exactly.

Positive hostile controls are the computable cycles `5,7,10,5` for `3n-1`
and `13,33,83,208,104,52,26,13` for `5n+1`. The latter has odd subsequence
`13,33,83,13`. Thus neither determinism, exact integer arithmetic, nor the
availability of arbitrarily many computed iterates supplies termination.

## 3. Exact stabilization criterion and its sharp height boundary

For compatible parity prefixes `w[0:k]`, let `r_k` be the unique least
representative in `[0,2^k)`. Set `r_0=0`. Prefix compatibility gives

```
r_(k+1)=r_k+e_k 2^k,   e_k in {0,1}.                (3.1)
```

Consequently the following are equivalent:

1. The whole branch has a nonnegative ordinary integer start.
2. The sequence `(r_k)` is bounded in the ordinary real order.
3. The sequence `(r_k)` is eventually constant.
4. Only finitely many lift bits `e_k` are one.
5. `limsup r_k/2^k < 1/2`.
6. `r_k=o(2^k)`.

Proof: (3.1) gives monotonicity and the equivalence of 2--4. Stabilization
at n gives `r_k=n mod 2^k` at every depth and hence the actual parity word
of n, proving 1; the converse is automatic. If infinitely many lift bits
are one, at those indices `r_(k+1)/2^(k+1)>=1/2`, contradicting 5 or 6.
If they are finite, the ratios tend to zero. To require a positive start,
add that some `r_k>0`; the all-zero branch represents zero.

The constant `1/2` is sharp. The computable 2-adic number

```
z=sum_(j>=1) 2^(j!)
```

has least residues unbounded but has `liminf r_k/2^k=0` and
`limsup r_k/2^k=1/2`. After each new factorial-position bit the ratio is
just above `1/2`; across the subsequent gap it tends to zero. Its parity
prefixes are computable by reducing z modulo a sufficiently high power
of two and iterating exact arithmetic. The limit is not an ordinary
integer. Computability and long quiet stretches of lift digits do not
force eventual stabilization. The simpler all-odd word represents `-1`,
with `r_k=2^k-1`; it is the inherited negative fixed-point shadow.

**What this fixes, and what it does not.** For one already fixed positive
n, stabilization is an identity, not an additional dynamical theorem:
`r_k=n` whenever `2^k>n`, regardless of what its forward orbit later does.
Thus one must prove that an infinite *bad* branch cannot stabilize, not
merely that an arbitrary compatible branch has a 2-adic limit. After the
input's binary length is exhausted, future lift bits are forced zero;
future parity letters remain determined by the changing integer state.
Treating those letters as fresh random choices loses the fixed integer.

## 4. Bounded-height compactness and an exact pruning rule

Suppose a finitely branching prefix tree describes a prefix-closed
property of integer orbit words. Fix one finite integer H. If at every
depth there is a word in the tree realized by some `1<=n<=H`, then some
single `n<=H` has all its prefixes in the tree. Indeed the finite sets
of eligible starts in `[1,H]` are nonempty and nested, so their
intersection is nonempty. No topological theorem is needed. Allowing H
to depend on depth destroys this conclusion, as the all-odd prefixes
with least starts `2^k-1` show.

The correct first-descent target is `forall n>1 exists k(n): T^k(n)<n`.
It is equivalent to Collatz by strong induction and does not ask for a
horizon independent of n. Equivalently, `forall H exists k(H)` such that
every `2<=n<=H` has descended by k(H). For a fixed H, the negation
`forall k exists n<=H` surviving through k is equivalent to `exists n<=H
forall k` surviving, by the finite nested-set argument. Removing the
common H invalidates the interchange. Mersenne starts refute a constant
lookahead horizon, not the possibility of a proof with n-dependent
certificates. A terminating search would itself compute k(n) if the
universal first-descent statement were already proved.

There is an exact height-aware carrier for first descent. For a parity
prefix write

```
T^j(n)=(A_j n+C_j)/2^j,  A_j=3^(number of odd steps),
C_0=0;
even extension: C_(j+1)=C_j;
odd extension:  C_(j+1)=3C_j+2^j.                   (4.1)
```

The starts that realize the word and do not descend below themselves
through its length k are exactly

```
n=r_k mod 2^k,  1<=n<=H_w,
H_w=min_(j<=k: A_j<2^j) floor(C_j/(2^j-A_j)),       (4.2)
```

where an empty minimum is infinity. Expanding prefixes need no upper
cutoff because `C_j>=0`; contracting prefixes impose the displayed
upper cutoff. This is an exact necessary-and-sufficient test, not a
drift approximation. It retains both the residue and the affine carry.
The script checks this carrier against actual trajectories in a finite
box. A proof that every bounded-start bad tree becomes empty would prove
first descent by strong induction; no bound establishing that emptiness
is supplied here.

Separately, if one fixed orbit is known to stay in `[1,H]`, repetition
occurs within H+1 visited states. This reduces to a finite cycle check.
For all starts, existence of a computable all-time height bound is
equivalent to absence of divergent orbits: one direction is immediate;
for the other, enumerate an orbit until its first repeat and report its
maximum. The assumption proves that this algorithm is total. It does
not classify the possible cycles and is not a weaker way to prove the
no-divergence half.

### 4A. A full-support weighted carrier replaces the null-set loophole

Fix any real `0<z<1` and define the exact actual no-descent set and its
discrete weighted mass by

```
S_k={n>=2: T^j(n)>=n for every 1<=j<=k},
Z_k(z)=sum_(n in S_k) z^n.                         (4.3)
```

These sets decrease with k. Every positive integer address has positive
mass, unlike a singleton under the Haar distribution on 2-adic inputs.
The uniform tail estimate

```
0<=Z_k(z)-sum_(n in S_k,n<=H) z^n<=z^(H+1)/(1-z)   (4.4)
```

holds for every k. Therefore

```
lim_k Z_k(z)=sum_(n in intersection_k S_k) z^n,
Collatz holds iff lim_k Z_k(z)=0 for one 0<z<1.     (4.5)
```

The same equivalence then holds for every such z. The limit identity
follows by first truncating at H, where all limits are finite sums, and
then using (4.4). Positivity of each mass proves that zero limit excludes
every bad fixed integer. This is an exact reduction, not progress toward
the missing decay estimate. A zero Haar density of bad addresses cannot
replace its premise: a surviving singleton contributes `z^n` forever.

At every finite depth, (4.2) computes Z_k exactly. For a word residue r
modulo `M=2^k`, let a be its least representative at least 2. With height
cap H_w, its contribution is zero if a>H_w; otherwise it is

```
z^a/(1-z^M)                                      if H_w=infinity,
z^a (1-z^(q M))/(1-z^M),
q=1+floor((H_w-a)/M)                              if H_w<infinity.
```

Distinct words give disjoint residue classes. Thus Z_k is a rational
function of z, exactly rational at rational z. The script checks this
formula against direct bounded-start orbit computations and (4.4).

The exact map back to the old density is

```
lim_(z increases to1) (1-z) Z_k(z)=W_k/2^k,         (4.6)
```

where W_k is the strictly positive prefix-multiplier count from
THM-4495. Each uncapped arithmetic progression contributes `1/2^k`;
each finite height-capped progression contributes zero. Thus the old
density is a boundary normalization of this carrier, and that operation
erases each fixed integer atom. Its known convergence to zero as k
increases does not establish (4.5). The missing information is retained
by fixing z strictly below one and keeping the mass unnormalized.

The cheapest hostile already rules out one-step strict contraction:
`Z_2(z)=Z_3(z)=z^3/(1-z^4)` for every `0<z<1`; at z=1/2 both equal 2/15.
After two steps precisely starts 3 modulo 4 remain. At the third step
their possible words are 110 and 111, neither of which descends. The
finite-height mass through address 27 also has a long exact plateau:
after all smaller starts have descended, its sole mass is `z^27` until
27's first descent at step 59. Hence an adaptive or block decrease
argument must retain the actual dormant mass; it cannot insist that
every clock tick consumes it.

A useful quantitative target would be an independently proved bound
`Z_k(z)<=E(k)` with E(k) tending to zero. Even an ineffective proof of
this limit would suffice. An exponential estimate is much stronger: it
would bound first-descent time by a constant times n plus a constant,
because a surviving start n forces `Z_k(z)>=z^n`. No such estimate, nor
a contraction for fixed blocks, is asserted here.

### 4B. The next 27, defined by first-descent records

For this finite census, a record start n has a larger first-descent time
`tau(n)=min{j>=1:T^j(n)<n}` than every smaller start at least 2. These are
records for the shortcut map, not records for maximum height or time to
reach 1. The script directly tested every `2<=n<=1,000,000`; no start
hit the explicit 10,000-step failure cap. The complete strict record list is

| start | first descent steps | shortcut steps to first 1 |
|---:|---:|---:|
| 2 | 1 | 1 |
| 3 | 4 | 5 |
| 7 | 7 | 11 |
| 27 | 59 | 70 |
| 703 | 81 | 108 |
| 10,087 | 105 | 142 |
| 35,655 | 135 | 204 |
| 270,271 | 164 | 256 |
| 362,343 | 165 | 228 |
| 381,727 | 173 | 236 |
| 626,331 | 176 | 319 |

Thus **703 is the next 27 under this explicit criterion**. None of the
record starts after 27 visits 27 or 41 before first reaching 1. They are
therefore distinct from the certified motif family whose members enter
41 after a prescribed initial block. No recurrence, limiting density,
or eventual pattern for these record starts is inferred from the census.
All displayed times use `(3n+1)/2` as one odd step; they must not be
confused with counts for the unshortened map.

The exact connection to (4.3) uses the smallest surviving integer
`m_k=min S_k`, with `m_k=infinity` if the set is empty and `z^infinity=0`:

```
z^(m_k)<=Z_k(z)<=z^(m_k)/(1-z).                    (4.7)
```

For the full integer universe S_k is in fact nonempty at every finite
depth: a sufficiently long all-odd Mersenne prefix supplies a start.
The infinity convention is useful for finite cutoffs. At z=1/2 the
weighted mass is within a factor two of `2^(-m_k)`, so this full-support
carrier cannot hide the smallest unresolved integer.

The first-descent records describe the exact plateaus of m_k: each new
record start is m_k from the previous record's descent time through one
less than its own. For example, `m_k=27` for `7<=k<=58`, `m_k=703` for
`59<=k<=80`, and `m_k=10087` for `81<=k<=104`. The census determines m_k
exactly through depth 175; at depth 176 it establishes only
`m_176>1,000,000`. The monotone escape `m_k->infinity` is equivalent to
the zero-mass limit, hence to Collatz. This is an exact restatement of
the outstanding obligation, not an asymptotic consequence of the finite
record table.

## 5. Actual segments refute endpoint-free average landing bounds

This section addresses the literal finite-segment setting of
[HYP-9161, landing multiplicity](../hypotheses/HYP-9161-landing-multiplicity-polylog.md)
as read on 2026-09-26. Its useful endpoint-stable repair is stated below.

For every integer `m>=6`, set

```
D=ceil(log_2 m)+6,  K=m+D.
```

Let the first m parity letters have prefix odd counts
`o_j=ceil(j log_3 2)`; append D zero letters. These bits are computable
without logarithmic rounding by comparing `3^o` and `2^j`. Their prefix
multipliers during the hover obey `1<=a_j=3^(o_j)/2^j<3`.
Let `r_K` be the unique start residue of the full word and choose

```
n=2^(4K)+r_K,   X=2^(4K+3),   L=log_2 X=4K+3,
theta L=D-3.                                     (5.1)
```

This is one exact positive integer, and all K specified steps occur in
its actual Collatz orbit. During the hover the normalized additive
carry `h_j` in `x_j=a_j(n+h_j)` grows by at most `1/3` at each odd step,
since every previous prefix multiplier is at least one. Hence

```
n<=x_j<=3n+j<4n   for 0<=j<=m.                     (5.2)
```

All following D steps halve. Thus all segment values lie below X and
are positive. They are also pairwise distinct. Indeed every segment
value is at least `n/2^D>3^K`. A repeat after p<=K steps would give
`(2^p-3^e)x=C`, where `0<=C<3^p` by the general carry estimate. If the
coefficient is nonpositive a positive cycle is impossible; otherwise
its integer value gives `x<3^p<=3^K`, a contradiction.

Call i a dipper if some later value within the segment and the horizon
`k=floor L=L` is smaller than `x_i/2^(D-3)`; charge it to the first such
position. For a hover source i<m, no hover value can be the landing,
because every hover value is at least n and its threshold is less than
n. At the end of the drop the value is at most `4n/2^D`, while its
threshold is at least `8n/2^D`; so every hover source is a dipper. At
drop times t<=D-5 the value is at least `32n/2^D`, at least its threshold.
Consequently every hover landing lies in the final five drop positions.

A source at drop time q has its first strict factor-`2^(D-3)` descent
exactly D-2 steps later. It lands within the segment iff q<=2. Thus there
are exactly three additional dipper sources, and their landing positions
also belong to the final five. We have proved

```
number of dippers=m+3,
number of landing positions<=5,
average landing multiplicity >=(m+3)/5=Omega(L),
theta L=ceil(log_2 m)+3=O(log L).                  (5.3)
```

Therefore no uniform bound `average<=C L^beta`, beta<1, and no uniform
bound `average<=C theta L`, holds for all such finite orbit segments.
These are actual distinct positive orbit values, not merely a family
of formally admissible parity words. In the retained m=96 example,
99 dippers land at only four positions.

**Essential endpoint qualification.** The entire constructed segment
has K steps, whereas its horizon is L=4K+3. The endpoint allowance
`O(k)` in the thin-divergence recursion can absorb the whole example.
The construction therefore does not refute a repaired bound

```
|dippers| <= C polylog(L) |landing points| + C' k,   (5.4)
```

nor a bound counting all applicable indices in an entire nonperiodic
orbit. Appending a full future horizon is not free: it can create many
additional landing points and change the average. The least failed
assumption in the endpoint-free hypothesis is uniformity under arbitrary
truncation. Its strongest survivor must retain boundary cost or impose
an explicitly endpoint-stable domain. The established thinness results
remain unchanged.

## 6. Repair of the incoming record proof

The incoming [multiplicity reassessment, Section 8](collatz_landing_20260926_multiplicity_reassessment.md)
claims ballot-thin leaders and running peaks. The peak proof's division
by the old window start z does not establish that its carry error is
small: a large peak does not imply that z exceeds the carry threshold.
The conclusion has a direct stronger repair.

Let `v=y_i` be a running peak, `k=floor(log_2 X)`, and suppose
`v>=2|b|(3/2)^k`. For each `1<=s<=k` with a predecessor available, apply
the affine formula directly to the suffix from `u=y_(i-s)`:

```
v=P_s u+beta_s,   |beta_s|<=|b|((3/2)^s-1)<v/2.
```

Since `0<u<v`, it follows that `P_s>(v/2)/u>1/2`. These P_s are exactly
the prefix products of the reversed k-word. Thus the reversed word
stays above barrier -1 bit; the claimed -1.6-bit barrier is unnecessary.
There are at most k peaks lacking a full predecessor window. Leaders
similarly require an O(k) allowance for the final truncated window.
The same ballot bound applies, and the asymptotic theorem is unchanged.
This repair uses the large endpoint where the hypothesis supplies it.

## 7. Reproduction and remaining constructive target

Incoming one-bit-band refinement (same-session integration): its proof
holds for b>0 or positive windows above5|b|. The weaker y>|b| premise
does not exclude three band points: b=-3 gives7,9,12 in(6,12]. Without
the high premise the landing need not even be a halving: b=-7 gives
6,3,1, with both earlier sources first dropping by factor2 at1, while6
lies outside the proposed band(2,4]. A distinct positive orbit visits
the finite core at most5|b| times, affecting at most5|b|(k+1) windows. Charge
these O_b(k) windows separately; the established thinness exponent survives.
Both witnesses are retained in the script.

Run `python3 04-computation/experiments/crossroads_poset_20260926_integer.py`
or use `-O`. The [retained output](crossroads_poset_20260926_integer.out)
is raw stdout. Controls include both hostile positive cycles, exact
least-residue recovery, sparse-address jumps, the affine no-descent
carrier, and actual hover/drop segments with their exact landing counts.

Two substantive targets remain: decay of the full-support mass (4.3),
and an endpoint-stable oscillation/transport inequality such as (5.4),
proved for actual indices of one orbit. A
fixed-integer proof cannot replace this with fresh uniform residue
sampling. Conversely, the finite-support criterion (3.1) and the
exact height cutoff (4.2) provide a lawful carrier on which an alternative
pruning or rank argument could work. What remains missing is a reason
that an infinite bad branch cannot stabilize at a positive start.
