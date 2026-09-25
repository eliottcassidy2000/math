# Divisor balance has an exact finite profile controller, but its quartic filter is not a Collatz quotient

**Status.** INHERITED PROVED, independently re-derived here: the precise
`F=S+U` classification and its values `{0,2,10}`. PROVED: the safe
prime-fourth-power exclusion, the eight live profiles with an absorbing
reject state, and the support-label requirement for multiplication.
FINITE-EXACT: direct divisor checks through `1000`, 12,869 exponent-profile
controls, and the listed minimal examples. No tournament decoder, Collatz
convergence claim, or general permission to discard degree-four information
is asserted. No novelty claim for the inherited classification.

Session: collatz seams, 2026-09-25, divisor-balance lane. The current
question asks what can validly be ignored above the cubic level and how
the earlier three divisor shapes fit the odd/doubled/higher-layer picture.

## Inheritance and scope recovery

Closest proved mechanism: the retained multiplicative fibre in
[THM-2422, operation fibres and twin-center ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md),
especially the divisor count and repeated-factor diagonal. The exact
classification already lives in
[arithmetic braids, divisor note, DB1--DB3](arithmetic_braids_20260917_divisors.md),
with its family extension in
[divisor balance family](collatz_mod6_20260917_divisor_balance_family.md).
These are proved result notes, not newly allocated canon IDs.

Canonical hostile: `N=4`, with `F=1`, `S=U=1`, so the equality is not a
conservation law. The correction is maintained in
[MISTAKES, 2026-09-21 Odd-square geometry](../../01-canon/MISTAKES.md),
the item beginning "FALSE pasted conservation". The least-used sidecar
is the named prime exponent profile: support size, total degree, squareclass,
and parity each forget relevant data. The earlier
[squarefree symmetry note](arithmetic_braids2_20260917_squarefree_symmetry.md)
already records that `2,8,32` share a squareclass but disagree on balance.

Anchor: recover the exact meaning before making a decoder. Niche: replace
an informal truncation by a sound finite classifier. Wildcard: distinguish
the unrelated meanings of "odd". The live board is divisor profiles,
properness boundaries, exponent truncation, support overlap, and parity.
The classification changes the board by supplying a finite reject region;
the multiplication and Collatz hostiles delimit the operations under which
that region remains closed.

## 1. What F, S, U, and the numbers 0, 2, 10 actually count

For `N>=2`, the inherited definitions are

    F(N) = #{d : d|N, 1<d<N},
    S(N) = #{d : d|N, 1<d<N, d is squarefree},
    U(N) = #{d : d|N, 1<d<N, d is prime}.

Thus `U` is neither total prime-factor multiplicity `Omega` nor always
the support size `omega`: at a prime, `U=0` because the prime itself is
not a proper divisor. For `N=prod p_i^a_i` with `r` distinct primes,

    F=prod(a_i+1)-2,
    S=2^r-1-[N squarefree],
    U=r-[N prime].                                         (1)

The endpoint corrections in (1) cannot be dropped. If the unit `N=1` is
admitted, its three defining sets are empty, giving `F=S=U=0`; handle
this separately from (1).

The exact classification is

    N>=2: F=S+U  iff  N=p, p^3, or p^2 q r,                (2)

where `p,q,r` are primes and are pairwise distinct in the last shape.
The corresponding equations are

    p:       0=0+0,
    p^3:     2=1+1,
    p^2qr:  10=7+3.

Therefore `{0,2,10}` is the set of **values of F at balanced numbers**.
It is not a set of defect values, locations of a parameter, or tournament
orders. Every balanced number has defect `D=F-S-U=0`.

The first mixed example is `60=2^2*3*5`. It has proper nontrivial divisors
`2,3,4,5,6,10,12,15,20,30`; seven of these are squarefree and three prime.
Minimality follows because three distinct primes must be at least `2,3,5`,
and squaring the smallest gives the least product, `2^2*3*5=60`.

### Complete elementary proof of the classification

Primes satisfy the equality by properness. A squarefree composite has
`F=S` and `U=r>0`, so does not. In the nonsquarefree case, (1) gives

    prod(a_i+1)=2^r+r+1.                                   (3)

For `r=1`, (3) forces `a_1=3`. For `r=2`, its right side is `7`, which
cannot be a product of two integers at least two. For `r=3`, the right
side is `12`, whose only factorization into three integers at least two
is `3*2*2`; hence the profile is `(2,1,1)` up to permutation.

For `r>=4`, at least one exponent is at least two, so

    prod(a_i+1)>=3*2^(r-1)>2^r+r+1.

The strict inequality is equivalent to `2^(r-1)>r+1`, true at `r=4`
and preserved when `r` increases. This excludes every remaining profile.
The three listed profiles verify the reverse implication directly.

## 2. Four meanings of "quartic" give different decisions

| Proposed condition | Meaning | Consequence for balance |
|---|---|---|
| Some `p^4` divides `N` | one prime exponent is at least four | safe exclusion |
| `Omega(N)>=4` | total prime multiplicity at least four | unsafe: `60` is balanced |
| `v_2(N)>=2` | divisible by four, on deeper dyadic layers | unsafe: `8` and `60` are balanced |
| A polynomial contains a degree-four term | degree of a function | no implication from (2) |

In particular, `p^2qr` has total prime degree four while its largest
individual exponent is only two. Its `F=10` does not identify it with the
ten-vertex tournament studied in the previous decoder experiment.

**Safe prime-fourth-power filter.** If some `a_i>=4`, then

    D=prod(a_i+1)-2^r-r-1
       >=5*2^(r-1)-2^r-r-1
       =3*2^(r-1)-r-1 > 0.                                (4)

The last expression is `1` at `r=1` and strictly increases afterward.
Thus a number divisible by any prime fourth power cannot be balanced.
This is rejection by a proved bound, not replacement by a lower power.

**Unsafe truncations.** Capping prime exponents at three sends `16=2^4`
(defect one) to `8=2^3` (balanced). Reducing exponents modulo four sends
`32=2^5` (defect two) to `2` (balanced). If the unit is admitted, the
modulo-four failure already occurs at `16 -> 1`. The overflow information
cannot simply be erased while retaining the balance predicate.

The safe filter is also **not closed under Collatz**. Ordinary halving
sends rejected `16` to balanced `8`. On odd integers, the accelerated
Collatz step sends rejected `81=3^4` to balanced `61`, since
`3*81+1=4*61`. Consequently, removing these states from a Collatz route
graph discards legal transitions rather than producing a proved quotient.

## 3. A positive construction: eight live profiles and one reject state

There is nevertheless a complete finite controller for balance under
**multiplication by a stream of certified primes**. The controller retains
one of the following sorted exponent profiles:

    (), (1), (2), (3), (1,1), (1,2), (1,1,1), (1,1,2),

or a ninth state `REJECT`. In a live state, it also retains the names of
its at most three distinct primes and which exponent belongs to each.
The empty profile represents the unit. This is finite control with a
small support register; it is not a finite automaton that recognizes
arbitrary unlabeled integer inputs without arithmetic information.

The accepting states are `(1)`, `(3)`, `(1,1,2)`, plus `()` if the unit
is admitted. The remaining live states `(2)`, `(1,1)`, `(1,2)`, and
`(1,1,1)` are transient: they are not balanced, but a suitable future
prime multiplication can make them balanced. Accepting states need not
be absorbing: a prime can become a square, which is not balanced.

On a prime input, the register tells whether it is new or which existing
exponent it increments. The exact transitions are:

| State | New distinct prime | Raise exponent 1 | Raise exponent 2 | Raise exponent 3 |
|---|---|---|---|---|
| `()` | `(1)` | -- | -- | -- |
| `(1)` | `(1,1)` | `(2)` | -- | -- |
| `(2)` | `(1,2)` | -- | `(3)` | -- |
| `(3)` | REJECT | -- | -- | REJECT |
| `(1,1)` | `(1,1,1)` | `(1,2)` | -- | -- |
| `(1,2)` | `(1,1,2)` | REJECT | REJECT | -- |
| `(1,1,1)` | REJECT | `(1,1,2)` | -- | -- |
| `(1,1,2)` | REJECT | REJECT | REJECT | -- |
| REJECT | REJECT | REJECT | REJECT | REJECT |

A dash means that exponent class is absent, so that transition is not
an available input description in the given live state.

### Soundness and the absorbing-reject proof

The eight live profiles are exactly the profiles of positive integers
that divide some balanced integer. To see this, every balanced integer
has one of the shapes in (2); enumerating its divisor profiles gives
precisely the displayed eight. Conversely each displayed profile embeds
in one of those shapes, so suitable extra prime factors complete it to a
balanced number. This statement includes the unit.

If `N` is rejected and a positive integer `M*N` could be balanced, then
`N` would divide a balanced integer, contrary to rejection. Thus rejected
numbers form an upward-closed set under divisibility, and `REJECT` is
absorbing under multiplication. Once this state is entered, support labels
can be discarded for this predicate and this operation.

For a live state, multiplication by one named prime either appends a new
exponent one or increases one existing exponent by one. The table performs
exactly those updates and sends any profile outside the eight to `REJECT`.
Starting at the unit, induction on the number of prime inputs shows that
the controller gives the true profile whenever live and a sound absorbing
rejection otherwise. Its final accepting status is exactly (2), with the
chosen unit convention. Reading prime factors in another order changes
no final result because the underlying named exponents are added exactly.

### Why names, rather than only sorted profiles, are necessary

The numbers `6` and `15` both have profile `(1,1)`. Under the same
multiplier four they become

    4*6=24=2^3*3,          D=1,
    4*15=60=2^2*3*5,       D=0.

Thus the profile alone does not define an operation-compatible quotient
for multiplication by specified integers. The missing information is
whether the multiplying prime overlaps the retained support and, when it
does, which exponent it raises. The profile controller becomes exact
when this support-register information is retained. It assumes certified
prime inputs or a factorization of each composite multiplier; it does not
replace primality testing or factorization.

## 4. The parity words refer to different predicates

The relevant distinctions are:

* An **odd integer** has `v_2=0`. Its oddness says nothing about the parity
  of its exponent multiplicities: `9` is odd but a prime square.
* An **odd prime exponent** contributes its prime to the multiplicative
  squareclass; even exponents cancel in that quotient. `2,8,32` share a
  squareclass while only `2,8` are balanced.
* An **odd function** satisfies `f(-x)=-f(x)`. It need not take odd values
  on odd integers: `f(x)=x^3+x` is an odd function whose integer values
  are all even. An odd integer polynomial contains only odd-degree terms,
  but an odd total degree alone is insufficient (`x^3+1` is not odd).
* An **odd-length directed cycle** concerns its number of edges. It is
  distinct from having an odd number of cycles or an odd number of paths.
* [THM-001, Redei](../../01-canon/theorems/THM-001-redei.md) states that a
  tournament has an odd **number of directed Hamiltonian paths**. It does
  not say that all directed cycles have odd length. A tournament can
  contain a directed four-cycle.

The inherited odd-function floor identity in
[floor sums and odd functions, S1--S2](collatz_mod6_20260917_floor_sums_odd_functions.md)
has a precise common mechanism: the involution `k -> n-k` sends residues
to their negatives, and the zero/fixed layer is retained. For an odd
function represented modulo `n`, the identity is

    sum floor(f(k)/n)
      = (sum f(k))/n - (n-1)/2 + Z_f(n)/2,    1<=k<n,

where `Z_f(n)` counts nonzero residue classes on which `f` vanishes.
The formula's correction is not determined by the word "odd" alone.
For `f(x)=x^e`, the exact correction uses
`Z_e(n)=n/prod p^ceil(a_p/e)-1`; exponent heights re-enter even after the
antipodal pairing has simplified the main term. These are inherited
identities, not a new bridge to Collatz.

## 5. Reproduction and boundary

    python 04-computation/experiments/seam_prime_balance_20260925.py

Output: [seam_prime_balance_20260925.out](seam_prime_balance_20260925.out).
The checker uses no external packages, floating point, or removable
assertions. It verifies (1) against direct divisor definitions for every
`1<=N<=1000`; tests all 12,869 sorted profiles of support size `1..8`
and exponents `1..8`; checks the classification, safe bound, live-divisor
criterion, and every reject transition in that finite universe; prints
the controller table and explicit hostiles. Completeness outside the
tested universe is supplied by the proofs above, not by extrapolation.

No inherited false statement was left in force: the old conservation
claim was already corrected. The positive result is an exact controller
for multiplicative growth with named-prime information. The stopping
boundary is the operation: it does not commute with Collatz's addition
and division, and it gives no rule for quotienting polynomial degrees,
tournament vertices, or odd cycles.
