# Powers, doubling, occurrence births, and the Paley-seven bridge

**Session69 power lane — PROVED elementary statements plus FINITE-EXACT controls.**
No priority claim. Source checkout read-only at
`61d7cca5d8df80169f939268269239c2f67ea724`.
The user's picture admits several distinct exact carriers; this note keeps
all-powers, repeated squaring, labelled occurrences, and integer-value support
separate. It does not settle the unrestricted prime-tournament atlas or full
Hamiltonian-path-complement metagraph composition.

## 0. Useful answer first

* In the **canonical odd-root doubling forest**, six really is the first
  positive integer whose square is a positive doubling of a later root:
  `6²=2²·9`, and `9>6`. Precisely, for `n=2^a u`, `u` odd, the predicate is
  `a>=1 and u>2^a`.
* Without “canonical odd root”, six is not first: `4²=2·8` is earlier.
* The later root nine is **not a new global seed** of the two-operation atlas:
  nine already occurs as `3²`; the earliest seed of 36 is three.
* All-powers rows and repeated-square rows have exact canonical origins and
  pairwise intersections. These are classification results for the proposed
  arithmetic carrier, not classification results for tournaments.
* A genuine, already-recovered bridge is `×2 mod7` on `{1,2,4}`: this is a
  Paley-seven orientation-preserving multiplier and Frobenius on exponent
  coordinates of `F8*`. It preserves actual arcs. It does not preserve integer
  order, integer square-later births, or a forbidden-Hamiltonian-count predicate.

## 1. Inheritance, source corrections, and concept board

Closest directly relevant source:
`05-knowledge/results/arithmetic_seams_20260921_operations.md`, §§1–2.
It proves the two deleted diagonal forests:

* distinct-parent addition omits `x→2x` **inside the ascending relation**;
* distinct-parent multiplication omits `x→x²` **inside proper divisibility**.

The roots are respectively odd positive integers and nonsquares at least two
(with one isolated after the square loop is suppressed). Thus the native
multiplication-diagonal interpretation is **repeated squaring**, not all
exponents. The sources behind this are:

* `THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`,
  swap-fixed operation fibres and their diagonal sidecar;
* `THM-2433-operation-fibre-deletion-incidence-and-startup-scar.md`,
  distinction between internal shadow, transitive closure, and induced order.

The older power/doubling phrasing was not located verbatim in bounded targeted
searches of the pinned tree. The above is recovered mechanism, not an assertion
that the precise user picture was already canonized. Sparse worktree contents
were insufficient; the recovery used `git grep ... HEAD` and `git show`.

Recovered Paley source:
`05-knowledge/results/collatz_paley_bridge_20261001.md`, Propositions 2–5,
especially Proposition 4's `F8` exponent/Fano dictionary. Its MISTAKE-554 repair
labels the Collatz overlap an **explained period-three coincidence**, not an
unrestricted dynamical equivalence. This report uses only the directly checked
finite-field/Paley mechanism, not the wider Collatz claims.

Current metagraph correction read: MISTAKE-674 says module shattering concerns
**labelled vertex subsets**, not unlabelled tree shape; such a correction must
not be weakened when discussing a prospective arithmetic/composition bridge.

| Lane/concept | Coordinate retained | Hostile or decisive probe |
|---|---|---|
| Anchor: arithmetic occurrence atlas | base, doubling exponent, power exponent | 36 has five nondegenerate addresses, not one |
| Odd-root square-later onset | odd part and comparison with source n | arbitrary root8 refutes first-six |
| Niche: repeated-square forest | maximal common power of two in exponents | 27 is not reached from3 by squaring/doubling |
| Primitive power components | primitive odd root, exponent gcd, v2 | 15 is a new seed but not prime |
| Wildcard: Paley7/Frobenius | labelled residues and orientation difference set | modulo7 destroys positive integer chronology |
| Full metagraph composition | path marking/ports and full decorated quotient | no supplied bridge from the arithmetic row join |

Meta-patterns used: **correct the object before sharpening the technique**;
**type every analogy and every implication**; **search the statement before the
method**. No new meta-pattern is proposed from this single lane.

## 2. Precisely define the two-operation carriers

Use positive bases, exponents `r>=1`, and doubling exponents `k>=0`:

`F_(k,r)(n)=2^k n^r`, and `S_n={F_(k,r)(n): k>=0,r>=1}`.

For repeated squaring use only `r=2^j`, `j>=0`:

`Q_n={2^k n^(2^j): k,j>=0}`.

These sets are value supports, not indexed cell multisets. Exponent zero would
introduce the value one in every row, and base zero would produce a degenerate
zero component; neither is included. Base one is retained as the root of the
power-of-two component, but its repeated identical power positions are removed
when counting nondegenerate addresses below.

**PROVED operation law.**

`F_(k,r) o F_(ell,s)=F_(k+r ell,rs)`.

Equivalently `P_r D_k = D_(rk) P_r`. Squaring and doubling do not commute:
`(2n)²=4n²`, not `2n²`. Thus a genuinely rectangular picture must specify which
operation is performed first. Every finite interleaving has the displayed
normal form; for the square-only alphabet its power exponent is a power of two.

## 3. All-powers canonical form and first-occurrence classification

Every integer with odd part greater than one has a unique expression

`n=2^a b^e`,

where `a>=0`, `e>=1`, and `b>1` is odd and not a proper perfect power. Explicitly,
if the odd part is `product p^v_p`, let `e=gcd_p(v_p)` and
`b=product p^(v_p/e)`. The pure powers of two form the separate component with
canonical seed one.

On the component indexed by b, use coordinates `(a,e)`. Then

`D_k(a,e)=(a+k,e)`, `P_r(a,e)=(ra,re)`.

**PROVED earliest-seed theorem.** The least numerical row base containing n is
b (or one for powers of two). Indeed n is in `S_b`; conversely every row that
contains n has the same primitive odd root and therefore base at least b.

Consequently, on admitting rows by increasing base n:

* n contributes an entirely new component iff n is an odd non-perfect-power
  integer, together with the initial n=1;
* every other row lies entirely in an already present earlier row;
* “new component” is not “prime”: 15 is already the smallest odd composite that
  is not a perfect power, so it is a new component too.

This is an exact first-occurrence classifier, provided all exponents in every
admitted row are visible. A finite-height display has a different chronology
and must retain that height bound.

**PROVED exact collision test.** Write `n=2^a b^e`, `m=2^c d^f` with odd parts
nontrivial. Then

`2^k n^r=2^ell m^s`

iff `b=d`, `er=fs`, and `ar+k=cs+ell`. Unique prime factorization proves both
directions. Support alone forgets k,r and the number of cell addresses.

**PROVED row-intersection law.** If `b!=d`, the rows are disjoint. If `b=d`, put

`g=gcd(e,f)`, `E=lcm(e,f)`, `A=max(af/g,ce/g)`, `L=2^A b^E`.

Then exactly

`S_n intersection S_m = S_L`.

To prove this, `er=fs` forces `r=(f/g)t`, `s=(e/g)t`, `t>=1`. The least allowed
2-adic exponent is `At`, and every larger exponent is attainable. L is the
least numerical common occurrence, and the entire intersection is its row.
For powers of two the analogous L is `2^max(a,c)`. No mixed component intersects
that pure-two component.

Examples: `S_6 intersection S_9=S_36` and
`S_6 intersection S_27=S_216`. This join is associative, commutative and
idempotent within each component because it represents set intersection
uniquely; it is not ordinary integer lcm (`lcm(6,9)=18`, not36).

## 4. Repeated-squaring canonical form

Every integer whose odd part exceeds one uniquely has form

`n=2^a v^(2^t)`, with v odd and nonsquare, `t>=0`.

In the preceding notation, write `e=2^t h` with h odd and set `v=b^h`.
The earliest square/doubling seed is v; initial novel bases are precisely odd
nonsquares, plus one. Thus 27 is novel here, but was old in the all-powers atlas.
This is the **first numerical difference between earliest-seed classifications**;
it is not the first difference between the models in every sense.

For `n=2^a v^(2^t)` and `m=2^c w^(2^s)`, rows are disjoint if `v!=w`. Otherwise
put `T=max(t,s)`, `A=max(a 2^(T-t),c 2^(T-s))` and `L=2^A v^(2^T)`. Exactly

`Q_n intersection Q_m=Q_L`.

The proof aligns the square heights at T, then doubles until the dyadic
coordinates match. Later common heights scale both required exponents by two.
For pure powers of two use again `L=2^max(a,c)`.

In particular `Q_6 intersection Q_9=Q_36`, but `Q_6 intersection Q_27` is empty.
The row for base three already distinguishes the models: it contains27 in
S_3 but not Q_3. Address multiplicities differ already at8; origins first
differ at27. These distinctions prevent a misleading unqualified “first”.

## 5. Exact occurrence multiplicities

Count labelled addresses `(n,k,r)` with `n>=2`, `k>=0`, `r>=1` producing a fixed
integer x; this is finite. The base-one addresses are excluded because all
r give the same power-of-two output, producing infinite raw multiplicity.

For `x=2^A b^E`, b odd primitive, the all-powers count is

`M(x)=sum_(r divides E) (floor(A/r)+1)`.

Indeed r must divide every odd prime exponent, and the dyadic exponent a of
the base may be any integer `0<=a<=floor(A/r)`; k is then forced to `A-ar`.
For `x=2^A`, `A>=1`, the count instead is `sum_(r=1)^A floor(A/r)` because
bases are `2^a` with a positive.

For `x=2^A v^(2^T)`, v odd nonsquare, the repeated-square address count is

`M_sq(x)=sum_(j=0)^T (floor(A/2^j)+1)`.

For pure `2^A` it is `sum_(j>=0,2^j<=A) floor(A/2^j)`.

Example:36 has five addresses in either model:
`36`, `2·18`, `4·9`, `6²`, `4·3²`. It has only one value and one earliest seed3.
At8 the counts are five versus four: the all-power-only extra address is `2³`.
At27 they are two versus one: `27` and `3³` versus just27.

These formulas classify exact repeats. Alone they do not classify additive
relations among distinct integer values. Section 6a supplies a separate exact
additive bridge by retaining a positive difference; no abc bound is established.

## 6. The exact first-six theorem and its limits

The canonical doubling forest has one ray for every odd u:
`u,2u,4u,...`. If `n=2^a u`, then

`n²=2^(2a) u²`.

The unique odd root of the square's doubling ray is `u²`. Therefore n² is a
**positive doubling** of a **strictly later odd root** iff

`a>=1 and u²>2^a u`, equivalently `a>=1 and u>2^a`.

The least possibility is a=1,u=3, yielding n=6. No n<6 qualifies. More
precisely the first occurrence at fixed a is

`n_min(a)=2^a(2^a+1)=2^(2a)+2^a`.

The list begins `6,10,14,18,20,22,26,28,30,34,36,...`.
For any fixed r>=2, the analogous positive-doubling/later-odd-root condition
for n^r is `a>=1 and u^(r-1)>2^a`; six is again the least possible n.

**Not rare in density.** The square-later set has natural density one half.
It is contained in the evens. Its exceptions among evens up to X satisfy
`u<=2^a`: summing over `a<=floor((log_2 X)/2)` gives `O(sqrt X)` possible odd
u, while the tail `a>floor((log_2 X)/2)` contributes
`O(sqrt X+log X)` by the bound `u<=X/2^a`. Thus only o(X) even integers fail.
The finite count is4859 qualifying n through10000. “First at six” is an onset
statement, not evidence of rarity matching a rare aliquot fixed-point property.

**Hostile domain changes.** If an arbitrary root m>n is allowed in
`n²=2^k m`, k>=1, the answer is every even n>=4: choose k=1 and m=n²/2.
This makes four first. If k=0 is admitted, every n>=2 trivially has later
root n², making two first. If “later” refers to the full power/doubling birth
order, nine is already present at three and is not a later-new component.

## 6a. An exact additive bridge to a binary-term abc family

Retain the positive square-root gap, not just the collision equality. For the
square-later n=2^a u,

`u²-n=u(u-2^a)`.

After dividing by u, this gives the positive additive triple

`(u-2^a)+2^a=u`.

All three terms are pairwise coprime, because u is odd. This is a bijection
between the square-later integers and primitive positive abc triples `(A,B,C)`
with `B=2^a>=2` and C odd: reconstruct n as BC. Here abc means the ordinary
coprime additive carrier `A+B=C`, not an assertion of the abc conjecture.
The source predicate `u²>n` is exactly positivity of A; no coordinate is lost
once the distinguished power-of-two term B is retained.

Six corresponds to `1+2=3`, not a radical-excess example. Under the conventional
quality test `C>rad(ABC)` (equivalently quality greater than one), the first
square-later n giving such a triple is72, from `1+8=9`, radical6. This minimality
is a finite-exact statement checked at every n below72, not an asymptotic abc
claim. The n<=10000 run checks this bijection and coprimality directly.

In these coordinates, an abc-type inequality would read

`u <= K_epsilon [2 rad(u(u-2^a))]^(1+epsilon)`.

The atlas formulas do not establish this inequality: they retain valuations
and exact equality, while the missing coordinate is the radical of the
*difference* `u-2^a`. That sidecar can change substantially along one component.
For the higher-power later-root predicate, the same normalization gives
`(u^(r-1)-2^a)+2^a=u^(r-1)`, again primitive when positive. This is a precise
way to pose the abc comparison without claiming a bound from multiplicative
collisions alone.

## 7. The genuine but restricted bridge to the prime tournament atlas

Let `Q={1,2,4}` in `Z/7`, and orient `i→j` iff `j-i in Q`. Since Q and -Q
partition the nonzero residues, this is a tournament. It is the Paley-seven
quotient P7. The following statements directly retain its native arc predicate:

1. `Q=<2>` and `2Q=Q`; so doubling the **residue labels** is an automorphism.
2. Translation also preserves differences. The21 affine maps
   `i→a i+b`, `a in Q`, `b in Z/7`, preserve all arcs.
3. With `F8=F2[alpha]/(alpha³+alpha+1)`, the bijection
   `i↦alpha^i` carries `i→2i mod7` to `x↦x²`. The **transported** arc condition
   is `alpha^i→alpha^j iff alpha^(j-i) in {alpha,alpha²,alpha⁴}`.
4. `{1,2,4}` is a `(7,3,1)` difference set: its six ordered nonzero differences
   exhaust `1,...,6` exactly once. Also `alpha+alpha²+alpha⁴=0`.

Thus the same three labelled exponents give the integer sum7, a modular
doubling orbit, and a field Frobenius orbit; the arc bridge requires the
specified difference predicate, not the integer sum alone.

**FINITE-EXACT independent controls:** all5040 vertex permutations were tested:
P7 has exactly21 automorphisms, all the displayed affine maps; all proper
non-singleton subsets were tested, and none is a module (so P7 is a genuine
prime quotient). Hamiltonian paths were independently counted as189, not7
or21. These finite checks support the explicit order-seven carrier only.

**Connection contract.** Source: the mod7 residue quotient with difference
set Q. Target: P7, or its transported `F8*` labels. Map: identity on residues,
then `i↦alpha^i`. Preserved predicate: every oriented edge, hence module and
Hamiltonian-path predicates under relabelling. Destroyed information when
passing from positive integers to residues: magnitude, order, integer
factorization, and the earliest-row base. Sidecar needed to restore those:
the full integer lift and original occurrence address. Cheapest hostile:
3 and10 are identical residues but have different parity and square-later
status. No theorem transports the square-later property through this quotient.

For decorated prime quotients, the affine maps can transport an attached
child-colouring along their actual vertex permutation; this is a legitimate
isomorphism test for that one quotient, not an unrestricted prime atlas.
For Hamiltonian-path-complement metagraph composition, a path/port sidecar is
still required. The arithmetic row intersection operation has no supplied map
to that path operation, and integer power-root primitiveness is not tournament
module-primality (15 already prevents confusing arithmetic seeds with primes).

The number21 here is the symmetry-group order. This does **not** prove that
21 is a forbidden Hamiltonian-path count; that separate global predicate
requires its own source theorem and hypotheses. Nor does it identify the
aliquot map's missing targets with tournament gaps.

## 8. Reproduction and scope

Run `python3 check.py` (stdlib only). It uses independent direct value
enumeration to test the origin and address-count formulas for every x<=5000;
checks both row-intersection formulas for every pair of bases1..120 with all
row values<=1000000; checks the square-later predicate for n<=10000; tests the
semigroup law on the printed finite parameter box; and fully enumerates
Paley-seven automorphisms, paths, and modules. The formulas themselves have
proofs above, rather than being inferred from this bounded computation.

Result: **57332 checks passed**. Normal and `python3 -O` outputs byte-match.
The checker uses an explicit `require` function, so optimization cannot erase
its validation gates. Input snapshots and SHA256 manifest retain source
provenance. No source repository files were edited by this lane.

### Exact disposition

* **Completed/high certainty:** carrier definitions, primitive-root origins,
  row intersections, address counts, repaired first-six theorem, density
  boundary, primitive binary-term abc correspondence, and finite Paley7 controls.
* **Recovered/high certainty:** operation-fibre diagonal-forest mechanism and
  already-audited period-three Paley/Frobenius dictionary.
* **Open:** a target-preserving bridge from integer occurrence births to
  unrestricted tournament prime classification or full complement-line
  composition. Missing action/target: supply a map preserving the actual
  quotient/module/path-port predicates, not just a shared small integer.
* **Operational blocker:** none. The inherited message-peek command reported
  no `.machine-id`; this lane did not fabricate registration or mutate the
  read-only worktree. Research and finite controls were unaffected.
