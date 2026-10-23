# Why 7/15 and 8/15 occur: native conductance, retained masses, and a ternary spectral address

Status: **PROVED** on the stated finite incidence and two-state Markov
carriers; **CITED** inherited rational-tree results credited below; **FINITE-EXACT**
controls in `check.py`. The ternary address is an arithmetic representation
of a spectral parameter. No assertion that its three branches are native
tournament moves, or preserve a native electrical guard, is made.

## 1. Inheritance and the map to test

The closest proved mechanism is the exact native order-four matrix in
[the Hamiltonian restriction packet](../../hamiltonian_observer_restrictions_20261023/coupling/RESTRICTIONS.md),
whose nontransitive component has eigenvalue `8/15` and gap `7/15`.
Complement-pairing compressed sources instead gives gap `8/15`.
[The all-size connectedness proof](../../hamiltonian_observer_restrictions_20261023/coupling/NATIVE_CONNECTEDNESS.md)
does not preserve arbitrary source-specific guards.

The arithmetic inheritance is
[THM-4057: Stern--Brocot depth and rational-edge gauges](../../../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md),
[THM-3333: the Gaussian-square Farey light cone](../../../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
and [the ternary Berggren decoder](../../ternary_berggren_20260925.md).
Their corrections matter: a reduced rational is not itself a tournament
edge; Stern--Brocot and Calkin--Wilf addresses differ; a Pythagorean lift
needs coprimality and the parity/content condition; branch counts alone
do not transfer a dynamical clock.

The live objects are actual incidence words, electrical conductance,
stationary mass, rational addresses, ternary recursion, and guard support.
The canonical hostile will be a single removed actual presentation that
disconnects the order-four incidence. The least-used sidecar is the mass
left by target leaves after electrical elimination.

Several older exact appearances have different types and should not be
collapsed into the new spectral calculation:

| recovered source | role of the fraction | retained boundary |
|---|---|---|
| [digits/repunits, section 2](../../procgen_repunit_20260924_digits_rotation_repunits.md) | `8/15` is the Haar one-step odd-Collatz last-digit probability `9 -> 9` | a probability law on random 2-adic starts, not an every-integer clock bound |
| [rank/Christoffel strategy note, section 4.2](../../procgen_tension_20260926_ranks_christoffel_duality.md) | `8/15` belongs to the exact level-six set of maximal odd-step cycle densities | that realizability set is NOT a Farey interval; a rational address does not prove a legal strategy exists |
| [gold-partition/balance note, section 4](../../one_twenty_eight_atlas_balance_gold_20261014.md) | `7/15` occurs among small-order minimum arc-precedence probabilities | its proposed universal `1/3` arc lower bound fails at order eight with `13/45` |
| [odd-square edge note, section 2](../../odd_square_20260921_edges.md) | `8/15` bounds a smaller-angle tangent in two guarded Collatz edge regimes | the all-regime sharp positive bound there is `65/72`, not `8/15` |

These inherited occurrences motivate tests. They are not evidence that a
transition entry, an eigenvalue, an angle tangent and a cycle density are
the same invariant.

## 2. Reconstruct the actual eight-word incidence

Fix the native path `0 -> 1 -> 2 -> 3`. Free bits, in order, reverse the
chords `(0,2),(0,3),(1,3)`. Compress them to graph pairs
`(0,1),(0,2),(1,2)`. On three vertices the compressed graph class is
determined by its edge count `e=0,1,2,3`.

The following integer incidence table is obtained directly from all eight
words. `D_-` and `D_+` are the two diamond tournaments; `S` is the strong
four-vertex tournament.

| compressed graph | transitive | D_- | D_+ | S | row mass |
|---|---:|---:|---:|---:|---:|
| empty | 1 | 0 | 0 | 0 | 1 |
| one edge | 0 | 1 | 1 | 1 | 3 |
| two-edge path | 0 | 0 | 0 | 3 | 3 |
| triangle | 0 | 0 | 0 | 1 | 1 |

The column masses are `1,1,1,5`. This is a count of normalized path-marked
presentations, not uniform tournament classes. The checker independently
enumerates actual arcs, scores, and Hamiltonian orders for each word.

For any nonnegative incidence matrix `W` with positive row/column masses
`r_i,c_t`, the source coupling and its conductance Laplacian are

    P = D_r^(-1) W D_c^(-1) W^T,
    L = D_r - W D_c^(-1) W^T.                        (1)

These formulas follow by choosing a target using a row's incidence weights,
then choosing a source using that target's incidence weights. They are also
the Schur complement of the bipartite electrical Laplacian after target
elimination. The stationary source measure is proportional to `r`.

## 3. The same 4/5 conductance explains both fractions

First delete the isolated empty/transitive component. The path and triangle
rows have identical normalized target profiles (both see only `S`), hence
their full source-transition rows are identical. Lumping them together is
exact; the discarded mode has eigenvalue zero. The two remaining source
masses are

    a=3, b=4.

Only target `S` connects them, with incidence weights `1,4`. Eliminating
that target gives raw conductance

    c = 1*4/(1+4) = 4/5.                             (2)

The two diamond targets are leaves of the first source. Eliminating them
adds no conductance, but their mass remains in `a=3`. Thus the lumped chain
and its gap are

    P_0 = [[11/15,4/15],[1/5,4/5]],
    gamma_0 = c(1/a+1/b) = (4/5)(1/3+1/4) = 7/15.    (3)

Its other eigenvalue is `8/15`. With normalized stationary masses
`(3/7,4/7)`, normalized conductance is `4/35`; the effective resistance is
`35/4`. Mixing raw and normalized conductances would lose a factor seven.

Now retain all eight words and pair their compressed source graphs by
complement: empty with triangle, and
one edge with two-edge path. Their masses are `a=2,b=6`. The shared strong
column still has entries `1,4`, so its conductance is still exactly `4/5`.
The transitive target is now an exclusive leaf on the first source, and
the diamonds are exclusive leaves on the second. Hence

    P_1 = [[3/5,2/5],[2/15,13/15]],
    gamma_1 = (4/5)(1/2+1/6) = 8/15.                 (4)

The products `ab` are both 12, while the total masses are seven and eight.
This is the precise reason the two fractions exchange the roles of gap
and nonconstant eigenvalue. It is an incidence/mass calculation, not an
assumption that numerical complementation of a fraction implements graph
complementation in general. Converse-merging the two diamond targets
does not change the calculation, because both are exclusive to the same
source in each relevant table.

**Mass is essential.** In a general two-source incidence, padding either
source by exclusive target leaves leaves the effective conductance fixed
but changes the gap `c(1/a+1/b)`. Padding both without bound drives the gap
to zero. Conversely, a fixed gap with varying normalized masses can have
unbounded source-to-target hitting time. These padding statements concern
general weighted incidence; no realization of every padded table by actual
tournaments is asserted.

## 4. A one-presentation guard deletes the whole bridge

The one-edge word meeting `S` is exactly code `2`: the only reversed native
chord is `(0,3)`. Among nonzero words, delete that single presentation and
retain `1,3,4,5,6,7`. The one-edge source then sees only the two diamonds;
the path/triangle sources see only `S`. The new conductance between these
source blocks is zero and the native coupling disconnects.

This is an exact literal-word guard, not the repository's electrical mod-8
predicate. It shows what a guard-specific proof must preserve: actual
bridge support, not merely nonempty endpoint populations. The unguided
native connectedness theorem cannot be applied after this deletion.

For two source blocks, positive surviving column overlaps give the exact
conductance

    c_guard = sum_t W_(1t) W_(2t)/(W_(1t)+W_(2t)),     (5)

omitting zero columns. The gap is then `c_guard(1/a_guard+1/b_guard)`.
This is a usable source-specific certificate when the actual guarded
carrier really has two states, or when identical-row elimination proves
that all remaining nonconstant modes are zero.

It is not enough to have a two-block quotient or even strong lumpability:
two disjoint copies of the same positive-gap two-state chain have that
same two-block quotient, but their full gap is zero. Internal modes and
actual connectivity must be retained. This is the cheapest hostile to
using a coarse Farey or conductance summary as a complete guard proof.

## 5. The exact Farey and Pythagorean address of the pair

The fractions `7/15` and `8/15` are NOT Farey neighbors: their determinant
has magnitude 15. Each is, however, Farey-adjacent to `1/2`:

    |7*2-15|=1,       |8*2-15|=1.

Their unreduced mediant `15/30` is `1/2`; the factor 15 is exactly what a
false Farey-neighbor claim would discard. In general the pair

    k/(2k+1), (k+1)/(2k+1)

straddles `1/2`, and each member is adjacent to it.

Normalize their odds as the Euclid pair `(m,n)=(k+1,k)`. It is coprime and
of opposite parity, giving the primitive Pythagorean triple

    (m^2-n^2,2mn,m^2+n^2)
      = (2k+1,2k(k+1),2k^2+2k+1).                  (6)

At `k=7` this is `(15,112,113)`. The odd leg equals the probability
denominator because `m-n=1`; that identity does not hold for arbitrary
Euclid parameters. The consecutive-parameter family is the inherited
parabolic Berggren spine, beginning at `(m,n)=(2,1)`.

### A genuine ternary tree of spectral coordinates

Write a reduced `g=n/(m+n)` with `0<g<1/2` and odd denominator. This is
equivalent to `m>n>0`, `gcd(m,n)=1`, and opposite parity. The Berggren
parameter branches in this ordered convention are

    A(m,n)=(2m-n,m),
    B(m,n)=(2m+n,m),
    C(m,n)=(m+2n,n).

They induce the exact rational maps

    A(g)=(1-g)/(3-4g),
    B(g)=(1-g)/(3-2g),
    C(g)=g/(1+2g).                                  (7)

Their ranges are respectively `(1/3,1/2)`, `(1/4,1/3)`, and `(0,1/4)`.
The boundary `1/4` has even denominator and is excluded; `1/3` is the
unique root. To decode, put `r=m/n`. For `1<r<2` the parent is
`(n,2n-m)`; for `2<r<3` it is `(n,m-2n)`; for `r>3` it is `(m-2n,n)`.
The omitted ratio `3` contradicts opposite parity after reduction, and
ratio `2` is precisely root `(2,1)`. Each parent decreases `m+n`.
Unimodularity preserves coprimality and the displayed branches preserve
opposite parity. Thus every reduced odd-denominator coordinate in `(0,1/2)`
appears exactly once, with a terminating constructive inverse.

Here `A^6(1/3)=7/15`; the complementary coordinate is `8/15`.
More generally `A^j(1/3)=(j+1)/(2j+3)`. This is a ternary-tree address
for the spectral pair, not merely a shared integer label.

For ANY retained positive stationary probability vector `pi=(pi_1,pi_2)`,
define the actual two-state chain

    P_g=(1-g)I+g J_pi,
    J_pi(i,j)=pi_j.                                 (8)

It is reversible, stochastic and PSD, with gap `g`, other eigenvalue
`1-g`, and normalized conductance `g pi_1 pi_2`. Holding `pi` fixed, (7)
therefore gives a genuine ternary family of weighted Markov chains with
the same stationary law. The reflected coordinates `1-g` give the other
half of the odd-denominator parameter range. A bit recording which
coordinate is the gap is necessary: the two actual native chains (3),(4)
have the same unordered pair `{gap, nonconstant eigenvalue}` and different
gaps; their full spectra are not equal.

Allowing `0<g<=1` in (8), the representation `(g,pi)` recovers every
irreducible two-state reversible PSD chain. The scalar `g` alone does not:
at `g=7/15`,
the native `pi=(3/7,4/7)` gives mean hitting time `15/4` from the first
state to the second, whereas balanced `pi=(1/2,1/2)` gives `30/7`.
In general this hitting time is `1/(g pi_2)`, unbounded as `pi_2` tends
to zero even with the spectral coordinate fixed. Neither branch address
nor gap pays an actual source clock without the mass sidecar.

The maps (7) change conductance. No theorem here realizes them by allowed
coupled moves in an arbitrary tournament, preserves an actual port guard,
or transfers the arithmetic parent descent to Collatz.

## 6. Calkin--Wilf addresses are exact electrical operations

For an effective conductance `c=a/b>0`, adding a unit conductance in
parallel sends `c` to `c+1`; adding a unit resistance in series sends it
to `c/(1+c)`. These are exactly the right and left Calkin--Wilf children.
Starting from one unit edge, the subtractive Euclidean decoder constructs
every positive rational conductance using these operations.

| conductance | Calkin--Wilf word | binary heap index |
|---|---|---:|
| `4/5` | `RRRL` | 30 |
| `7/15` | `RRRRRRLL` | 508 |
| `8/15` | `LLLLLLRL` | 258 |
| `8/7` | `LLLLLLR` | 129 |

The first row is the actual shared incidence reduction: four unit parallel
channels have conductance four; placing the one-unit channel in series
gives `4/5`. The other rows are exact network realizations of the listed
rationals, not replacements for the native stationary mass data.

The Stern--Brocot paths to `7/15` and `8/15` are instead `LLRRRRRR` and
`LRLLLLLL`. Reciprocal and complement must also be distinguished:
`x -> 1/x` is the global rational-tree reciprocal mirror, whereas
`x -> 1-x` reflects the `(0,1)` Stern--Brocot subtree. A Calkin--Wilf binary
heap index is not an actual native update clock or an electrical mixing
time. The edge count in the constructed unit network is yet another
quantity (one plus word length).

### Recovering the older dyadic denominator 15 without losing time

The one-step Haar valuation law is `Pr(v=j)=2^(-j)`, `j>=1`. Writing
`v=4K+r`, `r in {1,2,3,4}`, gives

    Pr(r) = 2^(4-r)/15,
    Pr(K=k | r) = (15/16) 16^(-k).                  (9)

The first coordinate probability is `8/15`, and the other three sum to
`7/15`. For an odd start ending in 9, the next odd Collatz value also ends
in 9 exactly when `v=1 mod4`: modulo5, `3n+1=3` and the required equation
is `2^v=2 mod5`. Thus (9) recovers that historical `8/15` directly.

There is an exact fair-bit sampler for (9): inspect successive four-bit
blocks until the first block containing a zero; retain the position `r` of
its first zero and the number `K` of preceding all-one blocks. Discarding
`K` keeps the four-symbol output law but loses the actual valuation.
For example `E[v | r=1]=1+4/15=19/15`. A finite nonzero-block receipt is
not a finite deterministic bound on every source's waiting time.

This is an exact shared representation technique: collapse repeated
blocks by summing their mass and retain their integer count to recover the
clock. It does not make the native incidence chain conjugate to Collatz.
In the former `8/15` is a spectral quantity; in this older instance it is
a particular transition probability under the specified Haar measure.

## 7. What this leaves for actual guard completion

The new incoming
[cyclotomic affine compiler](../../session43_observers_phase_guards_20261023/cyclotomic-lattice/NOTE.md)
gives a parallel warning: its entire universal affine congruence space
contains positive nonrealizable forward-polynomial candidates. It expressly
does not test incidence, spectral gaps or port compatibility. Here the
analogous missing data are visible directly as a deleted bridge word or
discarded stationary mass. There is no identified map equating its
polynomial lattice with these incidence conductances.

The constructive next test is therefore to extract the actual guarded
incidence table, certify surviving bridge support and the internal modes,
and retain `L` AND `D_r` when eliminating targets. A Farey/ternary address
can compress an exact rational parameter after these obligations are
settled; it cannot supply them. The positive all-size result of this packet
is the explicit spectral-coordinate tree and its inverse, while the exact
native connection is the shared `4/5` conductor behind the two order-four
fractions and its one-word guard failure.

## Reproduction

```text
python -B 05-knowledge/results/native_guard_moves_farey_20261023/farey/check.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/farey/check.py
```

The checker reconstructs all eight actual native words, verifies both
weighted couplings and the deleted-word hostile, checks every reduced
odd-denominator coordinate with denominator at most 151 by its terminating
parent and replay, enumerates the ternary tree through depth seven, and
verifies rational-network and address controls. The proof supplies the
unbounded claims; this finite universe does not assert arbitrary native
guard coverage.
