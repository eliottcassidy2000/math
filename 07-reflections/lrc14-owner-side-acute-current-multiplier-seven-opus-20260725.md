---
source: codex-2026-07-25-owner-side-acute-current
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On a scalar
  literal THM-2305 word, fix one deepest-bank boundary orientation,
  address/residue cell, and one side of the selected source owner's
  centered danger coordinate. At an aligned frequency K=m c_j, every
  nonempty such cell has nonzero scalar endpoint current for 1<=m<=7.
  The m=7 endpoint is exact and qualitative: strict comb interiors put
  all phases in one open half-plane. The corresponding THM-2355
  autocorrelation cone is acute for m<=3, and two one-sided cells obey
  the strict cross-cone support law when m+n<=6. Thus in the semi-hard
  primitive cases m=4,5,6, a vanishing orientation/address current must use
  and balance both owner sides, with an active endpoint beyond the exact
  inner threshold (7-m)/(14m) on each side. The phase-cone argument is sharp at m=8,
  but the supplied m=8 control is geometric rather than a canonical full
  THM-2305 word. No LRC profile is excluded.
depends_on:
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2304-deepest-boundary-cyclotomic-current-separation
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2355-quadratic-grouped-current-repair
related:
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2300-small-owner-multipliers-force-same-character-relation-multiples
  - THM-2358-real-endpoint-sign-cone-and-galois-colour-hostile
  - lrc14-two-frequency-orientation-current-leakage-opus-20260725
  - lrc14-quadratic-phase-tree-minimal-calibration-opus-20260725
script: 04-computation/lrc14_owner_side_acute_current_probe.py
output: 05-knowledge/results/lrc14_owner_side_acute_current_probe.out
script_sha256: 7806a56b9095bdb69b58d64d2bdc2b3901ecc526d0c939c09a2ff81079d36e5b
output_sha256: a5763b0c17054a3091ca5b4ffde887299047045ce767eded59a5fb89bbcf9886
hash_basis: working-tree bytes (LF)
---

# The owner-side current stays nonzero through multiplier seven

THM-2355 leaves two honest phase repairs: lawful quadratic pair probes or
an acute root/target cone. The existing root-character cone is not enough,
because it lives before integration over the continuous terminal coordinate.
There is nevertheless a scalar acute cone already present in the literal
word: the danger arc of the selected source owner.

Classifying the existing endpoint terms by the two geometric sides of that
arc doubles the labelled state without introducing a new gate. It yields an
exact noncancellation theorem through the semi-hard primitive range singled
out by THM-2300. This is a provenance partition, not multiplication by a
half-tooth indicator.

## 1. Scalar endpoint-current setup

Let `F` be a scalar literal-source THM-2305 word for source owner `j`, and
write

```text
c_j=13^lambda u_j,

K=m c_j,                         m>=1.              (1)
```

The source word lies in the owner's strict danger comb:

```text
F subset D_(c_j).                                    (2)
```

Fix one deepest boundary bank `a_*`, one of its orientations
`epsilon in {+1,-1}`, and any further relative-address and prime-to-thirteen
residue cell. Its endpoints have the form

```text
x=(14r+epsilon)/(14a_*).                             (3)
```

Let `L_0` be the lower Boolean gate. With THM-2304's right-sided convention,
an active endpoint contributes

```text
Delta(L_0 1_(D_(a_*)^c))(x)=epsilon L_0(x+).         (4)
```

The opposite Boolean polarity changes all terms by one common minus sign.
Thus after fixing the orientation, the active weights

```text
w_x=L_0(x+)
```

are positive.

For each endpoint put

```text
theta_x = the centered residue of c_j x
          in (-1/14,1/14).                           (5)
```

Split the cell by the owner-side label

```text
eta=sign(theta_x) in {+,-},                          (6)
```

with the isolated center case `theta_x=0` kept separately. Up to one common
nonzero orientation/polarity sign, its scalar current is

```text
C_(epsilon,eta)
 =sum_x w_x exp(-2*pi*i Kx)
 =sum_x w_x exp(-2*pi*i m theta_x).                 (7)
```

The second equality uses `c_jx in Z+theta_x`.

## 2. Open-half-plane noncancellation

On the positive owner side,

```text
0<theta_x<1/14,
```

so the arguments in (7) lie in

```text
(-pi*m/7,0).                                        (8)
```

On the negative owner side they lie in the reflected interval

```text
(0,pi*m/7).                                         (9)
```

For `1<=m<=6`, either cone can be rotated to have half-width

```text
m*pi/14 < pi/2.
```

Every rotated phasor then has strictly positive real part. Since the weights
are positive, the current cannot vanish.

At `m=7`, (8) and (9) are the strict open half-planes `(-pi,0)` and
`(0,pi)`. Every term has the same strict sign of imaginary part, so the sum
again cannot vanish. There is no uniform quantitative margin as an endpoint
approaches the boundary.

The center cell `theta_x=0` has width zero and is separately nonzero whenever
active. A point with `theta_x=+/-1/14` would invalidate the equality endpoint,
but the strict danger convention and THM-2304's nontrivial conductor
separation exclude coincidence with that lower owner boundary in the current
scope.

We have proved:

> **Owner-side multiplier-seven lemma.** Every nonempty scalar
> deepest-orientation/address/residue/owner-side cell has nonzero endpoint
> current at `K=m c_j` for every `1<=m<=7`.

All scalar combs and their Boolean words are reflection invariant, and
`T(x)=13x` commutes with reflection. Hence a positive literal word has equal
positive total mass on its two owner sides. The theorem is not claiming that
each finer address cell meets both sides; precisely that failure is useful.

## 3. Exact consequence at the live multipliers

Without the side split, the owner cone has width

```text
2*pi*m/7.
```

It is contained in an open half-plane only for `m<=3`. This recovers the
known low-multiplier noncancellation directly:

```text
Re Fhat(m c_j)
 >=cos(pi*m/7) measure(F)>0,          1<=m<=3.       (10)
```

For a primitive shell multiplier in the semi-hard cases `m=4,5,6`, the
unsplit word may cancel, but the owner-side lemma shows the exact anatomy.
THM-2300 separately forces some same-character multiple in this range of
the analysis; it does not assert that the forced relation multiple is the
primitive `K=m c_j` itself.

```text
a vanishing nontrivial separated orientation/address current at m=4,5,6

  must contain active endpoints from both owner sides

  and those two nonzero side currents must balance exactly.           (11)
```

Consequently, one **private one-sided nontrivial address** closes that
current channel. This is a smaller target than a complete component phase
tree or a general chirped Gram reconstruction.

It is important to retain the orientation before applying (11). Opposite
boundary orientations have opposite real jump signs and can cancel even
when each orientation/owner-side cell is separately acute.

Here and below, “nontrivial separated” excludes the isolated
`theta=0` center term, which belongs to the lower-conductor center cell.
There is also a geometric tax on any such two-side cancellation. Let
`alpha_+ in (-m*pi/7,0)` and `alpha_- in (0,m*pi/7)` be the arguments of
the two nonzero side resultants. If their sum vanishes, then

```text
alpha_-=alpha_++pi.
```

The two cone inclusions force

```text
alpha_+<-(7-m)pi/7,

alpha_-> (7-m)pi/7.                                (11a)
```

A positive weighted sum has argument inside the angular hull of its
summands. Hence each owner side must contain an active endpoint in its
outer strip

```text
|theta|>(7-m)/(14m).                                (11b)
```

For the live multipliers, the exact thresholds are

```text
m=4:  3/56,

m=5:  1/35,

m=6:  1/84.                                        (11c)
```

Thus a channel also closes if either owner side is confined to the
corresponding inner strip, even when both sides are present. This
**outer-side tax** is weaker than a private-address theorem but is directly
checkable in the finite endpoint atlas.

## 4. Match to THM-2355

The one-sided phase-cone width is

```text
w_m=m*pi/7.                                         (12)
```

THM-2355's autocorrelation survivor requires `w_m<pi/2`, hence exactly

```text
m<=3.                                               (13)
```

For two one-sided positive cells with multipliers `m,n`, the sum of cone
widths is

```text
(m+n)pi/7.
```

The strict cross-correlation support theorem therefore applies uniformly
when

```text
m+n<=6.                                             (14)
```

At `m+n=7`, an open-cone refinement is possible when both strict endpoint
cones and positive weights are retained, but it is not needed here and is
not advertised across coincident boundary cells.

Thus THM-2355 and the present lemma divide the work cleanly:

```text
m<=3:
  unsplit scalar current and one-sided autocorrelation are acute;

m=4,5,6:
  each side is nonzero, but a lawful cross-side probe or a private side
  is needed;

m=7:
  each side remains qualitatively nonzero, but the current relation
  multiplier bank is already at its known sharp boundary;

m>=8:
  even one side can cancel.                          (15)
```

## 5. Sharp hostile controls

### First unsplit failure: `m=4`

Let `a=1/(4m)=1/16`, choose a sufficiently small `epsilon>0`, and take two
equal intervals centered at `+/-a`, both inside `D_1`. At either endpoint
orientation, paired positions differ by

```text
2a=1/(2m).
```

Their frequency-`m` phases are antipodal, so both orientation sums vanish.
This is the local geometry in THM-2303's labelled profile-`(1,3,5)`
handoff. It is the smallest unsplit failure because (10) covers `m<=3`.
It is local, not a global scalar-cover counterexample.

### First one-sided phase-cone failure: `m=8`

The equal intervals

```text
(15/4096,17/4096),

(271/4096,273/4096)                                 (16)
```

lie entirely in `(0,1/14)`. Corresponding endpoints differ by

```text
256/4096=1/16=1/(2m),
```

so both entering and exiting orientation currents cancel exactly.

The same antipodal geometry can lie on one actual deepest-comb side whenever
the scale ratio `A=a_*/c_j` is divisible by `16`: take one endpoint with

```text
theta_1=1/(14A),

theta_2=1/16+1/(14A)<1/14.                          (17)
```

This proves sharpness of the phase-cone argument. It does **not** prove that
a canonical lower Boolean gate activates those two endpoints with equal
weights, so (16)--(17) are not promoted as a THM-2305 word hostile.

## 6. Connection contract

```text
source:
  a scalar literal THM-2305 word and THM-2304's oriented deepest endpoint
  current;

target:
  one nonzero primitive-colour scalar current in the semi-hard m<=6 bank;

map:
  split each oriented address/residue cell by the sign of the centered
  source-owner coordinate theta;

preserved:
  scalar Boolean legality because no new Boolean cut is imposed, source
  owner, exact terminal word, prescribed time, deepest provenance,
  orientation, address, residue, and shell multiplier;

destroyed by recombination:
  which owner side supplied the current;

positive theorem:
  every nonempty side cell is nonzero through m=7;

remaining cancellation:
  exact balance of the two nonzero owner-side currents;

cheapest decisive test:
  prove that one nontrivial deepest address in the marked THM-2349/2305
  word is private to one owner side for one of m=4,5,6, or construct a
  lawful quadratic twist that separates the two sides while retaining the
  same word and shell.                                  (18)
```

Artificially cutting into still narrower phase bins does not solve the
problem: it either introduces non-scalar endpoints, breaking conductor
separation, or merely decomposes a grouped current whose bins may cancel.
The owner-side label is canonical because it is the sign of an existing
centered scalar phase coordinate; it is used only to classify existing
endpoint terms and is not inserted as a new scalar comb factor.

PROVED THM-2358 (promoted at `59c933aae`) supplies a complementary
stopping rule, not a counterexample to this lemma. Its inverse target array
has opposite real signs even though its base coefficient, bare endpoint,
and deep leg are nonzero. The operations differ twice: one active factor is
the safe complement `g(26x+t/13)`, not a danger-owner restriction, and the
inverse target transform recombines shift packets with roots of unity.
After dephasing `d(13x+s/13)` at frequency `26`, a fixed owner-side packet
does lie in the present `m=2` cone, but target character `a` leaves the
shift-dependent factor

```text
zeta^((2-a)s).
```

Those thirteen cones rotate before they are added, while
`ghat(n)=-dhat(n)` for `n!=0` introduces a second signed mixer. Therefore
one must not transport this acute physical cone through target inversion
without a cone-preserving intertwiner on the *same* fixed
orientation/address/owner-side cell. Positivity and evenness of the original
factors, word nonconstancy, a marked triangle, and nonzero endpoint/deep
legs do not supply that intertwiner. This is precisely why (18) retains the
side and orientation labels until after selection.

Indeed, the THM-2358 frequency constraint is

```text
13n+26k+169l=26,       equivalently n+2k+13l=2.
```

Its inverse target support is the affine line `a+2b=2 (mod 13)`.
The positive `U_0` channel `(a,b)=(2,0)` assigns the frequency to the
danger-owner axis, while the negative `U_1` channel `(0,1)` assigns it to
the safe-complement axis. Their opposite signs are thus an exact
allocation/complement sign, not cancellation inside one positive
owner-side endpoint packet.

An independent hostile audit accepts every load-bearing part of THM-2358:
the real-cone/support-difference lemma; the two primitive `N=91` colour
signs; the target DFT normalization and affine support line; the exact
`26,25,25,11,11` interval ledger; the `1/91` transported-word tax; and all
three marked Fourier legs. An independent piecewise integration gives

```text
U_0= 0.0967462347008684,    U_1=-0.0187905630647440,

L_(0,0)=0.0553776032918977,
F_(0,0)=0.0153661256411735,
deep=0.138109483615507.
```

The only editorial repair is to describe the physical factors as
**nonnegative** real/even indicators and complements, not strictly positive
factors: `g=1-d` has zeros. The correct promotion scope is a local
three-coordinate formal-carrier obstruction. It proves that target inversion
need not preserve one acute endpoint cone; it does not give a canonical
LRC row, a profile decrement, or a root--target cone-preserving intertwiner.

No private-side theorem is presently available, so no LRC profile is
excluded.

## 7. Exact verification

Run

```text
python3 04-computation/lrc14_owner_side_acute_current_probe.py
python3 -O 04-computation/lrc14_owner_side_acute_current_probe.py
```

The companion checks every integer threshold in (10)--(15), both exact
orientation antipodalities in the `m=4` and `m=8` controls, the strict
one-sided interval containment, and the deepest-comb geometry (17), using
only rational arithmetic. Normal, optimized, and stored transcripts match
after LF normalization (the Windows process streams use CRLF).

An independent hostile audit checked the right-sided jump convention,
positive weights, strict/open endpoint treatment at `m=7`, the center-cell
boundary, the autocorrelation and two-cone thresholds, and the noncanonical
scope of the `m=8` sharpness witness.
