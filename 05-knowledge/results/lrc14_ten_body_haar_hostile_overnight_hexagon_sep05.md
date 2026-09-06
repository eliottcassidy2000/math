# The ten-body Haar floor fails; retain the component and its first owner exit

**Status: REFUTED universal `6/77` ten-body floor + FINITE-EXACT minimal
height bank; constructive synchronization corollary of proved THM-4335.**
This note concerns the **actual** body of the new scale-three consumer,
not another triple-network bound. No LRC counterexample or new fixed-clock
family closure is claimed.

Let `G_C={y: ||cy||>=1/14 for every c in C}`. The candidate

```text
|C|=10 and G_C nonempty  ==>  mu(G_C)>=6/77
```

is false. For the ten-speed body

```text
C_*=(1,2,3,5,7,8,9,11,12,13),
mu(G_(C_*))=14249/252252 <6/77.                         (1)
```

Both successive exact safe-band intersections and an independent global
endpoint sweep reproduce the **entire closed safe set**, including its
zero-length components. For tails `(1,5,11)`, even the tail-adaptive Haar
threshold equals `6/77`, so the same body fails that sufficient gate.

Nevertheless `3C_* union {1,5,11}` is safe at `x=1/14`. In fact every
body in the small bank below, with arbitrary odd tails, is already safe at
that fixed clock: `c<=13` ensures `14` does not divide `3c`, and an odd tail
is never divisible by fourteen. This is a hostile to the measure implication,
not a novel family closure or an LRC counterexample.

## Inheritance and the changed concept board

The immediate source is the conditional
[universal triple-bound Haar consumer](lrc14_universal_haar_consumer_empty_core_certificate_sep06.md).
Its correct implication is `mu(G_C)>=min_i E_i(T) => G_(3C union T)`
nonempty. It does not prove the required body floor. The closest proved
replacement mechanism is
[THM-4335 — owner permutations and component addresses](../../01-canon/theorems/THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md),
not an additional estimate of THM-4409's global three-speed network.

The canonical synchronization hostile has `C={1,...,10}`, tails `(1,5,11)`,
and the good body phase `y=2/11` whose three lifts are all spoiled. The
corrected near miss is importing a sufficient Haar floor into an arbitrary
body. The least-used sidecar is a **closed component's actual left endpoint
and the first effective-tooth exit**, including singleton components and
equality-safe endpoints. The concept board is: body measure; closed component;
owner permutation; first exit; rational event coset; actual entry map.

The strict-margin theorem
[THM-780 — quantitative safe-measure floor](../../01-canon/theorems/THM-780-haar-compactness-safe-measure-floor.md)
still supplies a genuine universal positive floor. Cited LRC through thirteen
total runners gives a ten-speed body a `1/11`-deep phase; at target `1/14`,
its explicit phase-partition bound is `mu(G_C)>=52^(-10)`. The false
`6/77` floor is therefore not replaced by an assertion of vanishing measure.
THM-780's stronger anchored Bohr packet retains event incidence that its
scalar floor discards.

The operational firewall was checked against
[THM-4300 — size-preserving response staircase](../../01-canon/theorems/THM-4300-lrc14-size-preserving-response-staircase-and-index-297-ideal.md).
That theorem acts on a specific labelled thirty-speed pool, its nine-body
deletion obligations, threshold `4/63`, and a realized response carrier.
It supplies no entry map sending an arbitrary actual ten-body into that
pool. Likewise THM-4370's seven-wall owner masks belong to an anchored
`2+12` row and cannot be identified with these three-sheet body lifts.
The quick diameter argument `max T>=11 max C` is also already proved:
THM-4330's unit equality-wall safe-band lemma `(MC9b)--(MC9f)` recovers
the older scale-three cone. It is not a new consequence to promote here.

## 1. Minimal finite bank and an honest small-clock control

The complete bank is every ten-element subset of `{1,...,13}`, exactly
`286` bodies. Exactly twelve have measure below `6/77`. The first in
height-then-lexicographic order is

```text
C=(1,2,3,4,5,7,8,9,11,13),
mu(G_C)=21514/315315 <6/77.                             (2)
```

Consequently height thirteen is the smallest possible maximum speed for
a ten-element positive-integer body violating this floor. This conclusion
uses the complete bank, including all smaller heights; it is not a claim
that (2) minimizes the measure globally. The minimum **inside this bank**
is (1).

The explicit row in (1) is closed by a fixed clock, so it must not be sold
as progress on an entry-filtered residual. A cheap follow-up retains tails
`T=(1,5,11)` and all ten-bodies with maximum at most fourteen for which no
reduced phase of denominator at most fourteen is safe for `3C union T`.
There are exactly `209` such bodies. They all pass the Haar gate in this
bounded bank; their smallest measure is

```text
4163/51480 >6/77,
C=(1,2,3,8,9,10,11,12,13,14).                          (3)
```

For `q<=14`, a reduced phase `p/q` is safe exactly when `q` divides no row
speed, since every nonzero residue has distance at least `1/q>=1/14`.
Thus this small-clock filter is exact. Within the stated height bound it
forces body members `8,10,13,14` and a multiple of three. The `209` cases
are generated completely, not sampled. Equation (3) is a finite positive
control only, not a repaired universal Haar theorem.

Incoming commit `910dad3281880a9ec940d28a24fb784892b66c76` independently
recovers the same height-thirteen bank in
[the empty-core body report](lrc14_haar_body_empty_core_sep06.md), and its
Section 2 refutes even the global small-clock-filtered repair. For
`C=(1,3,4,10,11,13,14,16,17,18)` and `T=(1,5,11)`, no reduced clock of
denominator at most fourteen is safe, yet
`mu(G_C)=534689/7796880<6/77`. Our height-fourteen positive control is
therefore genuinely bounded. The incoming report retains a later exact
physical witness, not an LRC counterexample.

## 2. Exact component synchronization using one endpoint and one exit

Let `C` be any finite nonempty positive integer body, let `T` have exactly
three distinct positive units modulo three, and let `I=[L,R]` be a connected
component of its actual closed safe set `G_C`. Zero-length components are
allowed. Oddness of the tails is not needed for this owner calculation,
although it remains part of the universal network theorem's typed domain.

At the left endpoint define, whenever the tail is active,

```text
n_w=nint(wL),  active iff |wL-n_w|<3/14,
kappa_w=-n_w w^(-1) mod3.                              (4)
```

If the active owners fail to cover `{0,1,2}`, choose a missing `j`. Then
`x=(L+j)/3` is a physical safe phase for `3C union T`.

Otherwise all three owners form a permutation. Define the first right exit
from these **addressed** effective teeth:

```text
U=min_(w in T) (n_w+3/14)/w.                           (5)
```

Then the component has an unspoiled lift if and only if

```text
U<=R.                                                 (6)
```

If (6) holds, choose a tail attaining the minimum and let `j=kappa_w`.
The explicit safe phase is `x=(U+j)/3`. Until this first event every other
active owner is unchanged; at the event the selected owner disappears
because danger is strict. Tied first exits only free more sheets.
If `U>R`, every point of `I` stays strictly inside the same three addressed
teeth, whose owners remain a permutation, and the entire component is
spoiled. This proves both directions. A singleton component is resolved by
the left-endpoint test; an active first exit is strictly to its right.

This is an operational corollary of THM-4335's proved owner/constant-address
argument, not a fresh claim to that argument. Given the actual body
components, it is an **exact decision and witness compiler** using three
owner evaluations and one extremal event per component. The operation count
does not grow with tail height, though integer bit lengths do. It avoids
enumerating all the tail teeth, and is strictly more informative than the
global Haar mass alone.

Two controls expose the endpoint semantics. With `T=(1,5,11)`:

```text
C={6},  I=[5/28,9/28],       U=31/154 is interior;
C={99}, I=[89/462,31/154],    U=R=31/154.
```

Both give the physical witness `x=113/154`, where tail eleven is exactly
at clearance `1/14`. Replacing `<=` by `<` in (6) loses the second genuine
certificate. In the canonical `C={1,...,10}` hostile, the component
`[5/28,13/70]` containing `2/11` has `U=31/154>R` and is entirely spoiled.
Moving to another body component, rather than relabelling that witness,
is necessary. The body (1) also contains the isolated safe component
`{3/14}`; measure-only summaries discard it completely.

## 3. The smallest remaining native event obligation

For a tail `w`, a body phase on either event coset

```text
w y = +/-3/14 mod1                                   (7)
```

already has an inactive tail, so the remaining two singleton owners cannot
cover all three lifts. Any actual body component reaching (7) gives an
immediate physical witness. The unresolved task is to force one such event,
an owner collision, or an inactive owner on an **actual entry-produced**
body; the global triple bound does not do this.

The forgotten operation available from THM-780 is to retain a deep body
phase `y0` together with its simultaneous Bohr-return packet. Its event-coset
interface tests (7) as an inhomogeneous finite cyclic code with the augmented
relation lattice and target phase retained. It does not assert that every
chosen deep phase's packet meets a useful event. That is precisely the
selection/entry obligation, now expressed using native objects rather than
a false numerical body floor or an unrelated septimal wall carrier.

## Reproduction

```bash
python3 -B 04-computation/lrc14_ten_body_haar_hostile_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_ten_body_haar_hostile_overnight_hexagon_sep05.py
```

The [script](../../04-computation/lrc14_ten_body_haar_hostile_overnight_hexagon_sep05.py)
compares full closed safe sets by two independent constructions for all
286 bodies, checks the separate 209-row fixed-clock-filtered bank, and
compares the endpoint/first-exit compiler against a literal effective-tooth
endpoint sweep on named bodies and tails. Every returned physical witness
is checked directly against all speeds; the interior, equality, entirely
spoiled, and singleton controls remain explicit. The
[output](lrc14_ten_body_haar_hostile_overnight_hexagon_sep05.out)
records the exact components and finite semantic digest. No shared navigation
or scarce theorem namespace is modified by this note.

Normal and optimized runs agree byte-for-byte: `644` component decisions
and `1,579` explicit gates. Frozen raw-LF SHA256 values are

```text
script 39bd4707a55cb5de792e3fa07e7d02730c2cae81f5cc6a3cc4931d25fd862f01
output 9f16f6c203fd5f053a37efc96919a1573f997bb02a45b111b11a3c3f8cdd7d80
```

Independent `observer_collision` audit: **PASS**. The referee reconstructed
the hostile's closed safe set by rational residue cells, independently ran
all `1001` ten-subsets through height fourteen for the `209`-row filtered
bank, and verified the first-exit proof and interior/equality/blocked
physical-witness controls. This validates the stated finite bank and inherited
compiler, not an arbitrary entry-produced body's synchronization.
