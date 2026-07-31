# LRC14 half-step chains in the rooted Hasse chart

**FINITE-EXACT / PROOF-COMPLETE RESEARCH NOTE.**  This is an unnumbered
handoff note, not a new canonical theorem.  It combines the exact block data
of THM-2818 with the rooted normalized-Hasse translator isolated in
THM-2820/THM-2201.  The positive coefficient-side result is exact.  The
negative result applies to a bridge determined by the scalar chain
polynomial; it does not rule out a future decorated physical cospan and does
not prove LRC(14).

## 1. Exact rooted chain polynomials

Let `X` denote one physical half-step

```text
h=T/(2*13^5).
```

Root each exceptional target-twelve row at the left endpoint of its first
block.  If one block has length `L`, its raw word is

```text
[L]_X=(1-X^L)/(1-X),
```

and its live word, which starts live and alternates, is

```text
E_m(X)=(1-X^(2m))/(1-X^2),       m=ceil(L/2).          (1)
```

Consecutive block roots are separated by the previous length minus one
plus the proved interblock jump `50`.  The exact root offsets are therefore

```text
clock 1: (0,194),
clock 2: (0,192,530),
clock 3: (0,192,530).                                  (2)
```

The three positive live-copy polynomials are

```text
L_1=E_73+X^194 E_48,

L_2=E_72+X^192 E_145+X^530 E_48,

L_3=E_72+X^192 E_145+X^530 E_37.                       (3)
```

The dead polynomials are

```text
D_1=X E_72+X^195 E_48,

D_2=X E_71+X^193 E_144+X^531 E_48,

D_3=X E_71+X^193 E_144+X^531 E_37.                     (4)
```

Thus the raw word is `R_e=L_e+D_e`, while the signed alternating contrast
is `A_e=L_e-D_e=2L_e-R_e`.  Equations `(2)--(4)` encode all three finite
chain ends; replacing a chain by an infinite `1010...` word would lose the
essential endpoint data.

## 2. The formal rooted Hasse translator is genuinely live

Reduce in

```text
F_13[C_13]=F_13[X]/(X^13-1)
           =F_13[eta]/(eta^13),       eta=X-1.
```

For `P=sum c_n X^n`,

```text
J_0(P)=sum c_n,              J_1(P)=sum n c_n,
b(P)=J_1(P)/J_0(P).                                  (5)
```

The exact first-jet ledger is

| clock | object | exact `(J0,J1)` | mod-13 `(J0,J1)` | `b` |
|---:|---|---:|---:|---:|
| 1 | raw | `(241,33624)` | `(7,6)` | `12` |
| 1 | live | `(121,16824)` | `(4,2)` | `7` |
| 1 | dead | `(120,16800)` | `(3,4)` | `10` |
| 1 | contrast | `(1,24)` | `(1,11)` | `11` |
| 2 | raw | `(528,162697)` | `(8,2)` | `10` |
| 2 | live | `(265,81528)` | `(5,5)` | `1` |
| 2 | dead | `(263,81169)` | `(3,10)` | `12` |
| 2 | contrast | `(2,359)` | `(2,8)` | `4` |
| 3 | raw | `(506,149178)` | `(12,3)` | `10` |
| 3 | live | `(254,74774)` | `(7,11)` | `9` |
| 3 | dead | `(252,74404)` | `(5,5)` | `1` |
| 3 | contrast | `(2,370)` | `(2,6)` | `3` |

Every displayed augmentation is nonzero.  Hence every polynomial is a unit
of the truncated local algebra.  In particular each live polynomial has a
regular thirteen-element translation orbit, and

```text
b(X^a L_e)=b(L_e)+a.                                  (6)
```

Incoming THM-2839 gives the structural reason: on the prime-power group
`C_13`, nonzero augmentation makes each integral mask a group-algebra unit,
so all thirteen Fourier characters and the full regular translate span
survive.  The companion below independently constructs the truncated
inverse and checks the full orbit.  THM-2839 also marks the positivity
boundary: because every live mask is nonnegative and has more than one
support point, no convolution inverse can be nonnegative.  Algebraic
recoverability here is necessarily signed and is not yet a positive
physical transport.

This is the positive bridge: once a lawful rooted physical translation is
supplied, every exceptional chain is an exact affine Hasse translator.
There is no denominator failure and no hidden cancellation.  In the
live-head roots, the three normalized coordinates are

```text
(b(L_1),b(L_2),b(L_3))=(7,1,9),                       (7)
```

all nonzero.

The qualification is immediate from `(6)`.  Absolute nonzeroness is not
root-gauge invariant: translating the three roots by `6,12,4` respectively
makes the corresponding `b` equal to zero.  What is intrinsic is the
unit-slope response to a **named** translation, not the absolute values in
`(7)`.

The finite chain boundary can also be read exactly from `(3)--(4)`:

```text
X L_1-D_1 =X^145,

X L_2-D_2 =X^143+X^481,

X L_3-D_3 =X^143+X^481.                              (7a)
```

Thus a one-half-step translation exchanges live and dead at every interior
tooth.  Its only defects are the terminal teeth of the odd blocks.  This is
the precise finite, rather than infinite, form of the semantic parity flip.

## 3. The scalar chain does not decide carrier-gauge transversality

THM-2820's two-axis carrier-address quotient observes

```text
q_eff=(q_1-z.L,q_2+z.L).                              (8)
```

The one-variable chain polynomial sees only its first formal motion
`q_1=1`.  It does not retain the second carrier axis or the address
increment.  Fix `z.W=1`.  Two completions of exactly the same observed
chain motion are

```text
L=0,       (z.L,q_1,q_2)=(0,1,0),       q_eff=(1,0),  (9a)

L=W,       (z.L,q_1,q_2)=(1,1,-1),      q_eff=(0,0).  (9b)
```

The second is precisely the marked pure-gauge direction
`(L,q)=(W,(1,-1))`; merely observing `q_eff=0` without that full address
identity would not suffice.  Therefore no
function of `L_e(X)` alone can decide whether its formal unit jet is
gauge-transverse.  The first missing algebraic sidecar is the omitted
carrier axis together with the typed address increment `L`.

## 4. The half-step and target-convolution axes cannot be identified

Every polynomial `(3)` lives inside the one fixed physical cell

```text
target t=12.
```

Multiplication by `X` advances an interval left endpoint by one half-step
while leaving `t=12`.  THM-2771's target generator instead moves between
the thirteen target labels, and its decoder `K_beta` convolves on that
target-label axis.

This gives a sharp equivariance no-go.  The formal orbit

```text
{X^a L_e: a in F_13}
```

is regular by Section 2.  Any target-preserving map sends all thirteen
members to `12`.  An intertwiner with the regular target shift would have
to send consecutive members to consecutive target labels.  Already at one
step it would require

```text
12=12+1=0 mod13,
```

which is impossible.  The exact companion checks all `39` orbit steps.
Thus one cannot simultaneously:

1. preserve the native target-twelve label; and
2. identify half-step translation with target convolution.

The repair is not another scalar normalization.  It must retain an
independent target variable `U`; the chain requires at least a bivariate
object `L_e(X,U)` before target convolution can even be typed.

## 5. Why a common-interval map still lacks its physical decorations

THM-2818 already supplies the relevant hostile ledger.

- The common interval `I` has all six native factors, while the selected
  equal live cofiber copies have masks `011111` and `011101`.
- `I` has delta-zero source and target carrier support; both selected
  cofiber copies have empty source carrier support.
- No endpoint translation identifies `I` with either source copy; even the
  one identity match occurs only on one target mask.
- Equal length, weight, content, and local ancestry do not retain these
  distinctions.

The scalar polynomial `(3)` records only half-step positions and live
coefficients.  It cannot certify preservation of a native-factor mask,
carrier mask, or endpoint origin that it does not contain.  This does not
prove that no enriched map exists; it proves that the scalar-chain functor
has already forgotten data on which the desired predicate depends.

There are two further deck distinctions.  THM-2771's target label is an
endpoint-dipole coordinate, not the physical root deck.  Even after adding
the independent target variable, an affine target-to-root identification
has `13*12=156` choices without an origin or generator orientation, and
still has thirteen origin choices after the generator is fixed.  None is
selected by `(3)` or `(7)`.

## 6. Minimal surviving target

### Status of the nearest-collar proposal

Incoming THM-2825 currently has status

```text
RESERVED / UNPROVED EMPTY STUB.
```

It therefore changes none of the proved conclusions above.  Its proposed
nearest right-to-common collar would translate each interval by exactly one
half-step.  If that geometric statement and its semantic parity audit are
promoted, it would supply a valuable named physical motion and a nearest
geometric origin.  Equation `(7a)` shows the remaining test sharply: an odd
half-step reverses live/dead on the chain interior, with the listed terminal
defects.  Such a collar does not transport the positive live coefficient
unless its promoted sidecar explicitly repairs that parity and the finite
ends.

Consequently THM-2825 could remove the *unnamed-translation* part of the
invoice, but a geometric `+h` collar alone cannot remove the semantic,
carrier-gauge, target-axis, or root-deck parts.  A genuinely positive bridge
would need either:

1. a collision-free even-half-step common collar preserving semantic
   content; or
2. a proved decorated odd-step collar which carries the live/dead toggle
   and terminal defects as part of its map.

The smallest object not killed by the preceding hostiles is a
**decorated birooted chain cospan**.  Concretely, it must retain

```text
sum_(live J)
  [native-factor mask of J,
   source/target carrier masks of J,
   endpoint address and common-interval image,
   local ancestry]
  X^(half-step index of J) U^(target label of J),      (10)
```

together with:

1. a physical cospan from each decorated cofiber copy to a common interval;
2. a common endpoint basepoint;
3. the two-axis carrier/address path needed to evaluate `(8)`; and
4. an affine target-role-to-root-deck intertwiner preserving that cospan.

Each entry is forced independently:

- omit the decorations and the two selected equal copies collapse despite
  different native masks;
- omit `U` and Section 4 forbids target equivariance;
- omit the address/second axis and `(9a)--(9b)` give opposite transversality
  verdicts;
- omit the affine deck map and the target/root numerals remain untyped.

The new positive fact is therefore useful but conditional:

```text
exceptional chain + lawful decorated transverse transport
  ==> exact affine Hasse label,

exceptional chain polynomial alone
  =/=> common physical atom, target convolution, or root-deck jet.       (11)
```

The lightweight exact companion is

```text
04-computation/lrc14_halfstep_chain_hasse_transport_probe_20260728.py
05-knowledge/results/lrc14_halfstep_chain_hasse_transport_probe_20260728.out
```

It reconstructs `(2)--(7)`, inverts every chain polynomial in the local
algebra, verifies all translation orbits, exhibits the gauge-hostile pair,
and checks the target-axis equivariance failure.  It uses no large physical
build.
