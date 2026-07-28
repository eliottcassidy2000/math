# Delayed-tail SCC diagnostic and an intrinsic guard-debt cycle

**Status: VERIFIED-EXACT ADDENDUM / NEW CARRIER SCOUT.**  This is not a canon
theorem, a row exclusion, or an LRC(14) proof.  It takes the newly proved
[THM-2693 delayed depth-four nilpotence](../01-canon/theorems/THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence.md)
as the current truth and asks which factor must change before an odometer lift
is worth attaching.

The closest proved mechanism is THM-2693: the fixed target-`a` delayed word is
positive through three states and empty at four, uniformly in every integer
lift.  The canonical hostile is the nonresonant near-carry period-four tail in
the affine periodic-tail scout: rail/present data survive, but every state
misses target `a`.  The least-used sidecar is the guard/source edge debt that
remembers how unshifted source danger propagates to the next guard phase.
Keeping that sidecar produces the new exact signal below.  Throughout, `a`
and `b` are delayed-tooth/cospan labels.  They are not THM-2305 terminal-word
labels.

## 1. Exact universe and automaton

On the canonical row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5),
R=13^6,
```

write `Rx=N+y`.  Every integer affine handoff has the skew-product law

```text
y'={13y},
N'=13N+k+floor(13y) mod R.                               (1)
```

The tail map is therefore `B(y)={13y}`, independently of `k`.  Four delayed
base words were rebuilt by exact interval arithmetic on the canonical grid

```text
T=297836897838480:

target a: D_(13^3),       safe at 2*13^5;
target b: D_(2*13^5),     safe at 13^3;

guard safe or guard danger;
safe at 14,27,40,53,66 in every case.                    (2)
```

The source-speed-`13` factor is then restored in each of its seven clock
positions.  The seven complements have union exactly equal to the uncut base
word: at every strict point one clock tooth is dangerous and the other six
are safe.  Thus projecting away the clock label loses a six-element fibre but
does not create base support.

The exact single-word censuses are:

| target | guard | intervals | grid mass |
|---|---|---:|---:|
| `a` | safe | 32,578 | 10,786,297,480,556 |
| `a` | danger | 14,906 | 4,936,568,463,918 |
| `b` | safe | 188,634 | 10,807,572,370,920 |
| `b` | danger | 86,950 | 4,981,755,322,960 |

The two guard rows for either target are disjoint and their union is exactly
the word with the guard factor deleted.

## 2. Recurrent SCCs of one repeated target tooth

This part has a two-line analytic reduction.  For target multiplicity
`m in {1,2}`, consider a recurrent tail satisfying

```text
||m B^j y|| < 1/14                    for every j>=0.     (3)
```

Write `y=c+e`, where `c` is a tooth centre and `|me|<1/14`.  Because
`13c=c mod 1` for the relevant centres and

```text
|13me|<13/14,
```

the next copy of `(3)` cannot use a nonzero wrap: either wrap would require
the excluded endpoint `13/14`.  Hence `|me|<1/(14*13)`.  Iteration gives
`|me|<1/(14*13^j)` for every `j`, so `e=0`.  The raw recurrent graphs are

```text
target a:  0 -> 0;
target b:  0 -> 0,    1/2 -> 1/2.                         (4)
```

The included half-open left endpoints map to excluded right endpoints, so
there is no endpoint-only cycle hidden behind the strict argument.  Speed 14
is dangerous at all vertices in `(4)`:

```text
14*0=0,                 14*(1/2)=7.                       (5)
```

Therefore its inherited safe factor deletes every recurrent vertex.  Guard
complementation cannot open an SCC.

This also makes the one-factor-deletion statement exact within universe
`(2)`.  Deleting any non-target factor leaves a centre-killing safe factor.
If speed 14 itself is deleted, speed 27 still kills the target-`a` centre and
the even speed 40 still kills both target-`b` centres.  In contrast, deleting
the repeated target factor opens exact fixed tails:

```text
guard safe:      y=5/12, 7/12;
guard danger:    y=1/12, 11/12.                           (6)
```

Each value in `(6)` is fixed by `B`, passes every remaining factor strictly,
and has six legal owner-clock labels for either target choice.

The strongest one-factor-deletion controls already survive the present
sidecar and an odometer lift:

```text
guard safe, y=5/12:
  N=(744626,769107),  k=(742582,399848),
  clocks 0->1,1->0;

guard danger, y=1/12:
  N=(748707,826230),  k=(746656,4488143),
  clocks 0->2,2->0.                                      (7)
```

All four lifts are nonzero modulo 13.  The selected delayed clock equals the
current shallow clock, the dynamically typed THM-2640 present factor is
strict at both events, and the two handoffs return to the initial full
`(N,y)` state.  These are positive controls for the diagnostic; deleting a
load-bearing target is not a lawful closure strategy.

## 3. Sharp homogeneous transient lengths

THM-2693 proves the target-`a` depth-four zero.  The factor-minimal bit words
make both target scales transparent:

```text
a: D_(13^3) at events 0,1,2 + safe_14 at event 3 -> empty;
b: D_(2*13^5) at events 0,1 + safe_14 at event 5 -> empty. (8)
```

For `a`, the first three target hits shrink the state-three phase to radius
`1/(14*13^2)`, so multiplication by 14 remains in the danger tooth.  For
`b`, two hits shrink the deviation from `0` or `1/2` to radius
`1/(28*13)`; multiplication by 14 gives radius below `1/14`.  Neither proof
uses the guard bit, clock label, carry, or any other ordinary speed.

Exact strict controls show:

```text
target a, either guard sector:  3 states survive, 4 are impossible;
target b, guard-safe sector:    5 states survive, 6 are impossible;
target b, guard-danger sector:  4 states survive, 6 are impossible.
```

The last bracket is deliberately not claimed sharp.  Thus swapping `a` for
`b` postpones homogeneous nilpotence but does not create recurrence.

## 4. A tooth schedule closes intrinsically after adjoining edge debt

The homogeneous obstruction patterns in `(8)` are `aaa` and `bb`.  This
suggests initially retaining the delayed-tooth label and searching the finite
symbolic skeleton that avoids both.  The first successful pattern is
alternating:

```text
delayed teeth: b, a, b, a;
guard sectors: danger, safe, danger, safe.                (9)
```

It has the exact period-four base orbit

```text
7/170 -> 91/170 -> 163/170 -> 79/170 -> 7/170,            (10)
```

because multiplication by 13 sends the numerators

```text
7 -> 91 -> 1183=6*170+163
  -> 2119=12*170+79 -> 1027=6*170+7.                     (11)
```

The four fractions are distinct, so the period is exactly four.

Positivity can be read without the large interval builder.  At the two
`b` states the phase-distance vector for
`(H,q1,...,q5,c1,c2,c3)` is

```text
(7/170, 36/85, 19/170, 6/17, 31/170, 24/85,
 79/170,79/170,6/85).                                    (12)
```

Thus the guard is dangerous, all five unit speeds are safe, `c2` is safe,
and `c3` is dangerous strictly (`6/85<1/14`).  At the two `a` states it is

```text
(79/170,42/85,77/170,7/17,63/170,28/85,
 7/170,7/170,7/85).                                      (13)
```

Thus the guard and units are safe, `c2` is dangerous, and `c3` is safe.
Choose owner clocks `0,1,0,1`.  Their shifted source distances are

```text
79/170, 219/1190, 79/170, 121/1190,
```

all strictly above `1/14`.  This proves every delayed factor in `(9)`--`(10)`
directly; the exact interval construction independently gives normalized
strict margins

```text
(1/883677340, 1/67975180, 1/883677340, 1/67975180).       (14)
```

There is an essential type distinction here.  At the two `b` vertices the
speed-one guard is dangerous, so those points are disjoint from the terminal
`A0/Q_b` word.  At the two `a` vertices the *unshifted* `c1` factor is
dangerous (`7/170<1/14`); only a shifted owner-clock complement is safe, so
those points are disjoint from the terminal `Q_a` word.  Thus `a,b` in `(9)`
cannot be silently reinterpreted as terminal endpoints.

This failure is complete at period four, not a defect of the displayed lift.
Every period-four tail is `j/(13^4-1)=j/28560`.  Exact enumeration imposing
the physical BABA teeth, all five ordinary safeties, and the opposite deep
safety leaves exactly

```text
j=1176,27384,       y=7/170,163/170.
```

These are the two `b`-phase rotations of `(10)`.  (The target teeth alone
also contain `j=0`, which speed 14 kills.)  Adding guard safety at the `b`
vertices makes this two-point set empty.  Adding unshifted-`c1` safety at the
`a` vertices also makes it empty.  Thus no alternative period-four fibre can
turn this delayed schedule into either literal terminal word.

The two failures are one edge variable, not independent vertex defects.  If
`y_{i+1}={13y_i}`, then

```text
||c1*y_i|| = ||13y_i|| = ||y_(i+1)||.                    (edge debt)
```

Consequently unshifted-`c1` danger at an `a` vertex is exactly *narrow*
radius-`1/14` guard danger at the following `b` vertex (and hence lies inside
the broad radius-`1/7` guard-danger sector).  The orbit carries an alternating
edge debt: danger on each `a->b` edge and safe on each `b->a` edge.  This is
the structural coordinate erased when vertexwise tooth labels are treated as
terminal words.

This observation removes the external schedule from the *new* object.  On a
tail/clock state `(y,ell)`, define

```text
A_debt = target-a delayed core + broad guard-safe
         + unshifted-c1 danger + ell-shifted c1 complement safe;

B_debt = target-b delayed core + narrow guard-danger
         + unshifted-c1 safe + ell-shifted c1 complement safe.           (debt sectors)
```

Here the delayed core includes all five ordinary safeties and the opposite
deep safety; narrow guard danger has radius `1/14` and implies membership in
the broad radius-`1/7` danger sector.  These predicates intrinsically label
the four orbit points

```text
B_debt,A_debt,B_debt,A_debt,
```

and the edge-debt word is `(0,1,0,1)`, with `1` exactly on `a->b`.
Consequently the four full states form a genuine recurrent SCC of the newly
defined guard-cospan edge automaton.  This is not a claim that they form an
inherited terminal-word SCC.

The guard refinement also gives a precise next experiment.  Split broad guard
danger into the disjoint narrow tooth and annulus

```text
D_1^(7) = D_1^(14) disjoint-union
          {1/14 <= ||y|| < 1/7}.
```

Because `c1=13` and `c2=13^3`, the following are global identities, not
numerical coincidences:

```text
1_(D_c1)(y) = 1_(D_1^(14))(B y),
1_(D_c2)(y) = 1_(D_1^(14))(B^3 y).                       (debt queue)
```

On the period-four SCC, the `+1` and `+3` debts at each `A_debt` state land
in the same `B_debt` sector class; the exact queue word is
`((0,0),(1,1),(0,0),(1,1))` in BABA order.  This suggests a signed
mapping-cone with a three-step debt queue, or equivalently a small Čech
complex for the narrow/annular guard cover.  Only the identities and this
orbit evaluation are proved here; exactness or cancellation of that proposed
complex is open.

There is an integral multiplicity hidden by the Boolean edge label.  On the
two-class period-four quotient, each `A_debt` vertex has both a `+1` and a
`+3` narrow-guard invoice landing in the same `B_debt` class.  Hence formally

```text
(B_* + B_*^3)[A_debt] = 2[B_debt],
```

not one copy of `[B_debt]`.  The Boolean word records only that this
coefficient is nonzero.  Any lawful integral mapping cone must either retain
or cancel coefficient `2`; reduction modulo two would cancel it for a reason
that is unavailable over the nonnegative integers.  This is a sharper
no-chain-map invoice, not a proved construction of the cone.

Before quotienting the two `A` and two `B` vertices, label the orbit
`B0,A1,B2,A3`.  The maps are

```text
A1 --B--> B2,       A1 --B^3--> B0,
A3 --B--> B0,       A3 --B^3--> B2.
```

Thus the unsigned debt matrix (rows `B0,B2`, columns `A1,A3`) and its
oriented analogue are exactly

```text
B_*+B_*^3 = [[ 1,1],[1, 1]],
B_*-B_*^3 = [[-1,1],[1,-1]].
```

Each `A` emits two units and each `B` receives two.  On the reflection
coinvariant/invariant class the unsigned map is multiplication by `2`, with
`Z/2` cokernel; Boolean OR erases that multiplicity and mod two the quotient
boundary vanishes.  The oriented matrix is anti-invariant and has total zero,
but that signed cancellation does not discharge the original nonnegative
two-unit current.

Now attach the odometer fibre *after* finding the base cycle.  The first
present-bearing lift found in the scout lay in the central rail-free gap and
was rejected by a hostile audit.  Retuning only the fibre coordinate gives
the rail-bearing lift

```text
N=(2572844,2254279,2572531,2254279),
k=(2594970,2227752,2599027,2228065).                      (15)
```

With future digits `(0,6,12,6)`, equation `(1)` holds cyclically at all four
steps.  The lift residues are `(1,7,2,8)` modulo 13, hence every handoff is a
lawful nonzero quotient lift.  The physical states are

```text
(437383487,383227521,437330433,383227509)/820557530.      (16)
```

Their intrinsic clock edges are exactly

```text
0->1, 1->0, 0->1, 1->0,                                  (17)
```

so every owner equals the next shallow clock and the last returns to the
first.  The debt-sector predicate reconstructs the delayed-tooth label from
the state itself, and the full `(N,y,clock,debt-sector)` state therefore
repeats after four handoffs.  Dynamically typing the present factor by the
actual future digit and the shallow clock gives strict margins

```text
(67/2871951355,
 361/2871951355,
 1177/11487805420,
 319/2871951355).                                        (18)
```

This retuned fibre also lies strictly in one THM-2584 route-two rail at every
event, with common source label `1`:

```text
rail keys:
 (6,1,1,0), (6,1,0,12), (6,1,1,0), (6,1,0,12);

weights:
 27581135604,27582102210,27581135604,27582102210;

strict rail margins:
 186812/2871951355, 361/2871951355,
 1123/2871951355, 319/2871951355.                         (19)
```

The rail owner is exactly the target of each intrinsic clock edge in `(17)`.
All four keys nevertheless come from the same fixed `PAT_QB` rail bank.  The
owner clock alternates, but the rail's categorical target does not.

Equations `(14)` and `(18)` also rule out an endpoint-only artifact.  None of
the four tails is speed-14 resonant.

This changes exactly one hypothesis of THM-2693: that theorem repeats the
same target-`a` tooth, so its analytic certificate has three consecutive
`D_(13^3)` occurrences.  At the `b` events in `(9)`, the roles of `c2` and
`c3` are swapped; the `a` occurrences are isolated, as are the `b`
occurrences.  The tail map, ordinary safe factors, guard sectors, owner-clock
cuts, dynamically typed present factors, and integer-lift law are unchanged.
Consequently there is no contradiction with THM-2693: `(9)` is precisely a
different delayed grammar, one of that theorem's stated open exits.  The debt
sectors make its labels intrinsic, but they change the state category.  This
is an exact recurrent SCC of that new category, not an inherited physical or
THM-2305 semantic recurrent SCC.

## 5. The primitive-unit sidecar is positive

At the four exact states `(16)`, retain

```text
carry c=N mod13=(1,1,0,1),
h=(0,6,12,6),
kappa=(1,1,0,0),
floor((2h+kappa)/13)=(0,1,1,0),                           (20)
```

and the closing half/root choices

```text
(right,2),(right,3),(left,2),(left,3).
```

At the dilated intermediate states the roots are `(1,1,12,12)`, while the
lift translations `2k mod13` are `(2,1,4,3)`.  Adding them gives the next
private root cyclically.  Thus the pointwise private-root handoff closes.

The complete target-`b` coefficient sweep changes only the two deep roles in
the inherited THM-2623 universe.  It tests

```text
162 rails x 2 guards x 2 halves x 2 kappas
          x 13 future digits x 7 clocks x 12 roots,
```

with all 26 future half-digit coefficients behind each fine check.  The exact
census is

```text
target-b positives:       2,986,852;
fine checks:                 176,904;
subdigit values:          18,398,016;
target-b global content:         520;
nonzero values not divisible by 26: 0.                  (21)
```

The inherited target-`a` bank has proved content `26`, so the enlarged broad
labelled bank has exact content `gcd(26,520)=26`.  But `B_debt` uses the
*narrow* guard tooth, so broad-danger coefficients are not the final typed
rows.  The exact word split is

```text
broad danger:  86,950 intervals, mass 4,981,755,322,960;
narrow tooth:  44,222 intervals, mass 2,533,638,252,600;
annulus:        42,728 intervals, mass 2,448,117,070,360. (22)
```

The last two supports are disjoint and merge exactly to the first.  On each
selected `b` row, the raw narrow and annular coefficients sum componentwise
to the independently computed broad row.  The complete refined bank tests
`27,597,024` coefficients and has

```text
refined positives:              3,167,868;
refined target-b content:             520;
refined nonzero values not /26:         0;
joint target-a/refined-b content:       26.             (23)
```

Thus `/26` remains the global labelled normalization after the split, rather
than only an integrality accident of the selected rows.  After that division
and reduction modulo 13, the four actual narrow-debt coefficient vectors are

```text
(9,0,0,0,0,12,0),
(0,0,0,2,0,0,0),
(9,0,0,0,2,4,0),
(0,0,0,11,0,0,0).                                      (24)
```

After normalization by private roots `(2,3,2,3)`, their multiplication
determinants in `F_13[z]/(Phi_7)` are

```text
(12,12,8,12),                                            (25)
```

hence all four rows are primitive units.  An independent hostile calculation
using polynomial gcd/resultant against `Phi_7=1+z+...+z^6` gives raw
resultants `(1,12,5,12)` and the same normalized values `(25)`.

There is also a useful conditional holonomy control.  The four normalized
unit polynomials are

```text
11+6z^5,       5z^3,       11+z^4+2z^5,       8z^3.
```

Their ordered product in `F_13[z]/(Phi_7)` is

```text
9+2z+8z^2+7z^3+6z^4+9z^5,                                (26)
```

with determinant `5` and exact multiplicative order `168`.  A separate SymPy
remainder/resultant computation reproduces all three values.  This is **not**
an actual endpoint obstruction: no unit transport between the four rows has
been constructed.  It says only that naive multiplicative gluing, if used as
a chain map, would carry a nontrivial order-168 sidecar that must be cancelled
rather than identified with the identity.

## 6. Exact semantic no-go and comparison with the half-odometer lane

The debt-sector label is now an intrinsic predicate, so the new automaton has
no externally chosen tooth transition.  Its lawful edge retains a positive
fixed-`PAT_QB` THM-2584 rail, delayed-tooth membership, guard side, future
digit, intrinsic clock interface, present membership, a closing private-root
chain, and a nonzero primitive unit.  It still does **not** produce

```text
a THM-2305 terminal word at these vertices,
endpoint-current transport,
or transport identifying the four primitive units.                     (27)
```

The literal `B_0(y)={13y}` repair is in fact ruled out, rather than merely
missing from this witness.  The independently audited exact proof packaged as
[THM-2701 literal singleton-word nilpotence](../01-canon/theorems/THM-2701-literal-singleton-word-one-step-dilation-nilpotence.md)
shows that the literal `Q_a union Q_b` language is positive through five
states and empty at six.  Its mechanism is factor-sparse: `Q_a` at state `j`
forces guard danger at `j+3`, while two initial `Q_b` hits force speed-14
danger at state five.  Thus no longer positive terminal-word schedule under
the same one-step base can turn the debt SCC into a literal recurrent SCC.

This should be compared with the proved
[THM-2698 central half-odometer cycle](../01-canon/theorems/THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary.md):

| lane | exact recurrent object | remaining boundary |
|---|---|---|
| literal `B_0(y)={13y}` | none; sharp nilpotence index six | change chronology or state category |
| `B_0` guard-debt automaton | this period-four primitive-unit SCC | signed guard-debt/Čech transport; order-168 naive unit holonomy |
| `B_1(y)={13y+1/2}` | THM-2698 full-word two-cycle | semantic `C_2` bibundle or direct `D^2` endpoint-current cospan |

THM-2698 is therefore the main positive endpoint frontier; the present `B_0`
cycle is its complementary cospan lane.  Its honest next use is as a signed
mapping-cone/Čech guard-debt control, where the alternating edge debt may
cancel, not as another attempt to attach positive literal endpoints
vertexwise.  The 165 scalar row obligations remain downstream in every lane.

## 7. Reproduction and scope

```bash
python3 04-computation/lrc14_delayed_tail_scc_factor_diagnostic_20260728.py
python3 -O 04-computation/lrc14_delayed_tail_scc_factor_diagnostic_20260728.py
python3 04-computation/lrc14_target_b_successor_content_mixed_unit_20260728.py
python3 -O 04-computation/lrc14_target_b_successor_content_mixed_unit_20260728.py
```

Both primary runs byte-match
`05-knowledge/results/lrc14_delayed_tail_scc_factor_diagnostic_20260728.out`;
the final normal run took `16.08s` and the final two-run byte check took
`30.38s` total.  The full refined companion run tested `27,597,024`
coefficients and produced its frozen output in `567.01s`.  The ensuing script
edit only froze the returned per-shard and aggregate tuples as mandatory
`require` gates.  An independent final normal/optimized companion replay is
routed separately because running two copies concurrently materially slows
the four-worker census.

The LF-byte SHA256 values are

```text
diagnostic script: c610ae5d0edef05813df2adf96a16ffae0e29bb2ca23de462715723ef2c3e9a6
diagnostic output: 9188b84440607b38b974666c60dc567bd233817b54d7091d1d4a3df0635b1547
coefficient script: fdd0c70865056d86d3e2dcd853ca12d4249b54749f7496743dc321b986f95d96
coefficient output: 3b5be77278ecc5af1e4759de7a55c6def351215ada56dbf8c5acf480cf9cfbc1
```

The universe is exactly `(2)` on the canonical typed row, plus the displayed
target-deleted controls, the complete period-four delayed-core lattice, and
the displayed intrinsic debt SCC.  The broad and narrow/annular coefficient
banks are exhaustive in the stated THM-2623 universe.  No exhaustive search
over longer mixed debt words is claimed.  The four narrow rows are primitive
units after exact global `/26` normalization, but no transport between those
units is supplied.  No terminal-word SCC, endpoint transport, row exclusion,
or LRC(14) conclusion follows.
