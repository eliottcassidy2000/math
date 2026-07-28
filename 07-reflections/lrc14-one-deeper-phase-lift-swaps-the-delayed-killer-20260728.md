# One deeper phase lift repairs the minimal delayed carrier and exposes the next killer

**Status: PROVED ELEMENTARY SCALE LEMMA + VERIFIED-EXACT SIDECAR.  Not a row
exclusion and not an LRC(14) conclusion.**

THM-2693 proves that the inherited delayed word is positive for three event
states and empty at four.  Its smallest certificate retains the target tooth
at speed `13^3` and speed-`14` safety.  The closest physical escape is to
prolong the THM-2657 odometer by one base-thirteen digit.  That prolongation
really does repair the two-factor carrier, but it cannot repair the full
word: the unchanged high-speed safety becomes the next depth-four killer.

This note separates those two statements, gives the sharp general interval
lemma behind them, and records slope `14` as the first wrap-side control.

## 1. The sharp nested-danger lemma

Let `b>=2`, let `0<rho<1/2`, and use the centered half-open arc

```text
I_rho=[-rho,rho) in R/Z,             E_b(z)={bz}.
```

If

```text
rho <= 1/(b+1),
```

then for every `n>=0` one has the exact half-open identity

```text
C_n := intersection_(j=0)^n E_b^(-j) I_rho
     =[-rho/b^n,rho/b^n).                              (1)
```

Indeed, take `z in I_rho` with `E_bz in I_rho` and write
`bz-k in I_rho`.  If `k>=1`, then

```text
z >= (1-rho)/b >= rho,
```

contrary to the excluded right endpoint.  If `k<=-1`, the strict upper
inequality in `bz-k<rho` gives `z<-rho`.  Hence `k=0`, and iteration proves
(1).  At the equality `rho=1/(b+1)`, the two apparent noncentral branches
land exactly on the excluded right endpoint, so (1) is literal rather than
an a.e. statement.

The radius threshold is sharp.  If `rho>1/(b+1)`, the positive interval

```text
[(1-rho)/b, min(rho,(1+rho)/b))                         (2)
```

lies in `I_rho` and returns to it through wrap branch `k=1`.  Thus the
single central tooth in (1) acquires a genuine side tooth immediately above
the threshold.

Now let `A,q` be positive integers, let `h,n>=0`, and suppose

```text
A divides q b^h,                 c=q b^h/A.               (3)
```

Let `D_(v,theta)` denote the event

```text
v y-theta in I_rho modulo one.
```

Multiplication by `A` sends the sparse factor intersection

```text
intersection_(j=0)^n D_(A b^j,0)
             intersection E_b^(-h) D_(q,theta)^c          (4)
```

onto

```text
C_n intersection {z : cz-theta notin I_rho}.             (5)
```

Put `delta=||theta||` and `a=c rho/b^n`.  Since multiplication by `A`
preserves Haar measure, (4) is empty exactly when the centered image arc of
radius `a` is contained in the danger arc of radius `rho`, namely

```text
delta+a <= rho.                                           (6)
```

Equality in (6) is included for positive `c`: the two half-open orientations
match.  If (6) fails, a positive endpoint sliver survives.  For zero phase,
the exact coefficient threshold and the least number of target contractions
are therefore

```text
c <= b^n,                  n_min=ceil(log_b c).           (7)
```

The least temporal span using a safe factor at state `h` is

```text
1+max(h,n_min).                                           (8)
```

If the proof is required to use only target factors strictly before the
safe event, then `n_min<h` is additionally required.  If (3) fails, the
safe value retains a nontrivial inverse-fibre phase; a branch label is then
load-bearing and the central scalar implication is unavailable.  Negative
slopes or speeds have the same positive-measure threshold after reflection,
but half-open equality can leave one null endpoint, so (1)--(8) deliberately
state the positive-speed form used by LRC.

## 2. The LRC parameters occupy the equality wall

For the THM-2693 carrier,

```text
b=13,       rho=1/14=1/(b+1),
A=13^3,     q=14,       h=3.
```

Here `c=14` and `n_min=2<h`.  Thus three target factors give

```text
z,13z,13^2z in I       =>       z in I/13^2,
```

and the state-three speed-`14` safety is impossible.  This is the sharp
factor-sparse depth-four proof: it uses no state-three target factor.

If the speed-`14` danger tooth is shifted by `theta`, its image radius on
the central cylinder is

```text
a=1/169=14/2366,
rho-a=155/2366.                                          (9)
```

Every nonzero thirteenth phase has distance at least

```text
1/13=182/2366 > 155/2366,                                (10)
```

so every such phase leaves positive support.  For `theta=+/-1/13`, the
surviving mass in the minimal cylinder is exactly `27/33124`; for the other
ten nonzero phases it is the full cylinder mass `1/1183`.  Zero phase is
empty.

## 3. A genuine one-level-deeper lift realizes the phase

Put `R=13^6` and define

```text
T_s^(7)(x)={13x+s/13^7},            s in F_13.
```

For the delayed base `y={Rx}`, this gives the exact affine map

```text
F_s(y)={13y+s/13}.                                      (11)
```

This is not a cosmetic replacement of one factor.  It is a physical circle
map obtained by adjoining one finer base-thirteen digit.  For each `s!=0`,

```text
y_s=s/13,                 F_s(y_s)=y_s,                  (12)
```

and

```text
13^3 y_s=0 mod 1,         14y_s=s/13 mod 1.              (13)
```

The first value is strictly target-dangerous and the second is strictly
speed-`14` safe, with minimum margin `1/182`.  Thus the minimal two-factor
carrier has twelve genuine recurrent interior fixed points.  The exact
replay also solves the fibre equation and produces a fixed point of the
original circle map for every `s`; at `s=1` it is

```text
x=5229043/62748517,              {Rx}=1/13.              (14)
```

All `13^2` two-edge phase words have positive three-state support.  Their
exact extrema are

```text
minimum: phase word (3,0),  6732 components, mass 18513/439239619;
maximum: phase word (1,8),  8294 components, mass 3509/67575326. (15)
```

At this minimal level, the one-deeper lift has supplied exactly the recurrent
base object which the old integer odometer could not supply.

## 4. The full delayed word still dies uniformly at four

The full raw guard-sector union also requires safety at

```text
c3=2*13^5=338*13^3.                                    (16)
```

For an arbitrary time-dependent phase word `(s_0,s_1,...)`, put

```text
z_j={13^3 y_j},             y_(j+1)=F_(s_j)(y_j).
```

Since `13^3(s_j/13)=169s_j` is integral,

```text
z_(j+1)={13z_j}                                             (17)
```

for every phase label.  Four target teeth therefore force

```text
z_0 in I/13^3.                                             (18)
```

But the state-zero high-speed phase is exactly

```text
c3 y_0=338z_0 mod 1,
```

and

```text
338<13^3=2197.                                            (19)
```

Equations (18)--(19) put the high-speed phase back in `I`, contradicting
its required safety.  This proof is uniform in all `13^3` three-edge phase
words.  Exact controls `000`, `123`, and `12,7,4` are empty; the analytic
argument proves the remaining `2194` words.  Hence the full phase language
has counts

```text
depth 1: 1/1,     depth 2: 13/13,     depth 3: 169/169,
depth 4: 0/2197.                                            (20)
```

There are two useful robustness checks.

First, shifting only speed `14` is weaker still: unchanged speed `27` also
satisfies `27<13^2`, so the original three-target prefix repeats with `27`
in place of `14`.  Second, any common base translation which preserves the
`13^3` target phase automatically preserves the `c3` phase, because (16)
makes `c3` an integral multiple of the target speed.  A common translation
therefore cannot create the needed target-versus-high-speed relative phase.

The one-deeper lift has not failed; it has exposed the next exact factor.
It swaps the speed-`14` killer for the `c3` killer.

## 5. Physical typing and what the map forgets

The source of (11) is the `C_(13^7)` affine stalk; its target is the delayed
base circle.  The quotient map is

```text
x -> y={13^6x}.                                           (21)
```

It preserves chronological circle dynamics and the target/high-speed phase
relation, and it shifts every ordinary speed `13k+1` by the same `s/13`.
It does **not** preserve the old THM-2657 integer-fibre recurrence: the old
lift had translation `m/R`, whereas (11) needs `s/(13R)`.  Projecting to `y`
forgets the new digit `s`; retaining it has not yet been identified with the
existing predecessor carry, private root, clock action, present packet, or a
semantic endpoint.

THM-2640's physical deck/gauge translation `y->y+s/13` is the closest
existing coefficient-side operation.  It shifts the ordinary unit/guard
family and fixes all speeds divisible by thirteen, but a fixed rail does not
descend under it; the whole carrier and its gauge label must be translated
equivariantly.  THM-2365's lawful present target action is not a substitute:
its delayed word is exactly neutral.  THM-2657's integer odometer is also not
a substitute, because it is invisible on `y`.

Thus the cheapest physical next test is a labelled `C_(13^7)` prolongation,
but the full-word escape requires one additional operation: a factor-wise
relative phase between the `13^3` danger and `c3` safety, deletion/change of
one of those factors, or a different base slope with genuine wrap branches.

## 6. Slope fourteen is the first wrap-side control

At slope `13`, `rho=1/14` is exactly the last no-wrap radius.  At slope `14`,

```text
1/14 > 1/(14+1),
```

so noncentral branches appear.  The target-only system already has the exact
two-cycle

```text
1/15 -> -1/15 -> 1/15,                                  (22)
```

and both points lie strictly in `I`.  Moreover `338/15` is high-speed safe.
This is a genuine wrap-side positive control, not an analogy.

The complete raw word nevertheless has no recurrent component under the
unshifted slope-`14` map.  Exact interval recursion remains positive through
fourteen states, with `20` components and mass

```text
109/4125810350281912434688,
```

then becomes empty at state fifteen.  Thus changing the slope moves the
nilpotence horizon from four to fifteen but does not by itself close the
carrier.

## 7. Replay and scope

Run

```bash
python3 04-computation/lrc14_c13_7_phase_lift_slope14_boundary_20260728.py
python3 -O 04-computation/lrc14_c13_7_phase_lift_slope14_boundary_20260728.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_c13_7_phase_lift_slope14_boundary_20260728.out.
```

The byte-matched script and output SHA256 values are respectively

```text
f7f682820e564f78b9b9600d3e7e083f49576e75f50ab5f1dc58691d431c5dda
e37d8c44d697dd41d632289e958805a0d85039b4ed455daed3767ad678efd221
```

The replay builds the canonical raw word from the THM-2693 modules,
enumerates the exact `13^2` depth-three phase atlas on one refined integer
grid, checks three depth-four hostile words independently, verifies all
twelve fixed lifts, and runs the complete slope-`14` interval recursion to
its first zero.

No inherited `C_(13^7)` carry/root law, present or clock gluing, endpoint
transport, scalar row closure, or LRC(14) conclusion is asserted.
