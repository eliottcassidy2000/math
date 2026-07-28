# The missing delayed coordinate is a reflection-symmetric `C_221` stalk

**Status: PROVED-ELEMENTARY projection, residue, and stalk laws +
FINITE-EXACT locally packet-typed design target with independent normal/optimized
replay.  This is not promoted canon, not a semantic endpoint transition, not
a row exclusion, and not an LRC(14) conclusion.**

The first foreign-base recurrent point in the inherited delayed word was the
old hostile `4/17`.  Auditing the physical handoff operations shows both why
the existing odometer could never see it and what the smallest missing
coordinate is.  The right state is not a bare `C_17` colour and not another
carry digit.  It is

```text
nu = N mod 221,                  221=13*17,
```

on the enlarged stalk `x=N/(17*13^6)`.  This one residue simultaneously
remembers the delayed seventeenth, predecessor carry, and exact high-probe
microphase.  On it, slope seven and quarter torsion force a unique reflected
two-cycle.  Suitable `C_(13^5)` lifts of that cycle meet every currently
implemented local packet factor for arbitrarily long prescribed finite words.

## 1. Operation typing: seven was never a delayed circle slope

Put `R=13^6` and `pi(x)={Rx}`.  The operations proved or tested in the
THM-2640--2694 handoff lane have three distinct types.

| operation | actual circle map | delayed map under `pi` | role |
|---|---|---|---|
| formal clutch | `c -> c+7 delta`, `r -> r+delta` | none | finite-label equation, not a circle map |
| THM-2657/2672 deck lift `G_k` | `x -> {x+k/R}` | `y -> y` | carry/root/configuration sidecar |
| THM-2680 dilation | `x -> {13x}` | `y -> {13y}` | chronological event map |
| THM-2684 signed dilation | `x -> {-13x}` | `y -> {-13y}` | reflected chronological event map |
| physical THM-2689 specialization | `x -> {eps*13x+k/R}` | `y -> {eps*13y}` | chronological event plus deck sidecar |

Indeed `R(k/R)=k` is integral.  Therefore for any word containing `d`
chronological letters and any number of deck letters,

```text
pi(word(x)) = {epsilon*13^d*pi(x)}.                       (1)
```

The literal statement “every affine handoff has delayed linear part
`+/-13`” is type-sensitive.  It is true for one chronological event; it is
false if configuration/deck correspondences are also called handoffs, since
they have delayed linear part `1`.  Reflection alone has part `-1`, but is a
symmetry rather than a new event.  The corrected statement is (1).

Consequently the proved mixed scout
[THM-2694](../01-canon/theorems/THM-2694-mixed-dilation-slope-seven-present-unit-long-word-and-first-gap.md)
has no hidden delayed slope: one `D` edge followed by any slope-seven orbit
still projects to multiplier `13`.  Its slope-seven leg is a genuine physical
configuration translation but only a sidecar on the delayed base.  The new
typing statement is compatible with THM-2694's proved 95-edge finite word and
its explicit pumping no-go; it does not upgrade that word to an infinite or
configuration-switching chronology.

## 2. The complete denominator-seventeen chamber

Let `W` be the raw THM-2693 word

```text
D_(13^3) intersect D_14^c intersect D_27^c intersect D_40^c
             intersect D_53^c intersect D_66^c intersect D_(2*13^5)^c.
```

At denominator `17`, the danger arc `[-1/14,1/14)` contains exactly residues
`0,+/-1`.  Since

```text
13^3 = 4                                      (mod 17),
```

the target tooth leaves candidates `0,+/-13`, namely `0,4,13`.  Zero fails
every safety factor.  At `4/17`, the centered target/guard residues are

```text
(-1, 5, 6, 7, 8, -8, 2)/17,                            (2)
```

and reflection gives the `13/17` row.  Hence exactly

```text
W intersect (1/17)Z/Z = {4/17,13/17}.                    (3)
```

Modulo `17`, `13` has order four:

```text
<13> = {1,13,16,4}.
```

The chamber `(3)` is the odd-power part `{13,13^3}`.  Either chronological
multiplier `+13` or `-13=13^3` sends it to the even-power part `{1,16}`,
where the **target tooth itself** fails.  Thus there is no one-event internal
edge.  This remains true although

```text
(+13)(-13)=1,               13^2=-1,               13^4=1 mod17.   (4)
```

The two- and four-step products return only after an illegal intermediate
event.  Deck sidecars cannot repair that event because they are delayed
identity maps.

For a denominator-seventeen affine phase

```text
F_(eps,t)(y)={eps*13y+t/17},
```

an exact chamber edge exists precisely for

```text
t in {3,5,12,14}.                                        (5)
```

The eight signed arrows are

```text
eps=+1:  t=3:4->4,  t=5:13->4,  t=12:4->13, t=14:13->13;
eps=-1:  t=3:13->4, t=5:4->4,  t=12:13->13,t=14:4->13.
```

In particular `y->{13y+3/17}` fixes `4/17`.  Its strict `W` component is

```text
(2446165/10396204, 2446177/10396204),
```

with left/right margins

```text
L_-=11/176735468,          L_+=193/176735468.             (6)
```

For every `H>=1`, the open interval

```text
(4/17-L_-/13^(H-1), 4/17+L_+/13^(H-1))                  (7)
```

stays in `W` for `H` states.  This raw fixed-point control is not globally
carry typed: its circle translation is `3/(17R)`, so

```text
R beta=3/17,                  c3 beta=6/221,              (8)
```

and a pure translation has predecessor jump `0` below one wall and `1` above
it.  Thus `(8)` lies outside THM-2657's integer-`R beta` classification and
outside the old `F_13` root-phase grid.

Among odd `q<100` coprime to `13`, the quarter-torsion phase
`(-1/4 mod q)/q` lies in `W` first at `q=17`; the exact accepted list is

```text
17,43,47,51,63,71,95.                                    (9)
```

So `17` is not only the first foreign fixed-point denominator: it is the
minimal transverse modulus accepted by the symmetric quarter-torsion
mechanism below.

## 3. Enlarging the stalk and finding the right state

Put

```text
L=17R,                         x=N/L,
T_m(x)={13x+m/L},              m=17k+a, 0<=a<17.
```

Write `Rx=n+y`.  Exact division gives the skew product

```text
y'={13y+a/17},
n'=13n+k+floor(13y+a/17).                                (10)
```

Thus the new `C_17` colour supplies the delayed phase, while the old fibre
lift `k mod13` steers the next carry.  A `C_17` colour alone does not transport
the packet globally; the pair is load-bearing.

There is a cleaner compression.  Set

```text
nu=N mod221.
```

Then

```text
a=nu mod17,
c=((nu-a)/17) mod13,
high-probe phase = 2nu/221.                              (11)
```

The affine map is simply

```text
nu -> 13nu+m                              (mod221).       (12)
```

Modulo `17`, multiplier `13` is invertible; modulo `13`, it resets the old
coordinate and `m` supplies the next one.  This is exactly the coordinate
combination separately lost by the delayed, carry, and coarse-root
projections.

The stalk fits the exact cyclic extension

```text
0 -> C_(13^5) -> C_(17*13^6) -> C_221 -> 0.              (13)
```

CRT separates `C_17` from `C_(13^6)`, but `(13)` does **not** split as
`C_221 x C_(13^5)`: the unique order-`221` subgroup maps to a subgroup of
order only `17`, because the quotient and kernel retain the shared factor
`13`.  This is the same nonsplit guardrail as the THM-2657 odometer, now one
transverse colour higher.

## 4. Slope seven plus quarter torsion uniquely forces the cycle

Demand a reflection-symmetric stalk cycle with opposite translations and
literal slope-seven motion in the refined high coordinate `u=2nu`:

```text
13nu+m=-nu,                  2m=7                 (mod221). (14)
```

Both `2` and `14` are units modulo `221`, so `(14)` has one solution:

```text
m=7/2=114,                  nu=-m/14=-1/4=55,
-m=107,                     -nu=166.                      (15)
```

Projection modulo `17` gives exactly

```text
nu: 4 <-> 13,               m:12 <-> 5,                  (16)
```

while `(11)` gives carries `3 <-> 9`.  Both fibre lifts have `k=6 mod13`.
The two-step `C_17` maps are

```text
r -> -r+8, whose unique fixed point is 4;
r -> -r+9, whose unique fixed point is 13.                (17)
```

On the refined high coordinate, the states and unshifted `D` images are

```text
u=(110,111)=+/-110,
13u=(104,117),
2m=(7,214)=+/-7.                                         (18)
```

Hence the root microphase transports literally:

```text
104+7=111,                  117-7=110             (mod221). (19)
```

This is a `C_17`-refined slope-seven toothpick, not the old `+/-7/R` deck
translation.  The congruence is structural; identifying these transverse
translations with a lawful continuation or switch of THM-2694's proved
finite operation remains open.

## 5. Exact physical lift and complete local packet audit

The `C_(13^5)` kernel in `(13)` can now be used for a physical landing.  The
exact selected numerators and lifts are

```text
N0=39123022,                 N1=41305508,
m0=25040740,                 m1=76541689,

N1=13N0+m0,                 N0=13N1+m1             (mod L),
m0=114 mod221,              m1=107 mod221.               (20)
```

Thus `x_i=N_i/L` is an exact two-cycle.  Full reconstruction of the canonical
THM-2584/2640 carrier gives

| state | clock edge | `(deep,c,h,kappa)` | right root | `D`-right root | source-1 right unit |
|---|---|---|---|---|---|
| `x0` | `1->4` | `(12,3,3,0)` | `6` | `6` | vector `(0,8,9,12,0,0,0)`, determinant `4` |
| `x1` | `4->1` | `(0,9,9,1)` | `6` | `6` | vector `(3,0,0,0,12,2,0)`, determinant `3` |

Both points lie strictly in sector zero, their selected present packets, a
source-one rail, both private half-edges, and the displayed intermediate
`D`-right half.  Each current owner is the following shallow clock, both
stored edges are nonconstant, and the right root glues literally on both
edges.  The common primitive-unit source set is

```text
{1,3,4,5,8,9,10,12}.                                    (21)
```

Packet support by itself admits all four affine interpolants between the two
states; their phase matrix is

```text
((3,12),(5,14)).                                         (22)
```

Clock chronology deletes both loops: state zero has owner `4`, so it must go
to state one with shallow `4`; state one has owner `1`, so it must return to
state zero.  The alternating cross-cycle is therefore selected by an intrinsic
label, not by choosing a convenient phase word after the fact.

The exact targeted fibre audit makes the selection boundary reproducible.
Its universe fixes `y=4/17,c=3` or `y=13/17,c=9`—the carries forced by literal
right-root gluing—and requires a nonconstant intrinsic edge, sector zero,
the corresponding present packet, some actual rail, the physical right root,
and at least one primitive right-edge unit.  It finds

```text
pre-unit states:              6260, 6252;
primitive-unit states:        6260, 5966;
initial states with a reverse-clock/common-source mate: 3441;
ordered admissible pairs:     902300.                    (23)
```

The pair `(20)` is the first exact witness in that declared ordering.  Thus
literal root gluing is not a unique numerical accident; it occupies a large,
explicit finite fibre family.

## 6. A positive all-retained-local-factor cylinder for every finite horizon

At each state, intersect the strict local components for

```text
rail, present packet, predecessor carry, half-digit, delayed sector word,
current right root, intermediate D-right root, shallow clock, owner clock.
```

The minimum two-sided radius over both states is exactly

```text
eta=11/853068347561612,                                  (24)
```

and the only limiting factors are the two delayed-sector components.  If the
initial point is `x0+d`, alternating `(m0,m1)` gives exactly

```text
x_j = x_(j mod2)+13^j d.                                 (25)
```

Therefore, for every prescribed `H>=1`,

```text
|d| < eta/13^(H-1)                                      (26)
```

keeps all `H` event states in the same locally packet-typed components.  The
starting cylinder has exact width

```text
(11/426534173780806)/13^(H-1)>0.                         (27)
```

This scout exhibits a delayed-handoff object that retains
rail, present, clock gluing, raw sector, predecessor carry, half-digit,
current and intermediate private root, and primitive unit for every finite
horizon.  It evades THM-2693 because `(10)` changes the delayed phase
transversely; it does not contradict the theorem's uniform zero for integer
`R beta` lifts.

## 7. Connections and exact missing consequence

- [THM-789](../01-canon/theorems/THM-789-two-sheet-erosion-and-symmetric-return-packet.md)
  already used `4/17` as its heavy-phase trapped return.  Here the same strict
  chamber is transported rather than merely revisited.
- [THM-2558, section 5](../01-canon/theorems/THM-2558-sparse-owner-fibre-all-slope-target-role-forcing-boundary.md)
  has the physical six-comb fibre `z=1/17`, missed root `3`, and
  `x_3=4/17`, where `k_a=17` is the unique clean first failure and all twelve
  selected heads miss it.  This confirms that `C_17` is also a sharp
  owner/blind-root modulus.  The physical `x` points and operations are not
  identified here.
- [THM-2114, uniform quarter torsion](../01-canon/theorems/THM-2114-finite-ring-needles-and-rank-eight-tree-equality-rigidity.md)
  uses a uniform order-four transverse point for a toothpick escape.  Equation
  `(15)` is the same quarter-torsion operation on the `C_13 x C_17` stalk, not
  a transfer of THM-2114's cover conclusion.
- [THM-2693](../01-canon/theorems/THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence.md)
  proved that deeper common `13`-adic phases only swap the delayed killer.
  The present mechanism explains the missing coordinate: the target/high
  relative microphase lives on `C_221`, not on another common `C_(13^j)` digit.

What is still missing is semantic, not local-support arithmetic.  An arbitrary
pair of rational states can always be connected by the interpolant

```text
m=N_1-13N_0 mod L.                                       (28)
```

Thus the existence of `(20)` does not prove that LRC chronology supplies
these transverse lifts.  No current theorem identifies `m0,m1` with a
target-active endpoint/owner operation, scalar return, or one of the `165`
row obligations.  The highest-leverage next computation is to insert the
semantic endpoint current into the exact cylinder `(26)`, or prove that every
lawful endpoint map is confined to the old integer-`R beta` subgroup.

## Reproduction

```bash
python3 04-computation/lrc14_mod17_typed_handoff_semigroup_scout_20260728.py
python3 -O 04-computation/lrc14_mod17_typed_handoff_semigroup_scout_20260728.py
python3 04-computation/lrc14_mod17_transverse_phase_typed_cycle_probe_20260728.py
python3 -O 04-computation/lrc14_mod17_transverse_phase_typed_cycle_probe_20260728.py
```

Normal and optimized transcripts are byte-identical.  SHA-256:

```text
projection script  5b43e3017b1dc0e5e82ac0e8d006f32537e70d42dab31ef6126c787f5e864bfb
projection output  6e84dc49f29b088652b55489234dde609bbead8a655b51a15f9165c5530b7a44
typed script       33530aa52760136e33de7ab200913652de2892ec9593377fd17e081209b99bc4
typed output       3a507ef2eb47e21091950ca10afddf9520684070442f2a5f25792313992bbd19
```
