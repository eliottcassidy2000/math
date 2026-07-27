---
id: THM-2464
title: "Linked-blocker depth threshold and clock-cell universality"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Literal
  lambda=1 is a root-multiplicity lock: clean
  two-root ordinary and four-root guard charts force the linked
  blocker safe. Lambda=2 is universal after phase refinement. For a
  fixed tuple, delayed linked factors are realizable iff their joint
  clock cell has interior. The THM-2462 tuple passes an exact
  same-depth lambda=10 control and a later positive-word filter.
  Semantic owner, scalar-cover, and endpoint-current debts remain;
  no scalar row is excluded.
source: codex-2026-07-26-linked-blocker-depth-threshold
depends_on:
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2458-clean-root-guard-danger-thirteen-chart-uniform-offset-hostile
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
  - THM-2462-mixed-radix-root-phase-orbit-universality
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2397-clean-root-same-parent-charged-role-partition
script: 04-computation/lrc14_linked_blocker_clock_cell_thm2464.py
output: 05-knowledge/results/lrc14_linked_blocker_clock_cell_thm2464.out
script_sha256: e4c66419cd49fa77187cd4f2f7cbd82ce4737d007dec6ee9600d5650414f5a3e
output_sha256: 8612b1ddd8f6be0d343efb2e07545a017c8aa52171811e88ebee4a36a22174b5
hash_basis: working-tree bytes (LF)
---

# THM-2464 -- linked-blocker depth threshold and clock-cell universality

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2462 showed that arbitrary finite families of strict open
root-phase charts can be realized by one fixed mixed-radix speed
tuple. Its labelled-bit corollary deliberately used independent
physical factors `c=13v`. An actual blocker is more rigid:

```text
c=13^lambda u
```

shares the unit `u` with its root-phase coordinate. This theorem
classifies that linkage.

There are two sharply different operations:

```text
refine phase cells, then rebuild the fixed speeds;

keep one fixed speed tuple, then shrink its realized parent cells.  (1)
```

The first has an exact depth threshold. The second has an exact joint
clock-cell criterion. Their common obstruction is not vague phase
compatibility: it is either the literal shallow multiplicity bit or
an empty finite cylinder.

## 1. The exact linked-blocker identity

Use the common oriented root section

```text
x_r=(y+r)/13,                    r in C_13.                       (2)
```

For an integer `13`-unit `u`, put

```text
c=13^lambda u,                  lambda>=1.
```

Because the root term is integral,

```text
c x_r
 =13^(lambda-1)u(y+r)
 =13^(lambda-1)uy+integer.
```

Therefore the blocker bit is independent of `r` and obeys the exact
identity

```text
1_(D_c)(x_r)
 =1_(D_u)(13^(lambda-1)y),                               (3)

D_a={x in T:||ax||<1/14}.
```

This is the linkage that an independent factor `c=13v` forgets.

## 2. Literal depth one is the root-multiplicity bit

Write

```text
N=round(uy),

delta=uy-N in [-1/2,1/2).                              (4)
```

Let `rho=u mod 13`. After the affine root relabelling which makes
`N+rho r` divisible by thirteen, the ordinary danger test is

```text
|delta+j|<13/14,             j in Z.                   (5)
```

Only `j=0,1,-1` can occur. Away from the two ordinary boundary
phases,

```text
# {r:x_r in D_u}
 =1,                       |delta|<1/14;

 =2,                       1/14<|delta|<1/2.            (6)
```

At the open-boundary phases

```text
|delta|=1/14,
```

the root count in (6) is still one, but the shallow blocker is safe.
Thus the precise statement is

```text
lambda=1 blocker dangerous
  -> ordinary root count one;

away from |delta|=1/14:
lambda=1 blocker dangerous
  <-> ordinary root count one.                         (7)
```

For a guard danger of radius `1/7`, the scaled inequality is

```text
|delta+j|<13/7.
```

It gives

```text
# guard roots
 =3,                       |delta|<=1/7;

 =4,                       |delta|>1/7,                 (8)
```

with the equality interpreted using the strict open danger.
Consequently

```text
four-root guard -> lambda=1 linked blocker safe.        (9)
```

The converse of (9) is false: a three-root guard with
`1/14<=|delta|<=1/7` also has a safe shallow blocker.

Equations (7)--(9) give the first exact compatibility boundary:

```text
clean two-root ordinary mask
or clean four-root guard
  -> its linked lambda=1 blocker is safe.              (10)
```

In particular, a shallow linked blocker cannot be assigned a
dangerous truth value on any of the six clean THM-2458 root roles.
This conclusion is conditional on the blocker quotient `u` actually
being the unit of that displayed root-moving role; it does not
identify an arbitrary blocker quotient with a guard/unit label.

### 2.1 Conditional complete-bank improvement

Suppose, in addition to the current canon, that:

1. the normalized `lambda=1` blocker is one of THM-2452's two
   nondeep Boolean coordinates; and
2. on the same canonical root section its quotient `u` is identified
   with one of the clean THM-2458 root-moving roles.

Then (10) fixes that blocker coordinate safe and halves the complete
bank:

```text
128 -> 64.                                               (10a)
```

Under precisely these extra hypotheses, THM-2459's nonloop
four-atom coarsening constants improve to

```text
service edge >=M_0/64^2=M_0/4096;

drift union >=D_0/123^2=D_0/15129;

root coefficient >=M_0/(4096*2028)
                 =M_0/8306688,                         (10b)
```

with the selected-loop drift branch at least `D_0/4096`.
THM-2460's three-word bank likewise drops from `384` to `192`.

This is a conditional invoice, not a current canonical improvement.
THM-2452 explicitly leaves the complete mask's identification with a
canonical `C_13` root section open, and no proved theorem identifies
the shallow blocker quotient with a displayed THM-2458 role. Thus
(10a)--(10b) do not change the live ledger.

## 3. Depth two is already universal before fixing the speeds

The obstruction in (10) disappears immediately at depth two.

For an ordinary two-root mask, either oriented branch of the
fractional phase coordinate

```text
t=uy mod 1
```

has length

```text
3/7.                                                    (11)
```

The chart separately retains the nearest-integer residue modulo
thirteen which fixes the root translation; refining `t` does not
change it.

For a four-root guard, either phase branch has length

```text
5/14.                                                   (12)
```

Choose strict target components for the linked blocker:

```text
B_1=(-1/14,1/14),                 dangerous;

B_0=(1/7,3/14),                   safe.                (13)
```

Their lengths are `1/7` and `1/14`. At `lambda=2`, the
blocker condition inside the phase coordinate is

```text
13t mod 1 in B_epsilon.                                  (14)
```

The starts of the inverse branches in (14) have mesh `1/13`.
An interval `J` of length `L` contains a complete inverse branch of
an interval `B` of length `ell` whenever

```text
13L>1+ell.                                               (15)
```

For (11)--(13), the worst exact margins are

```text
13(3/7)-(1+1/7)=31/7;

13(5/14)-(1+1/7)=7/2.                                    (16)
```

They are strictly positive. Hence every clean phase branch contains
a rational open subinterval on which either requested linked blocker
truth value holds. The retained subinterval has length `1/91` in
the danger case and `1/182` in the chosen safe case.

This proves:

> **Phase-before-speed linked-bit theorem.** Fix finitely many clean
> root charts and one independently supplied linked blocker factor
> for each chosen root-phase coordinate. For every chart prescribe
> the blocker safe or dangerous. If every blocker depth satisfies
> `lambda>=2`, refine the corresponding phase intervals as above and
> apply the mixed-radix finite-word lemma of THM-2462. There are
> pairwise distinct positive fixed `13`-unit speeds `u_i` which
> realize all old root masks and all prescribed physical blockers
>
> ```text
> c_i=13^(lambda_i)u_i                                  (17)
> ```
>
> on disjoint positive rational parent intervals. The intervals may
> be shrunk to equal weight.

The same proof works for every `lambda>=2`, because the inverse mesh
only becomes finer. Equation (10) shows that `lambda=2` is the exact
universal threshold for one linked bit on a clean root role.

This operation rebuilds the integer speeds after refining their
phase cells. It is stronger than what can be claimed for a previously
fixed tuple. It also treats each appended bit as a distinct physical
factor. Several transported requirements on one physical factor must
be tested jointly; Section 4 gives the exact rule.

## 4. Fixed-tuple relative-clock packets

Now keep the speeds fixed. Let

```text
u_1,...,u_s
```

be positive `13`-units, with physical collisions separately excluded,
and let

```text
Y_w subset T,                 w in W,                  (18)
```

be finitely many nonempty rational open intervals on which the old
strict root word is constant.

A **relative-clock packet** consists of finitely many tests at
relative depths `d>=0`. A test may be:

```text
1_(D_u)(13^d z)=epsilon,                              (19)
```

or a prescribed constant cell of a fixed rational finite-step
Boolean word,

```text
Q(13^d z)=epsilon.                                    (20)
```

For chart `w`, let `C_w` be the intersection of every requested set
in (19)--(20), with all requirements using the same `z`. Thus
`C_w` is the **joint** clock cell, not a product of marginal
nonemptiness statements.

> **Fixed-tuple clock-cell theorem.** The following are equivalent:
>
> 1. after arbitrarily large common base delay `n`, every `Y_w`
>    contains a positive rational open subinterval realizing its
>    entire relative-clock packet;
> 2. every joint cell `C_w` has nonempty interior.               (21)

### Proof

For sufficiency choose a rational open component

```text
J_w subset int(C_w)
```

of length `ell_w`. If the current interval has length `L_w`, choose
one common `n`, as large as desired, with

```text
13^n L_w>1+ell_w                     for every w.      (22)
```

The possible integers `k` for which

```text
(J_w+k)/13^n subset Y_w
```

form an open real interval of length

```text
13^n L_w-ell_w>1.
```

It contains an integer, proving the inclusion. Every blocker test at
relative depth `d` becomes the actual linked blocker

```text
c=13^(n+d+1)u,                                      (23)
```

by (3). Every delayed word becomes `Q(13^(n+d)y)`.

For necessity, the inverse image under `z=13^n y` of a finite
rational step set with empty interior also has empty interior. It
cannot contain a positive open subinterval. This proves (21).

Several packets can be installed recursively. After retaining one
set of intervals, choose the next common base delay above every
previous delay and above the new finite mesh threshold (22).
Consequently the base delays may obey arbitrary lower bounds or be
chosen in any prescribed cofinal subsequence.

If two roles share an exponent, they belong to the same `C_w`.
Individual feasibility is insufficient. If one physical factor is
used several times, every transported requirement on it must likewise
appear in the same joint cylinder. Asking one factor to be both safe
and dangerous at one clock gives the minimal empty-cell hostile.

For distinctness, it is enough to assume that base units occurring at
one exponent are pairwise distinct `13`-units. In general one must
exclude the finite collisions

```text
13^lambda u=13^lambda' u'.                              (24)
```

Shrinking the final intervals to a common rational length restores
equal chart weights.

## 5. Explicit same-depth linked realization of the THM-2462 atlas

The abstract criterion survives a hostile fixed-tuple test. Retain
the thirteen exact THM-2462 parent intervals, each of width

```text
L=1/6664627200,                                       (25)
```

and the first two fixed speeds

```text
u_K=14,

u_q=644=46u_K.                                        (26)
```

Prescribe the two chart words

```text
epsilon_K(g)=1_(g even),

epsilon_q
 =(0,0,1,1,0,1,1,1,1,0,1,1,0).                      (27)
```

The first word is a labelled parity control; the second is the exact
THM-2458 overlap word. They are not asserted to be canonical blocker
semantics.

Put `t=14z`. Then

```text
D_14(z)=D_1(t),

D_644(z)=D_46(t).
```

All four same-clock cells have the following explicit rational open
components in the original `z` coordinate (multiplying their
endpoints by fourteen gives the corresponding `t` components):

```text
danger,danger: (-1/9016,  1/9016),   width 1/4508;

danger,safe:   ( 1/9016, 13/9016),   width 3/2254;

safe,danger:   (321/9016,323/9016),  width 1/4508;

safe,safe:     (309/9016,321/9016),  width 3/2254.      (28)
```

Choose one common base delay

```text
n=9,                      lambda=n+1=10.              (29)
```

The two worst mesh margins are

```text
13^9 L-(1+3/2254)
 =3931001773/6664627200>0;

13^9 L-(1+1/4508)
 =3938393773/6664627200>0.                             (30)
```

Thus the appropriate component in (28) has a full pullback inside
every one of the thirteen old parent intervals. The linked physical
blocker speeds are

```text
c_K=13^10*14
   =1930018885886;

c_q=13^10*644
   =88780868750756.                                   (31)
```

Every old six-role root mask remains unchanged. The two new
root-constant bits equal (27) on all thirteen roots of every chart.
The retained interval widths are

```text
1/47805083173484

or

3/23902541586742.                                    (32)
```

Shrinking the larger cells gives thirteen pairwise disjoint
equal-weight intervals of width

```text
L_b=1/47805083173484.                                 (33)
```

This is a fixed-tuple, same-depth realization. It is strictly stronger
than the independent labels in THM-2462 Section 4.1. Its limit is
equally important: the chosen labels are controls, not a proof that
the physical factors are the source owner, expiration word, or
endpoint current of THM-2305.

## 6. Finitely many delayed words

Let `Q` be any fixed positive-mass rational Boolean circle word after all
physical speeds in the preceding construction have been fixed.
Choose a rational open interval

```text
J subset {Q=1},                  |J|<1.               (34)
```

If `Q=1` identically, choose any strict rational subinterval of the
circle. From (33),

```text
13^13 L_b=28561/4508>2>1+|J|.                          (35)
```

The mesh lemma therefore retains

```text
Q(13^13 y)=1                                          (36)
```

inside every one of the thirteen linked intervals. In the root
notation of THM-2397, this is the actual root-constant delayed word

```text
Q(Rx_r)=Q(13^13y),

R=13^14,                   clock exponent k=14.        (37)
```

The exact companion uses `J=(2/7,3/7)` and obtains thirteen
intervals of width

```text
1/2120125746145771.                                  (38)
```

Any finite list of fixed positive rational words can be added
recursively at larger chosen clocks. If a word must occur at the
same clock as a blocker packet, its support must instead be included
in the joint cell `C_w` of Section 4. A word depending on the clock,
on a changing root/tail gauge, or on an uninstalled semantic owner is
not a fixed observable and lies outside this theorem.

## 7. Sharp hostiles and strongest surviving conclusion

The exact failure boundaries are:

1. **Shallow multiplicity lock.** On

   ```text
   Y=(0,1/28),
   ```

   the depth-one unit blocker is dangerous everywhere, so a fixed
   shallow safe request has no positive realization. On every clean
   two-root/four-root phase cell the opposite lock (10) forces it
   safe.

2. **Empty joint cell.** One factor requested both safe and dangerous
   at the same clock has empty joint cell. More generally, if a
   same-clock family of danger combs covers the circle, its all-safe
   joint cell is empty. No amount of common delay repairs either
   obstruction.

3. **Fixed small clocks.** A globally nonempty target cell may miss a
   given parent at a prescribed small delay. The theorem gives
   arbitrarily large common delays, not every preassigned clock.

4. **Repeated transported uses.** Marginal nonemptiness for several
   copies of one physical factor is not enough. Their complete
   relative-clock cylinder must have interior.

5. **Changing semantic words.** The delayed-word corollary applies
   only after `Q` is one fixed rational step function in one common
   gauge. It does not manufacture the THM-2305 owner/word
   identification.

The source, target, map, and loss ledger are therefore:

```text
source:
  strict open clean root-phase charts plus compatible finite
  linked-factor/clock cylinders;

map:
  phase refinement before mixed-radix realization, or full inverse
  branches of the expanding map y -> 13^n y inside fixed parents;

preserved:
  literal old root masks, actual equations c=13^lambda u,
  finitely many fixed root-constant truth values, rational positive
  interval mass, and equal chart weights;

not supplied:
  a prescribed normalized shallow valuation profile, a full
  fourteen-speed scalar cover, semantic owner/expiration typing,
  repeated-factor compatibility outside one joint cylinder,
  endpoint phase, or a signed exact relation current.             (39)
```

Thus finite open phase compatibility, linked blockers of depth at
least two, and finitely many fixed delayed words still cannot exclude
the THM-2458 uniform-offset mechanism. A successful physical closure
must force one of the data deliberately left out of (39): most
sharply an anchored `lambda=1` dangerous bit, an empty semantic joint
cylinder, a bounded valuation profile which prevents delay, or a
canonical endpoint/current constraint.

No full scalar row is constructed or excluded. The ledger remains
`165`, and LRC(14) remains open.

## 8. Exact companion

Run

```text
python 04-computation/lrc14_linked_blocker_clock_cell_thm2464.py
python -O 04-computation/lrc14_linked_blocker_clock_cell_thm2464.py
```

The dependency-free `Fraction` companion:

- verifies the ordinary and guard shallow multiplicity laws for every
  nonzero residue and every root start, including all ordinary
  endpoint exceptions;
- checks all `1,248` orientation/residue/start/bit instances of the
  exact `lambda=2` phase refinement;
- verifies the four joint cells (28), both mesh margins (30), all
  twenty-six linked root-constant truth values, and all retained old
  root masks;
- checks the exact equal-width linked intervals and physical speeds;
- installs a positive delayed word at clock exponent fourteen; and
- includes the fixed shallow-clock and opposite-bit joint-cell
  hostiles, together with the explicitly conditional bank-halving
  arithmetic.

Both normal and optimized transcripts must reproduce

```text
05-knowledge/results/lrc14_linked_blocker_clock_cell_thm2464.out
```

byte-for-byte after LF normalization. Every truth-bearing executable
check uses explicit `require`; there are no `assert` or floating-point
truth tests.

QED.
