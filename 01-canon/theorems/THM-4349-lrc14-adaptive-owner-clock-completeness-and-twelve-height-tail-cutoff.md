---
id: THM-4349
title: "LRC(14) adaptive owner-clock completeness and twelve-height tail cutoff"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13 AND THM-2061/2066/4070. For every
  eleven-element positive body H, the dyadic seam 2H union {a,b} is
  1/14-lonely for every two distinct positive odd tails if and only if one
  finite adaptive owner clock N has empty complementary-owner relation
  R_N(H). THM-2061's inherited metric box says it is enough to test the
  finite tail box a,b<12 max(H).
  On the safe side one explicit, deliberately coarse clock is
  28(42h+1)^2 lcm(1,...,12h), where h=max(H).
  Equivalently, failure of every clock produces an actual bounded odd-tail
  counterexample to this seam, not merely a profinite ghost. This proves
  completeness of the owner-clock certificate on the dyadic seam; it does
  not prove that the certificate always occurs and is not LRC(14).
source: codex-2026-09-02-LRC-adaptive-owner-clock
depends_on:
  - LRCUpTo13
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-4070-lrc14-d2-mod14-two-bank-affine-ray-firewall
related:
  - THM-4075-lrc14-divisor-complete-dyadic-owner-word-closure-through-30
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-4087-lrc-literal-open-two-comb-compiler
  - THM-4347-lrc14-dyadic-owner-word-closure-through-forty
  - THM-4371-lrc14-minimal-wall-lcm-adaptive-owner-clock
---

# THM-4349 -- adaptive owner clocks are complete on the dyadic seam

**PROVED RELATIVE TO CITED `LRCUpTo13` AND THM-2061/2066/4070.**  Fix an
eleven-element set `H` of distinct positive integers and put

```text
h=max H,                 S(H;a,b)=2H union {a,b},       (1)
```

where `a,b` are distinct positive odd integers.  Then the following are
equivalent.

1. There is a clock `N` for which the THM-2066 complementary-owner relation
   `R_N(H)` is empty.
2. `S(H;a,b)` is `1/14`-lonely for every distinct positive odd pair `a,b`.
3. The same assertion holds just for

   ```text
   1<=a<b<12h,                  a,b odd.                (2)
   ```

4. No distinct positive odd pair has the continuum spoiled-set containment

   ```text
   G(H) subset Sigma_2({a,b}).                          (3)
   ```

When these conditions hold, condition 1 may use the explicit clock

```text
N_H=28(42h+1)^2 lcm(1,...,12h).                        (3a)
```

**Subsequent sharpening.** THM-4371 proves that the same equivalence already
holds at the much smaller exact wall-LCM clock
`7 lcm({2c:c in H} union {a<12h:a odd})`.  The bound `(3a)` remains a valid
deliberately coarse clock, but is no longer the best proved explicit bound.

This bound is intentionally coarse but depends only on the body height, not
on the detailed geometry of `H`.  Thus the four conditions are also
equivalent to the single finite statement `R_(N_H)(H)=empty`.

Moreover, if every finite clock has `R_N(H)` nonempty, the compatible
clock obstructions define two odd profinite integers, but the open safe arc
forced by lower-dimensional LRC makes both profinite integers ordinary
integers of absolute value strictly below `12h`.  Thus a failure of adaptive
clocks is an actual unsafe seam row, not a residue-only or profinite false
positive.

The theorem proves completeness of a proof method.  It does **not** prove
that condition 1 holds for every `H`; doing so would close this dyadic seam.
It makes no claim about the other LRC(14) component types, and LRC(14)
remains open.

The inherited metric cutoff `(2)` and continuum criterion `(3)` are
THM-2061's mechanism, not new claims here; Section 2 repeats the short
strict-boundary proof so the completeness argument is auditable in one
place.  THM-4079 and THM-4087 manufacture adaptive clocks from special
antipodal intervals.  The new implication is the converse at fixed `H`:
if every actual tail pair is safe, then some single owner clock certifies all
of them simultaneously.  The max-`72` body in Section 6 is the canonical
hostile to replacing this adaptive clock by the previously used fixed bank.

**Finite decision corollary.**  For fixed `H`, universal safety against all
odd tail pairs is decidable: exactly test the finitely many distinct odd
pairs in `(2)`.  If one fails, it is an actual unsafe seam row.  If all pass,
the rational-witness construction in Section 4 gives the explicit clock
`(3a)` with empty relation.  There is no body-height-independent constant.

The connection contract is:

```text
source:    the divisibility inverse system of labelled relations R_N(H)
target:    continuum spoiled-set containment by bounded actual odd tails
map:       inverse-limit compatibility, then open-arc residue stabilization
preserves: closed core safety, strict tail danger, both lift labels, odd parity
loses:     a body-height-independent or practically small first-clock bound
sidecar:   the LRC(12) margin 1/12 and THM-2066 divisor transport
test:      odd pairs a<b<12 max(H), or one exact empty R_N(H).
```

## 1. Exact objects and boundary convention

Write

```text
||x||=dist(x,Z),
G(H)={y in R/Z: ||cy||>=1/14 for every c in H}.          (4)
```

Safety in `(4)` is closed and danger is strict.  For a positive integer
clock `N`, retain the labelled rational packet

```text
A_N(H)={0<=r<N: 14|cr|_N>=N for every c in H}.           (5)
```

For odd tail residues modulo `2N`, `R_N(H)` is exactly THM-2066's relation:
its pairs cover both labelled lifts

```text
x_j=(r+jN)/(2N),              j=0,1,                    (6)
```

by strict tail-danger at every `r in A_N(H)`.  Repeated residue classes are
retained, because distinct integer tails can coincide modulo `2N`.

An odd tail is dangerous on at most one lift in `(6)`: the two tail phases
differ by one half.  It is dangerous on one of them exactly when

```text
||a r/N||<1/7.                                         (7)
```

Consequently, if two odd tails cover both lifts, **each** tail satisfies
`(7)` and their owner bits are complementary.  This individual eligibility,
not merely a union bound for the pair, is load-bearing below.

After factoring the common odd gcd of `a,b`, THM-4070 gives the exact
continuum counterpart:

```text
S(H;a,b) is not 1/14-lonely
  iff G(H) subset Sigma_2({a,b}),                       (8)
```

where `Sigma_2` consists of the body phases for which the two tails spoil
both lifts.  Formula `(8)` includes equality correctly: a lift not in the
strict spoiled set has both tail distances at least `1/14`.

## 2. The inherited open arc and the `12h` tail cutoff

Cited `LRCUpTo13`, used only at its settled eleven-speed consequence,
supplies `y_0` with

```text
||c y_0||>=1/12                 for every c in H.       (9)
```

Each function `y -> ||cy||` is `c`-Lipschitz on the circle.  Since

```text
1/12-1/14=1/84,                                       (10)
```

the closed circle arc

```text
I={y: dist(y,y_0)<=1/(84h)}                            (11)
```

has length `1/(42h)` and lies in `G(H)`; its interior lies in the strict
interior of `G(H)`.

Suppose now that `S(H;a,b)` is unsafe.  By `(8)`, both lifts over every
`y in I` are spoiled.  An odd tail can kill at most one lift, so both tails
must individually be eligible there.  By `(7)`,

```text
I subset {y: ||ay||<1/7},
I subset {y: ||by||<1/7}.                              (12)
```

For a nonzero integer `m`, the set `{y:||my||<1/7}` is a disjoint union of
`|m|` open arcs, each of length `2/(7|m|)`.  The connected **closed** arc
`I` must lie in one open component.  Hence its length is strictly smaller
than the component length, so

```text
1/(42h)<2/(7a),        1/(42h)<2/(7b),
so                         a,b<12h.                    (13)
```

This proves `3 => 2`; the reverse implication is immediate.  It also
records why the constant is `12`: it is exactly the quotient of the doubled
tail-danger width `2/7` by the inherited body margin width `1/(42h)`.

## 3. A discrete long-block lemma

The finite-clock extraction uses the following elementary rigidity fact.

**Lemma.**  Let `M>=1`, let `q` be a residue modulo `M`, and choose its
centred representative `q_0`.  If

```text
|j q|_M<2M/7                 for j=1,...,L,             (14)
```

then

```text
|q_0|<2M/(7L).                                        (15)
```

**Proof.**  Put `theta=|q_0|/M`.  The case `theta=0` is immediate.  From
`j=1`, `0<theta<2/7`.  Suppose inductively that the ordinary, unwrapped
number `(j-1)theta` lies in `(0,2/7)`.  Then
`j theta<4/7<1`.  Condition `(14)` says its fractional part lies in
`[0,2/7) union (5/7,1)`.  The second arc is impossible because
`j theta<4/7`; hence `j theta<2/7`.  Induction through `j=L` gives
`L theta<2/7`, which is `(15)`.  The negative sign is symmetric.  QED.

Now let `J` be any closed subarc of `(11)` of length

```text
ell=1/(42h+1).                                         (16)
```

For a large clock `M`, its grid points in `J` include `K` consecutive
residues with

```text
K>=ell M-1.                                            (17)
```

If an odd tail residue `z mod 2M` is eligible on `A_M(H)`, then for those
grid points

```text
|zr|_M<M/7.                                            (18)
```

Subtracting the condition at the first point from the next `K-1` conditions
gives `(14)` for the centred reduction `q_0` of `z modulo M`, with
`L=K-1`.  Therefore

```text
|q_0| < 2M/(7(ell M-2)).                               (19)
```

The right side tends to

```text
2/(7ell)=12h+2/7.                                      (20)
```

Thus, for all sufficiently large `M`, every eligible residue has

```text
|q_0|<=12h-1.                                          (21)
```

There is no parity loss when `M` is even: an odd `z` reduced by a multiple
of even `M` still has an odd centred representative.  If an explicit size
condition is wanted, writing `A=42h+1`, `(19)` is below `12h+1` whenever

```text
M>(4A^2+10A)/5.                                        (22)
```

The integer representative is odd, while `12h` is even, so the bound
`|q_0|<12h+1` furnished by `(22)` gives `(21)`.

## 4. Rational witnesses and one finite adaptive clock

Assume condition 3.  For each distinct positive odd pair in the finite box
`(2)`, the safe set of `S(H;a,b)` is nonempty.  All of its walls have rational
endpoints: they are defined by finitely many inequalities

```text
||m x||>=1/14,                 m integer.               (23)
```

A nonempty safe set is not the whole circle, so one of its component
endpoints is a wall for a physical speed `m in S(H;a,b)`.  Every such speed
satisfies `m<12h`, and the wall has the form

```text
x=(14k+/-1)/(14m).
```

After doubling, its denominator divides `14m`.  Therefore the single
height-dependent denominator

```text
Q_h=14 lcm(1,...,12h)                                    (24)
```

captures a rational witness for every bounded pair:

```text
y_(a,b) in G(H) \ Sigma_2({a,b}),
y_(a,b)=r_(a,b)/Q_h.                                    (25)
```

This remains true when the only witness is a boundary point: closed safety
and strict danger make the wall itself admissible.  Since
`Q_h>=168h`, the `Q_h`-grid also meets the arc `I`; fix such a point `y_*`.

Write `A=42h+1` and choose the lower and upper clocks

```text
N=Q_h A^2,
M=2N=28 A^2 lcm(1,...,12h).                            (26)
```

Then `M> (4A^2+10A)/5`, so `(21)` applies.  Notice that `M` is exactly the
explicit clock `(3a)`.  It is astronomically large and is offered as a
completeness bound, not as a practical replacement for small adaptive
searches.

We claim `R_M(H)=empty`.  Otherwise take `(u,v) in R_M(H)`.  Both residues
are individually eligible by Section 1.  Reduce `u,v`, originally modulo
`2M=4N`, to centred residues `a_0,b_0` modulo `M=2N`.  By `(21)` they are
odd nonzero integers with

```text
|a_0|,|b_0|<12h.                                       (27)
```

THM-2066's divisor transport from clock `M=2N` to clock `N` sends
`(u,v)` to `(a_0,b_0) mod 2N` in `R_N(H)`.  If
`|a_0|!=|b_0|`, the corresponding witness `(25)` is an `N`-grid point
because `Q_h|N`; sign does not change a literal danger mask.  At that point
the two projected residues do not cover both lifts, contradicting membership
in `R_N(H)`.  If `|a_0|=|b_0|`, the two masks are identical and cannot cover
both lifts over `y_*`; again they are not in `R_N(H)`.  Notice that this also
handles the repeated residue classes deliberately allowed in `R_N`.

Therefore `R_M(H)=empty` at the explicit clock `(26)`.  This proves
`3 => 1`.  THM-2066 proves
`1 => 2`, and Sections 1--2 prove the remaining equivalences.

## 5. The profinite formulation and why there is no ghost gap

There is an exact compactness interpretation.  Define the clopen set

```text
Zhat_odd={alpha in Zhat: alpha=1 mod 2};
```

this means odd profinite integers, not the prime-to-`2` factor of `Zhat`.
A profinite integer acts canonically on the torsion circle `Q/Z`: to evaluate
`alpha x`, reduce `alpha` modulo the order of `x`.  Call
`(alpha,beta) in Zhat_odd^2` a **profinite spoiled pair** if,
for every rational `y in G(H)` and both solutions of `2x=y`, at least one of
`alpha x,beta x` has circle distance strictly below `1/14`.

Then

```text
R_N(H) nonempty for every N
  iff a profinite spoiled pair exists.                 (28)
```

For this inverse-limit argument, work with the symmetric ordered version of
`R_N(H)`; ordering changes neither emptiness nor the certificate, and makes
the two tail coordinates unambiguous.  For the forward direction, use the
cofinal divisibility chain `N_k=k!`.
Every finite `R_(N_k)(H)` is nonempty, and divisor transport maps every node
at level `k+1` to a node at level `k`.  Koenig compactness gives a compatible
infinite branch, hence two odd profinite integers.  Every rational body-safe
phase occurs in a sufficiently fine factorial packet, so the branch spoils
all its lifts.  Conversely, reduction of a profinite spoiled pair modulo
`2N` gives an element of every `R_N(H)`.

The word “profinite” adds no obstruction here.  Indeed, let `alpha` be either
tail of such a pair.  Section 1 gives

```text
||alpha y||<1/7               for every rational y in I. (29)
```

On a proper closed subarc `J subset I`, apply `(14)--(19)` to the centred
residues of `alpha modulo k!`.  They are uniformly bounded.  Compatibility
then makes them eventually equal: once two centred residues have absolute
value at most `B` and the smaller factorial exceeds `2B`, congruence modulo
that factorial forces equality.  Since factorial moduli are cofinal,
`alpha` is the diagonal image of an ordinary integer.  It is odd.  The
rational form of `(29)` then places all of `I` in one actual danger
component: a strict failure has a rational neighbourhood, while equality can
occur only on a rational danger wall and is itself forbidden by `(29)`.
Equation `(13)` gives
`|alpha|<12h`.  The same holds for `beta`.

Finally `|alpha|` and `|beta|` are distinct: equal absolute values have the
same two lift masks and cannot cover both lifts over the nonempty arc `I`.
Replacing signs by absolute values preserves all masks.  The rational-cell
argument of Section 4 upgrades spoilage on rational `G(H)` to the continuum
containment `(3)`.  Thus `(28)` collapses exactly to an actual bounded unsafe
tail pair.

## 6. Scope, sharp next problem, and fixed-bank hostile

The max-`72` divisor-complete body

```text
{1,5,11,13,17,19,23,37,41,70,72}                     (30)
```

is a finite-exact hostile to the fixed bank `15,...,43`, yet clock `47`
certifies it in the
[exact audit](../../07-reflections/lrc14-owner-escape-fixed-bank-boundary-owner_escape_arbitrary_entry_20260902.md).
It demonstrates why the existential clock in condition 1 must be adaptive.
The present theorem explains that this is the only kind of delay possible:
an endless delay would already be an actual seam counterexample with both
tails strictly below `864`.

The next sharp problem is therefore no longer a vague all-clock compactness
question.  It is one of the following equivalent concrete tasks for a given
body `H`:

```text
find one empty R_N(H),
or find one unsafe odd pair a<b<12 max(H),
or rule out all continuum containments (3) in that finite box.             (31)
```

The explicit clock `(3a)` grows exponentially with `h`; no
body-height-independent or practical first-clock bound is asserted.  No
assertion is made that every eleven-body satisfies the safe side of `(31)`.
**QED.**
