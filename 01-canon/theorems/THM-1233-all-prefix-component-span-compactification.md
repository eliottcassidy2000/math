---
id: THM-1233
title: All-prefix component-span closure of every adjacent ratio
status: PROVED analytic composition / COMPUTER-EXACT constant audit.  Every putative six-comb cover of a complete c-slow gap has d_1/c<13/6 and adjacent-ratio bounds d_2/d_1<273/29, d_3/d_2<84/5, d_4/d_3<343/15, d_5/d_4<189/8, d_6/d_5<77.  Individual-cap ceiling propagation sharpens the absolute projective box to d_6/c<2345 while adding the previously missing all-prefix adjacent bounds.  This does not bound the carrier c or its phase clock and does not prove universal six-comb noncoverage or LRC(14)
source: codex-2026-07-18 all-prefix suffix-component session
depends_on: [THM-1196, THM-1198, THM-1232]
related: [THM-1176, THM-1178, THM-1199]
script: 04-computation/lrc14_all_prefix_component_span_compactification_codex_20260718.py
output: 05-knowledge/results/lrc14_all_prefix_component_span_compactification_codex_20260718.out
---

# THM-1233 -- all-prefix adjacent-ratio closure

At radius `1/14`, put

```text
D_s={t in R: ||st||<1/14}
```

and let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)],       |G|=6/(7c)   (1)
```

be a complete closed safe gap of the integer carrier `c`.  Suppose

```text
c<d_1<d_2<d_3<d_4<d_5<d_6.             (2)
```

Assume that the six strict danger combs cover `G`.  THM-1196 bounded the
last two adjacent ratios by looking at the survivor of the first five and
first four combs.  The same idea does not stop at `j=4`: THM-1198 also
controls the connected components of every suffix union.  The resulting
recursion removes all projective escape directions at once.

## 1. Suffix components cannot cross a slow tooth

> **Lemma 1 (universal suffix-component span).**  Let
>
> ```text
> s_1<...<s_r,                         1<=r<=6,       (3)
> U=D_(s_1) union ... union D_(s_r).
> ```
>
> Every connected component of the real lift of `U` meets at most one tooth
> of `D_(s_1)`.  Consequently every component has span strictly less than
>
> ```text
> 13/(7s_1).                                           (4)

**Proof.**  If a connected component met two teeth of `D_(s_1)`, then the
interval joining suitable points of those teeth would contain a complete
safe gap of `D_(s_1)`.  On that safe gap the other `r-1<=5` combs would give
a cover.  Enlarging their strict teeth to closed teeth would contradict
THM-1198's universal five-comb noncoverage theorem; its dual-load proof works
verbatim for any number at most five.  Thus at most one slow tooth occurs.

The component is therefore trapped strictly between the preceding and the
following `s_1` teeth.  This carrier cell consists of one tooth of width
`1/(7s_1)` and its two adjacent safe gaps of width `6/(7s_1)`, giving

```text
6/(7s_1)+1/(7s_1)+6/(7s_1)=13/(7s_1).                (5)
```

This first gives a weak span bound.  If equality held, the component would be
the full open interval between the two neighbouring `s_1` teeth.  In
particular the other at-most-five open combs would cover the interior of one
of the adjacent complete `s_1` safe gaps.  Their finite union of **closed**
tooth enlargements is closed and contains that interior, hence also its
closure.  This contradicts THM-1198.  The bound is therefore strict, proving
(4).  Notice that no harmonic-pressure hypothesis, gcd hypothesis, or bound
on the faster ratios is used. ∎

For suffixes of size at most three, chronology gives sharper spans.

> **Lemma 2 (one-, two-, and three-comb spans).**  For `p<q<r`, connected
> components obey
>
> ```text
> span(D_p)=1/(7p),                                   (6)
> span(D_p union D_q)<1/(7p)+2/(7q)<3/(7p),          (7)
> span(D_p union D_q union D_r)
>   <1/(7p)+2/(7q)+4/(7r)<1/p.                       (8)

**Proof.**  Equation (6) is the tooth width.  (A compact survivor component
contained in that open tooth has strictly smaller length.)  One `q` tooth is shorter than
the safe gap between adjacent `p` teeth, so a two-comb component contains at
most one `p` tooth.  Every attached `q` tooth extends by at most one tooth
width on each side, proving (7).

For (8), a component containing no `p` tooth already obeys the two-comb bound
(7).  Otherwise it contains at most one `p` tooth: a component of
`D_q union D_r` has span less than `3/(7q)<6/(7p)`, so it cannot bridge the
safe gap between two adjacent `p` teeth.  Remove that unique `p` tooth.  On
each side, every remaining piece lies in a component of `D_q union D_r` which
meets the `p` tooth.  By (7) such a component has span less than
`1/(7q)+2/(7r)`.  At most one such span is paid on each side, so

```text
span<1/(7p)+2[1/(7q)+2/(7r)].                         (9)
```

Since `p<q<r`, (9) is strictly smaller than `7/(7p)=1/p`. ∎

The coefficients `1,2,4` in (8) are the toothpick self-similarity: deleting
the slow tooth leaves one faster component on each side.  At four suffix
levels the naive recursion would give `1+2+4+8=15`; Lemma 1's carrier-cell
constant `13` is then better.  This is the exact point at which gap
noncoverage beats blind binary branching.

## 2. A survivor after every prefix

For `1<=j<=5`, define the compact open-tooth survivor

```text
E_j=G minus union_(i=1)^j D_(d_i).                    (10)
```

Let `n_i` be the number of `d_i` teeth meeting `G`, and put

```text
C_j=1+sum_(i=1)^j n_i,
Cbar_j=1+sum_(i=1)^j ceil((6d_i+c)/(7c)).             (11)
```

Deleting finitely many relatively open tooth intervals from `G` gives

```text
number of connected components of E_j <=C_j<=Cbar_j. (12)
```

THM-1198 supplies the six-bin probability density `f` with

```text
0<=f<=7/6,
integral_(one closed normalized comb) f<=7/36.        (13)
```

The complement of the first `j` **closed** combs is contained in (10).
The union bound in `f`-mass and then (13) therefore give normalized
Lebesgue survivor length at least

```text
[1-j(7/36)]/(7/6)=(36-7j)/42.                        (14)
```

Returning through the slow-gap scale `6/(7c)` proves the all-prefix floor

```text
|E_j|>=(36-7j)/(49c),                  1<=j<=5.       (15)
```

For the present strict faster slopes the inequalities can be made strict,
but the weak form (15) is sufficient.

## 3. Every adjacent ratio is bounded

Under the six-cover, `E_j` is contained in the suffix union

```text
D_(d_(j+1)) union ... union D_(d_6).                  (16)
```

Each of its at most `C_j` connected components lies in one connected
component of (16).  Apply Lemma 1 for suffix sizes four and five, Lemma 2 for
sizes one through three, and combine with (15).  The five rows are

```text
j  survivor floor       suffix-component span       consequence
1  29/(49c)             <13/(7d_2)                  d_2/c<(91/29)C_1
2  22/(49c)             <13/(7d_3)                  d_3/c<(91/22)C_2
3  15/(49c)             <1/d_4                      d_4/c<(49/15)C_3
4   8/(49c)             <3/(7d_5)                   d_5/c<(21/8)C_4
5   1/(49c)             <1/(7d_6)                   d_6/c<7C_5.       (17)
```

The last two rows recover THM-1196's component-span tails.  The first three
are the new all-prefix closure.

Put `x=d_j/c>1`.  Since `d_i<=d_j` for `i<=j`,

```text
ceil((6d_i+c)/(7c))<=ceil(x)<x+1,                    (18)
Cbar_j/x<j+(j+1)/x<2j+1.                             (19)
```

Divide the five rows of (17) by `d_j/c`, use (12), and insert (19).  Together
with THM-1198's first-tooth law, this proves the scale-free chain.

> **Theorem 3 (all adjacent ratios).**  Every cover (1)--(2) obeys
>
> ```text
> d_1/c   <13/6,
> d_2/d_1 <273/29,
> d_3/d_2 <84/5,
> d_4/d_3 <343/15,
> d_5/d_4 <189/8,
> d_6/d_5 <77.                                       (20)

Here the component row itself gives `d_3/d_2<455/22`; THM-1232's independent
functional gate `d_3/c<84/5`, together with `d_2>c`, gives the sharper value
displayed in (20).  Thus there is no longer a prefix-escaping branch in the
projective six-ratio cone.

## 4. A sharpened absolute projective bound

Keeping the exact tooth-count ceiling of **each** earlier speed is
substantially sharper than multiplying (20), and sharper still than replacing
all earlier caps by the largest current cap.  Write

```text
T(x)=ceil((6x+1)/7).                                  (21)
```

The first-tooth law and the `j=1` row of (17) give

```text
d_1/c<13/6 => T(d_1/c)<=2,
d_2/c<(91/29)(1+2)=273/29.                           (22)
```

For the next rung, the functional envelope is much stronger than the blind
four-suffix component cell: THM-1232 proves `d_3/c<84/5`.  Combining that gate
with the **individual** caps in the all-prefix rows (17) gives

```text
d_3/c<84/5                         => T(d_3/c)<=15,
d_4/c<(49/15)(1+2+9+15)=441/5      => T(d_4/c)<=76,
d_5/c<(21/8)(1+2+9+15+76)=2163/8   => T(d_5/c)<=232,
d_6/c<7(1+2+9+15+76+232)=2345.                      (23)
```

> **Corollary 4 (projective compactification).**  Every putative six-comb
> cover lies in the fixed compact ratio box
>
> ```text
> 1<d_i/c<2345,                           i=1,...,6.  (24)

The adjacent bounds (20) are the scale-free form; (24) strictly improves
THM-1232's former `40747` box by retaining the nonuniform tooth-cap ledger.

This is not a finite enumeration of integer packets: the primitive carrier
`c` and the phase address `k mod c` remain unbounded.  A sequence with
`c->infinity` can make several normalized ratios coalesce, while the offset
speeds and the clock phase retain nontrivial first-order geometry.  The
remaining object is therefore a compact projective packet plus a tangent
offset/phase stalk, not a finite list of rational ratios.

## 5. Kakeya, toothpick, and tournament carrier audit

The Kakeya needle is the protected gap `G`.  Its faithful hierarchy is now

```text
prefix survivor components
  -> suffix-union components
  -> their slowest tooth
  -> one faster component attached on each side.     (25)
```

This rooted binary attachment is the ladder's toothpick self-similarity.
The carrier-cell cap `13` truncates the binary coefficients exactly when
blind recursion becomes worse.

On speed vertices, use the pairwise observable `log(d_j/d_i)` and orient by
positive sign; strict speed order is the tie gauge.  The tournament is
transitive, with score histogram `0,1,2,3,4,5`, no directed cycles, six
singleton SCCs, and one Hamiltonian path.  It preserves the order needed in
(18), but destroys prefix ownership, component chronology, the carrier gap,
and which suffix component attaches on each side.

The challenged alternatives were runners, teeth, gaps, endpoints, wall
events, residues, Farey cells, components, and proof obligations.  The
smallest faithful vertices for this theorem are **suffix-union components**;
teeth are rooted sidecars and prefix survivors are the occupants.  A naked
tournament cannot express the non-crossing predicate in Lemma 1.

## 6. Exact audit and honest frontier

The companion dependency-free exact replay verifies every rational constant,
the prefix floors, all ceiling recurrences, the adjacent-ratio chain, and
finite endpoint-merging instances of Lemmas 1--2.  Normal and optimized runs
are byte-identical to the stored output.  Frozen SHA-256 hashes are

```text
source  876c3b4c5cf917fbd1d97e47813a2d5f56cb7c11eff87d3e5e037aec5646b094
output  882423017b2c1f79c424fb90b2c0645e423c7d072942dec41b14551b8102130f
```

The theorem adds all-prefix adjacent-ratio control to THM-1232's absolute
compact box.  It does **not** exclude that box, control the
unbounded carrier/phase clock, prove the continuum three-path inequality,
extract an AP from an arbitrary tight twelve-set, or prove LRC(14).  The next
honest target is the tangent degeneration of (24): classify the offset/phase
stalks arising when normalized speeds coalesce, and show that none can lift
to a strict integer cover.
