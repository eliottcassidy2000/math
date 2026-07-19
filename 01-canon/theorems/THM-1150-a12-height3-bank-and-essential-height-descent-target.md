---
id: THM-1150
title: Height-three A12 rigidity and the all-height owner-stalk law for essential lowering
status: PROVED all-height local lowering equivalence and central-cylinder obstruction; FINITE-EXACT height-three and named high banks. Primitive essential-height descent PEHD13 remains OPEN
source: codex-2026-07-18-S74
depends_on:
  - THM-1142  # exact essential-region lowering criterion
  - THM-1143  # all-height A12 mechanical carrier
related:
  - THM-769   # shallow/deep split
  - THM-770   # bounded shallow classification through height 12
  - THM-795   # proper one-coordinate lifts of AP rays
scripts:
  - 04-computation/lrc13_a12_h3_essential_regeneration_bridge_codex_20260718.py
  - 04-computation/lrc13_a12_lowering_owner_stalk_codex_20260718.py
outputs:
  - 05-knowledge/results/lrc13_a12_h3_essential_regeneration_bridge_codex_20260718.out
  - 05-knowledge/results/lrc13_a12_lowering_owner_stalk_codex_20260718.out
---

# THM-1150 -- essential lowering is a local owner-stalk rule

Let a shallow full-residue packet be labelled as

```text
W=(w_r)_(r=1)^12,       w_r=r+13h_r,       h_r>=0.
```

Write

```text
D_w={t in R/Z : ||wt||<=1/13}
```

for the closed danger comb.  Call `W` **accepted** when these twelve combs
cover the circle, equivalently when the strict lonely set in THM-1143 is
empty.

For a coordinate `r` with `w_r>13`, put

```text
w_r^-=w_r-13,
E_r(W)=(R/Z) \ union_(s!=r) D_(w_s).                     (1)
```

`E_r` is the open essential region left by deleting colour `r`.

## 1. Exact all-height lowering criterion

THM-1142 gives the global statement

```text
lowering r is accepted  iff  E_r(W) subset D_(w_r^-).    (2)
```

The new result makes (2) completely local on the thirteen-sheet cylinder.
Write `t=(j+u)/13`.  At a generic `u`, put

```text
k=floor(w_r u),
Q_r(u)={A,B}={-k r^(-1),-(k+1)r^(-1)} subset F_13.       (3)
```

Let `deg_W(j,u)` count how many of the twelve sheet edges contain `j`, and
define the unique-owner stalk

```text
U_r(u)={j in Q_r(u): deg_W(j,u)=1}.                      (4)
```

> **Theorem A (owner-stalk lowering law).**  If `W` is accepted, then in
> every generic chamber the sheets uncovered by the eleven-runner core are
> exactly `U_r(u)`.  Consequently, lowering `r` is accepted if and only if
>
> ```text
> U_r(u) subset Q_r^-(u)                                 (5)
> ```
>
> in every chamber of the common refinement of the old and lowered walls.

**Proof.**  If the full packet covers a sheet and the core does not, that
sheet is covered by `r` alone, so it belongs to (4).  The converse is
immediate.  Thus (5) is exactly the sheet form of (2).  ∎

The acceptedness hypothesis is essential.  In a loose packet, a sheet may be
uncovered both by the core and by `r`; it remains an essential obligation but
does not appear in the unique-owner set.

## 2. Functional form of the lowering drift

For generic `u`, write

```text
13u=m+theta,       0<theta<1,
w_r u=k+alpha,     0<alpha<1.
```

> **Theorem B (exact drift law).**  With
>
> ```text
> a(u)=k-floor((w_r-13)u),
> bar_a(u)=a(u) mod 13,
> ```
>
> one has
>
> ```text
> a(u)=m+1_(alpha<theta),                                (6)
> Q_r^-(u)=Q_r(u)+bar_a(u)r^(-1).                        (7)
> ```

**Proof.**  Since `(w_r-13)u=k+alpha-m-theta`, its floor is
`k-m-1_(alpha<theta)`, which proves (6).  Substitution in the two endpoints
of (3) proves (7).  ∎

The raw drift can equal `13` on the final fringe; in `F_13` that is shift
zero. Order the old edge as `(A,B)=(-kr^(-1),-(k+1)r^(-1))`. A two-point needle
on the 13-cycle intersects its translate in the following exact way:

```text
bar_a=0:   {A,B},
bar_a=1:   {A},
bar_a=12:  {B},
else:      empty.                                         (8)
```

Combining Theorems A and B gives a chamberwise decision rule:

```text
bar_a=0    permits any U_r subset {A,B};
bar_a=1    permits only U_r subset {A};
bar_a=12   permits only U_r subset {B};
else       requires U_r=empty.                            (9)
```

This is the exact functional form of the `H`-drift.  Its state is not a
scalar height: it retains the strip index, a strict comparison bit
`alpha<theta`, the reduction `a mod 13`, and the oriented private endpoint.
Equivalently in raw values, both `a=0` and `a=13` are free cases.

## 3. Central-cylinder obstruction

> **Corollary C.**  On
>
> ```text
> 2/13<u<11/13,                                           (10)
> ```
>
> the old and lowered edges are disjoint.  Hence every lowerable coordinate
> of an accepted packet has no unique ownership anywhere in this central
> cylinder.

**Proof.**  In (10), `m` is one of `2,...,10`, so (6) gives
`a in {2,...,11}`, so `bar_a=a`. The final case of (8) applies. ∎

The remaining lowering obligations live on the four fringe strips.  There
the `a=1` and `a=12` cases preserve an oriented leading/trailing endpoint;
central inessentiality is necessary but not sufficient.

## 4. Endpoint completeness

The tests above use generic chambers without losing wall cases.  Indeed,
`E_r` is open and `D_(w_r^-)` is closed.  Therefore

```text
E_r \ D_(w_r^-)
```

is open.  If containment fails at a wall or a thirteenth boundary, it fails
on a neighbouring generic point.  Equal wall times are grouped rather than
arbitrarily ordered.  Thus (5) and (9) are exact, not almost-everywhere
surrogates.

## 5. Exact height-three classification

> **Theorem D (finite exact).**  If `0<=h_r<=3` for every `r`, then `W` is
> accepted if and only if, as an unlabelled set,
>
> ```text
> W=c{1,...,12},       c in {1,2,3,4}.                    (11)
> ```

The exhaustive checker forms the common generic-sheet atlas of 758 chambers
and 9,854 sheet bits, splits the twelve coordinates `6+6`, and checks all

```text
4^12=16,777,216
```

packets.  Exactly the four rows in (11) survive.  The decision digest is

```text
83c336bcda06319f13bb871537fa7966d4470988c2851a450cd80b9fedcbb322.
```

This is an independent finite control and is weaker in bounded scope than
THM-770, which reaches height twelve by a different owner CSP.  Its purpose
here is to validate the exact all-height carrier and the transition into the
new lowering law.

## 6. Named high-bank stress tests

The same exact chamber checker proves, without extrapolation:

1. for each `L in {13,25,50,100}`, among the 4,096 packets with
   `h_r in {0,L}`, only `[12]` is accepted;
2. for each `c in {14,17,97}`, the complete ternary coordinate cube
   `w_r in {w_r(c)-13,w_r(c),w_r(c)+13}` has `3^12=531,441` rows and only
   its central dilation `c[12]` is accepted.

These are exact finite banks, not an all-height theorem.

## 7. Sharp guardrails against naive descent

Every legal coordinate lowering on the accepted rays

```text
c=2,3,4,14,17,97
```

fails: 59 of 59 lowerings are nonaccepted, and all 59 carry a central
unique-owner obstruction.  The witness digest is

```text
eda2ee57ef8a472f70c29caba4be90dc2712131331db3954d442b16c0163066e.
```

Thus high dilations are lowering-sinks; any global descent lemma must include
a primitive/dilation alternative.

Acceptedness cannot be omitted either.  For

```text
W=(1,2,3,4,5,6,7,8,9,23,11,25),
```

lowering the residue-12 coordinate passes the unique-owner proxy but fails
the real essential containment at `u=45/506`: the core-uncovered set is
`{9}`, the unique-owner set is empty, and the lowered edge is `{1,2}`.

## 8. The sharpened open bridge

The smallest all-height shallow target is now:

> **PEHD13 (OPEN).**  Every accepted primitive full-residue packet with
> `max h_r>=13` has a coordinate `r` with `h_r>=13` satisfying (9) in every
> chamber, equivalently `E_r subset D_(w_r-13)`.

It separates into two still-open obligations:

1. **central-colour elimination:** some high colour has no private trace on
   the central cylinder (10);
2. **oriented fringe alignment:** that same colour obeys the one-sided
   endpoint rules on the four fringe strips.

If PEHD13 were proved, repeated lowering with primitive normalization would
enter THM-770's bounded box.  THM-770 classifies the terminal packet, while
THM-795 excludes a final accepted proper lift from an AP ray.  This would
close arbitrary winding in THM-769's **shallow** full-residue branch.

It would not close THM-769's deep `s>=2` branch, the all-loose crown of
THM-1149, compact `INVcov`, or LRC(14).

## 9. Verification and carrier audit

Both scripts replay byte-identically with normal Python and `python3 -O`.
Frozen hashes are

```text
h3 source     0074a5e5c92283f8150f7ff024f8f4f597ecf70d71ddecfe0c94bdd09e86ede1
h3 output     b115d7c54089e8783606959d1a3de9234753a8a76281fd465d773b2046634e21
stalk source  9cd28fb756a27f86e3bbbe5af89f446ffca18c22d472facf66e00ca03cc908c2
stalk output  daf2269a29f3e159a08a3ab3ef2f68f95df34753df6a252ca45511637bb8c697
```

Tournament Analysis challenges runners as the vertex set.  The useful
vertices are labelled moving two-sheet needles or private cut obligations;
the pair observable is private ownership in a chamber, lowering translation
`bar_a r^(-1)` is the switch, and chronological grouped walls give the tie
Hamiltonian path.  The carrier

```text
(labelled wall chronology, sheet degrees, unique-owner stalk,
 strip/endpoint orientation)
```

preserves exact lowering.  A runner tournament, the degree vector alone, or
a central-owner count destroys chronology, owner identity, or fringe
orientation and cannot decide PEHD13.
