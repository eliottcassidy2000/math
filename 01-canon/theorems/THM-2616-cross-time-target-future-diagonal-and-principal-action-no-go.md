---
id: THM-2616
title: "Cross-time target/future diagonal and principal-action no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Enlarge THM-2600's complete constant-six carrier by retaining every
  delayed physical future digit h before imposing the numerical diagonal
  q=h with the present absolute target section q.  Across all 2,299,752
  nonzero-deep-root slots the 649,968 positive pre-route numerators have one
  global content 26; the route-two physical content is 338.  On the lawful
  q=h diagonal there are 49,872 positive fine entries, 1,703 positive
  rail/section pairs, and 1,483 unit Bockstein slices.  Every one of the 84
  base cells has 10 or 11 unit diagonal sections, every difference occurs,
  and eight constant sections work globally.  On the canonical q1=14 role,
  THM-2613's local shift s=-h and THM-2585's section convention q'=-s make
  this diagonal the framed partial label map q'=h.  Nevertheless the delayed
  digit carrier is exactly {1,...,11}; digits 0 and 12 are absent.  Hence it
  has no nonempty C13-invariant fibre, and every unit-section set has trivial
  translation stabilizer.  The result is a positive common-x cross-time
  correspondence, not a principal ancestry action, same-event state map,
  adjacent transition, row exclusion, or proof of LRC(14).
source: carry-transition-cell-2026-07-28-cross-time-diagonal
depends_on:
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - THM-2592-fallback-rail-digit-diagonal-pullback-and-primitive-bockstein
  - THM-2600-constant-six-middle-rail-common-x-atlas-and-uniform-bockstein-section
  - THM-2613-canonical-root-diagonal-opposite-shift-section
related:
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
  - THM-2609-external-target-section-itinerary-saturation-and-root-state-no-go
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2614-punctured-target-root-cosupport-and-principal-deck-no-go
script: 04-computation/lrc14_cross_time_target_future_diagonal_thm2616.py
output: 05-knowledge/results/lrc14_cross_time_target_future_diagonal_thm2616.out
script_sha256: 1f2de3f35873ed37a55098ff198f968cccbf3feed4a3cab2bdd016090c290117
output_sha256: 03e848358e7ec3ebaab22106253dc9e7e22c51a57771d0dbeae16005b4fb3366
hash_basis: LF-normalized bytes
---

# THM-2616 -- a cross-time diagonal exists, but it is not a deck action

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2609 shows that THM-2600's absolute target-section labels already have
every additive difference, but its physical arrival and delayed future digit
are both frozen at six.  One must not infer a state map from that coincidence.
There is, however, a lawful way to ask a stronger question: retain **every**
physical future digit first, form one globally primitive carrier, and only
then impose the numerical diagonal between the target section and future
digit.

That diagonal is unexpectedly positive and globally sectioned.  It still is
not a principal `C_13` transition.  The obstruction is visible before any
Fourier argument: the delayed word has only eleven physical digit states, so
the regular thirteen-cycle cannot act on its support.  This separates ordinary
support-theoretic sections from equivariant ancestry sections exactly.

## 1. The enlarged common-`x` carrier

Use THM-2600's canonical typed row, its `162` positive constant-six middle
rails, and THM-2592's notation.  A rail is

```text
e=(s,ell4,6,t),                    t in {12,0}.            (1)
```

For a present absolute target section `q`, present owner clock `ell5`, deep
probe `r`, and delayed physical future digit `h`, let

```text
A_(e,q;ell5,r,h)                                           (2)
```

be the common-denominator pre-route numerator of

```text
13 integral U(x)V(x-s/13) 1_e(x)
       F_(ell5,-q)(x) Delta_r(x) Q^+_(ell5,h)(Rx) dx.      (3)
```

Here

```text
q: external present target-section intervention;
h: physical digit floor(13 frac(Rx)) in the delayed word;
r: present translated deep-danger probe.                  (4)
```

All factors in (3) are evaluated on one physical `x` before any marginal or
transform.  Thus a positive entry is a genuine common-`x` Boolean pullback
component after expanding the two nonnegative Perron profiles.  But (4) is
load-bearing: `q` and `h` are coordinates at different roles/times.  Equality
of their printed residues is not yet a same-event state identification.

For each `ell5`, the thirteen half-open digit pieces form an exact disjoint
partition of the unconditioned delayed word:

```text
sum_(h in F13) Q^+_(ell5,h)=Q^+_ell5.                     (5)
```

The companion verifies (5) after intersection with every one of the
`162*13*7*12` nonzero-`r` fine fibres.  Thus the enlarged tensor is formed
before selecting a diagonal; it is not an artificial coupling of two
marginals.

## 2. One global primitive content is exactly `26`

Retain all

```text
162*13*7*12*13=2,299,752                                  (6)
```

nonzero-deep-root slots of (2).  Also retain the `191,646` `r=0` slots as
literal zeros, inherited from the present deepest-safe factor.  Exact integer
reconstruction gives

```text
#{A!=0}=649,968,

gcd{A_(e,q;ell5,r,h):A!=0}=26=2*13.                       (7)
```

The route-two factor `13` in THM-2592 makes the corresponding physical raw
content

```text
13*26=338=2*13^2.                                         (8)
```

Define the single globally primitive enlarged tensor by

```text
a=A/26.                                                    (9)
```

No future digit, rail, target section, or diagonal is reprimitivized
separately.  In particular, the diagonal subbank below independently has gcd
`26`, so it witnesses rather than changes the full-carrier content.

The value `26` is much smaller than THM-2600's constant-six content
`4,244,240`.  There is no contradiction: (7) takes the gcd over the strictly
larger carrier containing all future digits.  The Bockstein changes because
the primitive integral lattice has genuinely been enlarged, not because one
old slice was renormalized.

## 3. The target/future numerical diagonal

Now impose

```text
h=q.                                                       (10)
```

For a rail `e` and target section `q`, define its diagonal seven-clock
Bockstein from the globally primitive tensor (9):

```text
Y_ell5(e,q)
 =sum_(r=1)^12 a_(e,q;ell5,r,q) r^(-1) mod 13,

Y_(e,q)(z)=sum_(ell5=0)^6 Y_ell5(e,q)z^ell5
              in F13[z]/(Phi_7).                          (11)
```

Call `(e,q)` a unit diagonal slice when (11) is a unit.  The exact census is

```text
positive diagonal fine entries:          49,872 / 176,904;
positive diagonal rail/section pairs:     1,703 / 2,106;
unit diagonal rail/section slices:         1,483 / 2,106.  (12)
```

Every unit slice in (12) has a positive fine entry.  As in THM-2600, unitness
belongs to the complete `(ell5,r)` aggregate; it is never assigned to an
individual hidden Perron sheet.

For a base cell `(s,ell4)`, unite the two possible rails and put

```text
Q^diag_(s,ell4)
 ={q: some rail e above (s,ell4) has unit Y_(e,q)}.        (13)
```

Every one of the `84` sets is nonempty.  Their exact patterns are

| diagonal unit-section set | cells |
|:---|---:|
| `{1,2,3,4,5,6,7,8,9,10,11}` | 76 |
| `{1,2,3,4,5,6,7,8,9,10}` | 4 |
| `{1,2,3,4,5,6,7,8,9,11}` | 2 |
| `{1,3,4,5,6,7,8,9,10,11}` | 2 |

Equivalently,

```text
|Q^diag|=11 on 76 cells,
|Q^diag|=10 on  8 cells.                                 (14)
```

The eight exceptional cells are

| `(s,ell4)` | `Q^diag_(s,ell4)` |
|:---:|:---|
| `(2,2),(2,3),(7,2),(7,3)` | `{1,2,3,4,5,6,7,8,9,10}` |
| `(6,5),(11,5)` | `{1,3,4,5,6,7,8,9,10,11}` |
| `(6,6),(11,6)` | `{1,2,3,4,5,6,7,8,9,11}` |

## 4. Eight ordinary global sections and every difference

The intersection of all `84` section sets is

```text
intersection_(s,ell4) Q^diag_(s,ell4)
 ={1,3,4,5,6,7,8,9}.                                    (15)
```

Thus each of the eight constant choices in (15) gives an ordinary global
section of the **set-valued unit-support atlas**: in every base cell, some
positive rail realizes a unit slice with the same numerical equality `q=h`.
This is stronger than mere nonemptiness cell by cell.

Moreover every section set has every additive difference:

```text
Q^diag_(s,ell4)-Q^diag_(s,ell4)=F13                     (16)
```

for all `84` cells.  This follows already from `|Q^diag|>=10`, and the
companion checks it directly.  Hence neither finite support, global constant
sections, nor correction differences remain missing.

Equations (15)--(16) are still aggregate support statements.  The rail,
future clock, deep probe, and hidden positive Boolean component may vary with
the base cell.  They do not provide one sheetwise charged ancestry path.

## 5. Why the diagonal is not a principal `C_13` section

The exact delayed-digit phase support is

```text
h:                         0  1  2 3 4 5 6 7 8 9 10 11 12
positive owner phases:     0  6  7 7 7 7 7 7 7 7  7  6  0. (17)
```

In particular the physical future carrier available to (2) is exactly

```text
H_future={1,2,...,11}.                                    (18)
```

A nonempty free `C_13` orbit has thirteen points.  The eleven-point set (18)
therefore cannot be a principal `C_13` fibre.  In the displayed residue chart
its translation stabilizer is explicitly

```text
Stab_C13(H_future)={0}.                                   (19)
```

The same defect is visible after the unit test.  Every set in (13) is a
proper nonempty subset of `F13`, hence

```text
Stab_C13(Q^diag_(s,ell4))={0}                             (20)
```

on all `84` cells.  Translating the equality `q=h` by a nonzero residue exits
the available carrier at one of the missing future digits and fails to
preserve the unit section set.

This is the sharp holotopy distinction:

```text
ordinary support sections:       eight constant choices;
equivariant principal section:   impossible on this carrier. (21)
```

Abstractly, choosing one numerical diagonal would amount only to choosing
compatible coordinate origins in two printed `F13` charts.  On this typed row
there is a stronger fixed-gauge interpretation.  The ordinary all-root role is

```text
k=q1=14=1 mod 13.                                         (22)
```

THM-2613's canonical inverse-root/local-shift section is therefore
`s=-h`.  THM-2585 labels the literal target-shift slice `s` by `q'=-s`.
Consequently

```text
q'=h                                                       (23)
```

is the canonical **framed label orientation** on every physical future sheet
retained by (18), not an arbitrary sign choice.  Equations (18)--(20) still
show why this partial framed map is not a principal action: two future sheets
are absent.  Moreover THM-2613's gate is target-danger while the same-time
`F_(ell5,-q)` factor in (3) retains target safety.  They are later/complement
events, not one idempotent physical gate.  Thus (23) does not construct the
common action or thirteen-sheet bibundle classified by THM-2611.

## 6. Relation to the adjacent frontiers

The gain and loss ledger is exact.

1. THM-2609 keeps `h=6` and varies external `q`; this theorem first enlarges
   the physical delayed digit and then aligns `q=h`.  It therefore repairs a
   support correspondence, not the action law.
2. Equations (22)--(23) compose THM-2613's canonical local paired-shift label
   with THM-2585's absolute section convention on the eleven retained future
   digits.  This is a genuine partial `h -> q'` label connector.  It still
   does not identify the present deep probe `r` with `h`, and the safe/danger
   gate mismatch prevents a same-event physical composition.
3. THM-2614 retains the same-event pair `(q,r)` and finds a punctured product.
   Here `(q,h)` is cross-time and diagonal.  The former lacks the `r=0` sheet;
   the latter lacks future digits `h=0,12`.  Both fail principal closure for
   different physical reasons.
4. The arrival digit in (1) remains `v=6`, while `h=q` varies.  Except at
   `q=6`, the theorem deliberately gives up THM-2600's arrival/future
   diagonal.  The deep probe `r` and the owner clocks `ell4,ell5` also remain
   separate.

Thus no outgoing-deep-root-`r`-to-next-input map, adjacent-clock kernel,
semantic owner/repair endpoint, or relation-residue transport follows.  A
positive successor must add the two missing future sheets in one lawful word
cospan and connect `r` to the framed `h -> q'` map, or construct an enlarged
carrier with an honest `C_13` action and a selected equivariant graph.

## 7. Exact evidence and scope

Run

```text
python 04-computation/lrc14_cross_time_target_future_diagonal_thm2616.py
python -O 04-computation/lrc14_cross_time_target_future_diagonal_thm2616.py
```

The companion rebuilds the `162` rails and all thirteen delayed digit
prefixes.  Four deterministic exact shards exhaust the full universe (6),
checking (5) on every fine fibre and cross-checking selected vector sweeps
against thirteen independent scalar sweeps.  It then computes the one global
content, the complete positive and unit censuses, all `84` section sets,
their patterns, exceptional atlas, intersections, differences, and
translation stabilizers.  Every logical check is an explicit optimized-mode
guard.  Normal and optimized executions byte-match the stored transcript
after LF normalization.

An independent hostile audit replayed both executions, byte-matched the
stored transcript, and recomputed both declared LF hashes.  A second exact
evaluator replaced the optimized shared-endpoint digit vector by thirteen
separate scalar prefix sweeps on every fine fibre, checked each digit
partition against the whole delayed word, and recovered independently

```text
full positive census and global content:       649,968 and 26;
diagonal positive census and diagonal content:  49,872 and 26;
positive pairs and unit slices:                   1,703 and 1,483. (24)
```

It tested unitness by polynomial Euclid against `Phi_7` rather than by the
companion's multiplication determinant, and recovered all four section-set
patterns, the eight common sections, all differences, and every trivial
stabilizer.  The audit separately checked the sign composition: for
`k=14=1 mod 13`, THM-2613 gives the future-local slot `s=-h`, while THM-2585's
section `q'` selects `s=-q'`, hence `q'=h`.  This remains only the framed
cross-time label correspondence already scoped in Sections 5--6.  The
principal-action no-go is for the fixed translation-labelled digit carrier;
it does not exclude an enlarged carrier with a new physical action.

The theorem is confined to THM-2600's one canonical typed row and to the
enlarged carrier (2).  It does not claim that an aggregate unit Bockstein
charges one hidden sheet, that a constant support section has common owner or
root provenance, or that a different delayed word cannot supply the missing
digits.  No scalar row is excluded and LRC(14) remains open.

QED.
