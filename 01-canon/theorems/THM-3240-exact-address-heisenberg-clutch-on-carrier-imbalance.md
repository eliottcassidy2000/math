---
id: THM-3240
title: "Exact-address Heisenberg clutch on carrier imbalance"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/creative-reframes/2026-08-03
audit: >
  The assertion-independent exact companion pins the four proved inputs and
  exhaustively checks the group-law signs, faithful image, orbit and
  stabilizer data, raw exact-address constraint, section-change conjugacy,
  carrier-forgetting projection, and integral odometer hostile.  Ordinary
  and optimized runs byte-match the stored transcript.  An independent
  hostile audit rederived the quotient dimensions, action signs, raw
  constraint, orbit/core/centre claims, section conjugacy, physical loss
  ledger, and odometer boundary, and replayed every exact artifact.
depends_on:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-3228-four-jet-heisenberg-minimal-faithful-permutation-carrier-gate
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
script: 04-computation/lrc14_exact_address_heisenberg_clutch_thm3240.py
output: 05-knowledge/results/lrc14_exact_address_heisenberg_clutch_thm3240.out
script_sha256: d21c5b7d750fbfb06bec38159517e7654d0ac929a040045440cb1f50dcf806b8
output_sha256: 535548d9f8b6ed1223133bd09f3c11a7a11239abcbeb71109af7b485d0387a9e
hash_basis: LF-normalized bytes
---

# THM-3240 -- Exact-address Heisenberg clutch on carrier imbalance

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

For the fixed `p=13` exact-address row of THM-2763, a choice of one labelled
target axis and of the typed Bezout section gives the lawful standard
Heisenberg action on the target coordinate together with carrier imbalance.
It acts faithfully on the `2,197`-element imbalance quotient and lifts to the
`28,561`-element two-sided quotient after choosing which labelled endpoint
harmonic to hold fixed.

This is an exact finite permutation clutch on address residues.  It is not a
physical endpoint-current action and does not supply an LRC owner map.

## 1. Fixed row, quotients, and choices

Use THM-2763's fixed word

```text
p=13,
W=(1,14,27,40,53,66,13,13^3,2*13^5).              (1)
```

Write

```text
K_13=ker(W:F_13^9 -> F_13),
G=K_13/L_13,                                         (2)
```

where `L_13` is the six-dimensional owner-pivot relation packet.  Then
`dim G=2`.  Label the two target axes as in THM-2350 and fix the first one:

```text
K_13=L_13 direct_sum span(e_a,e_b),
[r_sharp]=s[e_a]+t[e_b] in G.                       (3)
```

Thus `s` is the active target coordinate and `t` is a spectator.  Choosing
`b` instead gives the symmetric construction with the labels exchanged; no
unlabelled preferred target axis is asserted.

The exact-address imbalance quotient is

```text
K_delta={(r,delta): r.W+delta=0},
G_delta=K_delta/(L_13,0),
delta=k-l.                                           (4)
```

Because `W_0=1`, the coordinate vector `e_0` is a typed Bezout section.  Put

```text
r_sharp=r+delta e_0 in K_13.                         (5)
```

Equations `(3)` and `(5)` identify

```text
G_delta = {(s,t,delta):s,t,delta in F_13}            (6)
```

as a set with labelled coordinates.  The identification uses both the
target-axis label and the section; those choices are part of the statement.

## 2. The Heisenberg action

Use the THM-2779 convention

```text
H_13=F_13^3,
(x,y,c)(x',y',c')
  =(x+x',y+y',c+c'-y x').                            (7)
```

Define the following permutation of `(6)`:

```text
rho(x,y,c)(s,t,delta)
  =(s+x, t, delta+c-y s).                            (8)
```

Equivalently, on a raw representative satisfying `(4)`, let

```text
delta'=delta+c-y s,
r'=r+x e_a+(y s-c)e_0.                               (9)
```

Then

```text
r'.W+delta'
 =r.W+x e_a.W+(y s-c)e_0.W+delta+c-y s
 =r.W+delta
 =0.                                                  (10)
```

Moreover,

```text
r'+delta'e_0
 =r+delta e_0+x e_a,                                 (11)
```

so `(9)` has precisely the quotient coordinates in `(8)`.  Replacing `r` by
`r+lambda`, with `lambda in L_13`, replaces `r'` by the same `lambda`.
Replacing either displayed lift `e_a` or `e_0` by that lift plus an
`L_13`-vector changes `(9)` only by an `L_13`-vector.  Hence `(8)` is
well-defined on `G_delta` for the fixed classes and labels.

## 3. Sign ledger and faithfulness

Apply `h'=(x',y',c')` first and `h=(x,y,c)` second.  Formula `(8)` gives

```text
s''=s+x'+x,
delta''
 =delta+c'-y's+c-y(s+x')
 =delta+(c+c'-y x')-(y+y')s.                         (12)
```

This is exactly `rho(hh')` for the product `(7)`.  Thus `rho` is a left
action; in particular the cross-term and its sign are forced.

If `rho(x,y,c)` is the identity, the first coordinate gives `x=0`.
Evaluating the last coordinate at `s=0` gives `c=0`, and varying `s` then
gives `y=0`.  Therefore

```text
rho:H_13 -> Sym(G_delta) is faithful.                (13)
```

The useful generators are

```text
X=(0,-1,0): (s,delta) |-> (s,delta+s),
Y=(1,0,0):  (s,delta) |-> (s+1,delta),
Z=(0,0,1):  (s,delta) |-> (s,delta+1),               (14)
```

with

```text
[X,Y]=XYX^(-1)Y^(-1)=Z.                              (15)
```

This agrees exactly with the standard affine action in THM-2779 after the
renaming `(r,w)=(s,delta)`.

## 4. Orbits, stabilizer, and the centre

The spectator `t` is invariant.  On each fixed-`t` fibre, `Y` and `Z` act
transitively on the `13^2` pairs `(s,delta)`.  Consequently

```text
G_delta has 13 orbits, each of size 169.             (16)
```

At `(s,delta)=(0,0)` the stabilizer is

```text
S={(0,y,0):y in F_13}.                               (17)
```

It is a noncentral, nonnormal subgroup of order `13`; conjugation by `Y`
already changes it.  The centre

```text
Z(H_13)={(0,0,c):c in F_13}                          (18)
```

acts by nonzero translations of `delta`, so every nonidentity central
element is fixed-point-free.  Each fibre in `(16)` is therefore the minimal
`169`-point faithful Heisenberg carrier of THM-3228, and `G_delta` is the
disjoint union of thirteen labelled copies.

## 5. Lift to the two-sided exact-address quotient

Retain both endpoint harmonics:

```text
K_full={(r,k,l):r.W+k-l=0},
G_full=K_full/(L_13,0,0).                             (19)
```

Put `delta=k-l` and use the same sharp coordinate `(5)`.  THM-2763 gives

```text
G_full = {(s,t,k,l):s,t,k,l in F_13}                 (20)
```

after the fixed section and target basis are chosen.  Hold the labelled
target harmonic `l` fixed and define

```text
s'=s+x,
t'=t,
k'=k+c-y s,
l'=l,
r'=r+x e_a+(y s-c)e_0.                               (21)
```

Then `k'-l'=delta+c-y s`, so `(21)` projects equivariantly to `(8)`, and

```text
r'.W+k'-l'=0.                                        (22)
```

The proof of the group law and well-definedness is unchanged.  This action
is faithful and has

```text
169 orbits of size 169, indexed by (t,l).            (23)
```

The fixed-`l` convention is a labelled lift, not a canonical splitting of
unlabelled endpoint data.  It commutes with changing both harmonics by the
same common-harmonic deck translation, but `G_full` retains that deck
coordinate rather than quotienting it away.

## 6. Choice ledger

There are three distinct kinds of choice, which must not be conflated.

1. Adding an `L_13`-vector to a chosen lift changes no permutation of the
   quotient.
2. Swapping the labelled axes `a` and `b` gives the conjugate construction
   with the other target coordinate active.
3. A genuinely different Bezout section can change the displayed action.
   For example `e_0-e_a` still pairs to one with `W`, and its coordinate
   chart is

   ```text
   C(s,t,delta)=(s-delta,t,delta).                    (24)
   ```

   Writing the input in the alternative chart as `(a,t,delta)`, direct
   transport gives

   ```text
   C rho(x,y,c) C^(-1)(a,t,delta)
    =(a+x-c+y(a+delta),
      t,
      delta+c-y(a+delta)).                            (24a)
   ```

   This is conjugate to `(8)`, but is not identical to `(8)` under the
   identity labelling of triples; the central generator at the origin
   already distinguishes them.

Thus the conjugacy class is stable under section transport, while the raw
formula is not section-free.

## 7. What survives after forgetting carrier imbalance

The projection

```text
pi:G_delta -> G,
(s,t,delta) |-> (s,t)                                (25)
```

is equivariant for the induced action

```text
(x,y,c).(s,t)=(s+x,t).                               (26)
```

Its permutation image has order `13`; the entire `y`-line and the centre are
killed.  Therefore the noncommutative Heisenberg clutch is invisible on the
old `169`-element target quotient.  In particular, the THM-2334 `A(q)` and
`B(q)` currents on that old quotient do not by themselves witness `(8)`.

This also explains why the `2,197`-element address is not a gratuitous
cardinality enlargement: the imbalance coordinate is exactly where the
commutator centre acts.

## 8. Exact boundary and hostile tests

The construction is purely modulo thirteen.  Formula `(9)` compensates in
the raw `e_0` guard coordinate.  THM-2763 shows that discarding carrier
provenance at precisely such a compensation can change atomic-factor and
all-`91`-unit masks.  No factor-labelled extended current is constructed
here.

Nor may `(14)` be read as the physical length-`13` odometer.  On
`F_13^2`, the carry-suppressed low shift `Y` has `Y^13=I`, whereas the
physical odometer has thirteen carry-wall states and satisfies `O^13=Z` in
the two-digit chart, as isolated by THM-2788.  An integral or physical lift
therefore needs an additional carry sidecar.

THM-3234 is a separate construction.  It concerns a Singer completion of a
cyclic `168`-owner set on a `169`-point affine plane.  The present theorem
concerns the fixed THM-2763 exact-address quotients of sizes `2,197` and
`28,561`; it uses no Singer cycle, supplies no owner/root map, and proves no
compatibility with that compactification.

More explicitly, this theorem does **not** construct an extended physical
factor current, preserve the all-`91`-unit predicate, prove a nonzero endpoint
determinant, survive target difference or Radon aggregation, exclude a
scalar row, or prove `LRC(14)`.  Its exact gain is a lawful, faithful finite
Heisenberg permutation action on the already-proved carrier-imbalance
address, together with the complete choice and information-loss ledger.

## 9. Exact companion

Run from the repository root:

```bash
python3 04-computation/lrc14_exact_address_heisenberg_clutch_thm3240.py
python3 -O 04-computation/lrc14_exact_address_heisenberg_clutch_thm3240.py
```

The assertion-independent companion uses exact arithmetic only.  It pins
THM-2350, THM-2763, THM-2779, and THM-3228 by LF-normalized SHA-256; checks
all `13^6=4,826,809` group-law coefficient pairs; enumerates the faithful
permutation image, stabilizer, centre, and both orbit decompositions; audits
`371,293` raw representative transitions in the explicit THM-2350 normal
form; and checks all `4,826,809` section transports.  It also verifies the
order-`13` carrier-forgetting image and the odometer hostile.  Ordinary and
optimized runs must byte-match the pinned transcript.

QED.
