---
id: THM-3084
title: "Cyclic root pointer, Kummer eigenline, and positive-cone separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the
  equal-weight norm on p cyclic coordinates, a physical root-permutation
  orbit lies in one diagonal Kummer fibre exactly when the coordinate vector
  is a single Fourier-character eigenline.  A nonconstant positive profile
  never does, while a same-atom charged root pointer identifies the two
  torsors exactly.  THM-2537 therefore has no internal predecessor-side
  mu_13 ambiguity after fixed-cylinder character projection; cross-packet
  semantic arrival remains open.  The p=2 analogue needs two independent
  same-atom contrasts to reconstruct a quartic V4 plane.
source: root-cyclic-pointer-kummer-2026-08-01
audit: >
  Two independent hostile audits rederived the equal-weight fibre criterion,
  cyclic eigenline and positive-cone separation, then checked the THM-2537
  same-cylinder pointer, point-before-sum boundary, residual semantic-arrival
  debt, and common-atom V4 comparison.  The first audit caught and the
  repaired theorem records the exact i_-a/P^-h orientation.  Immutable normal
  and optimized replay matched the 521-byte stored transcript and declared LF
  hashes; documentation passed.
depends_on:
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-3077-pointed-norm-relative-line-lift-and-relation-carry-obstruction
related:
  - THM-2548-seven-step-c91-transfer-and-full-norm-separation
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
  - THM-2613-canonical-root-diagonal-opposite-shift-section
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2962-planar-suspension-resolvent-and-v4-kummer-descent
  - THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan
script: 04-computation/cyclic_root_pointer_kummer_eigenline_thm3084.py
output: 05-knowledge/results/cyclic_root_pointer_kummer_eigenline_thm3084.out
script_sha256: c5e06c7bd2a5ec5a84e81be7d75a27526d68b7b06257a9cdaf534f687e5d8c71
output_sha256: e666ec5486061e42abc640ff087cac0f2805abad25c6942009c2e3c700bd74ed
hash_basis: LF-normalized bytes
---

# THM-3084 -- cyclic root pointer and Kummer eigenline

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3077 identifies the information lost by a weighted norm: its
norm-augmentation fibre is a diagonal Kummer torsor.  LRC already carries a
different-looking torsor, namely the physical choice of a root in a cyclic
root chart.  These torsors are not interchangeable merely because both have
thirteen elements.  Their exact meeting locus is a Fourier eigenline.

## 1. Cyclic shift versus the diagonal lift

Let `p` be prime, let `K` be a field with `char(K)` not dividing `p`, and
assume `K` contains a primitive root `zeta_p`.  Put

```text
T=(G_m)^p,
N(x)=product_(j in F_p)x_j,
A=T/Delta G_m,
F=(N,pi):T -> G_m times A.                                (1)
```

This is THM-3077 with all weights equal to one.  Its degree is `p` and every
geometric fibre is

```text
F^(-1)(F(x))=Delta mu_p x.                                (2)
```

Let cyclic coordinate shift act by

```text
(P x)_j=x_(j-1).                                          (3)
```

The norm in `(1)` is shift-invariant.  The following are equivalent:

```text
F(Px)=F(x),
pi(Px)=pi(x),
Px=lambda x for some lambda in K^*,
lambda in mu_p and x_j=c lambda^(-j) for some c in K^*.   (4)
```

The second and third conditions are equivalent by the definition of the
diagonal quotient.  Iterating `P^p=1` gives `lambda^p=1`.  Finally
`x_(j-1)=lambda x_j` recursively gives the last formula, and direct
substitution proves the converse.

If `lambda!=1`, primality of `p` makes it primitive.  Then

```text
P^h x=lambda^h x,                                         (5)
```

so `h->lambda^h` identifies the full physical `C_p` root orbit with the
full diagonal `mu_p` fibre `(2)`.  If `lambda=1`, the coordinate vector is
constant and its physical orbit is one point, while `(2)` still has `p`
points.  The collapsed root label therefore does not distinguish the Kummer
branches or provide a root-to-Kummer torsor conversion (although an externally
specified point `x` is of course already a pointed lift).

## 2. Positive-cone separation

Specialize to positive real coordinates.  If `(4)` holds, then
`lambda=x_(j-1)/x_j` is positive.  The only positive real `p`th root of unity
is one.  Hence

```text
x in R_(>0)^p and F(Px)=F(x)
 iff x_0=...=x_(p-1).                                     (6)
```

Thus no nonconstant positive real or positive rational root profile lets its
cyclic root pointer stand in for the diagonal Kummer pointer.  This is a
structural separation, not a shortage of estimates.  The literal positive
hostile `x=(1,2,...,p)` preserves `N` under shift but immediately leaves its
augmentation fibre.

At `p=13`, exact exhaustion of `{1,2}^13` gives precisely the two constant
aligned profiles and `8,190` failures.  The bound is sharp: both constant
profiles attain the collapsed case of `(6)`.

## 3. The charged eigenline survivor

For `a in F_p^*`, define

```text
v_a=(zeta_p^(a j))_(j in F_p),
i_a(c)=c v_a.                                             (7)
```

Then

```text
P^h i_a(c)=zeta_p^(-a h)i_a(c),
N(i_a(c))=c^p                    for odd p.                (8)
```

The norm identity follows because
`sum_(j=0)^(p-1)j=p(p-1)/2` is divisible by `p` for odd `p`.  Since `a!=0`,
the character is primitive and

```text
h -> zeta_p^(-a h)                                        (9)
```

is the canonical conversion from the physical root torsor to the diagonal
Kummer torsor on this line.  A physical origin on the **same nonzero charged
coefficient** therefore removes the entire internal `mu_p` ambiguity.

This survivor does not contradict `(6)`: Fourier projection is complex and
signed.  It leaves the positive Boolean cone before entering `(7)`.

## 4. Concrete LRC predecessor corollary

THM-2537 supplies the required same-atom pointer on its predecessor side.
On every fixed positive selector cylinder, equations `(48)--(51)` of that
theorem give

```text
S^sel_tau(iota_r(z))
 =g(z) 1_(r=s_tau),                                       (10)
```

with one occupied root `s_tau`, positive rational weight `g`, the THM-2349
terminal word and late owner, and a nonzero coefficient in every one of the
twelve charged root colours.

Fix `a!=0` before summing cylinders.  The normalized charged coefficient of
the singleton in `(10)` is

```text
c_a=(g/13)zeta_13^(-a s_tau)!=0.                          (11)
```

Across the thirteen physical root translates it becomes

```text
(c_a zeta_13^(-a h))_(h in F_13),                         (12)
```

which is exactly the eigenline `i_(-a)(c_a)` in the indexing of `(7)`.  With
the convention `(P x)_j=x_(j-1)`, physical translation by `+h` acts on this
orbit vector as `P^(-h)` and multiplies it by `zeta_13^(-a h)`.  Thus the
inversion is only an orientation convention, and `(9)` still gives the exact
torsor isomorphism.  The selected physical root `s_tau` fixes its phase and
hence points the `mu_13` lift.  Therefore:

> On each fixed THM-2537 selector cylinder, the selected predecessor-root
> current has no internal diagonal-Kummer ambiguity after nonzero character
> projection.

This is unconditional on the proved predecessor packet.  It does not prove
the chronological target in THM-2537 equation `(56)`.  The selected head is
still an empty predecessor, not a semantic arrival.  THM-2548 preserves
charged `C_91` modes and THM-2555/2613 provide later root/ancestry pointers,
but canon does not put those later pointers on the same semantic charged
coefficient.  That cross-object transport, not the internal root/Kummer
conversion, is the remaining pointer debt.

Pointing per fixed cylinder is load-bearing.  Integrating cylinders with
different roots first may cancel `(11)` to zero or erase the physical origin;
zero leaves the torus `(1)`.

## 5. The quartic `p=2` boundary

For `p=2`, the nontrivial eigenline is

```text
(c,-c),             P(c,-c)=-(c,-c).                       (13)
```

Thus one common odd contrast identifies one physical sheet swap with one
diagonal `mu_2` lift.  A quartic `V_4` plane needs two independent contrasts
on the same atom.  Each single character `V_4->mu_2` has a two-point kernel,
while two independent characters give

```text
V_4 -> mu_2 times mu_2                                    (14)
```

as an isomorphism; the third character is their product.  This is exactly
the rank-two information retained by THM-2655/2962.  THM-3072's finite
tomography shows why two separate views do not suffice: their origins must
be aligned on a common atom.  No new quartic physical intertwiner or JC(2)
closure is asserted here.

## 6. Sharp typing boundaries

The hypotheses in `(1)--(4)` cannot be dropped silently.

1. If `c=0`, `(7)` exits `G_m`; use the zero-carrier normal cones of
   THM-2983/2985 instead of a Kummer fibre.
2. For composite cyclic order, a nonprimitive character sees only its proper
   quotient; the orbit size is the order of `lambda`.
3. For unequal weights, cyclic shift need not preserve the weighted norm, so
   the physical orbit may leave the THM-3077 fibre before augmentation is
   tested.
4. A pointer on a separate packet carries an independent diagonal phase and
   does not select the lift of the target packet.
5. The singleton Boolean vector itself has zero coordinates and is not in
   `T`; it is its nonzero charged **orbit vector** `(12)` that lies on the
   Kummer eigenline.

## 7. Exact evidence and scope

Over `F_53` with primitive `zeta_13=16`, the exact companion verifies:

1. all `156` pairs `(a,h)` in `(8)`, including the equal-weight norm;
2. all `8,192` positive `{1,2}^13` profiles and the exact `2/8190` split;
3. one hundred deterministic non-eigenline torus hostiles which preserve the
   norm but leave the augmentation fibre;
4. all `156` fixed-root/nonzero-character singleton phase identities; and
5. the four-element `V_4` reconstruction and both two-point single-character
   fibres.

Run

```text
python 04-computation/cyclic_root_pointer_kummer_eigenline_thm3084.py
python -O 04-computation/cyclic_root_pointer_kummer_eigenline_thm3084.py
```

Both modes must equal the stored transcript after LF normalization.

This theorem identifies and, on THM-2537's predecessor cylinders, supplies
the exact common pointer needed by THM-3077.  It does not transport that
pointer to a distinct semantic packet, create a positive charged eigenline,
prove a target arrival, decrement an LRC row, align two separate JC views,
or prove LRC(14), JC(2), GMC(2), NC2, or DC(2).

**QED.**
