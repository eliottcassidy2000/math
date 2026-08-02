---
id: THM-3058
title: "K4 hafnian initial-face augmentation and unbounded cancellation jet"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  For every finite nonempty indexed fibre over a DVR with residue field
  F_q, the sum has the minimum input valuation exactly when the initial
  residues on the minimum face have nonzero sum.  This is also equivalent
  to nonvanishing for every lift with the fixed valuations and initial
  residues.  On the zero-initial-sum wall the same first-order data admit
  both an exact-zero lift and nonzero lifts at arbitrarily deeper order, so
  a fixed amplitude requires an unbounded filtered contraction jet.  For
  K4 this is an S4/V4-equivariant, V4-blind sidecar of the three selected
  matching monomials.  It applies only after THM-2290 endpoint-colour
  selection and pair aggregation and supplies no quartic origin or physical
  map.
source: codex-2026-08-01-k4-hafnian-initial-face-augmentation
depends_on:
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
  - THM-3049-k4-matching-monomial-tropical-root-extraction-clutch
  - THM-3046-quartic-resolvent-root-valuation-binary-ternary-clutch
related:
  - THM-3060-three-slot-physical-terminal-face-and-affine-tail-holotopy
script: 04-computation/k4_hafnian_initial_face_augmentation_thm3058.py
output: 05-knowledge/results/k4_hafnian_initial_face_augmentation_thm3058.out
script_sha256: b08f65ea08b99a43cd42cc53ad1a0d1890ff06597300d1f78339e52452eac4f8
output_sha256: ad069b02467e0cc41142586bb45619499dc0349495dd843e8cd00fbfac0f06bb
semantic_sha256: ed996d38ec1a9b6c4dc734ed5311af1e789d3a50c9b0d946f823df9dad4169a8
hash_basis: LF-normalized bytes
---

# THM-3058 -- K4 hafnian initial-face augmentation and unbounded cancellation jet

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## 1. Inheritance and exact finite-fibre theorem

The closest proved mechanism is
[THM-3049](THM-3049-k4-matching-monomial-tropical-root-extraction-clutch.md):
it retains the valuations of the three labelled `K4` matching monomials and
computes their binary/ternary clutch before contraction.  Its canonical
hostile `(1,1,-2)` shows that those valuations do not determine the sum.  The
corrected near miss is to append the three leading residues but still forget
which matching terms form the minimum face; the `lambda024/live` versus
`lambda222/zero` pair below defeats that quotient.  The least-used relevant
sidecar is the initial form of the already selected matching-fibre
contraction.  [THM-2290](THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete.md)
fixes the order: endpoint-colour selection and same-pair aggregation precede
this sidecar.

The underlying statement is not special to three terms.  Let `R` be a DVR
with fraction field `K`, uniformizer `pi`, normalized valuation

```text
v:K^* -> Z,
```

and finite residue field `k=R/(pi)=F_q`.  Use `v(0)=+infinity`.  Let `I` be a
finite nonempty index set and take nonzero amplitudes `A_i in K^*`.  Define

```text
H=sum_(i in I) A_i,
lambda_i=v(A_i),                 lambda=min_i lambda_i,
F={i in I:lambda_i=lambda},
alpha_i=bar(pi^(-lambda_i) A_i) in k^*,
sigma=sum_(i in F) alpha_i in k.                         (1)
```

For the fixed first-order datum

```text
D=((lambda_i,alpha_i))_(i in I),                         (2)
```

let `Lift(D)` be all labelled tuples `(B_i)` with those exact valuations and
initial residues.  Then the following are equivalent:

```text
(i)   sigma!=0;
(ii)  v(H)=lambda;
(iii) every (B_i) in Lift(D) has sum_i B_i !=0.          (3)
```

If `sigma=0`, the failure is sharp: `Lift(D)` contains one tuple whose sum is
exactly zero and, for every `N>=1`, another tuple with sum

```text
pi^(lambda+N).                                          (4)
```

The two tuples in `(4)` can be chosen to agree componentwise modulo relative
`pi^N`, hence through every depth `<N`.  Thus first-order data classify
**robust noncancellation across all lifts**, while the zero-`sigma` locus
contains both actual outcomes at every finite observation depth.

## 2. Associated-graded proof and lift boundary

Normalize the contraction by its minimum possible valuation:

```text
pi^(-lambda)H
 =sum_(i in F) pi^(-lambda)A_i
  +sum_(i notin F) pi^(lambda_i-lambda)pi^(-lambda_i)A_i. (5)
```

The second sum vanishes modulo `pi`; the first reduces to `sigma`.  Hence

```text
bar(pi^(-lambda)H)=sigma.                               (6)
```

If `sigma!=0`, the normalized sum is a unit and `v(H)=lambda`.  If
`sigma=0`, the normalized sum lies in `(pi)` or is zero, so
`v(H)>lambda` with the convention above.  This proves `(i) iff (ii)`, and the
same calculation applies to every tuple in `Lift(D)`.

It remains to prove that `(iii)` fails whenever `sigma=0`.  Choose a pivot
`j in F`.  For every `i!=j`, choose any lift

```text
B_i=pi^(lambda_i) u_i,             bar(u_i)=alpha_i.     (7)
```

Set

```text
B_j=-sum_(i!=j) B_i.                                    (8)
```

After division by `pi^lambda`, the residue of the right side of `(8)` is

```text
-sum_(i in F, i!=j) alpha_i=alpha_j,                    (9)
```

because `sigma=0`.  Thus `(8)` is nonzero, has valuation `lambda`, and has
the required initial residue.  It gives an exact-zero lift.  Replacing it by

```text
B_j+pi^(lambda+N)                                       (10)
```

preserves its valuation and leading residue, preserves its relative
`N`-jet, and changes the total sum from zero to `(4)`.  This proves the
remaining implications and the sharp lift boundary.

The decision bit is independent of the uniformizer.  If `pi'=u pi` for a
unit `u`, every minimum-face residue is multiplied by the same scalar
`bar(u)^(-lambda)`, so

```text
sigma'=bar(u)^(-lambda)sigma.                            (11)
```

The displayed scalar depends on the trivialization, but its vanishing does
not.

## 3. The unbounded filtered contraction jet

For a fixed amplitude tuple, define its depth-`N` contraction jet by

```text
J_N(H)=pi^(-lambda)H mod pi^N R in R/(pi^N),   N>=1.    (12)
```

The first jet is `J_1(H)=sigma`.  If `H!=0` and
`v(H)=lambda+d`, then `J_N(H)=0` for `N<=d` and
`J_(d+1)(H)!=0`.  If `H=0`, every `J_N(H)` vanishes.  Therefore

```text
H!=0 iff J_N(H)!=0 for some N,
H=0  iff J_N(H)=0 for every N.                           (13)
```

No uniform finite truncation of `(12)` decides `(13)` on the `sigma=0` wall.
For any proposed depth `N`, the two lifts `(8)` and `(10)` have identical
termwise relative jets modulo `pi^N`, and both contractions have zero image
in `R/(pi^N)`, yet one sum is zero and the other is `pi^(lambda+N)`.  Thus
exact fixed-amplitude nonvanishing requires the **unbounded filtered
contraction jet** when the first initial form cancels.  This is a no-go only
for bounded filtered-jet data; it is not a claim that an explicitly given
field element cannot be tested directly.

## 4. Exact finite-field cancellation count

Let `r=|F|`.  The number `Z_r(q)` of ordered words

```text
(alpha_1,...,alpha_r) in (F_q^*)^r
```

whose sum is zero is

```text
Z_r(q)=((q-1)^r+(-1)^r(q-1))/q.                         (14)
```

Indeed, `Z_1(q)=0`.  A nonzero `(r-1)`-word sum admits exactly one nonzero
last entry which makes the total zero, while a zero prefix admits none.
Therefore

```text
Z_r(q)=(q-1)^(r-1)-Z_(r-1)(q),                          (15)
```

and `(14)` follows by induction.  In particular,

```text
r=1: 0,              r=2: q-1,              r=3: (q-1)(q-2). (16)
```

These counts measure the additive residue wall.  They are not THM-3049's
multiplicative square/cube divisibility clutches, even though characteristics
two and three give distinctive small-face behavior.

## 5. Selected `K4` hafnian specialization

Fix a THM-2290 endpoint-colour query and first form its aggregated selected
scalar kernel `B`.  The three labelled matching amplitudes are

```text
A_1=B_01 B_23,        A_2=B_02 B_13,        A_3=B_03 B_12,
haf(B)=A_1+A_2+A_3.                                      (17)
```

If all three amplitudes are nonzero, apply `(1)--(13)` directly.  If some
are zero, apply the theorem to the nonzero indexed subfibre; the all-zero
case is already trivial.  A selected pair entry which is itself a sum of
parallel eligible edges must be aggregated **before** its valuation and
initial residue are taken.  Moving the construction before endpoint-colour
selection or pair aggregation repeats THM-2290's selector/phase-erasure
failure.

Vertex relabelling permutes the three matchings through

```text
S4/V4 ~= S3.                                             (18)
```

It carries `F` to the permuted minimum face and leaves the contracted sum
`sigma` unchanged.  Hence the initial-face augmentation is `S4/V4`-
equivariant, and the `V4` kernel acts trivially: it is **V4-blind**.  It does
not choose a labelled vertex, a quartic section, or an origin.  Under a
vertex-unit gauge `B_ij -> u_i u_j B_ij`, every term in `(17)` is multiplied
by `U=product_i u_i`, so `sigma` is multiplied by `bar(U)` and its
nonvanishing remains invariant.

If the full valuation vector `(lambda_1,lambda_2,lambda_3)` is retained, the
one new robust decision scalar is `sigma`.  If only THM-3049's 24-valued
clutch is retained, it does not determine `F`; the complete robust
augmentation is

```text
minimum matching face F + its initial-form sum sigma.     (19)
```

## 6. Sharp hostile controls

### 6.1 Unbounded depth

Over `Q_5`, for every `N>=1`,

```text
(1,1,-2)               has sum 0,
(1+5^N,1,-2)           has sum 5^N.                     (20)
```

The two triples have the same valuations, the same leading residues
`(1,1,3)`, and identical component jets modulo `5^N`.  This is the canonical
zero-versus-deep-live boundary in `(12)--(13)`.

### 6.2 Same-clutch face loss

The two literal `K4` matching triples

```text
(1,5^2,-2*5^4):        lambda=(0,2,4),  F={1},   sigma=1, sum!=0;
(5^2,5^2,-2*5^2):      lambda=(2,2,2),  F={1,2,3}, sigma=0, sum=0       (21)
```

have the same THM-3049 clutch

```text
(lambda mod 2, sum lambda mod 3)=((0,0,0),0)             (22)
```

and the same labelled leading residues `(1,1,3)`.  Put the displayed three
values on edges `01,02,03` and put `1` on edges `12,13,23`; then the three
matching products are exactly the displayed triples.  Thus neither the
24-valued clutch nor the clutch plus the detached residue word recovers the
minimum face or robust noncancellation.

### 6.3 THM-3046 oriented Pluecker wall

For THM-3046's quartic root-difference matching products,

```text
P_1-P_2+P_3=0.                                           (23)
```

Apply the theorem to the signed triple `(P_1,-P_2,P_3)`.  Its exact sum is
zero, so its minimum is attained at least twice and its initial-face sum is
forced to satisfy `sigma=0`.  This is an exact hostile **on** the cancellation
wall, not a nonvanishing mechanism.  An arbitrary THM-2290 selected `K4`
kernel does not inherit `(23)`.

Equation `(23)` is only the labelled quartic-root specialization already
proved in THM-3046.  Nothing here constructs a quartic origin, transports an
arbitrary matching fibre to a quartic, selects a physical owner, or supplies
an affine/Keller/LRC map.

### 6.4 Separate incoming positive physical specialization

**CONDITIONAL comparison; message-level payload, not a proved dependency.**
[MSG-2992](../../agents/broadcast/MSG-2992-from-kind-pasteur-2026-08-01-thm3060-physical-k3-face--exa.md)
starts from THM-3060's proved-candidate rank-two cancellation at four slots
and proposes a distinct first transverse augmentation for fixed lows `a<b`,
`S=x_0+x_1`, and moving variable `z`.  After the rank-two face cancels, the
proposed surviving packet retains a fixed-only quadratic `Q` together with

```text
U_(3,C)L_3^3,                    U_(4,C)L_4^4.           (24)
```

At the common projective point `P=[1:-1:0]`, the payload gives

```text
Q(P)=G_2=L((f_(n+a)-f_(n+b))^2)/L(f_n^2)>0,
det(L_3,L_4)^2=(4^h-3^h)^2,                             (25)
```

and hence proposes the positive first-augmentation asymptotic

```text
R_C ~ G_2^12(4^h-3^h)^24 U_(3,C)^8 U_(4,C)^6
    =K 12^(24C) C^(-17),                    K>0.         (26)
```

If canonicalized, `(26)` would be a concrete physical instance of the
abstract pattern here: the first face has `sigma=0`, while a later filtered
contraction class is positive.  It does not overwrite or prove the general
DVR theorem, and THM-3058 does not supply its physical map.  In particular,
one must not feed the naive termwise first face into `(1)`: the message
records **305 equal-base patterns which cancel**.  Those equal-base
contributions must first be contracted into the correct face amplitudes,
exactly as THM-2290 requires pair aggregation before the K4 application.

## 7. Exact verification and provenance

The theorem was recovered from the scratch packet committed as
`631ae23d4c76c7cbba40a6aef9a1b3c02d878c8d`:

```text
.scratch/k4_hafnian_initial_face_sidecar.md
.scratch/k4_hafnian_initial_face_sidecar.py
.scratch/k4_hafnian_initial_face_sidecar.out             (27)
```

The canonical companion strengthens that three-term packet with the general
count `(14)` and the arbitrary finite-fibre lift construction.  Run

```text
python3 04-computation/k4_hafnian_initial_face_augmentation_thm3058.py
python3 -O 04-computation/k4_hafnian_initial_face_augmentation_thm3058.py
```

Both executions LF-byte-match the stored transcript.  The script uses an
explicit `require` function and no truth-bearing assertions.  Its q=5 bank
exhausts:

```text
87,380 nonzero residue words of lengths 1,...,8 for (14);
39,216 zero-sigma lift signatures of lengths 1,...,5;
1,906,624 nonzero triples modulo 5^3, in 1,728 signatures;
5,040 ordered distinct-root Pluecker packets;
10,368 S4/V4 face-action checks;
1,048,576 vertex-unit gauge checks.                       (28)
```

It also retains the `lambda024/live` versus `lambda222/zero` hostile and
twelve successively deeper versions of `(20)`.  The exhaustive computation
is evidence for the general proof, not its source.

```text
PROVED CANDIDATE HERE: general finite-fibre initial-face equivalence;
                       exact zero/live lift boundary;
                       unbounded filtered-jet no-go;
                       finite-field zero-sum count;
                       selected K4 specialization and S4/V4 equivariance.

VERIFIED-EXACT HERE:   the q=5 universes and hostile controls in (28).

NOT PROVED HERE:       promotion through independent hostile audit;
                       a selector or pair aggregation before THM-2290;
                       recovery from the 24-valued clutch without F;
                       the conditional physical augmentation (24)--(26);
                       a quartic origin, physical owner/map, Keller, LRC,
                       or tournament consequence.
```
