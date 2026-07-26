---
id: THM-2439
title: "Cyclic marker replica degree and homometric Gram boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The
  translation-equivariant lexicographic marker of a nonconstant
  Boolean C_13 word is an exact squarefree polynomial in its
  translated truth bits. On the full cardinality class 1<=|S|<=10
  allowed by THM-2435's current information, its intrinsic
  independent-replica degree is exactly 10, with 2753 nonzero
  coefficients for one marked root. On all nonempty proper masks the
  intrinsic degree is 12, with 2816 nonzero coefficients in the
  zero-extension. Thus a same-gauge replica-current algebra through
  degree 10 is exactly sufficient for the THM-2435 marker, but current
  THM-2334/2410 probes do not provide that algebra. Two non-dihedral
  ambient (13,4,1) difference sets have identical complete cyclic
  self-Gram/Fourier-power banks but distinct marked words, and support
  four is the first such hostile. No canonical replica-to-Abel
  intertwiner, endpoint reference, terminal current, row exclusion,
  or proof of LRC(14) follows.
source: codex-2026-07-26-cyclic-marker-replica-boundary
depends_on:
  - THM-2435-top-blocker-essential-parent-and-punctured-stalk-carrier
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2410-full-coordinate-projector-local-gram-and-integrated-phase-boundary
  - THM-2421-all-clock-septimal-ancestry-endpoint-event-detector
  - THM-2424-coprime-common-root-crt-and-unit-residue-spectrum
script: 04-computation/lrc14_cyclic_marker_replica_boundary_thm2439.py
output: 05-knowledge/results/lrc14_cyclic_marker_replica_boundary_thm2439.out
script_sha256: adc6095961b8eb51b9c8cf55f60fd8d909769369c09ad5531724704cba97ac18
output_sha256: bb34248ab18e88909f2310f541be9ad451c59f14809d01de2d56dd76ec98ddc7
hash_basis: working-tree bytes (LF)
---

# THM-2439 -- the ambient cyclic marker has degree ten; energy does not mark

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2435 produces a translation-equivariant marked root with a flat
quotient spectrum. Its stopping boundary can now be made exact:

```text
complete cyclic self-Gram / Fourier power
  != marked root;

same-gauge squarefree replica algebra through degree ten
  => marked root exactly.                                      (1)
```

The distinction in (1) is not cosmetic. The thirteen translated root
truth values are coordinates of one cyclic profile, not thirteen
independent Boolean interventions.

## 1. The prime-cyclic lexicographic marker

Let

```text
C=Z/13Z
```

and let `S` be a nonempty proper subset. For `a in C`, define

```text
w_(S,a)=(1_S(a+r))_(r=0)^12.                                   (2)
```

The thirteen rotations are distinct. A nontrivial stabilizer would,
by primality, make `S` empty or all of `C`. With `1>0`, let

```text
sigma(S)=the unique a maximizing w_(S,a) lexicographically.       (3)
```

The maximal word starts with `1`, so `sigma(S) in S`, and

```text
sigma(S+t)=sigma(S)+t.                                          (4)
```

Put

```text
M_a(S)=1_(sigma(S)=a).                                          (5)
```

For the coefficient counts below only, use the explicit extension

```text
M_a(empty)=M_a(C)=0.                                            (6)
```

The intrinsic degree conclusions on nonempty physical domains will
not depend on this convenient extension.

## 2. Exact Boolean replica degree

Write `x_r=1_S(r)`. Every Boolean function has the unique squarefree
multilinear expansion

```text
M_a(x)
 =sum_(T subset C) c_(a,T) product_(r in T)x_r,                  (7)

c_(a,T)
 =sum_(U subset T)(-1)^(|T|-|U|) M_a(U).                        (8)
```

On THM-2435's inherited root masks,

```text
1<=|S|<=10.                                                     (9)
```

Every monomial of degree greater than ten vanishes there. Hence

```text
M_a(S)
 =sum_(T subset S, |T|<=10)c_(a,T)                              (10)
```

for all masks in (9). For `a=0`, the exact expansion through degree
ten has

```text
2753 nonzero coefficients,

max{|T|:c_(0,T)!=0, |T|<=10}=10.                                (11)
```

The degree is intrinsic to the **full cardinality class** in (9), not
an artefact of (6). Let

```text
T_1={0,1,...,9},

T_2={0,2,3,4,5,6,7,8,9,10}.                                    (12)
```

Exact Möbius inversion gives

```text
c_(0,T_1)=3,                    c_(0,T_2)=0.                     (13)
```

If a polynomial of degree at most nine agreed with `M_0` on every
mask in (9), give it arbitrary value `alpha` at the empty mask. Its
degree-ten Möbius coefficient on each `T` would be `c_(0,T)+alpha`
and would have to vanish. Equations (13) would force both
`alpha=-3` and `alpha=0`, a contradiction. Thus:

```text
the exact restricted replica degree is 10.                      (14)
```

This is the sharp universal invoice under the cardinality information
(9). THM-2435 does not assert that every mask in (9) is physically
realized; its smaller actual mask family could have lower degree.
Nor is (14) a lower bound for a different operation such as a lawful
cyclic shift or Abel probe. It is specifically the independent
squarefree-replica complexity of the coordinate indicator `M_a`.

On all nonempty proper masks, the zero-extension (6) has

```text
2816 nonzero coefficients,

exact degree 12,                    c_(0,C)=0.                   (15)
```

For two degree-twelve supports,

```text
c_(0,C minus {6})=-10,

c_(0,C minus {10})=2.                                          (16)
```

Changing the empty-mask value shifts every degree-twelve coefficient
by the same number; changing the full-mask value affects only degree
thirteen. Hence (16) rules out every degree-at-most-eleven extension.
The exact intrinsic degree on all nonempty proper masks is therefore
twelve.

Translation rotates (7), so the same degree and coefficient counts
hold for every marked root `a`.

## 3. The exact conditional current bridge

Let `Y` be a fixed parent chamber and let `f` be a Boolean root-profile
packet in one common gauge on `Y`. Write

```text
x_r=tau_r f
```

for its thirteen root translates, evaluated pointwise, and assume

```text
1<=sum_r x_r<=10                         a.e. on Y.              (17)
```

Define `M_a(f)(y)` fibrewise by (2)--(5). Let `J` be any linear
current functional on the span of the resulting same-gauge packets.
If the lawful packet algebra contains every chamber-restricted
squarefree replica product

```text
1_Y product_(r in T) tau_r f,              |T|<=10,              (18)
```

where the empty product is `1_Y`, then (10) gives the exact marked
current

```text
J(M_a(f))
 =sum_(|T|<=10)c_(a,T)
   J(1_Y product_(r in T)tau_r f).                              (19)
```

Thus a degree-ten independent-replica service is sufficient, and (14)
says no universal degree-nine service of this same algebraic type can
be.

This is a conditional algebra theorem, not a claim that (18) already
exists in LRC canon. THM-2334 and THM-2410 translate the nine base
factors inside one present packet; they do not supply products of ten
translated copies of the **derived whole packet**. Constructing those
copies is precisely the missing replica-to-current intertwiner.

## 4. Why complete-subcube tomography does not supply the replicas

THM-2383 starts with a labelled Boolean intervention cube and recovers
its Walsh components after polarized spanning references and two
complex quadratures are supplied. Here the thirteen bits in (2) are a
cyclic translation orbit in the **codomain** of one packet, not
independent toggles in the domain of an intervention cube.

Applying THM-2383 directly would therefore require a genuinely
labelled intervention cube whose correlated monomial bank realizes:

```text
the replica products (18),
a mean anchor,
supportwise spanning oriented references,
and the imaginary quadrature.                                 (19a)
```

Once (19a) exists, (19) already reconstructs the marker. Without it,
the application is a type error.

THM-2421 supplies a useful analogy: its complete signed endpoint-event
word retains its own ancestry chambers while scalar Dirichlet energy
does not. It does **not** construct the present `C_13` marker.
THM-2424 transfers the Fourier energy of an already constructed marker
to physical residue classes; it does not construct the marker current
or its endpoint phase.

## 5. A sharp homometric energy hostile

In `C`, put

```text
A={0,1,3,9},                    B={1,2,5,7}.                     (20)
```

Both are cyclic `(13,4,1)` difference sets:

```text
|A intersection (A+t)|
 =|B intersection (B+t)|
 =4                              if t=0,
 =1                              if t!=0.                       (21)
```

Consequently their complete cyclic self-Gram distances agree:

```text
|S triangle (S+t)|=0,6,6,...,6.                                (22)
```

With normalized `1/13` Fourier transform, their complete power banks
also agree:

```text
|Shat(0)|^2=16/169,

|Shat(k)|^2=3/169                 for every k!=0.                (23)
```

Nevertheless,

```text
sigma(A)=0,                       sigma(B)=1,                    (24)

w_(A,sigma(A))=1101000001000,

w_(B,sigma(B))=1100101000000.                                  (25)
```

The two words are not related by a cyclic translation or reflection;
their cyclic gap words are `(1,2,6,4)` and `(1,3,2,7)`. Thus even the
complete labelled cyclic self-Gram/power bank cannot recover the
marked profile.

Full complex Fourier phase, or a polarized cross-Gram bank with a
spanning oriented reference, does distinguish the masks. The hostile
targets energy-only self-data, not those stronger measurements.

An exact exhaustive scan proves that support four is the first such
non-dihedral homometric collision: none exists at support one, two,
or three, and at support four exactly one autocorrelation signature
contains two dihedral orbits.

The displayed pair is affine-unit equivalent, under exactly
`x -> ux+7` for `u=7,8,11`. Hence the hostile is scoped to the
translation-oriented
marker gauge and is not claimed to survive quotienting by the full
affine group. That caveat does not weaken (22)--(25): the labelled
shift banks themselves are identical.

The displayed masks are ambient `C_13` hostiles. They are not claimed
to occur in THM-2435's physically realized root-mask family.

## 6. Exact companion

Run

```text
python3 04-computation/lrc14_cyclic_marker_replica_boundary_thm2439.py
python3 -O 04-computation/lrc14_cyclic_marker_replica_boundary_thm2439.py
```

The dependency-free companion:

- checks all `8099` masks of sizes one through ten;
- computes the `2753` restricted and `2816` universal nonzero
  Möbius coefficients, their degreewise census, and the
  degree-ten/twelve witnesses;
- checks all `106470` marker-equivariance instances;
- verifies (21)--(23) exactly through integer autocorrelation and
  cyclotomic reduction, without floating point; and
- exhausts supports one through four to locate the first
  non-dihedral homometric collision and verifies the affine-unit
  caveat.

Normal and optimized runs must reproduce

```text
05-knowledge/results/lrc14_cyclic_marker_replica_boundary_thm2439.out
```

byte-for-byte.

## 7. Consequence and scope

Under only THM-2435's present cardinality information, complete cyclic
self-energy is not a universal marker service. One exact sufficient
service is:

```text
parent/chamber-resolved replica natural transformation through
degree ten

  + oriented complex endpoint reference

  +, at positive depth, a physical ancestry section.            (26)
```

Degree ten is sharp on the ambient cardinality class but may not be
minimal on the physically realized THM-2435 mask family. THM-2439
does not construct the service, identify a canonical Abel current,
preserve a terminal word, remove a valuation shape or scalar row, or
prove LRC(14).

## 8. Independent audit

Three independent audits reproduced the exact Möbius degrees, counts,
witnesses, equivariance bank, homometric census, affine-unit caveat,
normal/optimized transcript, and LF hashes. One audit caught the
missing chamber/cardinality hypothesis in the conditional current
bridge; another supplied the THM-2421 event-word hostile. Both scope
repairs are incorporated above.

QED.
