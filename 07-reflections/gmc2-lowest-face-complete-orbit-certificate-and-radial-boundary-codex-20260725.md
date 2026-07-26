# GMC2 lowest-face complete-orbit certificate and the radial boundary

**Status:** FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.

The canonical truth is already stronger than this note: THM-2022,
`THM-2022-gmc2-frobenius-lowest-balanced-face`, proves NC2 and GMC(2)
for arbitrary finite support and arbitrary complex coefficients, including
arbitrary radial polynomials. This note gives a new nonnegative diagnostic
for its one-variable face seed and explains exactly why the diagnostic
cannot replace the Frobenius proof.

## The exact face-seed residual

Let

```text
f(u)=sum_(q=-M)^N c_q u^q,

c_(-M)c_N!=0,                   M,N>=1,             (1)
```

and put `g=f^m`. For an integer `K>=1`, let

```text
R_a g(u)=g(zeta_K^a u),        a in C_K.            (2)
```

Normalized `L2(T)` makes the Laurent monomials orthonormal. THM-2375's
finite-unitary complete-orbit identity gives

```text
V_K(g)
 :=||g||_2^2
   -1/(2K) sum_(a in C_K)||R_a g-g||_2^2

 =sum_(n congruent 0 mod K)|g_n|^2.                 (3)
```

The exact zero-residue alias-freedom criterion is

```text
supp(g) intersect K Z is a subset of {0}.
```

For the width data in (1), its sharp monotone sufficient condition is

```text
K>m max(M,N).                                        (4)
```

Under (4), the support of `g` lies in `[-mM,mN]`, and zero is its only exponent
divisible by `K`. Hence

```text
V_K(f^m)=|CT(f^m)|^2.                                (5)
```

Thus the scalar constant-term seed has a finite, nonnegative
complete-orbit certificate once the power `m` is fixed.

This is weaker than injectivity of the full exponent-labelled bank. To
prevent any two exponents in `[-mM,mN]` from sharing a residue, the sharp
width-only condition is

```text
K>m(M+N).
```

The constant-term certificate needs only the smaller zero-residue
condition (4).

## One uniform finite certificate

THM-2111,
`THM-2111-effective-compound-root-bound-for-one-variable-constant-terms`,
sets

```text
C=binom(M+N,min(M,N))                                (6)
```

and proves

```text
CT(f^m)!=0
```

for some `1<=m<=C`. Choose one cyclic order satisfying

```text
K>C max(M,N).                                        (7)
```

Then (4) holds simultaneously for every possible witness power
`m<=C`, and (5) gives

```text
there exists 1<=m<=C with V_K(f^m)>0.                (8)
```

Equation (8) is a phase-safe finite certificate for the nonzero
whole-face seed used by THM-2022. It applies after lowest-face extraction
with arbitrary original radial coefficients because the extracted
one-variable Laurent polynomial retains the complete labelled face
coefficients: THM-2022 proves that equal charge and equal face height
cannot hide a projection collision. If that face already contains a
neutral charge, its nonzero coefficient is detected at `m=1`; THM-2111
supplies the two-sided case.

It does not prove that one of the residuals is positive. That existence
is exactly the imported content of THM-2111.

## The alias threshold is load-bearing

The strict inequality in (4) is uniformly sharp even within the
two-sided class (1). At `m=1`, take

```text
f(u)=u^(-1)+u^K.
```

Then `CT(f)=0`, while the exponent `K` lies in the trivial residue, so

```text
V_K(f)=1.                                            (9)
```

The nonzero exponent `K` has aliased into the trivial character. A
complete orbit is not enough without an exponent-range certificate.

## Why the bank must be applied after taking the power

Put

```text
f_+(u)=u+u^(-1)+u^2+u^(-2),

f_-(u)=u+u^(-1)+u^2-u^(-2).                         (10)
```

Their labelled one-body Fourier norm banks are identical. Nevertheless,

```text
CT(f_+^2)=4,

CT(f_-^2)=0.                                         (11)
```

The missing sign appears only in a cross-term after multiplication.
Therefore the norm bank of `f` cannot choose the face seed. Formula (5)
is useful precisely because it measures `f^m`, after the relevant phases
have interacted.

## Arbitrary radial coefficients: the exact no-go

For circular Gaussian measure, use the orthogonal basis

```text
phi_(q,n)
 =Z^(q_+) W^(q_-) L_n^(|q|)(s),

s=ZW,

||phi_(q,n)||^2=(n+|q|)!/n!.                       (12)
```

For `q>=1`,

```text
L_0^(q)(s)=1,

L_1^(q)(s)=q+1-s.
```

Consider

```text
P_-=Z^q s
    =Z^q((q+1)L_0^(q)-L_1^(q)),

P_+=Z^q(2(q+1)-s)
    =Z^q((q+1)L_0^(q)+L_1^(q)).                    (13)
```

The two polynomials have:

- the same nonzero charge `q`;
- the same labelled Laguerre-mode magnitudes;
- the same charge/Laguerre norm, exactly `(q+2)!`; and
- the same scalar moment sequence, identically zero.

But their radial orders are

```text
ord_s(P_-)=1,

ord_s(P_+)=0.                                        (14)
```

Thus even the full labelled charge-and-Laguerre norm bank cannot select
the lowest radial face. For a Gamma radial law of shape `alpha`, replace
`q+1` by `q+alpha`; the same obstruction remains.

## Structural diagnosis

The two useful decompositions are transverse:

```text
Laguerre/chaos grading:
  orthogonal and unitary,
  but not multiplicative;

s-adic/Newton filtration:
  multiplicative and Frobenius-stable,
  but not orthogonal.                                (15)
```

Their change of basis is triangular Laguerre/Pascal. The sign in (13)
is a phase cancellation inside that triangular transform. An independent
unitary phase on the Laguerre index is not an algebra automorphism:
products mix radial degrees through Laguerre linearization. It therefore
cannot be pushed from `P` to all powers `P^m`.

THM-2022 uses the second structure because the proof must propagate one
face through a large prime power. THM-2375 uses the orthogonal/unitary
charge grading; even the hypothetical refinement by labelled Laguerre
norms stays on the nonmultiplicative side. The complete-orbit certificate
(5) is a useful diagnostic **after** multiplication, but the hostiles
show that no map from the unpolarized one-body norm bank can supply the
powered constant-term data or the Frobenius face selector. Enriched
polarized reference data is a genuinely different sidecar.

The radial pair shows that degree one is the sharp first boundary for
recovering the `s`-adic/Newton order even from all unpolarized labelled
charge-and-Laguerre magnitudes. This is a boundary of that face selector,
not of NC2 itself; powered constant-term phase loss already occurs in the
constant-radial pair (10).

## Exact reproduction

Run

```bash
python3 04-computation/gmc2_lowest_face_complete_orbit_certificate.py
python3 -O 04-computation/gmc2_lowest_face_complete_orbit_certificate.py
```

and compare both transcripts, after LF normalization, with

```text
05-knowledge/results/gmc2_lowest_face_complete_orbit_certificate.out
```

The exact companion checks `24` zero-residue-alias-free powers, `24`
real and `8` Gaussian-rational residue identities, `49` zero-residue and
`49` full-bank width invoices, the in-scope alias and phase hostiles,
and `12` radial Laguerre hostiles with exact monomial expansions,
`s`-adic orders, norms, and charge-null typing. It explicitly records
that NC2 is already proved by THM-2022 and that this note adds no new
NC2 step.

LF-normalized SHA256:

```text
script c0eeb3a11f147c9a65945d3e452cd3642f8ba72637196b66f9541b04614450b7
output f5e679e57f660d6357628e6ab6b515502a61cbba11f099ba9c6c76e7ec11e373
```

Independent hostile audit accepted the zero-residue/full-bank
distinction, the sharp thresholds and in-scope alias witness, the
THM-2111/2022 scope, every radial expansion and normalization, both
phase-loss boundaries, the LF hashes, and normal/optimized/stored replay.
