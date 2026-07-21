---
id: THM-2041
title: "Frobenius stability of exact-period projectors"
status: >
  PROVED. In the cyclic group algebra over characteristic p, p coprime to q,
  every exponent-dilation-invariant integral twist is fixed by Frobenius.
  Consequently its twisted constant term at level p*m is the p-th power of
  the level-m constant term. Primitive-support and Ramanujan exact-period
  kernels satisfy the invariance because multiplication by p permutes the
  units modulo q. Nonvanishing therefore propagates through every p-power
  amplification. This preserves a whole tied exact-period layer, but does not
  produce the initial nonzero layer required by LRC(14).
source: codex-2026-07-21-NC2-transfer
related:
  - THM-2022
  - THM-1820
  - THM-1830
  - HYP-8800
  - HYP-2978
  - HYP-2979
  - HYP-3036
script: 04-computation/lrc_frobenius_exact_period_projector_codex_20260721.py
output: 05-knowledge/results/lrc_frobenius_exact_period_projector_codex_20260721.out
---

# THM-2041 -- Frobenius stability of exact-period projectors

## 1. Cyclic group-algebra statement

Fix integers `q>=1` and a rational prime `p` with `p` not dividing `q`. Let

```text
A = k[u,u^(-1)]/(u^q-1)
```

over a field `k` of characteristic `p`, and let `ct_q` denote the coefficient
of the identity class `u^0`. Exponent multiplication

```text
sigma_p(u^a)=u^(p*a)
```

is an automorphism because `p` is a unit modulo `q`. For a polynomial with
coefficients in `F_p`, the freshman's dream is the group-algebra identity

```text
X^p = sigma_p(X).                                      (1)
```

For coefficients in a finite extension of `F_p`, the right side also applies
the coefficient Frobenius.

Call an integral twist `Theta` **p-stable modulo q** when

```text
sigma_p(Theta)=Theta.                                  (2)
```

Since its reduced coefficients lie in `F_p`, (1)-(2) give `Theta^p=Theta`.
Thus, for every `Lambda in A` and `m>=0`,

```text
Theta * Lambda^(p*m)
  = Theta^p * (Lambda^m)^p
  = (Theta * Lambda^m)^p.                              (3)
```

The only exponent classes contributing to `ct_q(Y^p)` are the p-dilates of
the identity class, and `p` is invertible modulo `q`. Therefore

```text
ct_q(Theta * Lambda^(p*m))
  = ct_q(Theta * Lambda^m)^p.                          (4)
```

In particular, a nonzero twisted constant term stays nonzero at all levels
`p^r*m`.

## 2. Exact-period kernels are p-stable

Let `U_q=(Z/qZ)^*`. The primitive-support kernel

```text
E_q(u)=sum_(a in U_q) u^a                              (5)
```

is p-stable for every `p` coprime to `q`, since multiplication by `p`
permutes `U_q`.

The same statement in the character basis is the elementary Ramanujan
identity

```text
c_q(p*n)=c_q(n),                 gcd(p,q)=1,           (6)
c_q(n)=sum_(a in U_q) zeta_q^(a*n).
```

Indeed `a -> p*a` merely permutes the primitive characters. Hence the cyclic
kernel

```text
R_q(u)=sum_(x mod q) c_q(x) u^x                        (7)
```

is also p-stable, as is any rational linear combination of exact-period
projectors whose denominators are p-units. Equations (4)-(7) prove that
Frobenius retains the **whole primitive-period packet**, not one selected
phase.

## 3. Relation to THM-2022

THM-2022 preserves the whole lowest balanced Wick face after Kummer removes
non-dilated allocations and the factorial height removes dilated off-face
allocations. The surviving sum becomes `Q^p`. Formula (4) is the exact-period
analogue of the last step: once an LRC certificate has been organized as a
p-stable projector, internal cancellation in that tied packet is transported
as one Frobenius power.

There are two important differences.

1. THM-2022 has DvdK to produce a nonzero base face constant term. No theorem
   in this argument produces a nonzero LRC exact-period moment.
2. The Gaussian Wick factorial supplies a strict p-adic filtration separating
   off-face channels. A raw LRC sinc or safe-phase sum has no corresponding
   monotone factorial. Formula (4) preserves a layer already isolated by some
   other argument; it does not isolate it.

Thus THM-2041 is a genuine reusable transfer lemma, but not TTNC and not
LRC(14).

## 4. Composite period versus coefficient characteristic

The theorem works for composite `q`, including `q=14`, at every prime
`p` coprime to `q`. This does **not** repair the composite-14 obstruction in
the Sungkawichai--Trakulthongchai polynomial method. Their interpolation is
performed over the non-field `Z/14Z`, where nonunits and null polynomials are
load-bearing. Here `q` is the phase period while `p` is a separate coefficient
characteristic. Keeping those roles distinct is essential.

## 5. What it buys for finite LRC certificates

For an actual primitive safe-residue count `D_q`, one has
`0<=D_q<=phi(q)`. Hence a congruence `D_q != 0 mod p` with `p>phi(q)` is already
an honest positivity certificate. THM-2041 can propagate such a certificate
once it is represented by a p-stable twist. By contrast, a nonzero signed
Ramanujan trace need not imply `D_q>0`; endpoint/open-boundary labels and the
safe-phase inequality cannot be discarded. This exactly matches the quotient
guardrail in HYP-2978/2979 and HYP-3036.

## 6. Verification and methodology

The companion script checks (4) for primitive-support and Ramanujan kernels
over all tested `2<=q<=24`, all coprime primes through `23`, five moment
levels, and deterministic group-algebra inputs. It also checks (6) exactly on
a larger range and includes a non-invariant-twist negative control.

Tournament Analysis is intentionally not imposed on the arithmetic test.
The operative relation is equality under a cyclic-group automorphism, not an
antisymmetric pairwise comparison. Turning phase classes into tournament
vertices would require an extrinsic tie-breaker and would forget the whole
projector sum that (4) is designed to preserve. Candidate vertices considered
were runners, phases, primitive periods, Frobenius orbits, projector kernels,
and proof obligations; the faithful quotient is the projector packet.
