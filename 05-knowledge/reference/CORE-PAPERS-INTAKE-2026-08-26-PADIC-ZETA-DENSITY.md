# 2026-08-26 intake: p-adic-zeta density and the specialization/Cartier gap

**Audit date:** 2026-08-26--27. This is a source and dependency audit. It
separates statements proved in repo canon from claims made in an external
unrefereed manuscript. A broken proof step makes its dependants unsupported;
it does not by itself make their conclusions false.

## 1. Source identity and immutable pins

The two supplied links do not identify the same work.

- [`arXiv:2608.23661`](https://arxiv.org/abs/2608.23661) is Raphael Cerf,
  *theta(p_c,Z^d)=0?*, a 708-page percolation manuscript. It contains no
  p-adic-zeta or Cartier argument. Local PDF SHA-256:
  `37c5759165df93e14da20db1de9f0117d67e0d16a1d1dd385ce16d85c0422a97`.
- The relevant source is Christopher D. Long,
  [`octonion/p-adic-zeta-density`](https://github.com/octonion/p-adic-zeta-density),
  audited at immutable commit
  [`4c87bcdf4d7d62d0f1981f16e228901f02cd9f57`](https://github.com/octonion/p-adic-zeta-density/tree/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57).
  The 51-page manuscript is titled *Multiweight arithmetic holonomy and
  density-one irrationality of p-adic zeta values*. The repository has no
  tag and its four manuscript commits contain no repair of the two steps
  below. TeX/PDF SHA-256:
  `cd100ccf2103093d65a07e7ce77bda4f753594dba7cfb8e7f34a3ca20e4f6b87` /
  `eed1b53981e12967ee8db0781c7bc784fec2a443b27a57543fe549eae204d62a`.

**External status:** **AUTHOR-CLAIMED / UNREFEREED RESEARCH DRAFT; TWO
LOAD-BEARING PROPOSITIONS OPEN.** No density or irrationality headline is
promoted to proved canon by this intake.

## 2. Exact specialization failure

For

```text
S=k[u][[f]],                 ev:S -> k[[f]],
ev(u)=f,                     ev(f)=f,
```

the manuscript needs scalar vanishing below a pivot to imply coefficientwise
vanishing in `S`. This is false. The minimal witness is

```text
u-f != 0 in S,               ev(u-f)=0.
```

The exact replacement, proved in
[THM-4255](../../01-canon/theorems/THM-4255-specialization-kernel-and-transverse-hasse-jet-repair.md),
is

```text
ker(ev)=(u-f)S,
ev^(-1)(f^n k[[f]])=(u-f)S+f^nS.                    (1)
```

Therefore specialization order records only graph-ideal-plus-tail membership.
It cannot be promoted to coefficientwise source order unless the actual
restricted digit module has an additional injective short-jet theorem.

There is also a Cartier incompatibility. If `C_r^f` extracts residue classes
of the explicit source `f`-exponent while keeping `u` in the coefficient ring,
then for every `ell>=2`,

```text
C_0^f(u-f)=u,          C_1^f(u-f)=-1,
C_r(ev(u-f))=0.                                        (2)
```

Cartier after graph specialization is well-defined, but it is not this
coefficientwise source decomposition. Specialization moves `u^a f^i` from
channel `i mod ell` to scalar channel `i+a mod ell`.

## 3. Proposition 6.2

The proof is at
[TeX lines 1144--1164](https://github.com/octonion/p-adic-zeta-density/blob/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57/multiweight_arithmetic_holonomy_all_primes.tex#L1144-L1164),
PDF page 14.

After taking pole-grade symbols, the proof says that every coefficient below
the scalar Bost pivot vanishes and concludes that the unspecialized symbol
`Pi_b` lies in `f^(n+C) A[[f]]`, where `A` still contains the torsor
coordinate. This is exactly the invalid implication refuted by (1). Descending
induction with Lemma 6.1 begins only after that step and cannot restore it.

**Status ledger:**

- Proposition 6.2's bound `1<=a<=s_*`, assuming the stated integrality:
  **SURVIVES**.
- Proposition 6.2's representation `n=j+ell*r` in the stated narrow window:
  **OPEN / UNPROVED**.
- Lemma 6.1, the coefficient-ring Taylor no-backflow statement:
  **VALID IN ITS STATED ALGEBRA**; it does not prove specialization strictness.
- Proposition 6.4's Fitting-length bound: **NOT INVALIDATED**; it does not
  reconstruct the missing prime window.

An abstract hostile to the inference uses

```text
D=10, C=0, ell=23, u=1+f^15, z=f^23,
ell^(-1)(u-1)z -> ell^(-1)f^38.
```

The scalar pivot `38` has no representation `38=j+23r` with `r>=1` and
`0<=j<=10`. This is a countermodel to the algebraic inference, not an actual
modular-geometric counterexample.

## 4. Proposition 12.3

The proof is at
[TeX lines 1997--2017](https://github.com/octonion/p-adic-zeta-density/blob/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57/multiweight_arithmetic_holonomy_all_primes.tex#L1997-L2017),
PDF page 25.

The sentence passing from `Pi_b(q)=Gamma_b(q,q^ell)` to
`ord_q Gamma_a(q,q^ell)=n` has no consistent supporting type:

1. If `Gamma_b` is still coefficient-valued before torsor specialization, the
   scalar pivot does not determine its coefficient-ring order, by (1).
2. If `Gamma_b` is already a scalar specialized germ, Lemma 12.1's
   componentwise bound on the unspecialized global coefficient lattice does
   not apply to its coefficients.

Lemma 12.2's exact block no-backflow statement is valid once its hypotheses
are supplied. The missing statement is the hypothesis that the actual
specialized nonzero scalar digit has order at most `D+C<ell`.

An independent abstract hostile uses

```text
D=12, C=1, ell=29, u=1+f^17, z=f^29,
ell^(-1)(u-1)z -> ell^(-1)f^46.
```

There is no `46=j+29r` with `r>=1` and `-1<=j<=13`.

## 5. Downstream dependency audit

The following are **OPEN / unsupported by this manuscript**, not refuted.

Through Proposition 6.2:

- Proposition 6.3 and equation (59);
- Proposition 6.6 and equation (62);
- Theorems 1.1 and 1.2;
- the finite-count argument in Theorem 1.3, including the `p=13` bound;
- Theorem 1.5, the genus-zero density-one conclusion.

Through Proposition 12.3:

- Proposition 12.4 and equation (83);
- Proposition 13.5's large-prime hybrid bound and Theorem 13.6;
- the denominator/slopes assembly and slopes clauses of Theorems 13.7 and
  13.11; their separate geometric capacity identities are not invalidated;
- Proposition 14.2 and Theorems 1.4 and 1.6, including positive-genus finite
  counts and the all-prime density conclusion.

The manuscript's numerical terminal certificate remains only a conditional
calculation. It cannot validate either missing Cartier gate.

## 6. Strongest repair targets

For an actual restricted source `V subset S`, the exact first gate is

```text
J_(V,n): V/(V intersect f^nS) -> k[[f]]/(f^n).       (3)
```

The desired order inference holds exactly when (3) is injective, equivalently
when

```text
V intersect ((u-f)S+f^nS)=V intersect f^nS.          (4)
```

The quotient by `V intersect f^nS` must not be omitted. A separate gate must
show strictness for the `ell`-adic coefficient filtration; for example
`u-1 -> ell f^N` under `u=1+ell f^N` shows that graph normalization alone
does not supply it.

Two possible repairs remain:

1. **Bivariate route:** retain torsor and base variables, compute (4) on every
   genuine digit space, retain all shifted Cartier channels, and prove
   `ell`-adic associated-grade strictness.
2. **Scalar-first route:** after actual specialization, prove every nonzero
   scalar digit `w_r(q)` is `ell`-integral and satisfies
   `ord_q w_r<=D+C<ell`; Lemma 12.2 then gives the claimed block window.

Without either repair, the strongest evident denominator-pivot survivor is
only the coarse pair

```text
1<=a<=s_*,                   ell<=n+C,
```

which loses the short degree window needed by the stated density estimates.

## 7. Relation to prior repo work

- The older `p-adic-zeta-irrationality` draft's 22 singleton claims remain
  **AUTHOR-CLAIMED / UNREFEREED**. This audit does not automatically refute
  them, but it makes their torsor product-digit gate a targeted priority:
  compute (4) at the narrow `(p,s)=(5,5)` cell. THM-4089's terminal formula
  optimization and THM-4091's LCM-coordinate theorem remain safe.
- THM-3488 is a positive firewall: its Rule-30 calculation retains the
  transverse Hasse channel when the base Cartier value vanishes.
- For the planar Jacobian lane, the typed transfer is to retain a wall-normal
  strict transform or conormal Hasse jet before lifting a `W=0` identity.
  Current `W=0` exclusions were recomputed directly and are not invalidated.
- For LRC(14), the analogy is methodological only: local Hasse/Frobenius data
  can be height-blind, so owner height remains a required sidecar. No
  loneliness estimate follows.
- Incoming THM-4257 is a positive density-transfer control, not a p-adic
  theorem. It first proves an exact finite exponent-orbit compiler, monotone
  suffix lifts, an all-height union, and complete high-bit cylinders; only
  then does it pass to density one. A repaired p-adic route needs the analogous
  exact graph-normal source and lift law before any asymptotic density step.

## 8. Does not prove

This intake proves no p-adic-zeta irrationality or density theorem, no
counterexample to one, no LRC(14), no planar Jacobian or Dixmier conjecture,
and no Rule-30 prize. It proves the general graph-kernel firewall in THM-4255
and identifies exact additional obligations for the external proofs.
