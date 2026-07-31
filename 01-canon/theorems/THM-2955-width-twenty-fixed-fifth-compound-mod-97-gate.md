---
id: THM-2955
title: "Width-twenty fixed fifth-compound mod-97 gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  support (0,1,2,20), THM-2949's fixed rank-35
  Macaulay cofactor is nonzero at every integer depth n>=0.  After
  exact division by c300*q200^5 and an explicit degree-425
  negative-root factor, the primitive degree-506 quotient has no root
  modulo 97.  This is a discrete integer-depth certificate beyond the
  bounded Gregory--Newton atlas; no irreducibility or positive-real-ray
  claim is made.
source: codex-gmc-width-twenty-modular-depth-2026-07-29
audit: >
  Independent hostile audit accepted the cofactor typing, exact
  quotient and integral primitive normalization, coefficient order,
  complete mod-97 gate, outside-grid controls, rank-gap implication,
  replay hashes, and discrete-only scope.
depends_on:
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
related:
  - THM-2956-koszul-gale-fixed-fifth-compound-exchange
script: 04-computation/gmc_width_twenty_fixed_cofactor_mod97_gate_thm2955.py
output: 05-knowledge/results/gmc_width_twenty_fixed_cofactor_mod97_gate_thm2955.out
script_sha256: d5fc09fc0b7de7b4d9e364bdfb7c5c01dbcea96dc00377b753ecf197edb8feb0
output_sha256: ed68e586f378324dfbe8637191e8d06be3519e1b8443ddb654aa40dd1b2d1f7e
thm2949_dependency_sha256: 9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54
hash_basis: LF-normalized bytes
---

# THM-2955 -- width-twenty fixed fifth-compound mod-97 gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `P_20(n)` be the fixed `20Q+10C+5F` rank-`35` cofactor of
THM-2949 for normalized support

```text
(0,1,2,20).                                             (1)
```

Then

```text
P_20(n) != 0                         for every integer n>=0. (2)
```

By THM-2947's proved conjugate-pair corank gap, the degree-seven
Macaulay map therefore has full rank at every physical depth for this
support.  Thus first-window SFC(4) holds on
`(n,n+1,n+2,n+20)`.

This certificate is genuinely different from THM-2949's finite
Gregory--Newton atlas.  The old exact scanner finds no one-sign Newton
vector at any base `0<=r<=4M=80`.  Direct exact signs are

```text
n=0,...,61:      positive,
n=62,...,116:    negative,
n=117,...,200:   positive.                             (3)
```

In particular the cofactor has positive-real sign-changing seams, even
though it never vanishes at an integral depth.

## 1. Exact quotient

The cofactor has exact degree

```text
deg P_20 = 1065 = 55*20-35.                            (4)
```

For the pure coefficients

```text
q200=[x0^2]Q,                         c300=[x0^3]C,     (5)
```

the companion verifies

```text
deg q200=19,      every coefficient of q200 is positive,
deg c300=39,      every coefficient of -c300 is positive. (6)
```

Construct directly

```text
B_425(n)
 =(2n+1)^5
  (n+1)^6
  (n+2)^21
  [prod_(r=3)^6 (n+r)^26]
  [prod_(r=7)^10 (n+r)^24]
  [prod_(r=11)^13 (n+r)^20]
  [prod_(r=14)^20 (n+r)^19].                          (7)
```

Its exponent census is

```text
5+6+21+4*26+4*24+3*20+7*19=425.                       (8)
```

Exact polynomial division gives

```text
P_20/(c300 q200^5 B_425) = unit * R_506(n),            (9)
```

where `R_506` is primitive-normalized with positive leading
coefficient and degree `506`.  No factorization or irreducibility of
`R_506` is used.  Every factor removed in `(9)` is visibly nonzero for
real `n>=0`.

## 2. The mod-97 gate

The coefficient-vector SHA-256 of `R_506` is

```text
a0176d5c76c67c931de4217da581d1599e29d573514917f416638c1d2497355a. (10)
```

Exact Horner evaluation gives

```text
R_506(r) != 0 mod 97                     for every r in F_97. (11)
```

The complete table, in residue order `r=0,...,96`, is

```text
33,54,20,7,5,1,12,70,84,38,83,13,71,22,6,42,65,48,62,92,
62,23,1,26,95,9,79,25,1,88,7,11,22,58,31,92,54,76,53,2,
27,40,40,40,36,30,19,61,6,19,49,92,52,82,88,70,14,26,25,
79,13,52,12,12,95,48,44,43,32,75,29,29,26,68,2,73,63,19,
24,39,57,40,38,69,85,66,40,47,87,57,10,10,22,37,94,2,35.
                                                               (12)
```

Its table digest is

```text
81743f0f0a65a8c3d3b92e313b8ee73b3376cd656b9c9943d93be34f07181978. (13)
```

If `R_506(n)=0` for an integer `n`, reduction modulo `97` would make
the entry indexed by `n mod97` vanish, contradicting `(11)`.
Equations `(6)--(9)` therefore prove `(2)`.

## 3. Exact replay and scope

Run

```text
python -B 04-computation/gmc_width_twenty_fixed_cofactor_mod97_gate_thm2955.py
python -O -B 04-computation/gmc_width_twenty_fixed_cofactor_mod97_gate_thm2955.py
```

Both outputs byte-match the stored transcript.  The companion imports
the promoted THM-2949 chart, pins its LF hash, reconstructs all `1066`
interpolation values, checks three direct determinants outside the
grid, performs the exact divisions, primitive-normalizes `R_506`,
tests all `97` residues, and checks the sign runs `(3)`.

The theorem closes one width-twenty support.  It proves neither a
complete width-twenty atlas nor an arbitrary-width theorem.  The
positive-real seams show why the conclusion must remain a discrete
integer-depth statement.

## 4. Independent hostile audit

An independent reconstruction replayed both interpreter modes against
the stored transcript and obtained the declared LF hashes.  It also
rebuilt `P_20` and the quotient in `(9)` without using the frozen
residue table, checked that every coefficient of `R_506` has
denominator one, verified the ascending FLINT coefficient order used
by Horner evaluation, and compared direct primitive-quotient
evaluations with representative entries on both sides of the sign
seams.  The audit rechecked all `97` nonzero residues, all three
outside-grid determinants, and the precise THM-2947 rank-gap handoff.

The audit specifically rejects stronger interpretations: the
certificate neither proves irreducibility of `R_506` nor removes its
positive-real roots.  It only excludes integral roots, exactly as
required by physical factorial depth.

**QED.**
