---
id: THM-3443
title: "The weighted-lift infinity torsor carries a nonzero fundamental leading amplitude"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  On the cyclic Q=infinity branch torsor of THM-3440, the normalized source
  coordinate x has leading profile j -> const*zeta_n^j, and its mode-one
  Fourier component is a formal unit descending to the unramified parameter.
  At n=91, CRT identifies the leading residue with a nonzero (1,1) line after
  linked branch/root markings.  The full normalized Puiseux germ is mode-mixed,
  beginning with a nonzero (-1,-1) correction.  THM-3450 subsequently gives
  an explicit carrier isomorphism after all character and amplitude gauges are
  marked; it does not give an unmarked H1 map, bispectrum identity, physical
  current, full-germ factorization, or LRC(14) consequence.
source: codex2 Puiseux reconstruction-current derivation, 2026-08-15
audit: independent exact Puiseux/Fourier/CRT referee plus normal/-O/stored replay; full-germ purity corrected to leading-residue/mode-one scope
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-3440-weighted-lift-cyclic-infinity-torsor-and-7x13-character-grid
related:
  - THM-2512
  - THM-3431-valuative-persistence-multiset-and-lrc-jc-boundary
  - THM-3450-marked-d5-carrier-isomorphism-and-full-germ-margin-obstruction
script: 04-computation/jc_weighted_lift_infinity_character_thm3443.py
output: 05-knowledge/results/jc_weighted_lift_infinity_character_thm3443.out
script_sha256: 912b264ad66011f96c3b895825cae528abff54367e0fab261e34a47b3655a042
output_sha256: 812a719c476a72b5e663798d6eb78a396286173e15a0778329af7a8f6c987404
semantic_sha256: 2cba88397477fe124be6a39891d2a55ceb512caee11443b83476e0c0c2ccfee8
hash_basis: LF-normalized bytes
---

# THM-3443 -- the weighted-lift infinity torsor carries a nonzero fundamental leading amplitude

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Exact local expansion

Use a lawful weighted lift from THM-3438.  Let

```text
d=deg(p),             n=d+1,
p(w)=a_d w^d+lower,   a_d!=0,
R(w)=integral_0^w p=(a_d/n)w^n+lower.                 (1)
```

Fix a generic `P=P_0` and the target line `(A,B,C)=(Q,P_0,1)`.  Put

```text
t=Q^(-1),             s^n=t,
w=s^(-1)W.                                               (2)
```

The inverse equation `R(w)-P_0w+cQ=0`, after multiplication by `t`, becomes

```text
(a_d/n)W^n+c+terms divisible by s=0.                    (3)
```

Choose `rho` with

```text
(a_d/n)rho^n+c=0.                                       (4)
```

Over a coefficient field `K` containing `rho` and `zeta_n`, Hensel's lemma
gives `n` branches

```text
W_j(s)=rho zeta_n^j+O(s),
w_j(s)=s^(-1)W_j(s),             j in Z/n.              (5)
```

These are the marked cyclic branches of THM-3440.

The weighted reconstruction has

```text
gamma=(P_0-p(w))/c,          x=1/gamma                  (6)
```

on `C=1`.  Substitution of `(1)` and `(5)` gives

```text
gamma_j(s)=-(a_d/c)rho^(n-1)zeta_n^(j(n-1))s^(-(n-1))
           +O(s^(-(n-2))),                              (7)

X_j(s):=s^(-(n-1))x_j(s)
       =-(c/a_d)rho^(-(n-1))zeta_n^j+O(s).              (8)
```

The exponent in `(8)` is `+j` because `-(n-1)=1 mod n`.

## 2. Exact Fourier nonvanishing

With the marked branch ordering, define

```text
Xhat(k;s)=(1/n) sum_(j mod n) X_j(s) zeta_n^(-kj).      (9)
```

Orthogonality and `(8)` give

```text
Xhat(1;0)=-(c/a_d)rho^(-(n-1))!=0,
Xhat(k;0)=0                         for k!=1.            (10)
```

Therefore `Xhat(1;s)` is a unit of the formal power-series ring.  In
particular it is identically nonzero as an algebraic Puiseux germ.  This is a
structural consequence of the weighted reconstruction, not a numerical
sample and not a genericity guess.

The coherent deck law is stronger and also records the limitation:

```text
X_j(s)=zeta_n^j X_0(zeta_n^(-j)s),
Xhat(k;zeta_n s)=zeta_n^(1-k) Xhat(k;s).                (10a)
```

Thus the coefficient of `s^m` can occur only in mode `k=1-m mod n`.
In particular

```text
Xhat(1;s) in K[[s^n]]=K[[t]]                            (10b)
```

is a unit, whereas the other Fourier modes have zero constant term but need
not vanish identically.  After the unique ramification scaling
`s^(-(n-1))`, the source coordinate `x` therefore supplies a natural
Keller-side **leading residue profile** and a nonzero mode-one component.  It
does not make the whole normalized germ a pure character vector or a physical
current.

## 3. The degree-91 `(1,1)` line

For THM-3438's degree-91 seed, `a_d=-1`, `c=1`, and `(4)` is

```text
rho^91=91.                                               (11)
```

Let `tau` be the oriented `C_91` inertia generator.  In the CRT coordinates
of THM-3440, use the factor generators

```text
(1,0)=tau^78,             (0,1)=tau^14.                 (12)
```

Normalize the corresponding primitive roots by

```text
eta_7=zeta_91^78,         eta_13=zeta_91^14.            (13)
```

Since every branch label has the CRT reconstruction

```text
j=78(j mod 7)+14(j mod 13) mod 91,                      (14)
```

one has

```text
zeta_91^j=eta_7^(j mod 7) eta_13^(j mod 13).            (15)
```

Thus the fundamental character in `(8)--(10)` is exactly the `(1,1)` tensor
character at leading order relative to the explicit factor-root normalization
`(12)--(13)`, and

```text
Xhat_91((1,1);0)=rho^(-90)!=0.                          (16)
```

The audit computes the first terms on the branch `j=0`:

```text
W_0(s)=rho+s/90+s^2/(180rho)+O(s^3),
X_0(s)=rho/91-s^2/(16380rho)+O(s^3).                   (17)
```

The nonzero `s^2` term lies in mode `90`, whose normalized CRT address is
`(6,12)=(-1,-1)`.  Hence `(16)` upgrades THM-3440's common character
**carrier** to one explicit nonzero Keller-side leading amplitude and a
mode-one formal unit, while `(17)` refutes purity of the entire germ.

## 4. Typed bridge and missing sidecar

After linked choices of the Puiseux root `s`, inertia generator, branch zero,
`rho`, and factor roots `(13)`, the JC leading-residue line in `(16)` is
one-dimensional.  THM-2512's LRC `(1,1)` coefficient also lives in a
one-dimensional `7 tensor 13` eigenspace.  Once bases are selected there is a
unique linear map sending the marked LRC basis to the leading coefficient
vector—or to the mode-one component `(10b)`.  This is a marked carrier
identification, not yet a natural D5 map.

What is not canonical is the amplitude comparison: changing `rho` or changing
the inertia/root markings incoherently rescales or readdresses the JC basis,
while the LRC coefficient has its own ancestry, centering, and response
normalization.  Moreover `(17)` shows that a map of full response germs must
respect deck semivariance and mode mixing, not merely one residue line.  A
marked identification therefore still needs an amplitude calibration;
dimension matching and separate nonvanishing do not provide one by themselves.

THM-3450 supplies exactly that calibration on the marked mode-one line: it
multiplies the LRC primitive coefficient by the Keller mode-one unit.  The map
is explicit but belongs to an `R^times` torsor and changes under independent
gauge reversal.  THM-3450 also proves that the full Keller germ cannot factor
through the doubly-centred interaction alone: its first two nonprimitive terms
occur in the two ANOVA margin sectors.  Thus the remaining D5 object is a
gauge-linked Rees connection carrying both margins and ancestry semantics, not
another bare one-dimensional intertwiner.

| field | exact content |
|---|---|
| source | weighted inverse branches and reconstructed coordinate `x` |
| target | a nonzero leading eigenline and mode-one formal unit |
| map | ramification normalization `(2),(8)` followed by branch Fourier transform `(9)` |
| preserved | leading cyclic character, CRT address under linked gauges, valuation, and formal nonvanishing |
| destroyed | higher-mode data under residue projection, finite-target values, LRC ancestry/centering, positivity, physical time, and response normalization |
| successor / needed sidecar | THM-3450 gives the marked mode-one amplitude isomorphism; a gauge-linked Rees connection with both ANOVA margins is still missing |
| cheapest hostile | reverse `tau` while fixing roots: `(1,1)` becomes `(6,12)`; coherent root reversal restores the label |

## 5. Boundaries

- The result is a formal germ over target infinity, not a finite Jelonek
  component or a claim about every target value.
- The nonzero objects are the leading residue coefficient and the mode-one
  component of a reconstructed source coordinate.  The whole normalized germ
  is mode-mixed by `(17)`.  None of these is a bispectrum or proves THM-2512's
  LRC contraction nonzero.
- The `(1,1)` label uses the explicit CRT root normalization `(13)`.  With
  conventional roots `zeta_7=zeta_91^13`, `zeta_13=zeta_91^7`, the same
  leading character has exact label `(6,2)`, not `(1,1)`.  This is root-gauge
  readdressing, not scalar rescaling inside one fixed labelled line.
- No additive map to THM-3437's characteristic-zero Tor tower is asserted.
- THM-3450's successor map is characteristic-zero and marked; it is not a map
  from the mod-13 `H^1` class and does not extend to the full germ.
- Nothing here proves LRC(14), `JC(2)`, or a physical D5 current.

## 6. Exact companion and audit

The deterministic companion verifies the Newton/Hensel equations through
`s^2`, exponent conversion in `(7)--(10)`, the degree-91 constants, deck
semivariance, mode support, both CRT root conventions, all generator/root
gauges, and the full-germ-purity hostile `(17)`.  Normal and optimized runs
reproduce the stored transcript byte-for-byte after LF normalization:

```text
python -B 04-computation/jc_weighted_lift_infinity_character_thm3443.py
python -B -O 04-computation/jc_weighted_lift_infinity_character_thm3443.py
```

The script, transcript, and semantic hashes are recorded in the frontmatter.
