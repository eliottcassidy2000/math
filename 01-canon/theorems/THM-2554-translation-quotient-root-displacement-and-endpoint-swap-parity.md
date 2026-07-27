---
id: THM-2554
title: "Translation-quotient root displacement and endpoint-swap parity"
status: >
  PROVED + VERIFIED-EXACT. For a word stratum carrying a measure-preserving
  simultaneous C_13 root action, every equivariant source/later-root joint
  table is circulant and is classified exactly by the displacement d=b-h.
  The semantic diagonal is the free d=0 orbit, not a rotation-fixed locus;
  root marginals and C_13 localization alone cannot force it. Any permitted
  nonzero displacement supplies an off-diagonal perfect matching and a
  zero-arrival Hall coupling. If a sidecar-preserving endpoint swap is
  additionally supplied, it acts by d->-d. For a nonzero rational ledger in
  canonical primitive integer normalization, odd total quotient weight then
  forces positive zero-displacement weight. Both hypotheses are sharp: one
  nonzero offset is an odd zero-arrival hostile without endpoint swap, while
  the pair +/-1 is swap-invariant, Hall-perfect, zero-arrival, and has all
  twelve relative mixed characters nonzero. No current LRC theorem supplies
  both a genuine later target-active field and the required chronological
  endpoint-swap/oddness, so no row is removed.
source: codex-2026-07-27-translation-swap-arrival
depends_on:
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
related:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - HYP-9045-parity-rank-box-refutations-and-rotation-redirect
  - HYP-9055-shifted-word-palindrome-layer
script: 04-computation/lrc14_translation_swap_root_displacement_thm2554.py
output: 05-knowledge/results/lrc14_translation_swap_root_displacement_thm2554.out
script_sha256: d2ac6e35f992a757bd4bedc269be92f0c9efa30ac481e2c2c2f6efa464d0a0c6
output_sha256: 78cd17ce25f7ae9f86e423fd42af3db7514998d5536f389089dadeff2a962391
hash_basis: working-tree bytes (LF)
---

# THM-2554 -- translation quotient, root displacement, and endpoint-swap parity

**PROVED + VERIFIED-EXACT.**

THM-2538 identifies the remaining semantic statistic as a product of a
selected source/head field with a genuinely later target-active root field on
one ancestry base. THM-2545 then gives the exact weighted Hall criterion for
forcing their diagonal. The present theorem asks what the native root
translation symmetry does to that Hall problem.

It makes the missing coordinate one-dimensional, but it does not force it:

```text
simultaneous C_13 rotation of (head,later root)
  -> quotient coordinate d=later root-head;

semantic arrival
  -> the orbit d=0, not a rotation fixed point;

word/sidecar-preserving endpoint swap
  -> d maps to -d;

endpoint swap plus odd primitive quotient weight
  -> positive d=0 mass.                                      (1)
```

Thus the correctly primed rotation idea is a quotient operation, not by
itself a nonvanishing argument. A second, genuinely chronological symmetry or
a direct zero-displacement witness is still required.

## 1. The diagonal rotation quotient is displacement

Let `R=F_13`. On ordered pairs use the simultaneous translation action

```text
u.(h,b)=(h+u,b+u).                                          (2)
```

The difference

```text
d(h,b)=b-h                                                   (3)
```

is invariant. Conversely, two pairs have the same difference exactly when
they differ by a unique translation in (2). Hence the `169` ordered pairs
split into exactly thirteen free orbits

```text
O_d={(h,h+d):h in R},                 d in R,                (4)
```

each of size thirteen. In particular,

```text
diagonal={(h,h):h in R}=O_0.                                (5)
```

No nonidentity translation fixes any point of `O_0`; the semantic diagonal
is a whole free orbit. Moreover `(h,h)->(h,h+c)` is an equivariant bijection
from `O_0` to every `O_c`. Any invariant depending only on rotation-orbit
size, stabilizer, or fixed-point count therefore treats arrival and every
nonzero displacement identically.

This is the exact scope boundary for the HYP-9045/9050 rotation-localization
idea. A `C_13` congruence localizes onto nonfree sectors, but the action (2)
has no such sector. The diagonal is recovered only after taking the quotient,
where it becomes the distinguished value `d=0`.

## 2. Equivariant fields give circulant joint tables

Let `(Omega,nu)` carry a measure-preserving `R` action. Fix a terminal-word
stratum invariant under that action, and let

```text
h,b:Omega->R                                                (6)
```

be equivariant root maps:

```text
h(u.omega)=h(omega)+u,
b(u.omega)=b(omega)+u.                                      (7)
```

Here `h` is the selected head and `b` is stipulated to be a genuinely later
target-active root; existence of the latter is not assumed by current LRC
canon. Both labels must use the **same root torsor and the same translation
gauge**. Independent root charts do not define the difference (3).

The domain in this theorem is cemetery-free. If a partial selector has a
fixed cemetery value, first restrict to an invariant active subpacket on which
both roots lie in `R`; cemetery mass is not a fourteenth displacement and may
not be included in the parity count below. Their nonnegative joint table

```text
C(s,t)=nu({omega:h(omega)=s,b(omega)=t})                     (8)
```

obeys `C(s+u,t+u)=C(s,t)`. Put

```text
m_d=sum_s C(s,s+d).                                         (9)
```

Then invariance and (4) give the exact normal form

```text
C(s,t)=m_(t-s)/13.                                          (10)
```

If `M=sum_d m_d`, both one-root marginals are uniform:

```text
sum_t C(s,t)=sum_s C(s,t)=M/13.                             (11)
```

The desired semantic hit is simply

```text
H=sum_s C(s,s)=m_0.                                         (12)
```

For completeness, with `zeta=zeta_13`, the unnormalized two-dimensional
transform is

```text
Chat(alpha,beta)
 =0                                      if alpha+beta!=0,
 =sum_d m_d zeta^(beta d)                if alpha+beta=0.   (13)
```

Thus the complete rotation-covariant mixed spectrum is exactly the Fourier
transform of the displacement ledger. Separate source and later-root
marginals retain only `M`; they forget all twelve nonconstant displacement
coordinates, including the value `m_0` within that family.

## 3. Rotation-covariant Hall reduction

Suppose an ancestry sidecar permits offsets in `D subset R`. If `D` contains
one nonzero `c`, the table

```text
C_c(s,t)=(M/13) 1_(t=s+c)                                  (14)
```

has the uniform marginals (11), is supported on the allowed graph, and has
`H=0`. Its support is an off-diagonal perfect matching. Equivalently, after
deleting the semantic diagonal, every Hall inequality holds because
`S->S+c` is a bijection for every root subset `S`.

Consequently, within a translation-invariant word stratum:

- root marginals plus rotation covariance never force arrival;
- an allowed nonzero displacement already gives a zero-arrival coupling;
- a support-graph argument forces arrival for every invariant coupling only
  when its positive-mass support is restricted to `D={0}`.

This specializes THM-2545's general bipartite Hall search. It does not say
the physical compatibility graph is translation invariant; that action and
the word/ancestry preservation in (7) must first be proved.

## 4. Endpoint swap makes zero displacement the unique parity-fixed class

Now additionally suppose there is a measure-preserving involution `iota` on
the same word/owner/ancestry stratum with

```text
h(iota omega)=b(omega),
b(iota omega)=h(omega).                                    (15)
```

The endpoint swap commutes with every simultaneous translation in (2), so
the group generated on pair space is `C_13 x C_2` (equivalently `C_26`),
**not** a dihedral group. It descends under (3) to

```text
d -> -d.                                                    (16)
```

Thus `m_d=m_(-d)`. Because thirteen is odd, zero is the unique fixed
displacement under (16); the other twelve values form six pairs.

Assume the nonzero rational vector `(m_d)` is given its **canonical primitive
integer normalization**: multiply its rational ray by the unique positive
factor making

```text
n_d in Z_(>=0),              gcd{n_d:n_d!=0}=1.             (17)
```

This normalization is essential. Parity of an arbitrary denominator-cleared
multiple is not invariant under rescaling, whereas (17) is. Endpoint-swap
symmetry gives

```text
W:=sum_d n_d
  =n_0+2 sum_(d in {1,...,6}) n_d,

W congruent n_0 mod 2.                                     (18)
```

Therefore

```text
W odd  ->  n_0 odd  ->  m_0>0  ->  H>0.                    (19)
```

This is the promised rotation-then-parity criterion. It does not contradict
the mod-two no-go in HYP-9045: (18) is parity on a primitive **integer orbit
ledger after quotienting by the odd group**, not `C_2` cohomology with
odd-order `7 x 13` coefficients.

The same proof applies separately in each terminal-word stratum. A swap
which swaps word, owner, deep, or chronological roles does not license (18)
inside the stratum needed by THM-2545.

## 5. Both extra hypotheses are sharp

Without endpoint swap, one nonzero orbit `O_c` has primitive quotient weight one
and zero arrival. Thus oddness alone is useless.

With endpoint swap but even weight, take `c=1` and

```text
C_pm(s,t)=1_(t=s+1)+1_(t=s-1).                             (20)
```

This is endpoint-swap invariant, has primitive displacement weights
`n_1=n_(-1)=1`, uniform source and later-root marginals equal to two, and
zero diagonal. Its off-diagonal support contains two perfect matchings. It
also has every one of the twelve relative mixed characters nonzero: by (13)
the value at `beta!=0` is

```text
zeta^beta+zeta^(-beta)!=0,                                 (21)
```

since an odd-order root of unity cannot have square `-1`. Complete mixed
root colour therefore does not replace the parity sidecar. Conversely, its
individual root marginals are uniform and have no nontrivial colours; this
is exactly the mixed-Haar information loss in THM-2538.

There is an especially sharp equal-mass comparison. Count every atom with
unit weight. The aligned swap-invariant ledger `2 O_0` and the displaced
ledger `O_1 union O_(-1)` both have total mass `26`, uniform marginals equal
to two, simultaneous-translation symmetry, and endpoint-swap symmetry. Their
arrivals are respectively `26` and `0`. Thus even after fixing all those
statistics, only the displacement coordinate distinguishes the semantic
answer.

## 6. LRC consequence and exact remaining object

The theorem changes the description of the open semantic edge, not its
status. On a common invariant ancestry base the relevant object is not an
arbitrary `13 x 13` joint table but the thirteen-entry displacement current

```text
m_d=mass{later target-active root-selected head=d}.          (22)
```

There are now three sharply typed ways to finish that interface:

1. construct the actual later field and prove `m_0>0` directly;
2. prove its word-resolved compatibility support contains only `d=0`; or
3. construct a word/owner/ancestry/chronology-preserving endpoint swap
   and prove odd canonical primitive quotient weight.

Current canon supplies none of these. In particular, chronological source and
arrival roles are directed, so swapping their printed roots is not a lawful
swap. Choosing a tournament orientation of the pairs `+/-d` destroys,
rather than supplies, (15). THM-2545's two-root swap hostile is the smallest
unstructured version of the same obstruction; (20) is its fully
rotation-covariant and endpoint-symmetric form.

The exact HYP-9055 shifted-word palindrome probe lands on this boundary rather
than crossing it. Its singleton clock count is odd, but its pair count is
always even, and its clock mirror is not the chronological endpoint swap
(15). Its varying mod-`7`/`13` pair residues remain useful localization data,
but neither the singleton parity nor that mirror proves `m_0>0`.

No scalar row is excluded. The LRC(14) ledger remains `165`.

## Exact referee

Run

```bash
python3 04-computation/lrc14_translation_swap_root_displacement_thm2554.py
python3 -O 04-computation/lrc14_translation_swap_root_displacement_thm2554.py
```

Both executions byte-match the stored output. The referee classifies all
`169` ordered pairs into thirteen free orbits, reconstructs invariant tables
from arbitrary displacement ledgers, checks all `8,192` Hall subsets for the
`+/-1` hostile, verifies its twelve cyclotomic mixed modes exactly, and tests
the primitive endpoint-swap parity identity on `2,916` ledgers.

**QED.**
