---
id: THM-2445
title: "Twenty-four-cell graft owner conditioning"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  complete post-THM-2436 branch nu_7(c_3)<=M, the
  canonical locally drifting (c_3,q_*) graft has an exact
  six-by-four positive partition. Some cell has drift energy at
  least 1/576 of the isolated energy. Twenty-three cells retain a
  literal blocker-danger or first guard/unit-failure label. The sole
  all-safe ghost either has a drifting source-danger companion or,
  if that companion is circulant, admits an unfiltered C_7 source
  completion at one fixed eligible target/deep coefficient. Thus
  arbitrary bare owner-mask cancellation reduces to typed repair
  currents or a source-residue-refined current. THM-2442 subsequently
  restores any fixed positive canonical delayed word to both ghost
  subbranches while preserving the common partial bare endpoint.
  Semantic word/repair alignment for the other twenty-three cells and
  fully masked terminal endpoint transport remain open; no scalar row
  is removed.
source: codex-2026-07-26-post-deep-c3-graft-conditioning
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2409-unfiltered-septimal-source-completion-and-word-phase-boundary
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2426-compositional-thirteen-root-final-septimal-lane-exclusion
  - THM-2442-delayed-word-septimal-source-completion
  - THM-2448-right-endpoint-cospan-transition-atlas
script: 04-computation/lrc14_twenty_four_cell_graft_conditioning_thm2445.py
output: 05-knowledge/results/lrc14_twenty_four_cell_graft_conditioning_thm2445.out
script_sha256: fe989a84df2141bf0ebd4a7635a731900531091bdbc3eefcb9e6c3e7d008d1d6
output_sha256: 3f3ec3ed6aaedfcfd32d9f81d69fe497b90d66f6a54b397a12fc068f6c6bf138
hash_basis: working-tree bytes (LF)
---

# THM-2445 -- a drifting graft has only one ownerless bare cell

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

After THM-2436, the complete surviving valuation branch is:

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5)),

                              nu_7(c_3)<=M.             (1)
```

Choose a fixed label `q_*` among the guard and five ordinary roles
with `nu_7(q_*)=M`, and put

```text
L_*=2 if q_*=H,              L_*=1 otherwise.
```

Choose the THM-2309/2365 owner pivot so that `q_*` is the actual
label `k_c` grafted to the deepest target coordinate;
THM-2367 Section 8 permits any of the six such labels. After
common-gcd reduction, THM-2367 puts this isolated lawful
`(c_3,q_*)` graft on its noncirculant side. Merely choosing an
arbitrary isolated pair without this owner-pivot identification would
not supply the target quotient used below. The remaining question was
whether the seven other present factors could cancel that drift
without leaving a typed physical current.

They cannot do so invisibly. There is an exact finite alternative:

```text
literal blocker status or guard/unit first failure,

genuine source-owner drift,

unfiltered source-residue completion.                         (2)
```

The `24` cells share one common **partial bare** right endpoint. The
source completion co-shifts the inserted source factor in both
endpoint copies. The delayed word is deliberately absent.

## 1. The isolated table

Retain the notation and lawful target coordinates `(r,s,t)` of
THM-2365/2367. Let

```text
T_0:F_13^3 -> Q_(>=0)                                  (3)
```

be the isolated `(c_3,q_*)` overlap table, embedded with its unused
second target coordinate. It contains:

- the deepest probe `d(c_3x-r/13)`;
- the deepest safe complement `g(c_3x-t/13)`; and
- the lawful `q_*` safe graft `u_(L_*)`.

Put `P` for THM-2365's circulant projection, `Q=I-P`, and

```text
D_0=||Q T_0||_2^2>0.                                  (4)
```

The strict positivity is exactly THM-2367's criterion combined with
(1), for both `L_*=1` and `L_*=2`. The arithmetic reduced blocker
quotient

```text
C=c_3/gcd(c_3,q_*)
```

is a `7`-unit and remains present in the address typing of every table
below. It is not being reified as an additional physical factor.

## 2. The six-by-four partition

The six guard/unit roles consist of `q_*` and five other labels. In
any fixed order, let

```text
g_1,...,g_5
```

be the five remaining lawfully co-shifted safe indicators and put
`d_i=1-g_i`. Define

```text
S_0=product_(h=1)^5 g_h,

S_i=(product_(h<i)g_h)d_i,             1<=i<=5.       (5)
```

Pointwise, including after every lawful target shift,

```text
sum_(i=0)^5 S_i=1.                                     (6)
```

Let `j` be the selected source blocker and `a` the other nondeep
blocker in the THM-2365 owner pivot. For

```text
sigma=(epsilon_j,epsilon_a) in {0,1}^2,
```

put

```text
B_sigma
 =d_j^(epsilon_j)g_j^(1-epsilon_j)
  d_a^(epsilon_a)g_a^(1-epsilon_a).                   (7)
```

Then `sum_sigma B_sigma=1`. Insert (5) and (7) on the one moving
present endpoint of the isolated integrand, leaving its partial bare
right endpoint fixed, and call the resulting nonnegative rational table

```text
T_(i,sigma).                                           (8)
```

This accounts for all nine present roles:

```text
six guard/unit roles + {source blocker j, target blocker a, c_3}.
```

The isolated core keeps `q_*` and `c_3`; (5) partitions exactly the
other five guard/unit roles and (7) exactly the other two blockers.
This is the asymmetric endpoint deletion identity of THM-2370
Section 6. Missing right-endpoint factors are `1`; they are not being
silently restored.

Equations (6)--(8) give the exact identity

```text
T_0=sum_(i,sigma)T_(i,sigma),            6*4=24.       (9)
```

Every cell retains the same deepest probe and deepest safe complement,
so

```text
T_(i,sigma)(t,s,t)=0                                   (10)
```

almost everywhere. It also retains the same target quotient, the
`q_*` graft, the reduced arithmetic `7`-unit
`C=c_3/gcd(c_3,q_*)`, and the common partial bare endpoint.
Thus every cell has the full THM-2365 finite transform and relation
address typing.

Applying `Q` to (9) and using the triangle inequality gives a cell
`(i,sigma)` such that

```text
||Q T_(i,sigma)||_2>=sqrt(D_0)/24,

D(T_(i,sigma))>=D_0/576.                               (11)
```

Cauchy--Schwarz also gives the aggregate companion bound

```text
sum_(i,sigma)D(T_(i,sigma))>=D_0/24.                  (11a)
```

THM-2365 then gives eligible deep/target energy

```text
sum_(a!=0,(b,a+h)!=(0,0))
 |B_(i,sigma)(a,b,h)|^2
 >=D_0/7488.                                           (12)
```

Among the `2016` eligible colours, some amplitude is at least

```text
sqrt(D_0)/(24 sqrt(26208)).                            (13)
```

## 3. Twenty-three cells are already typed

Exactly one cell has neither a blocker-danger nor a first-failure
label:

```text
Z=T_(0,(0,0)).                                         (14)
```

For every other cell:

- if `i>0`, it is first a fixed literal guard/unit-failure layer,
  regardless of its blocker bits, and retains `d_i` and all preceding
  safe labels;
- if `i=0` and `sigma!=(0,0)`, all six guard/unit roles are safe and
  the cell retains a fixed nonempty local blocker-danger status.

Thus if the cell selected in (11) is not `Z`, equations (11)--(13)
already give a blocker-status-labelled or guard/unit-repair bare
current with an exact `91`-unit deep leg.

The status label lives on the indicator side. The nonzero target
character is not asserted to lie in the corresponding semantic
THM-2305 pure/fork support. That additional identification requires a
word/ANOVA intertwiner and is not part of this theorem.
Likewise, a first-failure cell is a labelled THM-2379-type repair
stratum, not automatically its canonical repair endpoint:
THM-2401's target/repair alignment and common-filter debt remains.

## 4. The unique ghost has a source companion

Suppose (11) selects `Z`. Let

```text
O=T_(0,(1,0)),

U=O+Z.                                                 (15)
```

Thus `O` replaces the stationary source-safe bit `g_j` by the genuine
source-danger bit `d_j`, while `U` omits that source factor. All three
tables retain (10), the same bare endpoint, and the same target
quotient.

If

```text
D(O)>0,                                                (16)
```

then `O` is a genuine source-owner bare current and THM-2365 applies
directly. No lower bound for `D(O)` in terms of `D_0` is claimed.

It remains to consider

```text
D(O)=0.                                                (17)
```

Since `U=O+Z`, (17) gives

```text
Q U=Q Z,

D(U)=D(Z)>=D_0/576.                                    (18)
```

Moreover `QO=0`, so every finite-transform coefficient with nonzero
target satisfies

```text
B_U(alpha,b,h)=B_Z(alpha,b,h).                         (18a)
```

Choose one eligible coefficient supplied by (12)--(13):

```text
B_U(alpha,b,h)!=0,

alpha!=0,                  (b,alpha+h)!=(0,0).         (19)
```

## 5. The same coefficient acquires every source colour

Work with the unfiltered word `Q_terminal=1`. For
`ell in F_7`, insert the translated source danger

```text
d_(j,ell)(x)=d(c_jx-ell/7)
```

in both bare/present copies and call the resulting table `O_ell`.
At the strict-open indicator boundary this is lawful because the two
source copies have the same unfiltered phase and `d^2=d` almost
everywhere; no independently transported word factor is being frozen.
Fourier coefficients are taken through one joint Poisson--Abel
parameter. No factorwise idempotence is asserted at finite Abel
radius.
The strict-open translates partition one almost everywhere, so

```text
sum_ell O_ell=U,                    O_0=O.             (20)
```

At the fixed eligible coefficient (19), put

```text
z_ell=B_(O_ell)(alpha,b,h) in Q(zeta_13).             (21)
```

Then

```text
sum_ell z_ell=B_U(alpha,b,h)!=0.                       (22)
```

Because (17) makes `O` circulant and (19) has nonzero target,

```text
z_0=0.                                                 (23)
```

For `kappa in F_7`, define the normalized source transform

```text
Z_kappa
 =1/7 sum_ell z_ell zeta_7^(kappa ell).                (24)
```

The cyclotomic fields of conductors `7` and `13` are linearly
disjoint. Hence `Phi_7` is irreducible over `Q(zeta_13)`. If one
`Z_kappa`, `kappa!=0`, vanished, the degree-at-most-six polynomial
`sum z_ell X^ell` would be a scalar multiple of `Phi_7`. All seven
coefficients would be equal, contradicting (22)--(23). Therefore

```text
Z_kappa!=0                     for every kappa!=0.     (25)
```

Parseval and the sharp anchored-zero Cauchy inequality from THM-2409
also give

```text
sum_(kappa!=0)|Z_kappa|^2
 >=|B_U(alpha,b,h)|^2/294
 >=D_0/4438167552.                                    (26)
```

Thus one fixed eligible deep/target coefficient survives in every
nonzero septimal source colour. The exact Abel expansion then contains
an atomic relation address with nonzero source residue and
`gcd(m,91)=1`, while retaining the labelled `(c_3,q_*)` graft. The
atomic exact address and ordinary frequency `X` may depend on
`kappa`; no common triangle, source residue modulo thirteen, or
all-coordinate-unit address is asserted.

## 6. Consequence and exact boundary

The complete post-THM-2436 branch (1) now has the bare-endpoint
trichotomy

```text
blocker-status / guard-unit-repair drift with (11),

genuine source-owner drift (16),

same-coefficient all-six source completion (25).        (27)
```

This is stronger than the untyped deletion conclusion of THM-2370,
but it respects that theorem's sharp clone hostile: no monotonicity of
drift under multiplication is used.

The operation in Section 5 by itself is lawful only for the unfiltered
terminal word. With a delayed word, shifting the source danger alone
fails to shift the word's source occurrence, exactly as in THM-2409.
THM-2442 Section 6.1 subsequently supplies the missing lawful skew:
in branch (17), it restores any fixed positive canonical word while
preserving the same coefficient `(alpha,b,h)`, labelled graft, and
common partial bare endpoint. In branch (16), THM-2365 Section 7 gives
the corresponding target-neutral BV restoration. Thus the ghost's
delayed-word phase is discharged.

The theorem plus that corollary still does not identify a
blocker-status label with a THM-2305 word support, retain a preselected
marked triangle, or transport the current to the canonical terminal
endpoint. The first-failure branches require a lawful
factor-repair/common-filter intertwiner. The no-failure blocker
statuses are local until the source, expiration, root, and clock data
are fixed. Restoring the fully masked endpoint can add cross-layer
terms not present at the partial bare endpoint. No scalar row is
removed, and LRC(14) remains open.

Reducing the clock modulo seven does not repair this. At `R=13`,
`y=1/7` gives a unique present source root but

```text
d(y-1/7)g(13y-13/7)=d(0)g(0)=0,
```

while the fixed word-safe factor has `g(13y)=g(6/7)=1`. Conversely,
at `y=1/26` the phase-zero product

```text
d(y)g(13y)=1
```

does not have the anchored zero used in Section 5. These are interior
points, not seam artifacts. Clock parity controls the coarse source
label but not the within-cell dilation/carry.

## 7. Exact companion

Run

```text
python 04-computation/lrc14_twenty_four_cell_graft_conditioning_thm2445.py
python -O 04-computation/lrc14_twenty_four_cell_graft_conditioning_thm2445.py
```

The dependency-free exact companion:

- exhausts all `32` five-bit states and verifies the six-cell
  lexicographic partition;
- exhausts the four blocker states and the full `24`-cell product,
  finding exactly one ownerless/no-failure ghost;
- checks every energy and coefficient denominator in (11)--(13) and
  (26);
- verifies over `Q` that, for each `kappa!=0`, reduction of the
  source-character polynomial modulo `Phi_7` has one-dimensional
  constant kernel; and
- checks the sharp `1/294` anchored-zero energy control.

Both transcripts must reproduce

```text
05-knowledge/results/lrc14_twenty_four_cell_graft_conditioning_thm2445.out
```

byte-for-byte.

## 8. Independent hostile audit

Two independent audits reconstructed the proof in complementary
directions. The first checked the `24`-cell partition, diagonal-zero
typing, all constants in (11)--(13) and (26), the coefficientwise
`Q_terminal=1` source completion, the two `R=13` hostiles, and exact
normal/optimized/stored replay. The second checked the guard-maximal
extension directly against THM-2309/2350/2367: the `L_*=2` graft has
the same noncirculance criterion, quotient, five-factor complement,
source typing, and constants as the ordinary `L_*=1` graft.

Both audits emphasized the same stopping boundary recorded above:
the result is a partial-bare-endpoint conditioning theorem. THM-2442
later restores the delayed word on the ghost, but neither theorem
supplies semantic word/repair alignment for the other cells or the
fully masked endpoint. THM-2448 subsequently expands every fixed
current through the omitted right factors into a finite
complete-mask/transition cospan; positive same-root physicalization
remains open.

QED.
