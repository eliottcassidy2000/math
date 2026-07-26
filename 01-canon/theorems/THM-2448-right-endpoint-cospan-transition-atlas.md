---
id: THM-2448
title: "Right-endpoint cospan transition atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Fix any
  nonzero THM-2445 marked current at the common
  partial bare endpoint, including its exact X,m, target/deep colour,
  and any terminal word already retained by THM-2365/2442. Splitting
  the missing bare endpoint through the
  twenty-four THM-2445 labels, then recursively splitting the omitted
  tail factors in the sole matched first-failure branch, leaves at
  most sixty-nine fixed-current cospan pieces. Some piece has modulus
  at least 1/69 of the original and is either a complete matching
  local factor mask or has a fixed-order earliest guard/unit or
  blocker mismatch. In the source-completed ghost only twelve pieces
  are needed, with modulus at least 1/12 and the same nonzero source
  colour. These are coefficient-level two-sided cospans, not
  nonnegative drift tables; matching labels do not by themselves give
  a positive same-root filtration or a canonical owner/terminal
  current. No scalar row is removed.
source: codex-2026-07-26-right-endpoint-cospan
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2442-delayed-word-septimal-source-completion
  - THM-2445-twenty-four-cell-graft-owner-conditioning
related:
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2408-endpoint-prony-resultant-clock-separation
  - THM-2419-affine-relation-shell-pushforward-and-observer-homogenization
script: 04-computation/lrc14_right_endpoint_cospan_transition_thm2448.py
output: 05-knowledge/results/lrc14_right_endpoint_cospan_transition_thm2448.out
script_sha256: 70b30775ca333424fe05265cb6686b27d2fd1ae5a3fa4afd2d6bc749a207ec13
output_sha256: 3ff1e437511de19c56a7bdafccaa4f49b0abf0d1500ac5a1371a0d76aeb9d359
hash_basis: working-tree bytes (LF)
---

# THM-2448 -- the missing endpoint has a finite transition cospan

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2445 turns every row in the surviving `nu_7(c_3)<=M` branch into a
nonzero marked current carried by one of twenty-four labelled present
cells, or by its source-completed ghost. For the generic drifting and
`D(O)>0` branches, THM-2365 gives the ordinary delayed-word BV
restoration. For the circulant ghost, THM-2442 restores a fixed
positive delayed word without changing the chosen coefficient. All
routes deliberately stop at a common **partial bare** right endpoint.
For a generic typed branch before semantic word/repair alignment, the
argument below remains unconditional with terminal word `1` and
conditional for any stronger word-retaining current.

This theorem inserts the missing local factors at that endpoint without
changing the already selected Fourier data. It proves the exact finite
alternative

```text
complete matching local factor-mask cospan
or
fixed-order earliest two-sided factor mismatch.                 (1)
```

The first line is still only a cospan statement: equality of branch
labels does not produce positive same-root descendants or identify
semantic owner labels. The second line is a complex fixed-frequency
transition current, not a
nonnegative deletion table. Those two distinctions are the remaining
endpoint debt.

## 1. Fix the current before splitting the endpoint

Use the notation and fixed order of THM-2445 Section 2. The five
guard/unit factors other than the grafted `q_*` are

```text
g_1,...,g_5,       d_h=1-g_h,
```

and their first-failure partition is

```text
S_0=product_(h=1)^5 g_h,
S_i=d_i product_(h<i)g_h,       1<=i<=5.             (2)
```

Let `j` be the source blocker, let `a` be the other nondeep blocker,
and let

```text
B_sigma
 =d_j^(epsilon_j)g_j^(1-epsilon_j)
  d_a^(epsilon_a)g_a^(1-epsilon_a),

sigma=(epsilon_j,epsilon_a) in {0,1}^2.              (3)
```

Thus

```text
sum_i S_i=1,       sum_sigma B_sigma=1               (4)
```

almost everywhere. The strict-open endpoints form a null set. At every
lawful target/deep twist, the present and bare copies of these factors
are co-shifted by the same THM-2350/2365 coordinate action; (4) remains
valid after each such shift.

Fix one left/present cell

```text
alpha=(i,sigma).                                      (5)
```

After the finite target/deep transforms of THM-2365, the lawful source
transform when present, the delayed-word restoration of THM-2442, and
the absolutely iterated `m`-then-`X` collapse, fix a nonzero marked
current

```text
J_alpha=J_alpha(X,m;chi)!=0.                           (6)
```

Here `chi` records the already chosen target/deep finite character.
When a lawful source transform is present, its character `kappa` is
recorded separately. Any terminal word already supplied by
THM-2365/2442 is also fixed. Equation (6) is a THM-2334
residue-current aggregate; no individual relation address is selected.

The order of operations is load-bearing. We first make the Boolean
partition of each whole endpoint layer, then Poisson--Abel smooth the
whole layer. We do not smooth factors separately and invoke
idempotence at positive Abel radius.

## 2. The twenty-four right labels

On the partial bare right endpoint insert (4). For

```text
beta=(i',sigma') in {0,...,5} x {0,1}^2              (7)
```

let `J_(alpha,beta)` be the marked current obtained by multiplying the
right endpoint by the correspondingly co-shifted
`S_(i')B_(sigma')` at every finite target twist, while keeping every
selector in (6) fixed. This is the same action on present and bare
packets as THM-2365 equation (5), inherited from THM-2350 equation
(16). Linearity at positive Abel radius and the boundary passage of
THM-2334 equations (29)--(33) give

```text
J_alpha=sum_beta J_(alpha,beta).                      (8)
```

Every summand in (8) has the same:

- exact integers `X,m`;
- target/deep finite character `chi`;
- any delayed terminal word already retained in `J_alpha`;
- reduced arithmetic and graft label;
- deepest probe and deepest safe complement.

Consequently each underlying unspecialized whole layer retains the
table-level deep diagonal zero used by THM-2365. A single scalar
`J_(alpha,beta)(X,m;chi)` is not itself said to satisfy a
table-valued diagonal identity. The source-completed ghost instead
uses the restricted twelve-cell bank of Section 4; every still-missing
mask there is source-neutral, so its selected septimal source
character is unchanged.

For a fixed `alpha`, twenty-three right labels are unequal to it. The
fixed-order transition type of an unequal pair is unique:

1. if `i!=i'`, treat `0` as infinity and take the smaller nonzero
   first-failure index; at that guard/unit factor one endpoint is safe
   and the other dangerous;
2. if `i=i'` and `sigma!=sigma'`, take the first differing blocker bit
   in the declared `(source,target)` blocker order.

"First" here is relative to the once-fixed THM-2445 order. It is not an
intrinsic ordering of the physical roles.

Across all twenty-four possible left labels, the exact ordered-pair
census is

```text
matched                                      24,
guard/unit mismatch depths       160,128,96,64,32,
blocker mismatch positions                    48,24.
                                                            (9)
```

The entries in (9) total `576`.

## 3. Completing the sole matched branch

If `i=0`, the matched right term already specifies all five
guard/unit factors and both blocker bits on both endpoints. Suppose
instead that `1<=i<=5`. The first-failure word `S_i` leaves the factors
with indices `h>i` unspecified.

Inside the sole matched term, split each omitted factor on both
endpoint copies:

```text
(g_h^L+d_h^L)(g_h^R+d_h^R)
 =g_h^L g_h^R+g_h^L d_h^R+d_h^L g_h^R+d_h^L d_h^R.  (10)
```

Stop a branch as soon as an off-diagonal middle term in (10) occurs.
Continue only the two diagonal terms. A stopped branch has a
fixed-order earliest later guard/unit mismatch; a branch that reaches
the end has a complete matching local factor mask.

If `n` tail factors remain, the number `T(n)` of terminal cospan
pieces satisfies

```text
T(0)=1,
T(n)=2+2T(n-1)=3*2^n-2.                               (11)
```

For left unit labels `i=0,1,2,3,4,5`, the matched-term terminal counts
are therefore

```text
1,46,22,10,4,1.                                       (12)
```

The other twenty-three right labels were already off-diagonal.
Combining them with (12), the complete fixed-left decompositions have

```text
24,69,45,33,27,24                                    (13)
```

terminal pieces. The triangle inequality applied once to this complete
decomposition proves that some terminal current `J_*` satisfies

```text
|J_*|>=|J_alpha|/69.                                  (14)
```

It has exactly one of the two types in (1), and retains all the data
listed after (8). The denominator `69` is uniform; the branch-specific
denominators are those in (13).

The phrase **complete local factor mask** includes the five split
guard/unit roles, the two blocker bits, and the already retained
`q_*`/`c_3` data. It does not mean that positive same-root co-support,
absolute roots, semantic THM-2305 owner labels, or terminal
orientations have been identified.

## 4. The source-completed ghost needs only twelve cells

In the THM-2442 completion of the THM-2445 ghost, the source factor has
already been inserted in both present and bare endpoint copies and the
source transform has already been taken. Only

```text
five guard/unit roles + the other nondeep blocker              (15)
```

remain missing. Their right atlas has `6*2=12` cells. The fixed ghost
left label is all-safe with the remaining blocker safe. Hence there is
one complete matching local-mask cell, two transitions at each of the
five fixed-order guard/unit depths, and one blocker transition:

```text
1 + 2+2+2+2+2 + 1 = 12.                              (16)
```

For the fixed nonzero ghost current `J_Z(X,m;chi,kappa)`,

```text
J_Z=sum_(beta=1)^12 J_(Z,beta),

max_beta |J_(Z,beta)|>=|J_Z|/12.                     (17)
```

Every term retains the same nonzero `kappa mod 7`, exact `X,m`,
target/deep character, and delayed word. No additional source split is
performed.

The `D(O)>0` source-owner alternative of THM-2445 has unit label
`i=0`, so its generic right expansion already has the sharper
denominator `24`; its unique matched branch, if it is the surviving
one, is already a complete local mask.

## 5. What the atlas does and does not prove

### 5.1 A matched label is not yet a positive common-filter service

The lawful target action co-shifts the present and bare copies; this
theorem does not freeze one copy while moving the other. Nevertheless,
a surviving `J_(alpha,alpha)(X,m;chi)` is one fixed-frequency product
of two separate endpoint Fourier coefficients. It is not an integral
of pointwise same-root co-support. That positive object appears only
after the appropriate full frequency recombination and physical root
identification.

Thus a complete matching local mask is a **cell-matched two-sided
cospan**, not yet the positive common filtration with surviving
descendants required by THM-2401 equations (6), (10)--(15). A
same-root physicalization/reference sidecar would upgrade it; this
theorem does not supply one.

### 5.2 A transition current is not a positive drift table

The terms in (8) and (10) are fixed-frequency cross-endpoint currents.
They can be complex and signed. There is no pointwise inequality
between them, no inherited drift-energy floor, and no direct
application of THM-2379 to an off-diagonal label without a lawful
target-action identification.

This is why the tempting `24`-by-`24` nonnegative-table proof is false.
Even in one coordinate gauge, identical Boolean partitions collapse
pointwise to diagonal cells only after the full physical-frequency
recombination. Individual fixed-frequency cross products can survive
off diagonal and cancel; they are precisely the complex cospan
currents treated here.

### 5.3 The matched term can vanish

Linearity cannot force the matched branch to carry any current. For a
positive rational control, take a right packet that is the disjoint
union of a half-circle matched layer and a quarter-circle transition
layer. At physical frequency `2`, the half-circle coefficient is zero
while the quarter-circle coefficient is nonzero. Tensoring with fixed
nonzero target/deep factors gives

```text
J_matched=0,       J_transition=J_alpha!=0.            (18)
```

So the alternative (1), rather than unconditional matched-endpoint
survival, is sharp. This is an abstract rational endpoint-partition
control, not a claim that the two intervals form a realized THM-2445
comb cell. Likewise the constants in (14) and (17) cannot be improved
from linear conservation alone: abstract complex components may share
the current equally.

### 5.4 The exact residual

THM-2365 restores the delayed word on the generic drifting and
`D(O)>0` routes, while THM-2442 removes delayed-word cancellation for
the source-completed ghost. For any resulting nonzero word-retaining
current, and unconditionally with word `1` in every typed cell, this
theorem reduces the missing-right-endpoint problem to:

```text
complete cell-matched cospan
or
one typed two-sided factor transition.                (19)
```

To close the graft, one still needs semantic word/repair alignment on
the generic typed cells and then either:

- a same-root positive physicalization that turns the first line of
  (19) into a canonical owner/terminal current; or
- a lawful target/repair action or polarized reference that converts
  the second line into such a current.

THM-2383, THM-2401, THM-2408, and THM-2419 describe four versions of
that missing sidecar. None follows merely from the atlas. No scalar
profile is removed and LRC(14) remains open.

## 6. Exact companion

Run

```text
python 04-computation/lrc14_right_endpoint_cospan_transition_thm2448.py
python -O 04-computation/lrc14_right_endpoint_cospan_transition_thm2448.py
```

The companion:

- enumerates all `24^2=576` ordered endpoint-label pairs;
- verifies the matched, guard/unit, and blocker censuses (9);
- proves the recurrence and all denominators in (11)--(14);
- enumerates the twelve-cell ghost atlas (16);
- checks the `1/24`, `1/69`, and `1/12` sharp triangle controls and a
  zero-matched-term hostile;
- uses explicit `require` checks which remain active under `-O`.

The normal, optimized, and stored transcripts are byte-identical.
LF-normalized SHA-256 hashes are recorded in the frontmatter.

## 7. Independent audit

Two independent hostile audits reconstructed the theorem from its
dependencies. They checked:

- the lawful target co-shift on both present and bare copies;
- whole-layer Poisson--Abel linearity and the iterated `m`-then-`X`
  boundary passage;
- preservation of `X,m,chi`, delayed word, and ghost source colour;
- the `24`-cell census, tail recurrence, `69`-piece bound, and
  twelve-cell ghost;
- the distinction between a table-level diagonal zero and a scalar
  fixed-current coefficient;
- the matched-label/common-filter boundary and zero-matched-term
  hostile.

Both audits replayed the normal, optimized, and stored companion and
verified the LF hashes. No unresolved proof defect remains.

QED.
