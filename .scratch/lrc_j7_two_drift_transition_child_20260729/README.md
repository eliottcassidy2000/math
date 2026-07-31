# LRC(14) critical seven-tail aligned/drift transition packet

**Status (2026-07-29): scratch audit and handoff.**  The exactly
five-aligned/two-drift face is now proved empty in THM-2941, with THM-2928
providing an independent address-mask proof.  The four-aligned/three-drift
face is not closed, but ordered typed recursion gives `z_2<=6515` and the
subsequent exact scalar screen sharpens this to the necessary bound
`z_2<=2163`.  Nothing here proves LRC(14).

## 1. Inheritance and live concept board

- Closest proved mechanisms: THM-2928's aligned Boolean tensorization and
  address masks; THM-2941's safe-surplus/projected-wall reduction; and
  THM-2893's literal-residual first-apex gate.
- Canonical hostile boundary: a measure-theoretic endpoint handoff does not
  cover the strict-open seam (MISTAKE-146 and THM-2954).
- Corrected near miss: singleton excess and Gram data erase the body-cell
  address.  They give finite caps but do not decide the cover.
- Least-used sidecars: the ordered forbidden aligned prefix in the apex tree,
  and the valuation/residual-unit block of an endpoint link.
- Current board: aligned multiplicity surplus; projected De Morgan residual;
  ordered typed apex; body-by-phase incidence; strict-open seams; endpoint
  link blocks; three-drift owner words.

The useful recurring move is to keep three representations of one cover at
once:

```text
scalar:    drift singleton excess pays aligned multiplicity surplus;
projected: non-full drift columns must fit inside an aligned comb union;
recursive: a low residual apex is typed as drift or least aligned label.
```

Each representation forgets information retained by the other two.

## 2. General aligned safe surplus and projected wall

Let `E` be six body speeds, let

```text
C_E = T \ union_(e in E) D_e,     h=mu(C_E),
L=14 lcm(E),
```

and suppose the seven tails contain `k` aligned speeds `La_i` and
`d=7-k` ordered drifts `z_1<...<z_d`.  If `u_k` is the proved lower bound
for the safe mass of the `k` normalized aligned multipliers, integration of
the aligned multiplicity gives

```text
Delta = integral_(C_E)(m-1) >= (u_k-d/7)h = eta_k h.   (2.1)
```

The exact constants currently used are

```text
k        2       3         4          5        6
u_k    66/91   55/91    558/1183   478/1365  61/273
eta_k   1/91    3/91     51/1183    88/1365  22/273.
```

For `k=5`, the same inequality has a Gram derivation:

```text
(m_A-1)_+ >= (2/5) binom(m_A,2),
sum_(i<j) mu(D_(a_i) cap D_(a_j)) >= 44/273,
```

so `Delta>=88h/1365`.  This identifies the safe-surplus and pair-Gram
pressures as two coordinate descriptions of the same phenomenon.

For a carrier with `r` interval components, the correct singleton
discrepancy estimate is

```text
delta_E(z):=mu(C_E cap D_z)-h/7 <= 6r/(49z).            (2.2)
```

The earlier scratch coefficient `99r/490` was superseded.  Since

```text
Delta=sum_q delta_E(z_q),
```

ordering the drifts gives the rootwise first-drift cap

```text
z_1 <= floor(d(6r/49)/(eta_k h)).                       (2.3)
```

There is also a stronger geometric wall.  Delete the aligned tails and use
settled smaller-`n` LRC on the body-plus-drift partial family.  Project the
resulting safe arc by `t -> Lt mod1`.  Compact containment in the proper
open aligned union gives

```text
max(E union Z) > alpha_k L,
alpha_k = k/[7(14-k)(1-u_k)].                           (2.4)
```

Explicitly,

```text
k       2        3          4            5          6
alpha  13/150   13/132   2366/21875   2275/18627  819/5936.
```

This projected statement keeps content even when a full-cell argument is
too coarse.

## 3. The projected residual and its seam semantics

For selected body cells `j` and drift set `Z`, put

```text
E_z(j)={u in T: ||z(j+u)/L||<1/14},
P_(E,Z)=union_j (T \ union_(z in Z) E_z(j)).
```

Equivalently,

```text
T\P_(E,Z)=intersection_j union_(z in Z) E_z(j).         (3.1)
```

If the original cover exists, every phase in `P_(E,Z)` must be covered by
the normalized aligned union `U_A`.  Hence

```text
P_(E,Z) subset U_A,
mu(P_(E,Z)) < mu(U_A) <= 1-u_k.                         (3.2)
```

The inequality is strict because `P` is compact, `U_A` is a proper open
set, and compact containment has positive distance from its complement.

The closure convention in the exact interval programs is safe here:
closures are used only for mass and component arithmetic.  Deleting a
closed representative produces a positive-mass subset of the literal
open-set residual, so every literal cover still covers the tested subset.
No endpoint handoff is promoted into pointwise coverage.  Body-grid
boundary points project to phase zero, which lies in every normalized
aligned danger comb.  This is the seam distinction required by
MISTAKE-146.

The incidence form exposes the duality with THM-2928.  Let

```text
I(j,u)=1  iff the drifts cover phase u in cell j.
```

THM-2928 fixes `u` in the aligned safe set and forbids a full column
`I(j,u)=1` for all `j`.  The projected argument uses

```text
P={u: exists j with I(j,u)=0}
```

and proves that the set of non-full columns is too large for `U_A`.

## 4. Exactly five aligned tails: proved empty

The canonical projected-suffix bank has `4,702` rows `(E,z_1)`.

### 4.1 High-excess bank

There are `4,084` rows with

```text
delta_E(z_1)>=88h/1365.
```

Delete `D_(z_1)`.  If the residual has mass `h_1` and `r_1`
components, THM-2893 forces one of the six remaining labels below

```text
A_1=floor(36r_1/(7h_1)).
```

The aligned labels are at least `L`, while the projected wall forces
`z_2>2275L/18627`.  These typed floors close `3,827` rows.  The remaining
`257` rows yield exactly `42,912` nonaligned integer `z_2` candidates.
Exact cell projection kills every candidate by

```text
mu(P_(E,{z_1,z_2})) >= 887/1365 = 1-u_5.
```

An independent aligned recursion closes `39,913` candidates at its first
gate and all remaining `2,999` before the one-label terminal; no aligned
multiplier above four is used.

### 4.2 Subcritical bank

For

```text
g=88h/1365-delta_E(z_1)>0,
```

equation `(2.2)` gives

```text
z_2 <= floor((6r/49)/g).
```

There are `2,290` analytic first rows, `618` of which survive the canonical
projected-suffix predicate.  The analytic intervals contain `7,218,110`
row-labelled `z_2` candidates.  Exact singleton integration and the suffix
predicate leave `194,073` pairs on `590` `(E,z_1)` rows.  Exact projected
residual measure kills all of them.  The closest margin is
`1/378105`.

Thus the projected proof closes both banks.  THM-2928 independently closes
the same face by body/divisor address masks:

```text
no six-body/seven-tail cover has five aligned tails and two drifts.
```

## 5. Four aligned tails: ordered typed-apex compression

THM-2941's canonical projected-suffix filter leaves `87,975` candidate
rows on all `3,003` body roots, with `15<=z_1<=182`.  After deleting
`z_1`, four aligned labels and the ordered drifts `z_2<z_3` remain.

At a state with `s` already deleted aligned labels, let `R_s` have mass
`h_s` and `r_s` components.  There are `p=6-s` remaining tails.  The
first-apex cap is

```text
A_s=floor(6p r_s/[7(7-p)h_s]).                          (5.1)
```

The exact floor uses a small but load-bearing non-strict integer variant of
THM-2893.  If

```text
c_R(w)<=h_R/7+gamma/w
```

and `p<=6` integer labels cover `R`, put

```text
T=7p gamma/[(7-p)h_R].
```

If every label exceeded `floor(T)`, integrality would give `w>T`, whence

```text
sum c_R(w)<p(h_R/7+gamma/T)=h_R,
```

a contradiction.  Thus some label is at most `floor(T)`.  Taking
`gamma=6r_s/49` proves `(5.1)` directly.  Replacing the coefficient by
`6r_s/49+1` to fit only the strict form of THM-2893 proves finiteness but
does **not** justify this exact numeric cap.

The coefficients for `p=6,5,4,3,2` are

```text
36/7, 15/7, 8/7, 9/14, 12/35.                          (5.2)
```

Maintain:

```text
B_s       current strict lower bound for z_2;
m_s       last deleted aligned multiplier (m_0=0).
```

The recursion is exhaustive and disjoint:

1. `B_s<=z_2<=A_s` is a bounded-drift branch.
2. On `z_2>A_s`, also `z_3>A_s`.  The low apex must therefore be the
   least remaining aligned label `Lm`.
3. The forbidden-prefix sidecar forces
   `m_s<m<=floor(A_s/L)`.  Delete `D_(Lm)`, set
   `B_(s+1)=max(B_s,A_s+1)`, and recurse.
4. After four aligned deletions, the complementary branch has only the two
   drifts, both above the apex, and is impossible.

Allowing a multiplier below `m_s` would be a safe over-approximation but
would lose the claimed ordered/disjoint tree.  The concurrent scout was
corrected to enforce this invariant before its result was accepted.

The exact ordered run reports

```text
suffix rows / roots                 87,975 / 3,003
initial aligned rows / states          706 / 758
recursive states by depth       758, 1052, 1427, 2131
bounded branches by depth         758, 990, 1293, 543
closed complements by depth           0, 54, 684, 2131
uniform necessary z_2 cap                         6,515
raw nonaligned (E,z_1,z_2) bank              50,285,016.
```

Every row's union of typed `z_2` intervals merges to one interval.  These
are exact **necessary-filter candidate** counts, not realized-cover counts:
merging deliberately forgets which aligned prefix generated a value.

On the direct `z_2<=A_1` branch, the scalar surplus splits `49,129,013`
candidates into

```text
scalar impossible         40,630,772
positive gap, finite z_3   2,032,220
no positive-gap z_3 cap    6,466,021.
```

Screening the **full merged recursive bank**, including candidates reached
after aligned-first branches, gives

```text
raw candidates                         50,285,016
scalar impossible                      41,770,842
positive gap, finite z_3 cap            2,042,669
no positive-gap z_3 cap                 6,471,505
necessary screened bank                 8,514,174
screened uniform z_2 cap                     2,163.
```

All `87,975` first rows retain at least one pair, so this is a compression,
not a row closure.  In the last class, the exact first-two-drift excess
satisfies

```text
delta_E(z_1)+delta_E(z_2)>=eta_4 h.
```

This means only that the universal positive-discrepancy estimate no longer
bounds `z_3`; a negative exact `delta_E(z_3)` may still violate the total
surplus condition.  The finite-gap class should next be expanded with its
explicit `z_3` caps, while the no-positive-gap class needs projected or
literal-residual structure.

The separate full Cartesian superset run uses every six-body root and every
nonaligned `15<=z_1<=182`.  Its first stage contains `504,503` rows:
`498,981` force `z_2<=1888`, while `5,522` rows on `91` roots permit an
aligned first apex, with multiplier-cap histogram

```text
((1,4416),(2,739),(3,278),(4,26),(5,30),(6,28),(7,5)).
```

Its complete ordered recursion was still running when this handoff was
written; it is an independent superset check, not a dependency of the
canonical `6515` bound.

## 6. Endpoint-link sidecar for three drifts

For two speeds `u,v`, exact exit/entry abutment of strict danger teeth
requires

```text
14 gcd(u,v) | u+v.                                     (6.1)
```

THM-2954 rewrites this as:

```text
v_2(u)=v_2(v),  v_7(u)=v_7(v),
unit(u) == -unit(v) mod14.                             (6.2)
```

Thus endpoint links form disjoint complete bipartite blocks indexed by the
valuation pair and one of

```text
{1,13}, {3,11}, {5,9}.
```

This exactly classifies the exceptional zero residue in the older
two-clock overlap formula.  It does **not** permit a zero-width transition
in a pointwise strict-open cover: the common endpoint is safe.  Every
actual owner switch must pay a positive overlap quantum.  For three drifts,
no three types can be pairwise zero-linked, so any owner pattern realizing
all three pair adjacencies must pay a positive transition on at least one
edge.  This is a useful typed sidecar for the post-`z_2` bank, not yet a
global exclusion.

## 7. Reproduction and negative results

Primary scratch programs in this directory:

```text
residual_first_apex_audit.py
high_excess_recursive_aligned_closure.py
subcritical_exact_pair_bank.py
subcritical_projected_pair_filter.py
k4_typed_first_apex_scout.py
k4_recursive_typed_apex_cap_scout.py
```

Run exact controls both normally and with assertions disabled:

```bash
python3 SCRIPT.py
python3 -O SCRIPT.py
```

The canonical `k=4` ordered run currently lives at

```text
../lrc_z130_survivor_attack_20260729/
  k4_exact_suffix_typed_stage_scout.py
  k4_exact_suffix_typed_stage_scout.discovery.out
```

`subcritical_projected_recursive_audit.py` is an abandoned broad
Fraction-based recursion.  It was stopped without a verdict after the
branching proved uninformative.  It is not evidence for closure.

## 8. Orthogonal completion of the four-aligned/three-drift face

The typed-apex recursion above remains useful as a finite description, but it
is no longer the shortest closure of this face.  The independently developed
THM-2928 divisor-fibre route now gives the exact necessary chain

```text
544,571 all-top occurrences
 -> 419,511 all-divisor fibre survivors
 -> 29,364 common-phase/parity survivors
 -> 19 eight-status threshold-transport survivors
 -> 0 local unit-needle survivors.
```

The main `29,364 -> 19` transport is sound because an actual joint high/low
status table is feasible for the real eight-cell polytope, while the real LP
is an upper relaxation.  At a target threshold `y`, a covered fibre with load
at least `y` must have raw status capacity **at least** `y`; changing this
inclusive comparison to a strict one rejects a literal positive control.

In a `q`-fibre, put `g=gcd(d,q)` and `n=d/g`.  A global
`ell=ceil(d/7)`-term unit AP restricts to at most

```text
ceil(ell/g)=ceil(n/7)
```

terms of a unit AP modulo `n`: solvability selects one AP-index class modulo
`g`, and division by `g` followed by inversion of `q/g mod n` preserves a unit
step.  Granting an independent maximal local AP to each drift in each fibre is
therefore an upper relaxation.

The scratch candidate

```text
lrc14_three_drift_terminal_independent_audit_candidate.py
```

contains no payload or Lorenz dependency.  It embeds and semantic-hashes the
`19` identities, reconstructs their literal `D=L` body slices, and uses an
incidence/pair-cover algorithm different from the primary grouped recursion.
It also exhausts the generic transport and trace laws:

```text
needle/quotient cases                         692,416
literal fibre traces                       14,279,768
monotone status words                              20
LP versus integral controls                    10,520
fractional LP optima                               17
maximum LP/integral gap                           1/2
positive literal threshold unions                252
terminal local kills                              19
terminal local survivors                           0.
```

Ordinary and optimized outputs are byte-identical.  At the audited primary
source state their hashes and runtimes were

```text
audit source  7e31fb0582f36a21872865acf78d59bdf4e5fc2067b06bd34ea62476bdd9397a
audit output  f934e30114a7a73d71db99a25ff35ad5848463c88217b4d1d92a90f22eadcbfb
ordinary      real 84.10, user 77.93, sys 0.60
optimized     real 83.66, user 77.68, sys 0.54
```

The old independent Lorenz ledger is only a sidecar.  Its earlier `29,221`
count missed two double-denominator-two rows because it compared the whole
support with twice a coarse capacity.  The repaired test compares one exact
parity load with the third mask's exact parity-fibre cap and gives `29,219`.
Neither Lorenz count is a dependency of the `29,364 -> 19 -> 0` proof.
