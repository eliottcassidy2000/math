---
id: THM-2565
title: "Target-active self-return and future-root overlap"
status: >
  PROVED + VERIFIED-EXACT.  On every one of the 165 live rows, resolve a
  positive THM-2559 target-informed head packet by its physical first
  base-thirteen digit and reuse that same resolved Boolean field after a
  sufficiently long forward iterate.  Exact BV mixing then gives an actual
  positive old-head/later-root diagonal.  If rho is the head mass, its
  diagonal mass is at least rho^2/26 for all sufficiently large delays; one
  delay works for all 165 rows.  The future copy contains literal k_a danger,
  and putting k_a first in the freely ordered five-role bank types it as the
  unique target-active first-failure label.  This is a cemetery-free literal
  THM-2545 table for the target-informed selector.  Its word is retained only
  as the source-sibling stratum label: the head does not emit or inherit that
  terminal word.  The selector is not THM-2537's canonical t_tau(e), the
  paired blocker is not target-co-shifted, and no THM-2334 current, scalar-row
  exclusion, or LRC(14) conclusion follows.
source: lrc-semantic-frontier-2026-07-27-target-active-self-return
depends_on:
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
  - THM-2559-target-informed-chord-and-universal-old-repair-packet
related:
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
  - THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary
script: 04-computation/lrc14_target_active_self_return_thm2565.py
output: 05-knowledge/results/lrc14_target_active_self_return_thm2565.out
script_sha256: 4e11fc2c59a720870518f2e64f16efb61c5cabe33726d3999f8cb217bb2fff27
output_sha256: 865a36517db738a922ca1373f9f770f3ad35d5da4a3bdc67c75fa66a39089bf5
hash_basis: working-tree bytes (LF)
---

# THM-2565 -- a target-active head eventually returns to its own root

**PROVED + VERIFIED-EXACT.**

THM-2555 isolates two missing inputs in the old arrival construction: a
genuinely later target-active Boolean field and overlap of its root marginal
with the selected-head marginal.  THM-2559 subsequently constructs a positive
target-active head field on every live row.  The cheapest future field is not
an independently labelled copy.  It is the **same physical field evaluated
later**.

That stationary choice removes the offset gauge.  If the root-resolved masses
are `p_h`, both endpoint marginals are the same vector, so their limiting
diagonal overlap is

```text
sum_h p_h^2>0.                                               (1)
```

Exact mixing turns (1) into a positive same-ancestry return at every
sufficiently large delay.  The result closes one real semantic seam, but only
in THM-2559's target-informed head chart.

## 1. The positive root-resolved field

Fix one of the `165` live rows.  Choose a THM-2559 nonzero slope stratum
`delta` with positive target-informed head packet

```text
A=T_delta,                   rho=integral_T A(x)dx>0.        (2)
```

Thus, with `z=Tx`, every point of `A` has the form

```text
x=iota_(t(z))(z)=(z+t(z))/13,                               (3)
```

where `t(z)` is the selected marker in the literal failure mask of `k_a`.
In particular

```text
A(x)=1  implies  d_(L_a)(k_a x)=1.                          (4)
```

Resolve `A` by the physical first digit

```text
A_h(x)=A(x)1_(floor(13x)=h),        h in F_13,

p_h=integral_T A_h(x)dx.                                    (5)
```

The fields `A_h` are disjoint rational Boolean BV functions and

```text
A=sum_h A_h,              rho=sum_h p_h.                    (6)
```

Equation (3) makes the root typing literal:

```text
A_h(x)=1  implies  h=t(Tx)=floor(13x).                       (7)
```

No ancestry-sheet quotient or abstract relabelling occurs in (7).

## 2. Exact self-mixing forces a diagonal

For `N>=1`, form the actual root-to-root table before integration,

```text
C_N(h,b)=integral_T A_h(x)A_b(T^N x)dx,                     (8)

H_N=sum_h C_N(h,h).                                         (9)
```

The second factor in (8) is the same field, not an independently gauged
future packet.  At a point counted by `C_N(h,b)`, its root is therefore

```text
b=floor(13T^N x),                                           (10)
```

the immediate physical digit at the later horizon.

Apply THM-2555 equation (28) with `G_h=A_h`.  Put

```text
c_h=min(p_h(1-p_h),Var(A_h)/12),

P=sum_h p_h^2,
E=sum_h c_h Var(A_h).                                       (11)
```

Then

```text
|C_N(h,h)-p_h^2| <= c_h Var(A_h)/13^N,

H_N >= P-E/13^N.                                            (12)
```

Since `rho>0`, Cauchy--Schwarz gives

```text
P>=rho^2/13>0.                                              (13)
```

If `E=0`, (12) already gives `H_N=P` for every `N`.  Otherwise, every
delay satisfying

```text
13^N>=2E/P                                                  (14)
```

obeys the explicit floor

```text
H_N>=P/2>=rho^2/26>0.                                      (15)
```

Increase `N`, if necessary, beyond every clock used in the definition of
`A`.  The second factor is then genuinely later.  Equations (7), (9), and
(10) show that every diagonal atom has

```text
old selected-head root h = future immediate root b.         (16)
```

There are finitely many rows and thirteen root pieces on each row.  Taking
the maximum of their finite BV and chronology thresholds gives one common
delay for all `165` rows.  The positive slope `delta`, field `A`, mass `rho`,
and numerical floor remain row-dependent.

## 3. Why the future copy is categorically target-active

THM-2445 allows the five roles other than the fixed graft `q_*=k_b` to be
placed in any fixed first-failure order.  Put `k_a` first.  Its first cell is
then simply

```text
S_1=d_(L_a)(k_a x),                                         (17)
```

with no earlier five-bank safety factor.  By (4), every occurrence of the
future field `A_b(T^N x)` lies in that `k_a`-labelled first-failure cell.
THM-2461 proves that `k_a` is the unique target-active member of this
five-role bank.  Thus the second endpoint in (8) is not a target-neutral
future owner with a root label attached after the fact: it contains a
literal, categorically target-active unit/guard failure at time `N`.

This use of (17) types the categorical role only.  It does not assert that
the separately fixed graft factor is safe, restore THM-2445's partial-bare
current, or co-shift the blocker paired with `k_a`.

The root gauge is equally concrete.  Replacing the future labels by
`b+c` would define a different family of fields.  In (8), both endpoint
families use the same formula (5), the same base-thirteen origin, the same
orientation convention, and the same role `k_a`.  A simultaneous change of
root gauge transports both labels and preserves (16); an offset applied only
to the future endpoint is not a symmetry of the constructed packet.

## 4. The terminal word is a source-side label, not a head word

Fix the terminal word `sigma` used in THM-2559.  On its active direct fibre,

```text
1_(Q_(j,sigma))(13^k iota_r(z))
 =1_(Q_(j,sigma))(13^(k-1)z)                                (18)
```

is root-constant.  The existence of the occupied sibling
`S_delta(iota_(s(z))(z))` therefore proves that the literal factor in (18)
equals one at the head root as well.  The sibling `S_delta<=W` is the genuine
`E_j` source carrying terminal word `sigma`; the head `A=T_delta` is disjoint
from `W` and is not itself an `E_j` source.

Accordingly, `sigma` may be retained on (8) as the fixed source-sibling word
stratum from which the head was selected.  The late owner factor `g(z)` is
also retained in the definition of `A`, and the future copy contains the
corresponding later factor.  But neither literal membership in (18) nor reuse
of the label permits the statement

```text
the selected head emits or inherits terminal word sigma.    (19)
```

This distinction is exactly THM-2537's source/head warning.

## 5. Exact Hall consequence

Restrict to the positive packet

```text
Omega_N={x:A(x)A(T^N x)=1}.                                 (20)
```

On `Omega_N`, use the constant source-side word stratum `sigma`, the old
head map

```text
h(x)=floor(13x),                                            (21)
```

and the genuinely later target-active root map

```text
b(x)=floor(13T^N x).                                        (22)
```

There is no cemetery value on (20), because the future target-active field
is present by definition.  Its THM-2545 table is literally (8), restricted
to (20), and (15) says

```text
sum_h nu({x in Omega_N:h(x)=b(x)=h})=H_N>0.                 (23)
```

Thus THM-2545's primitive known-table criterion
`C^sigma_(h,h)>0` is met: this is an actual positive Hall diagonal, not an
inference from separate margins.  The abstract theorem permits any selected
empty-head map, so (20)--(23) are a literal instance of it after replacing
the canonical selector by THM-2559's target-informed selector.  They do not
prove that the margins force every admissible coupling to have a diagonal;
no Hall-deficient set or overload is asserted.

THM-2538's cross-Kakeya reconstruction can therefore be applied to this
same-base packet if one wants its mixed root coefficients.  Positivity is
already visible in (23); reconstruction is not needed to prove it.

## 6. The offset-one hostile is finite-delay, not stationary

THM-2555's offset-one cylinder correctly shows that chronology, full sheet
retention, and overlapping marginals do not force a diagonal at a prescribed
short delay.  It does not contradict eventual mixing of fixed fields.

For an exact three-digit control, write `d_1,d_2,d_3 in F_13` and define

```text
U_h={d_1=h, d_2=h+1},              V_h={d_1=h}.              (24)
```

Both root marginals have common support.  At delay one,

```text
U_h(x)V_h(Tx)=0                    for every h,              (25)
```

because (24) would require `d_2=h+1` and `d_2=h`.  At delay two,

```text
sum_h mu(U_h(x)V_h(T^2x))
 =13/13^3=1/169>0.                                         (26)
```

The same fields have gone from a perfect offset coupling to a positive
diagonal once the inspected digit blocks separate.  General rational BV
fields need not separate after two digits, but (12) proves the same eventual
statement quantitatively.  The hostile forbids a zero-delay shortcut; it does
not provide a persistent obstruction to the stationary self-return.

## 7. Exact stopping boundary

The theorem changes the semantic ledger as follows:

```text
target-informed old head + independently named future field   obsolete;

target-informed old head + its stationary later copy
  -> positive same-root target-active return                   proved. (27)
```

It does **not** prove THM-2537 equation (56).  That equation uses the
canonical head `t_tau(e)` and the weight `g Psi_tau(e)`, whereas THM-2559's
marker is selected directly from the `k_a` failure mask and may have another
slope.  Nor does (23) make the head emit a terminal word, produce a lawful
THM-2334 target current, co-shift the blocker paired with `k_a`, restore the
full-X endpoint in THM-2563, or turn the nonnegative return into a scalar-cover
contradiction.

Consequently no row is excluded and LRC(14) remains open.  The next
high-leverage composition is to retain (23)'s stationary future-root
diagonal while completing THM-2563's missing fixed-head target residue and
paired blocker action.

## 8. Exact companion

Run

```bash
python3 04-computation/lrc14_target_active_self_return_thm2565.py
python3 -O 04-computation/lrc14_target_active_self_return_thm2565.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_target_active_self_return_thm2565.out.
```

The dependency-free referee checks the exact thirteen-root Cauchy floor, the
first-role truth table, the three-digit offset-one hostile and its delay-two
return, and the finite-delay covariance threshold algebra.  The all-row
conclusion is the proved BV argument above, not an extrapolation from those
finite controls.

**QED.**
