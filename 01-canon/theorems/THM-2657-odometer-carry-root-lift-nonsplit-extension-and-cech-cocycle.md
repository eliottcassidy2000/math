---
id: THM-2657
title: "Odometer carry/root lift, nonsplit extension, and Cech cocycle"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Every circle translation which globally lifts a nonzero THM-2640
  predecessor-carry/root increment delta has the form tau=k/13^6 with
  k=7 delta mod 13.  Every such k is a thirteenth-unit, so tau has order
  13^6, not 13.  Hence the physical translation extension
  0 -> C_(13^5) -> C_(13^6) -> C_13 -> 0 is nonsplit.  The minimal
  set-theoretic section k(delta)=7 delta has exact wrap cocycle
  7 floor((delta_1+delta_2)/13)/13^5, representing the nonzero class 7
  modulo 13.  This is an odometer/Čech obstruction before any present-factor
  mismatch; it supplies no positive common-component relation or LRC exit.
source: carry-transition-cell-2026-07-28-odometer-clutch
depends_on:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
related:
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2642-cyclic-difference-relation-saturation-and-thick-holotopy-no-go
script: 04-computation/lrc14_odometer_carry_root_lift_cech_thm2657.py
output: 05-knowledge/results/lrc14_odometer_carry_root_lift_cech_thm2657.out
script_sha256: a835ae3ab0516572ad4596471a1a01ccc59839b426fbc91b4caff89a67991e9f
output_sha256: 41398b3f621c29efe83b29aae5a5bae0732406aa203c2e795eb57f9758bc6f9a
hash_basis: LF-normalized bytes
---

# THM-2657 -- the physical slope-seven clutch is a nonsplit odometer

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2640 finds the unique formal clutch slope `c->c+7 delta` compatible with
`r->r+delta`, and exhibits the smallest physical translation realizing it
on the carry/root/delayed-prefix subpacket.  The failure is stronger than its
speed-one guard hostile: no physical translation lift can have order
thirteen.  The carry coordinate is the first digit of a length-six odometer,
and the desired `C_13` action meets a nonsplit cyclic extension.

## 1. Classification of every physical translation lift

Put

```text
p=13,              R=p^6,              S=R/p=p^5,
c3=2S.                                                     (1)
```

For `x` away from the finite carry boundaries, let

```text
c(x)=floor(Rx) mod p.                                      (2)
```

Fix a nonzero desired root increment `delta in F_p^*`.  Suppose a circle
translation `x->x+tau` moves `(2)` globally by the compatible carry increment
`7 delta`.  Write

```text
R tau=n+alpha,                 n in Z, 0<=alpha<1.          (3)
```

If `alpha!=0`, then as `{Rx}` crosses `1-alpha`, the difference

```text
floor(R(x+tau))-floor(Rx)
```

jumps from `n` to `n+1`.  It cannot be one constant residue modulo `p` on
the whole circle.  Therefore `alpha=0`; every global lift has

```text
tau=k/R,                       k in Z/RZ.                   (4)
```

Its carry increment is `k mod p`.  The high-speed probe moves by

```text
c3 tau=2k/p                    mod 1.                       (5)
```

Thus it sends the absolute root `r` to `r+delta` exactly when

```text
2k=delta mod p
 <==> k=7 delta mod p.                                    (6)
```

Equation `(6)` simultaneously gives the THM-2640 carry slope.  It classifies
all lifts, not just the minimal representative: for each `delta` there are
exactly

```text
S=p^5                                                       (7)
```

choices `k=7delta+p j`, `j in Z/SZ`.

## 2. The extension is nonsplit

For `delta!=0`, every `k` in `(6)` is a `p`-unit.  Hence the order of its
circle translation is

```text
ord(k/R)=R/gcd(k,R)=p^6.                                   (8)
```

In particular its thirteenth iterate is the nontrivial fine translation

```text
13 tau=k/p^5,                                              (9)
```

of order `p^5`.  No lift of a nonzero quotient element has order dividing
`p`.

Equivalently, consider the exact cyclic extension

```text
0 -> C_(p^5) -> C_(p^6) -> C_p -> 0,                      (10)
```

where the quotient map sends `k` to `2k mod p` and the kernel consists of
the multiples of `p`.  The `p`-torsion subgroup of `C_(p^6)` is the unique
subgroup

```text
{j p^5:0<=j<p}.                                           (11)
```

Every element of `(11)` is itself a multiple of `p`, so `(11)` lies inside
the kernel and maps trivially to the quotient.  A splitting would have to
map a quotient generator to a nontrivial element of `(11)`, which is
impossible.  This proves `(10)` is nonsplit.

## 3. The exact wrap cocycle

Choose representatives `delta in {0,...,12}` and the minimal set-theoretic
section

```text
s(delta)=7 delta in Z/RZ.                                  (12)
```

This is a section because `2s(delta)=delta mod13`, but it is not a
homomorphism.  Its defect is

```text
omega(delta_1,delta_2)
 =s(delta_1)+s(delta_2)-s((delta_1+delta_2) mod13)

 =91 floor((delta_1+delta_2)/13)             in Z/RZ.     (13)
```

As a physical circle translation, `(13)` is

```text
omega/R
 =7 floor((delta_1+delta_2)/13)/13^5.                     (14)
```

There are exactly `78` wrapping and `91` nonwrapping ordered pairs.  In the
kernel coordinate obtained by dividing `(13)` by `p`, a wrap has value `7`.
Modulo `p` this is the nonzero class

```text
[omega]=7 in H^2(C_p,C_(p^5))=C_(p^5)/p C_(p^5) ~= C_p.  (15)
```

Thus changing the set-theoretic section alters `(13)` by a coboundary but
cannot remove its class.  This is the exact Čech clutch behind the missing
physical `C_13` action.

## 4. Relation to the present-factor defect

For any lift `(4)` with `k` a `p`-unit, a present speed `v` acquires phase

```text
v tau=vk/R.                                                (16)
```

This lies on the lawful `1/p` phase grid exactly when

```text
S divides vk <==> S divides v.                            (17)
```

On the canonical speed tuple

```text
(1,14,27,40,53,66,13,13^3,2*13^5),                       (18)
```

only `c3=2*13^5` satisfies `(17)`.  THM-2640's guard defect is therefore the
first visible failure of a more general statement.  But `(8)`--`(15)` show
that even deleting every low-speed present factor would not create a
physical order-thirteen translation action: the odometer extension itself
does not split.

## 5. Holotopy boundary and next test

The theorem rules out a **translation action**, not a positive relation.
Distinct lifts may still have positive chart overlaps.  Those overlaps form
a Čech nerve over the delayed base.  Pairwise positivity fills only its
one-skeleton; lawful composition requires common higher simplices or a
balanced component/gain section.  THM-2642 cannot be invoked until those
same-base intersections are supplied, and THM-2637 still requires a fixed
branch rather than a thick relation.

Likewise, `(15)` does not exclude a signed two-clock cancellation.  The next
nonlinear test is the adjacent-clock `R,13R` response matrix suggested by
THM-2639: ask whether a determinant or resultant kills the `C_(13^5)` kernel
even though neither individual clock splits it.  THM-2624's signed-versus-
positive boundary remains load-bearing.

No common physical overlap component, positive composable relation,
fixed branch, owner/root endpoint, row exclusion, or LRC(14) conclusion
follows.

## 6. Exact companion

Run

```bash
python 04-computation/lrc14_odometer_carry_root_lift_cech_thm2657.py
python -O 04-computation/lrc14_odometer_carry_root_lift_cech_thm2657.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_odometer_carry_root_lift_cech_thm2657.out.
```

The companion exhausts all `13^6` grid translations, checks every carry/root
congruence and lift fibre, verifies the `p`-torsion nonsplitting boundary,
enumerates all `169` cocycle pairs, and rechecks `(17)` on `(18)`.

QED (candidate; independent hostile audit pending).
