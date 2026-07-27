---
id: THM-2591
title: "Theta-zero selector Cech coboundary and C91 holonomy no-go"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Every admissible choice of a positive THM-2586 theta-zero diagonal edge
  over the seven owner-clock cells is a vertex-valued root selector and hence
  changes a seven-chart transition law only by a Cech coboundary.  The exact
  THM-2586 zero sets leave 128 choices for each displacement other than 6,7
  and 16 choices for each of 6,7, for 1,312 selectors total.  All have zero
  selector holonomy.  By contrast THM-2542's constant nonzero root-deck step
  a has holonomy 7a != 0 in F_13, unchanged by every selector gauge and by
  every cyclic chart reidentification.  Thus the positive arrival/deep/later
  diagonal cannot by chartwise root selection trivialize the C91 local
  system.  A successful comparison needs a genuine transition-dependent
  mixed-square correction of cyclic sum -7a, or the known thirteen-fold
  trivializing cover.  The theorem does not identify the two clock carriers,
  construct that mixed 2-cell, produce a target/relation current, remove a
  row, or prove LRC(14).
source: codex-2026-07-28-deep-cech
depends_on:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2586-depth-five-arrival-to-future-root-diagonal
related:
  - THM-2551-horizontal-transfer-transverse-projector-bicomplex-boundary
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2587-deep-root-coshift-incidence-wall-and-theta-selector-no-go
script: 04-computation/lrc14_theta_selector_cech_holonomy_thm2591.py
output: 05-knowledge/results/lrc14_theta_selector_cech_holonomy_thm2591.out
script_sha256: 49c1ecca2de0e68ca96c2a72faf8bf11d35e345ab8492f55a25b1b4d2d927a77
output_sha256: 6fee5dee756786fb2e586362c3b7e688f0f3f0968e373c61ad91ec8eebfb468b
hash_basis: LF-normalized bytes
---

# THM-2591 -- a root selector cannot kill the seven-chart holonomy

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2586 puts a positive current-arrival/deep/later-root diagonal in every
one of the `12*7` nonzero-displacement/owner-clock cells of one canonical
typed packet.  THM-2542 puts a nontrivial root-deck local system over a
seven-cycle of transported charts.  It is tempting to use the root selected
in each positive THM-2586 cell as the missing chart origin.

That cannot work.  The reason is cohomological but elementary: a root chosen
at each chart is a zero-cochain, so it changes transition data only by an
exact one-cochain.  THM-2542's class has nonzero cyclic sum.  The exact
THM-2586 positivity restrictions do not create an exception.

The comparison below grants an identification of the two seven-cycles.  No
such physical identification is asserted by either dependency.  Proving
failure even after granting it makes the no-go stronger, not a construction
of the missing map.

## 1. The complete theta-zero selector bank

Fix a collision displacement

```text
s in F_13^*.                                               (1)
```

For `ell in F_7`, THM-2586 retains the two theta-zero rail edges

```text
e_0=(v,t)=(0,0),              e_6=(v,t)=(6,12),             (2)
```

for which the rescaled deepest root is

```text
w=7t=v.                                                     (3)
```

Its exact edge-zero sets are

```text
Z_0={(7,4),(7,5),(7,6)},
Z_6={(6,4),(6,5),(6,6)}.                                   (4)
```

Define the admissible root set

```text
A_(s,ell)
 ={0 : (s,ell) notin Z_0} union {6 : (s,ell) notin Z_6}.    (5)
```

A theta-zero selector at displacement `s` is any function

```text
h:F_7 -> F_13,                 h_ell in A_(s,ell).           (6)
```

The two sets in (4) are disjoint, so every fibre in (5) is nonempty.  More
precisely:

```text
s notin {6,7}:  A_(s,ell)={0,6} for all seven ell;
s=6:            h_4=h_5=h_6=0, four other values free;
s=7:            h_4=h_5=h_6=6, four other values free.      (7)
```

Thus the exact selector counts are

```text
N_s=2^7=128                    for ten values of s,
N_6=N_7=2^4=16,

sum_(s in F_13^*)N_s=10*128+2*16=1312.                     (8)
```

THM-2586's deterministic `81/3` priority selector is only one member of
this bank.  The theorem below treats every positive choice simultaneously.

## 2. The seven-chart class and its gauge law

Fix independently a nonzero THM-2542 marker root

```text
a in F_13^*.                                               (9)
```

The symbols `a` and `s` have different types and are not identified.  On
the oriented chart cycle, with edge `ell` directed from `ell-1` to `ell`,
THM-2542's transition one-cochain is

```text
g_ell=a.                                                    (10)
```

Changing the root origin in chart `ell` by `h_ell` changes it to

```text
g^h_ell=g_ell+h_ell-h_(ell-1)
       =a+h_ell-h_(ell-1).                                 (11)
```

Using the opposite sign convention for the selected root replaces `h` by
`-h` and changes nothing below.  Summing (11) around the cycle gives

```text
Hol(g^h)
 =sum_(ell in F_7)g^h_ell
 =7a+sum_ell(h_ell-h_(ell-1))
 =7a != 0                  in F_13.                         (12)
```

The selector contribution telescopes because it is the coboundary `dh` of
a globally defined vertex labelling.  Equation (12) holds for every
`h in F_13^7`, hence in particular for all 1,312 positive selectors in
(8).  Consequently

```text
there is no admissible h with g^h_ell=0 for every ell.       (13)
```

Rotation of the proposed clock identification merely reindexes (12).
Reflection reverses the cycle and changes the class to `-7a`, which is
still nonzero.  Thus no dihedral choice of chart origin or clock orientation
repairs (13).

This also explains why allowing the other THM-2584 rail does not by itself
solve the problem.  Any arrival root, deepest root, or rail displacement
chosen separately at each chart is still vertex data.  Its change along
chart edges is another coboundary and has cyclic sum zero.  What is missing
must be attached to transitions, not only to vertices.

## 3. The exact mixed-square boundary invoice

Suppose a proposed lawful comparison supplies an additional transition
correction

```text
c=(c_ell)_(ell in F_7) in C^1(C_7^graph;F_13).              (14)
```

Here `c_ell` is meant to be the boundary contribution of a mixed square
coupling chart transport to the semantic arrival map.  Flattening the
combined transition would require

```text
a+h_ell-h_(ell-1)+c_ell=0             for every ell.        (15)
```

Summing (15) yields the necessary invoice

```text
sum_ell c_ell=-7a.                                          (16)
```

This condition is gauge invariant.  In particular, another chartwise
factor, common multiplier, selected root, or vertex sidecar changes `c` only
by a coboundary and cannot alter its sum.  Equation (16) is the cheapest
decisive test for the next proposed bridge:

```text
source:    one common Boolean carrier with chart and arrival data;
map:       its transition-dependent mixed-square boundary c;
preserve:  positivity, current-arrival=later-root typing, chart covariance;
must hit:  Hol(c)=-7a;
may lose:  a preselected edge/frequency, but not the physical common carrier.
                                                                    (17)
```

A candidate failing the fourth line of (17) is another gauge and cannot
close the seam.  Passing it is necessary, not sufficient: the entries of
`c` must still come from one lawful physical construction rather than seven
independently spliced coefficients.

THM-2542 gives the other exact option.  On an `n`-fold clock cover the pulled
class has holonomy

```text
n*7a,                                                        (18)
```

which vanishes exactly when `13` divides `n`.  The minimal trivializing
cover therefore still has degree thirteen.  A selector coboundary on that
cover does not lower the degree.

## 4. Relation to the bicomplex boundary

THM-2551 proves on its abstract product carrier that horizontal C91 transfer
commutes with the all-unit target projector and preserves its kernel
exactly.  Equation (12) is the finite chart analogue of the same stopping
phenomenon: a product/horizontal operation cannot manufacture the missing
transverse class.

This is a comparison of mechanisms, not an additional dependency.  The
next positive object must be genuinely mixed.  Algebraically it must have a
nonzero square commutator; geometrically it must pay (16).

## 5. Scope and exact gain

The proved implication is

```text
positive THM-2586 diagonal in every owner-clock cell
  + any admissible chartwise root selection
  + even a granted cyclic identification with THM-2542 charts
  -> only a Cech coboundary;

THM-2542 nonzero root-deck class
  + that coboundary
  -> the same nonzero holonomy 7a.                           (19)
```

This closes one tempting shortcut at the post-THM-2586 frontier.  It does
not:

1. identify THM-2586's owner-clock with THM-2542's scheduler clock;
2. identify either carrier with THM-2569's selected old head;
3. rule out a direct old-head map outside the THM-2542 local system;
4. construct the correction `c`, a THM-2334 relation residue, or a
   THM-2365 target co-shift;
5. exclude the canonical typed row or prove LRC(14).

The remaining question is narrower than "choose a better root."  Construct
one lawful transition-dependent mixed square on the same positive ancestry
and test its seven-entry boundary against (16).

## 6. Exact companion

Run

```text
python3 04-computation/lrc14_theta_selector_cech_holonomy_thm2591.py
python3 -O 04-computation/lrc14_theta_selector_cech_holonomy_thm2591.py
```

The dependency-free companion uses only the two exact zero sets (4).  It
enumerates all 1,312 admissible selectors, checks every selector coboundary,
checks all `1,312*12=15,744` marker/selector pairs and boundary invoices,
checks all 220,416 dihedral relabellings, prints the permutation
`a -> 7a` of `F_13^*`, and verifies that the minimal trivializing cover has
degree thirteen for every nonzero marker.

The stored transcript is

```text
05-knowledge/results/lrc14_theta_selector_cech_holonomy_thm2591.out.
```

Independent hostile audit is still required before promotion.

QED (candidate).
