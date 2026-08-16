# Three-times-p affine covers, residue automata, and harmonic rank spectra

**Research reflection / provenance, not a truth source.**  Exact claims are
routed to the audit-pending THM-3469 and its proved dependencies THM-3405 and
THM-3455.  LRC(14) remains open.

## Inheritance board

The closest proved mechanism is THM-3455's cap-seven atom sieve.  The canonical
hostile is the strict boundary `p=35`, where phase `x=85` survives all eight
owners.  The corrected near miss is THM-3464: its primitive mixed-order
`q=123` packet does not exclude an inherited rank-eight realization.  The
least-used decisive sidecar is the combination of odd-phase parity with the
order-three backbone; either one alone leaves the exceptional septimal gap.

THM-3469 turns the isolated `q=123` packet into

```text
q=3p,
U_p=(1,p-1,p+1,2p-1,2p,2p+1,3p-6,3p-1).
```

Its exact boundary is not mysterious numerical coverage.  After an odd phase
is written as `x=mp+y`, six affine channels erase the nearest multiple, owner
`2p` erases the residue class `3|x`, and owner `3p-6` folds the remaining
septimal gap.  The only failure has `p=7s`, `|y|=3s`, `s` odd, and `3` not
dividing `s`; equivalently `p mod 42` is `7` or `35`.

The lower bound needs a separate type bridge.  On every odd divisor `Q`, the
sheet permutation `phi(ell)=2ell+1` conjugates a fixed-zero mask of owner `s`
to the half-twist mask of canonical owner `2s`.  This is what licenses the
passage from THM-3455's literal half-layer atom sieve to the full
zero-cochain/divisor rank.  Without this coordinate, the eight-owner upper
bound would not by itself exclude a seven-owner fixed-zero escape.

## The finite-state object is not a tournament

The useful state has coordinates

```text
(x mod 3, nearest-multiple index m mod 6, gap depth |y|-h, parity).
```

The channel selected at one state is an owner-labelled hyperedge certificate.
There is no intrinsic binary orientation between two owners, so forcing the
eight owners into a tournament discards the union predicate.  A tournament of
size four or six can still be used as a diagnostic on channel competition,
but it must carry ties, missing edges, the orientation gauge, and the full
incidence sidecar.  The proof itself is closer to a small deterministic
automaton with one default channel, one residue-class override, and one
boundary repair.

This gives a precise version of the recurring `1+2k` motif.  The `1` is the
order-three fixed/backbone class, while the paired `2k` channels are the
`+/-` affine choices around nonzero residue classes.  Here `k=3`, producing
six nearest-multiple channels.  The extra septimal repair is not another
tournament vertex: it is a boundary event attached to those channels.

## Ternary trees and recurrence indices

The phase classifier branches first into `x mod 3=0,+1,-1`; the two nonzero
branches then carry six cyclic nearest-multiple states.  This is a genuine
ternary decision root, but iteration does not by itself produce Berggren's
ternary tree.  The Berggren connection enters only after intersecting the
affine family with the distinguished parabolic U-spine

```text
q_t=(2t+1)^2+2.
```

That intersection is exactly

```text
t mod 21 in {5,8,12,15}.
```

Thus the branch transplant is arithmetic: a finite residue automaton selects
four of every twenty-one U-spine indices.  It does not assert a Fibonacci
recurrence for the cover ranks.  The next positive target is to determine
whether other Berggren branches make the nearest-multiple channel table
periodic after a finite state lift, or whether new divisor coordinates are
required.

## Subsets of the natural numbers and the harmonic series

Every rank condition defines a subset of the family index `k`, and hence a
subseries of the harmonic series.  That observation alone gives no asymptotic
coefficient.  THM-3469 earns one because the exact rank classifier is periodic
with minimal period

```text
3*5*11*13*17*23*29=24,322,155.
```

For exact rank eight the residue census gives

```text
delta_8=57344/147407,
sum_(k<=N, rank(k)=8) 1/k = delta_8 log N+O(1).
```

Because the affine labels `q=42k-3` are linear, the same subset also has

```text
sum_(k<=N, rank(k)=8) 1/(42k-3)
  =(4096/442221) log N+O(1).
```

This separates three notions that are easy to blur: arbitrary subsets of
natural numbers, periodic subsets with logarithmic density, and reciprocal
sums over the quadratic labels `q_t`, which converge.

On the U-spine intersection, the lane indicator and the inherited rank word
combine to a minimal annotated period `11781`.  Its rank-4/6/7/8 counts are
`748,272,144,1080`; the ambient `t`-harmonic coefficients are respectively
`4/63,16/693,16/1309,120/1309`.  This is stronger than the bare lane
coefficient `4/21`: it records how the recurrence branch redistributes the
four exact grades.

## New frontiers opened by the template

1. **Affine-cover synthesis.**  Replace `3p` by `ap` and search for a minimal
   residue backbone plus paired nearest-multiple channels and boundary repairs.
   The invariant is a strict integer gap budget, not pairwise dominance.
2. **Ancestry classification within a grade.**  At `q=123`, inherited and
   primitive rank-eight realizations coexist.  Classify realization types as a
   monoid under divisor pullback and primitive affine insertion, rather than
   assigning one ancestry label to the rank.
3. **Current transplant.**  Feed the explicit `q=123` and `q=291` packets into
   the relation-residue lift.  Coverage survives residue-preserving owner
   lifts, but grouped endpoint coefficients may cancel; the cheapest decisive
   test retains the owner-by-sheet incidence matrix and evaluates the grouped
   target current exactly.
4. **Spectral closure.**  The affine template supplies a periodic positive
   family of zero-cochain covers, not the LRC `7 tensor 13` bispectrum.  A real
   bridge must map owner labels to current classes while preserving ancestry,
   target activity, and the nonvanishing spectral contraction.

The broader lesson is that residue automata can prove infinite cover families
and exact harmonic spectra, while tournaments and ternary words remain useful
only as typed shadows of the labelled incidence/current object.
