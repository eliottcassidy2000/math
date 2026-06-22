---
id: HYP-2920
title: LRC14 AP/Goddyn-Wong common-condition stack and top-petal audit
status: PROOF-TARGET / computational atlas; not a proof of LRC14
source: codex-2026-06-22-S124
related:
  - HYP-2913
  - HYP-2914
  - HYP-2917
  - HYP-2918
  - HYP-2919
  - HYP-2909
  - HYP-2910
  - HYP-2893
  - THM-523
  - THM-568
  - THM-569
---

# HYP-2920: AP/GW common-condition stack

This hypothesis records a necessary-condition atlas for the two known
primitive LRC14 tight rows:

- AP: `{1,...,13}`.
- Goddyn-Wong: `{1,...,11,13,24}`.

It is a companion to the S54 "kind battery" and the S38
divisibility-threshold / doubling-operad notes.  The point is not that the
last filters are proved.  The point is to split the AP/GW resemblance into
layers, so that every future non-AP/GW tight candidate must say exactly which
layer it violates.

## Layered necessary-condition program

The theorem-level conditions are:

1. **Punctured divisor cover:** a primitive tight row must contain a multiple
   of every `q=2,...,13` and no multiple of `14`.  Otherwise a `q`-grid witness
   gives `M(S)>1/14`, or the denominator-14 apex is blocked.  This is the
   q-threshold reread of THM-523 plus the apex-denominator work around
   THM-568/THM-569.
2. **Apex unit cover:** every unit residue modulo `14` must be represented.
   Otherwise some `a/14` with `a in (Z/14Z)^*` is strictly lonely.

The AP/GW fingerprint layers are stronger proof targets:

3. **Exact odd skeleton and shell profile:** the residues
   `1,3,5,7,9,11,13` occur exactly once, there are six unit residues, six even
   residues, one half-step residue, and no zero residue.
4. **Co-finite residue support and maximal zsum:** the nonzero residue support
   misses at most one class, and the unordered antipodal pair count
   `#{i<j : s_i+s_j=0 mod 14}` is maximal (`6`).  GW loses the AP pair
   `(2,12)` but regains the count via `(4,24)`.
5. **Literal AP complement binders:** the unit-shell binders are not merely
   congruent to zero modulo `14`; the representatives in residues
   `(1,13)`, `(3,11)`, `(5,9)` literally sum to `14`.
6. **Single even dipole:** the residue defect relative to AP is either zero or
   one missing even class plus one doubled even class.
7. **Minimal one-petal acceleration:** all representatives are the AP minimal
   ones except possibly one replacement `m -> 2m`.
8. **Top-petal gate:** the only tight petal seen is the top composite petal
   `12 -> 24`.  The other minimal petals `8 -> 16` and `10 -> 20` survive many
   coarse filters but are loose off apex.

The last two bullets are deliberately strong.  They should be read as the
sharpest current single-swap proof target, not as a theorem.

## Exact audit

The script
`04-computation/lrc14_ap_gw_common_conditions_codex_s124.py` stores its output
at `05-knowledge/results/lrc14_ap_gw_common_conditions_codex_s124.out`.

Known rows:

- AP and GW pass every filter and have exact `M=1/14`, denominator set `{14}`.
- The residue liar `{1,...,11,13,26}` has the same mod-14 residue multiset as
  AP and passes the residue-shell filters, but it fails the divisor threshold:
  `q(S)=12`, `M=1/12`, and the first failed stack condition is
  `punctured_q_cover`.  This is the executable version of the rule
  "residues lie; divisibility tells the truth."
- The near miss `{1,...,11,13,36}` passes through the single-even-dipole layer
  but fails the one-petal layer and escapes off apex with exact `M=3/41`;
  the output records the Farey check `det[[1,3],[14,41]]=-1`.

Bank results:

- AP single replacements `v<=300`: `3732` primitive rows.  The theorem-level
  q/unit filters leave `963` rows; the AP/GW fingerprint through zsum and
  literal complement sums leaves `395`; the loose one-petal layer leaves `99`;
  `minimal_one_petal` leaves `4`; the `top_petal_or_ap` gate leaves exactly
  AP and GW, with `tight=2`, `loose=0`, `below=0`.
- AP two replacements with values `<=40`: `27730` primitive rows.  The same
  stack leaves exactly AP and GW, again with `tight=2`, `loose=0`, `below=0`.
- Terminal shrink inside the minimal-petal premise: before the final
  `top_petal_or_ap` gate, both broad local banks already have the same
  four-row exact core, namely AP, GW, `8 -> 16`, and `10 -> 20`.  This core
  has maximum speed `24`; the replacement-ceiling sweep stabilizes at `24`
  (`16` sees two terminal rows, `20` sees three, and `24,30,40` all see four).
  Exact `M` on those four rows then leaves AP and GW tight while `8 -> 16`
  has `M=2/23` and `10 -> 20` has `M=2/27`.
- Primitive bounded `13`-subsets of `[1,19]`: `27132` rows.  Since `24` is
  outside the bank, the final stack leaves AP alone, with `tight=1`,
  `loose=0`, `below=0`.
- Minimal AP-doubling petals `v -> 2v` for `v=7,...,13`: only `v=12`
  has no coprime blocker in the Goddyn-Wong window `[14-v,27-2v]`, and it is
  the only row in the ledger with `M=1/14`.  The full exact list is:
  `7->14` has `q=15`, gate false, `M=1/11`; `8->16` has `q=14`, gate false,
  `M=2/23`; `9->18` has `q=14`, gate false, `M=2/23`; `10->20` has `q=14`,
  gate false, `M=2/27`; `11->22` has `q=14`, gate false, `M=2/25`;
  `12->24` has `q=14`, gate true, `M=1/14`; `13->26` has `q=14`, gate false,
  `M=2/27`.

Tournament Analysis on the condition carrier used the conditions as vertices,
pass-set inclusion over the single-swap bank as the pairwise observable, and a
listed Hamiltonian path to orient ties.  The resulting carrier is acyclic in
directed triples (`directed_3cycles=0`) but has many incomparability/tie flips,
which is exactly the warning: these filters are a hierarchy of lenses, not a
single monotone theorem.

## Interpretation

The useful abstraction is now:

```text
tight non-covering row
  => q=14 divisor-transversal and unit-cover
  => AP/GW residue skeleton plus maximal antipodal zsum
  => even dipole
  => top Goddyn-Wong petal 12 -> 24, or else off-apex escape.
```

The computed evidence says that the tempting broad condition "one petal" is
still too weak: it leaves many loose off-apex rows.  The sharper missing lemma
is a rigidity theorem for the petal height and site.  In the tested banks, the
only tight acceleration is the high-composite, divisor-compensating top petal
`12 -> 24`; `8 -> 16`, `10 -> 20`, and higher lifts such as `12 -> 36` all
escape to denominators other than `14`.

This completes the local AP-petal census inside the audited banks, including
the two-replacement bank `<=40`.  More strongly, once the proof target has
entered the minimal-one-petal branch, the exact terminal census is bounded by
maximum speed `24`; the larger banks only stress-test earlier necessary
filters and near-miss escape modes.  It does not prove the full LRC14
tight-locus census.  It converts the remaining proof into a stack of named
failure modes: a hypothetical third tight row must either break cofinite/zsum
rigidity, break literal complement binding, realize a multi-dipole/two-petal
structure, or produce a new top-balanced unbounded phenomenon not seen by
these AP/GW banks.
