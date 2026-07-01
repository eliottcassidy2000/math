---
id: HYP-3813
title: THE COVERING-MIN PHASE CLOUD = AP(step n) + ANTIPODAL KILLER -- one geometric picture unifying the Chebyshev 2-point dual, the three-gap theorem, the phase-residue, the runner-cloud tournament, and the Phi6-irreducibility. In phase-residue coordinates (S68, p(v)=n*v mod Phi6), the covering-min construction's runner cloud at the binding t* is EXACTLY the ARITHMETIC PROGRESSION {n*k mod Phi6 : k=1..n-2} (step n, the small runners 1..n-2) TOGETHER WITH the KILLER at p(n(n-1)) = -n mod Phi6. The observer (0) sits in the gap between the killer (-n) and runner 1 (+n): clearance = n/Phi6 = M_C on BOTH sides. VERIFIED (n=14): cloud = {14,28,...,168} + {169}; three-gap sizes {1, n, 2n} = {1,14,28} (three-gap theorem, HYP-3762); observer clearance 14/183 = M_C; FLANKING cloud points {14,169}={+n,-n} = runners {1, killer} = the Chebyshev 2-point alternation (S73); the killer is the iota-antipode (-n) of the slowest runner (+n), symmetrizing the gap. The cloud's ROTATIONAL tournament is NEAR-REGULAR (scores mostly (n-2)/2..., one +1 one -1; c3=90 at n=14; H odd by Redei). BEATER-OBSTRUCTION in cloud terms: covering forces the small speeds {1..n-2}, whose t*-phases are the AP of step n, leaving the size-n clearance at 0; so the observer's gap cannot shrink below n at t* => M >= n/Phi6. This is the S79 Phi6-metric-irreducibility as a CLOUD fact (the AP tiles Z/Phi6 with step n). UNIFIES S68(phase-residue)+S73(Chebyshev 2-pt)+S79(Phi6-irreducible)+S70(runner-cloud tournament)+HYP-3762(three-gap) into ONE picture: covering-min = an AP phase cloud with an antipodal killer, observer in the size-2n gap, Chebyshev-flanked at +-n
status: MIXED (exact n=14 + established-structure synthesis; general-n via S68). VERIFIED (n=14, exact): cloud = AP(step 14) {14..168} + killer 169=-14; three-gap {1,14,28}={1,n,2n}; clearance 14/183=M_C; flanking {14,169}={+-n}={runner 1, killer}=Chebyshev 2-point (S73); rotational tournament scores near-regular (one 7, one 5, rest 6), c3=90; H odd (Redei, not brute-forced at 13 vertices). GENERAL-n: p(k)=nk mod Phi6 makes the small speeds an AP of step n (S68), so the cloud = AP(step n)+killer(-n), clearance n/Phi6, three-gap {1,n,2n}, for all n. HONEST: a UNIFYING geometric reframe (the cloud picture) integrating five prior threads + a clean 'why n/Phi6' (AP-step n); not a new bound (the no-beater is S79/OPEN-Q-108). Only n=14 tournament invariants computed.
source: klein-2026-07-01-S80
depends_on:
  - HYP-3812   # S79: Phi6-metric-irreducibility (the cloud fact: AP tiles Z/Phi6 with step n)
  - HYP-3800   # S68: phase-residue p(v)=nv mod Phi6 (the cloud coordinates)
related:
  - HYP-3806   # S73: Chebyshev 2-point dual {1, killer} = the flanking cloud points at +-n
  - HYP-3762   # three-gap theorem: the cloud gaps {1, n, 2n}
  - HYP-3802   # S70: runner-cloud tournament (the rotational tournament of the cloud)
  - HYP-3715   # t*=n/Phi6 binding; the killer identity n(n-1)=-1 mod Phi6 (=> p(killer)=-n)
  - THM-523    # covering condition (forces the small speeds => the AP)
results:
  - 04-computation/phase_cloud_tournament_integration_klein.py
  - 05-knowledge/results/phase_cloud_tournament_integration_klein.out
---

# HYP-3813 — the covering-min phase cloud = AP(step n) + antipodal killer (one unifying picture)

## The picture
In phase-residue coordinates (S68: `p(v) = n*v mod Phi6`), the covering-min construction's runner cloud at
the binding `t*` is
> **the arithmetic progression `{n*k mod Phi6 : k = 1..n-2}` (step `n`, the small runners) + the KILLER at
> `p(n(n-1)) = -n mod Phi6`.**
The observer `0` sits in the gap between the killer (`-n`) and runner `1` (`+n`), with **clearance `n/Phi6 =
M_C` on both sides**. VERIFIED (n=14): cloud `= {14,28,...,168} + {169}`, clearance `14/183 = M_C`.

## Five threads, one picture
- **Chebyshev 2-point dual (S73/HYP-3806)**: the two cloud points **flanking** the observer's gap are
  `{+n, -n} = {runner 1, killer}` — exactly the length-2 alternation. The killer is the `iota`-antipode
  (`-n`) of the slowest runner (`+n`), symmetrizing the gap.
- **Three-gap theorem (HYP-3762)**: the cloud's gaps are the three sizes **`{1, n, 2n}`** (`= {1,14,28}`):
  `n` between AP points, `1` between runner `n-2` and the killer, `2n` for the observer's (doubled) gap.
- **Phase-residue (S68/HYP-3800)**: the coordinates; `p(k)=nk` makes the small speeds an AP of step `n`.
- **Phi6-irreducibility (S79/HYP-3812)**: the AP of the small speeds **tiles `Z/Phi6` with step `n`**, so
  the observer's gap is `2n` and the clearance `n` — irreducibly at the composite `Phi6 = n^2-n+1`.
- **Runner-cloud tournament (S70/HYP-3802)**: the cloud's rotational tournament is **near-regular**
  (n=14 scores mostly `6`, one `7`, one `5`; `c3 = 90`; `H` odd by Rédei), the AP giving a circulant-like
  order perturbed by the killer.

## Why `n/Phi6` (the clean mechanism, in cloud terms)
The covering condition (THM-523) forces the small speeds `{1..n-2}` (multiples of `q = 2..n-2`). Their
`t*`-phases are `{n, 2n, ..., (n-2)n} mod Phi6` — an **arithmetic progression of step `n`**. This AP tiles
`Z/Phi6`, and the gap containing the observer `0` has clearance exactly `n` (the nearest AP point is runner
`1` at `+n`). So `M(S) >= min-clearance = n/Phi6` at `t*`. The killer `n(n-1) = -n mod Phi6` adds the
antipodal point that **symmetrizes** the gap (making the two-sided Chebyshev alternation) and supplies the
covering multiples of `q = n-1, n`. So: **covering => the AP(step n) => clearance `n` => `M >= n/Phi6`**, and
the extremizer is the AP + antipodal killer.

## Net
The covering-min extremizer, in phase coordinates, is an **arithmetic-progression phase cloud (step `n`)
with an antipodal killer**, the observer in a size-`2n` gap, flanked by the Chebyshev 2-point alternation
`{+n, -n} = {runner 1, killer}`, with three-gap sizes `{1, n, 2n}` and a near-regular cloud tournament. This
single geometric picture integrates the phase-residue (S68), the Chebyshev dual (S73), the
`Phi6`-irreducibility (S79), the runner-cloud tournament (S70), and the three-gap theorem (HYP-3762) — and
makes the `n/Phi6` value transparent (the AP step is `n`). A unifying reframe, not a new bound.
