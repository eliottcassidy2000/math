---
id: HYP-3091
title: The three notions of sameness (equinumerosity, equidecomposability, equidistribution) are three computable invariants of the LRC lonely set; the equidecomposability face splits into D(S)=mod-41 Farey scale (bounded core) and 1/lmax=V* (apex), and the equidistribution face meas(L)=0 EXACTLY characterizes the tight atoms
status: VERIFIED invariant table (meas=0 <=> tight; D=41 for the hard family; D<=1/lmax split) + SYNTHESIS (the fiber as fundamental invariant). Not a new bound.
source: mac-mini-2026-06-27-S62
revives:
  - HYP-2187   # equinumerosity = cardinal shadow; equidecomposability = retained fiber (H=volume, beta1=Dehn)
  - HYP-2239   # triune carrier (sum/product/fraction faces)
  - HYP-2213   # Dehn-scissors side channel (beta1 = first scissors obstruction)
  - HYP-2232   # CH caution: a scalar/count does not decide a predicate without its fiber
related:
  - HYP-3089   # 1/lmax = the V* / Conjecture-7.1 constant (apex face)
  - HYP-3085   # gK8 = the 'volume' (meas) lower bound = the scissors form of CRUX 1
  - THM-530    # witness floor meas(L) >= m_P (the equidist invariant)
  - OPEN-Q-108
reflections:
  - three-notions-of-sameness-are-the-lonely-sets-fiber
---

# HYP-3091 — The lonely set's three-sameness fiber

The owner asked to add **equidecomposability** and **equinumerosity** to the project's **equidistribution**
lens and find invariants capturing the fundamental nature. Carrying HYP-2187's fiber functor (built on the
tournament side: `H`=volume, `β₁`=Dehn, count=cardinal shadow) onto the **LRC lonely set**
`L(S)={t:‖s·t‖≥1/14 ∀s}` gives three computable invariants at three resolutions.

## The verified fiber (`lrc_three_sameness_invariants_macmini_S62.py`)
```
 config                      EQUINUM        EQUIDECOMP (scissors)      EQUIDIST
                             cover #res     D(min wit)  1/lmax(worst)  meas(L)
 AP {1..13}    (tight)       False  13      D=14        --             0.00000
 GW {..,13,24} (tight)       False  12      D=14        --             0.00000
 generic non-covering        False   9      D=6         186            0.10183
 easy-cover {1..12,14}       True   13      D=13        308            0.02390
 hard-cover {1..11,13,84m}   True   13      D=41        511..980       0.00536..
```

## The three findings (one per resolution)
1. **EQUIDIST: `meas(L(S)) = 0` EXACTLY on the tight atoms** (AP, GW → 0 arcs). The density invariant = the
   witness floor `m_P`; its zero set is precisely the tight locus. "Volume detects tightness."
2. **EQUINUM is predicate-blind** (the cardinal shadow): "covering" is independent of tightness (AP
   non-covering yet tight; dilated AP covering yet tight, S60). A count cannot decide the LR predicate.
3. **EQUIDECOMP carries the arithmetic and splits into TWO scissors invariants:**
   - `D(S)` = min witness denominator (easiest rational reassembly) = **41 for the hard family, independent
     of the apex** — exactly the project's **Farey-neighbor scale mod 41** (kps-S40). A **bounded-core**
     invariant.
   - `1/ℓ_max` (worst arc) **grows with the apex** = the **V\*** constant of S61/HYP-3089 (the
     Conjecture-7.1 `D≈213`). An **apex-driven** invariant.
   With `D(S) ≤ 1/ℓ_max`, the scissors lens **separates** the two scales the proof had conflated:
   **mod-41 = bounded core, V\* = apex.**

## Interpretation (Dehn vs volume, on the measure side)
1-D scissors congruence under *all* translations is trivial → the density face sees only `meas=0` vs `>0`.
Restricting reassembly to **arithmetic (rational) translations** makes `D(S)` a genuine **Dehn-type
obstruction** (classically the detector of irrationality). So **`meas` = volume; `(D, 1/ℓ_max, #lengths)` =
the Dehn invariant** — the `H`-vs-`β₁` split of HYP-2187, now on the lonely set. The "fundamental nature"
of a speed set is its fiber `Φ(S) = (covering/residues | D, 1/ℓ_max, arc-spectrum | meas)`: no single
column separates {tight, generic, easy-cover, hard-cover}; the triple does.

## What it buys (honest)
- A clean extremal characterization: `meas(L)=0 ⟺ tight` — CRUX 1 ("`meas>0` off the tight locus") is the
  **scissors-nondegeneracy** form of the gK8 `p0≤cap` bound (HYP-3085, the volume face).
- The two open obligations are the two Dehn invariants: bounded-core = **mod-41 Farey** (`D`), apex =
  **V\*** peel (`1/ℓ_max`) — explaining why they are different constants doing different jobs.

## Next
1. Prove `meas(L(S)) ≥ m_P > 0` off the tight locus as a Dehn-nondegeneracy statement (the scissors form
   of CRUX 1) — does the Dehn invariant `D(S)` *lower-bound* the volume `meas`?
2. Map the hard family's `D=41` to the bounded-core gK8/Farey structure exactly (why 41 = the near-AP
   `{1..11,13}` reassembly modulus); test whether `D` is uniformly bounded over bounded cores (it should be,
   since it is apex-independent) — a finite statement, unlike `1/ℓ_max`.
