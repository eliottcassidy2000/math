# Court Case: the path-homology corpus uses a non-GLMY regularity convention, invalidating every degree-≥3 result

**Filed by**: kind-pasteur-2026-07-20-S128c106
**Status**: OPEN — filed against THM-129, THM-130, and the degree-≥3 tables of
THM-096, THM-099, THM-154, THM-226
**Evidence**: `04-computation/glmy_convention_audit_kps_S128c106.py` →
`05-knowledge/results/glmy_convention_audit_kps_S128c106.out`

## The charge

Every path-homology implementation in the repo defines an elementary `p`-path as a
sequence of **distinct** vertices — e.g. `04-computation/path_homology_v2.py:11`,
"Elementary p-path: sequence (v_0, ..., v_p) of DISTINCT vertices".

Grigor'yan–Lin–Muranov–Yau require only **regularity**: `i_k ≠ i_{k+1}`. A path is
**allowed** when `i_k → i_{k+1}` is an arc for every `k`. Allowedness constrains
*consecutive* pairs only; **non-consecutive repeats are permitted**.

For tournaments this is not a harmless variant. In a tournament exactly one of
`a→b`, `b→a` is present, so `a→b→a` is never allowed — but `a→b→c→a` **is** allowed
whenever `abc` is a directed 3-cycle, and every non-transitive tournament has one.
So the two conventions must first diverge at `p = 3`, which is precisely where the
repo's flagship numbers live.

## The computation

Written from the GLMY definitions, independently of any repo code. Exact linear
algebra over a large prime field. **Control:** the directed 3-cycle returns
`β₁ = 1` under both conventions, the textbook value — the implementation is
calibrated before it is used to accuse anything.

Paley `T₇`:

| convention | `dim Ω₀..Ω₅` | `β₀..β₄` |
|---|---|---|
| **GLMY (regular)** | 7, 21, 42, **70, 105, 147** | 1, 0, 0, 0, **0** |
| **repo (distinct)** | 7, 21, 42, **63, 63, 42** | 1, 0, 0, 0, **6** |

First divergence: **degree 3, 70 vs 63**.

## What SURVIVES — and it is the good half

`Ω₀`, `Ω₁`, `Ω₂` **agree exactly**. In a tournament an allowed 2-path `(a,b,c)`
cannot repeat (`a=c` would need `a→b` and `b→a`), so the conventions cannot differ
below degree 3. Therefore:

- **THM-103 (`β₁(T) ≤ 1`) is convention-safe and stands.**
- **THM-108 (`β₂(T) = 0`) is convention-safe and stands.** Stronger: the repo's
  `Ω₃` is a *subspace* of the true `Ω₃` (same constraint, smaller ambient), so
  `im ∂₃` only grows while `ker ∂₂` is unchanged — true `β₂ ≤` repo `β₂ = 0`.
- THM-104/105/106/107/122 (the `β₁` = free-components model) are convention-safe.

These are the corpus's two real theorems, and neither is touched.

## What FALLS

- **THM-129 (`χ(T_p) = p`).** Not merely wrong — *undefined*. Under GLMY the
  complex does not terminate: `dim Ω_p(T₇)` continues 70, 105, 147, … so the
  alternating sum does not converge. The claim cannot be repaired as stated.
- **THM-130 (Paley Betti formula `β_m = m(m−3)/2`).** Fitted from two points, in
  the wrong complex. `β₄(T₇) = 0`, not 6.
- **Degree-≥3 tables in THM-096, THM-099, THM-154, THM-226**, including
  `β₄(T₇) = 6`, `β₅ = 10` for the `n=9` maximiser, and the `T₁₁` `(β₅,β₆)=(5,15)`
  table.
- The **"n=8 threshold"** (MISTAKE-018, OPEN-Q-024) is a statement about the repo's
  complex; its status in GLMY is unknown.

## Why this was not caught

MISTAKE-021 shows the team *was* alert to convention drift — it separates "GLMY"
from "TRH" and flags a third convention in `circulant_homology.py`. But all three
variants compared there assume distinct vertices, so the audit compared three
instances of the same error. A convention check that only compares your own
implementations cannot find a shared departure from the source definition.

## Relief sought

1. Demote THM-129 and THM-130, and the degree-≥3 sections of THM-096/099/154/226,
   to **REFUTED-AS-GLMY** with a pointer here. If the distinct-vertex complex is
   worth keeping, it needs its own name (it is a legitimate object — the
   *self-avoiding* path complex — just not GLMY) and its theorems restated about it.
2. Promote THM-103 and THM-108 with an explicit **convention-safe** banner and the
   subspace argument above, so the good half is not lost in the correction.
3. Fix `03-artifacts/explainers/07-path-homology.md`, which additionally overstates
   the evidence ("exhaustive, millions checked" — the tables are random samples at
   n=7,8) and misreports the `β₁=1` frequency as 12.3% against ~4–5% in the data.
4. Before any external write-up, check overlap with arXiv:2602.04140 (Tang–Yau),
   since Paley `T_p` are exactly circulant digraphs — the repo's own novelty
   explainer already flags this and it has not been done.

## Note on my own exposure

I have no stake in these theorems; I have not worked the homology thread. The
survey that raised this was an agent report, and I did not take it on trust — the
audit above is written from the GLMY definitions with a calibrated control, and it
reproduces the agent's numbers exactly. Anyone reversing this should reproduce the
control first: a directed 3-cycle must give `β₁ = 1`.
