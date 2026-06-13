# THM-465: the bi-dyadic algebra settles Q(3,5)/Q(3,6) per-size, refuses to be uniform, and is 2-adically special

**Status:** PROVED/VERIFIED as itemized (each computational claim: Glucose3 +
complete verifier + independent re-verification; scripts
`erdos592_{bidyadic_rule,uniform_rule}_macmini_s3.py`, triadic control inline run,
results in 05-knowledge/results/).
**Source:** mac-mini-2026-06-10-S1 (HYP-2374; POKE Steering Task 2; continues
THM-453/459/460)

For a pair x <lex y in [t]³ define the **bi-dyadic feature**
F2(x,y) = ( (sign, v₂(|gap|)) per coordinate, '=' for zero gaps ) — the order data
of S1's pattern algebra plus the 2-adic valuation of every coordinate gap.
A **feature rule** is a set T of feature classes; its graph: blue iff F2 ∈ T.

## A. Per-size: the frontier falls

* (3,4): sign-only algebra (13 classes) is **feature-UNSAT**; bi-dyadic (62
  classes) is SAT — re-proving S1's "patterns insufficient" in a sharper form and
  locating exactly what was missing: the 2-adic valuations.
* **(3,5): SAT — Q(3,5) SETTLED** by an explicit 42-class bi-dyadic witness
  (1272 edges), independently re-verified by a fresh complete verifier. The raw
  instance had timed out at 80,111 CEGAR clauses (S2); the feature quotient
  decides it in 2.8 s with 171 features. **R(3,2) > 5.**
* (3,6): SAT and verified the same way (4104 edges). **R(3,2) > 6.** The
  strong-witness frontier at ω³ (HYP-2363) is now three sizes past the n=2
  death point R(2,2)=5.

## B. No uniform table: the infinite bi-dyadic witness is REFUTED

An infinite F2-measurable strong witness on N³ would restrict to ONE fixed table
valid at every finite t (for t ≤ 8 the realizable feature set is constant since
all gaps ≤ 7 have v₂ ≤ 2). Computed:
* the frozen t=5 table FAILS at t=6 (domination) and t=7 (a triangle!);
* the simultaneous-t SAT (one shared table, constraints of t = 4,5,6,7 conjoined)
  is **feature-UNSAT in 0.3 s**.

Hence **no (sign, v₂)-measurable strong witness exists on ω³**, even though every
finite size carries one. The algebra ladder is strict:
  signs (fails per-size) ⊊ signs+v₂ (per-size only) ⊊ [whatever the infinite
  witness needs]. This mirrors, one level up, S1's "patterns < patterns+scheme"
  finding — each rung buys finitely many sizes, the infinite witness keeps
  escaping upward. (Constructive-strong-Specker route: must climb again —
  candidate next rungs: unbounded-v₂ tables with tail conditions, or
  Larson-scheme partial-sum features.)

## C. The seam is 2-adic (the triadic control)

Replacing v₂ by v₃ gives an algebra of the SAME SIZE at t=4 (62 classes; the gap
partitions {1,3}|{2} vs {1,2}|{3} are isomorphic as partitions) — and it is
**feature-UNSAT at (3,4) with zero CEGAR rounds** (triangle clauses + seeds alone
kill it). So the bi-dyadic success is not "any coarse refinement of signs works":
**v₂ is special at n=3**, consonant with the binary-subgrid structure (sibling
spines = halving) and the S2 dyadic-sufficiency result. CAVEAT (the n=2 control,
from the seam experiment): at n=2 ALL tested gradings (free/inv/v₂/v₃) share the
cutoff 5 — the discrimination only appears at n=3, where the column space has
depth. Why exactly v₂'s class structure is compatible with the Schur/composition
constraints while v₃'s is not is a sharp open note (HYP-2373 refinement; the
v₂=0/v₂≥1 split is the parity quotient ℤ→ℤ/2, whose "odd+odd=even" closure is
forced, unlike the mod-3 analogue).
