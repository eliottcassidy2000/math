# Starting to formalize DvdK: the positive-coefficient case, and where cancellation is the crux

*boxeph-2026-07-21-S228. Owner: aim earnestly at formalizing DvdK; make it simpler / circumvent it; let the
thoughts spill over to LRC. Builds on S226/S227 (the two-charge DvdK Lean), codex THM-2067 (the general
DvdK1, Galois orbit-product), and adopts codex MISTAKE-238 (my S227 LRC descent retracted). New Lean file
`04-computation/lean/TournamentH7/TournamentH7/GMC2DvdKPositive.lean` builds kernel-pure.*

## What DvdK1 actually is, and where the difficulty lives

`DvdK1` (the one-variable input to GMC(2)): for a two-sided Laurent polynomial `f=Σc_iz^{q_i}` (a `+`
charge and a `−` charge, distinct exponents, nonzero coefficients), `CT(f^m)≠0` for some `m≥1`. The exact
Lean statement expands the constant term as a sum over **balanced compositions**:
`CT(f^m)=Σ_{Σr=m, Σq_ir_i=0} multinomial(r)·∏c_i^{r_i}`. Two observations locate the whole difficulty:

1. **The multinomial coefficients are positive.** As a *polynomial in the coefficients `c`*, `CT(f^m)` has
   all-nonnegative coefficients; it is nonzero (as a polynomial) exactly when a balanced composition of size
   `m` exists — an elementary *feasibility* fact (two-sided ⟹ a balanced composition exists).
2. **The only hard part is cancellation.** For *specific complex* `c`, a nonzero polynomial can evaluate to
   `0` (cancellation). This — not feasibility — is DvdK's content, and it is exactly what killed my S222/S223
   saddle bypass (the counterexample `f=u²+u+u⁻¹−u⁻²`: aperiodic, cofinite feasibility, yet `CT(f^m)=0` for
   every odd `m`, only `CT(f⁴)=−12≠0`). Codex's **THM-2067** (Galois orbit-product) resolves the cancellation
   in general.

## Formalized: the positive-coefficient case (no cancellation), any support

The corollary of (1)+(2): with **positive real coefficients** there is *no cancellation* — the constant term
is a sum of positive terms — so `DvdK1` is elementary for *any* support. I formalized this, kernel-pure
(`#print axioms = [propext, Classical.choice, Quot.sound]`):

- `ct_pos_of_balanced` — for `c_i>0` and *any* given balanced composition `r₀` of size `m`, `CT(f^m)>0`
  (the sum of nonnegative terms containing one strictly-positive term);
- `exists_balanced_of_twosided` — a two-sided support (`q_i>0`, `q_j<0`) always admits a balanced
  composition: `|q_j|` copies of the `+` charge, `|q_i|` of the `−`, at `m=|q_i|+|q_j|`;
- `dvdk1_positive` — **the positive-coefficient `DvdK1`**: two-sided support + `c_i>0` ⟹ `∃m≥1, CT(f^m)>0`.

This is a genuine start on formalizing DvdK: the entire **no-cancellation regime** of `DvdK1` is now proved
inside Lean, for arbitrary support and arbitrary number of charges — DvdK-premise-free. It complements the
S226/S227 two-charge case (which is DvdK-premise-free even for *complex* coefficients, because with only two
charges the balanced composition is unique, so cancellation cannot occur there either).

## Making it simpler / circumventing it — the honest map

The clean structural picture the formalization exposes:

```
   DvdK1  =  "two-sided => CT(f^m) != 0 for some m"
        │
        ├── FEASIBILITY (a balanced composition exists)      -- elementary (exists_balanced_of_twosided)
        ├── POSITIVE coefficients (no cancellation)          -- elementary (dvdk1_positive)   [FORMALIZED]
        ├── TWO charges (unique composition, complex ok)     -- elementary (S226/S227)         [FORMALIZED]
        └── >=3 charges, complex coefficients (CANCELLATION) -- the crux = codex THM-2067 (Galois orbit-product)
```

So "making DvdK simpler" = isolating that **cancellation is the sole difficulty**, and everything else is
elementary and now formalized. A genuine *circumvention* would reduce the complex case to the positive case
— but the retracted S222/S223 shows feasibility/Zariski-density does **not** do this (the vanishing locus of
`{CT(f^m)}` can meet the complex torus even when the positive orthant is safe); the Galois-orbit structure
(THM-2067) is what actually rules it out. The right next Lean target is therefore THM-2067's orbit-product,
sitting on top of the elementary base now in hand.

## Spillover to LRC: the cancellation parallel

The same "positive vs. cancellation" split is the honest read on why LRC is hard. The LRC covering measure
`|G_δ|=Σ_{k·v=0}∏ĝ(k_j)` (THM-515) is a lattice sum with **signed** (sinc) weights `ĝ` — it is *all*
cancellation, with **no positive regime**: unlike GMC's angular part, there is no coefficient sign choice
that makes loneliness a sum of positive terms (covering *is* the sign cancellation). That is exactly why the
`χ`/topological criterion (codex THM-2047, my S212) is needed — a measure/volume argument is cancellation-
blind (my S211 volume-ceiling), and the sign structure must be controlled by *symmetry*: the mirror `ι`
(S212) and the doubling homeomorphism (codex THM-2075) are the LRC analogue of THM-2067's Galois orbit for
DvdK — structure that tames sign cancellation. (Honest: my S227 attempt to compose these into a `χ=0`
terminal-core descent is **retracted** — codex MISTAKE-238: the tower transports the nonempty *core* sets,
not the full-set emptiness, and THM-2077 gives the terminal core a safe *interval*, so `χ>0` there. The
doubling identity, the homeomorphism, and mirror-parity survive individually; they do not compose to that
reduction.)

## Scope

Concrete progress on formalizing DvdK: the positive-coefficient `DvdK1` (any support, no cancellation) is now
kernel-pure in Lean, isolating cancellation as the sole remaining difficulty (handled on paper by THM-2067).
Not the general complex `DvdK1` (that is THM-2067's Galois orbit-product, the next Lean target). The LRC
spillover is a structural parallel (cancellation is the crux in both; symmetry tames it), with the S227
over-reach honestly retracted.

Links: HYP-8925, THM-2067, THM-2022, THM-2075,
[[a-kernel-pure-lean-proof-of-the-two-charge-dvdk-seed-boxeph-S226]],
[[where-gmc2-reaches-lrc14-the-ct-functional-and-where-it-stops-the-volume-ceiling-boxeph-S211]].
