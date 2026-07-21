# The moment–nullcone template, and the complexity hierarchy it lives in

**kind-pasteur-2026-07-20-S128c128.** Companion to THM-1750. The owner asked for creative
unifying frames — to think abstractly and let past concepts connect. This is the frame that
kept surfacing once I looked for it: **the last several sessions of GMC/TNC "detection depth"
work are not a detour from the tournament project; they are the tournament project's own
`λ = 0 ⟺ transitive` seen one and two rungs up an arithmetic-complexity ladder.**

## The one sentence

> An invariant worth studying here is almost always **a projection of a power onto the trivial
> component of a symmetry** — a *moment* — and the interesting question is almost always about
> its **nullcone**, the objects all of whose moments vanish. Because the generating function of
> a moment sequence is recurrence-governed, the nullcone is cut out at a **finite depth**, and
> it is always the **most degenerate** object in sight.

Write `φ(Y) = ⟨Y⟩_G` for the projection onto `G`-invariants, `X` an object, and study
`F(t) = Σ_m φ(X^m) t^m`. Three of this project's pillars are this, verbatim.

## The ladder

| rung | object | `φ` | `X` | `F(t)` | nullcone | depth |
|---|---|---|---|---|---|---|
| **rational** | tournament | trace | adjacency `A` | `−t (log det(I−tA))'` | transitive | `n` |
| **algebraic** | TNC | circle average `CT` | Laurent `Λ` | diagonal of `1/(1−tΛ)` | one-sided | `D` |
| **holonomic** | GMC(2) | Gaussian `E` | `P(Z,Z̄)` | Laplace of the diagonal | charge-one-sided | `K` |

The three rungs are the three floors of `rational ⊂ algebraic ⊂ D-finite`. Nothing about the
ladder is metaphorical: **Cayley–Hamilton is the rational special case of the holonomic
recurrence.** A matrix is a period you can sum in rank-one pieces, `Σ tr(A^m)t^m =
Σ_i λ_i t/(1−λ_i t)`; a diagonal of a rational function is algebraic; a period integral is
holonomic. The recurrence order — `n`, then `D`, then `K` — is just how high `F` sits on the
ladder. The transitive tournament and the one-sided Gaussian are the *same theorem* evaluated
at two different arithmetic complexities.

I find this clarifying about my own recent work. "TNC has detection depth `D`" (THM-1710) and
"GMC on bounded charge is a finite Gröbner test" (THM-1740) felt like new machinery. They are
the machinery the project has had since `λ = 0 ⟺ transitive` (THM-895) — Newton's identities
and Cayley–Hamilton — pushed off the rational floor.

## The nullcone is the degenerate object

Every one of these nullcones is the "boundary" case, the object with the least structure:

- **transitive** — the unique acyclic tournament order, `F ≡ 0`, spectrum `{0}`;
- **one-sided** — a Laurent polynomial that doesn't straddle `0`, `F ≡ 1`;
- **charge-one-sided** — a Gaussian polynomial whose charges don't straddle `0`.

"All moments vanish" is a *rigidity*: it forces `X` to the extreme of its family. This is the
same shape as the project's other rigidity results — the tight AP in LRC is the extreme
configuration, the regular tournament is the spectral extreme (`ρ = (n−1)/2`, THM-1555), the
Paley/BIBD extremizers. The frame suggests reading each of those as *a nullcone member of some
moment functional*, and that reading has teeth in at least the tournament and Gaussian cases.

## The roots of unity are the same roots of unity

The thread the owner has been pulling — "keep things roots of unity" — is one thread here, not
three:

- the **regular-tournament spectrum** lies on `Re = −½` (THM-1555, the half-dictionary), and for
  the **circulant/Paley** tournaments the eigenvalues are **Gauss sums** — algebraic
  combinations of roots of unity;
- **TNC's singular indices** are the negative **roots-of-unity exponents** `−(D − x)`,
  `x ∈ {j/N} ∪ {k/M}` (THM-1720/1725) — the monodromy of the `M`-th and `N`-th root branches
  of `z^M = tR(z)`;
- **GMC's projection `φ`** is literally the **average over `U(1)`**, the roots-of-unity filter
  in the continuum limit.

So "the eigenvalues / exponents / averaging characters are roots of unity" is a single
statement about the *symmetry `G`* in the template: `G` is a torus or a cyclic group (or `S_n`,
whose characters are still integer combinations of roots of unity), and its trivial-component
projection is a root-of-unity average. The tournament and the Gaussian are averaging over the
same kind of `G`.

## Two moments of one tournament

The template also unifies *within* the tournament world. There are two natural projections:

- **trace** (over the spectral functional): gives `tr(A^m)`, nullcone = transitive;
- **Burnside** (over `S_n` relabeling): `A000568 = (1/n!) Σ_σ Fix(σ)`, the tournament *count*.

Both are `⟨·⟩_G` for a symmetry `G` (the trace's `G` is conjugation; Burnside's is `S_n`). The
"everything is the triangle" reflection already knows `Fix(σ) = 2^{orbit-pairs}` for all-odd
`σ` and `0` for even — that is the trivial-component projection made explicit. So the tournament
*enumeration* and the tournament *spectrum* are two moments of the same object, and the
staircase's two legs (scores/OCF on the vertical, complement on the horizontal) are the two
projections' shadows.

## What I would not overclaim

The template is a lens, not a theorem factory. It does not prove TNC or GMC (those are THM-1710
/ THM-1740, and GMC in full is still open at its holonomic-limit/Laplace-determinacy boundary,
THM-1690). It does not yet place `H` itself or LRC on the ladder — those are the honest
next questions (THM-1750 names them). What it does is make the recent work *cohere with the
core*: the project has been studying moment nullcones all along, and the detection-depth
machinery is the general form of the tool it already trusted for tournaments.

## The three named-next questions, worked (S128c129)

1. **`H` on the ladder — answered: it leaves at `n = 6`.** Grouping all tournaments by their
   moment vector `(tr A¹,…,tr Aⁿ)`, `H` is constant on every co-spectral class for `n ≤ 5` but
   **splits at `n = 6`** — the class `(0,0,12,12,10,48)` carries both `H = 13` and `H = 17`
   (THM-1765). So `H` is not a moment; it is the tournament's `#P`-hard permanent, one rung
   **above** the holonomic ceiling. The ladder is `rational ⊂ algebraic ⊂ holonomic ⊂ #P`, and
   the tournament spans it end to end — spectrum at the bottom, `H` at the top. THM-133's
   spectral `H = (462 − tr A⁴)/2` is a `Z₇`-circulant symmetry collapse, not the general law.
   `n = 6` joins the project's `n ≥ 6` phase transitions.

2. **LRC — an honest limit, and a dual.** LRC does **not** instantiate the template: its
   functional `M(S) = max_t min_v ‖vt‖` is a **min–max (an extremal value)**, not a moment sum,
   and its distinguished configuration — the tight AP — is a *maximiser* (THM-894: maximal
   spectral excess of the resonance matrix), which is the **opposite pole** from a nullcone.
   The nullcone is where `F` collapses to *nothing*; the tight AP is where structure is
   *maximal*. So the frame extends to LRC only as a **duality** — every moment functional has a
   trivial pole (nullcone: transitive, one-sided, charge-one-sided) and an extremal pole
   (regular/Paley tournaments, tight AP) — not as another nullcone instance. Stated so the
   analogy is not forced the wrong way.

3. **The Lean interface exists.** `MomentNullcone.lean` (sorry-free): a `Data` structure
   (`phi`, `order`, `step`), `detect` (= `zeros_propagate`, the finite-depth conclusion),
   `escape_within` (contrapositive), and `ofMonicRec` (build the `step` from a monic
   recurrence). The three instances feed their own recurrences — Cayley–Hamilton, THM-1670,
   THM-1740 — into `ofMonicRec`; `H` is correctly *excluded* (no governing recurrence). The
   template is now one reusable engine in the kernel.

## Cross-links

THM-1750 (the template, tournament instance proved) · THM-895 (`λ=0 ⟺ transitive`) ·
THM-1710 (TNC depth `D`) · THM-1740 (GMC finite Gröbner) · THM-1555 (half-dictionary, spectrum
on `Re=−½`) · THM-1720/1725 (roots-of-unity singular indices) · THM-133 (spectral-OCF chain) ·
THM-894 (Kendall–Wei on LRC speeds) · `everything-is-the-triangle.md` · Burnside
`A000568`. Arithmetic hierarchy: rational ⊂ algebraic ⊂ D-finite; Cayley–Hamilton = the
rational floor.
