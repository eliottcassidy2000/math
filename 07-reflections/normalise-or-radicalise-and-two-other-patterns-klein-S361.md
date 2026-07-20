# Normalise or radicalise — and two other patterns that kept recurring

*klein-2026-07-20-S361. The owner asked me to look back through past work for what becomes
reminiscent of the recent patterns, and to think abstractly. Three patterns recur hard enough
to be worth naming. The first is the one with a concrete payoff and it predicts its own
applicability; the second reframes a classical argument; the third is about me.*

## 1. Normalise the degenerate locus away, or pay for radical membership

Almost every question in this cluster has the shape *"is the solution variety contained in the
degenerate locus `D`?"* — `D` being "the span collapses", "the polynomial is constant", "`P` is
one-sided". That question is **not** ideal-inconsistency, because points of `D` are genuine
solutions. It is **radical membership**, which is expensive and, in my hands, error-prone.

But sometimes a group `G` rescues it:

> **If `G` acts on `V ∖ D` transitively enough to fix a normal form, then normalising turns
> "is `V ⊆ D`?" into "is the normalised system INCONSISTENT?" — a Gröbner basis of `⟨1⟩`,
> cheap and unambiguous.**

And crucially, **the pattern predicts when it is available**:

| | `D` | `G` | normalisable? | test used |
|---|---|---|---|---|
| **TNC** (THM-1600) | extreme coefficient vanishes | `Λ↦μΛ`, `u↦λu` | **yes** — `D` is away from the origin | inconsistency, `⟨1⟩` |
| **S357 system** (THM-1590) | `α = β = 0`, the **origin** | scaling | **no** — every scaling fixes the origin | radical membership |

Checked both ways this session: the normalised TNC systems at `(1,1)` and `(2,2)` return
`⟨1⟩`; the S357 system returns `[a₀, b₀]`, never `⟨1⟩`, exactly as predicted. **Scaling can
normalise "some coordinate is nonzero"; it can never normalise "all coordinates vanish."**
That is decidable in advance, and had I asked it in S357 I would not have used a
radical test with the depth set wrong.

**The repo has been doing this move for a long time without naming it.** Three instances, all
the same abstraction:

- **Tilings are a transversal for switching** (THM-1430 §3): `cut(K_n)` acts on arc space; the
  base path is a spanning tree, so each coset has exactly one representative with the path
  arcs standard. The tiling model *is* a group normalisation.
- **The half tiling is a coordinate fundamental domain** (THM-549): fold by the grid
  reflection `σ`, and the half-region is the normal form.
- **The merged metagraph** is the same fold by complement.

So the tiling machinery and the nullcone elimination are the same manoeuvre in different
clothes — which is worth knowing, because it means the repo's instinct for fundamental domains
transfers directly to the algebraic-geometry side of the GMC work.

## 2. One relation buys a disjunction; a family buys rigidity

Recorded in S353 and it keeps paying. The classical Vieta/`e`–`π` argument concludes only
*"at least one of `e+π`, `eπ` is transcendental"* — because there is exactly **one** relation
on the pair (sum, product). Schanuel's conjecture is precisely the upgrade from disjunction to
conjunction.

Everything in this cluster that concluded **rigidity** did so by having a *one-parameter
family* of relations:

| result | the family | conclusion |
|---|---|---|
| EMP (THM-1510) | `E[h^m] = 0` for every `m` | `h ≡ 0` |
| the master theorem (THM-1500) | vanishing for every `m` | `f = f(0)/(1+s)`, forced |
| THM-1550 | `Π(t) = ct` for every `t` | `R` constant |
| LRC pair-sum ruler | one symmetric function only | a **set** of candidate pairs |

The diagnostic: **when a result stalls at "at least one of these fails", ask whether a family
of relations is available or only a single one.** That is a question about the *setup*, not
about effort, and it is answerable before starting.

## 3. My errors are in the instrument, not the argument

Uncomfortable but consistent, across seven consecutive sessions:

| session | the defect | the mathematics |
|---|---|---|
| S339 | a **hardcoded verdict** in a print statement, printed against its own data | right |
| S343 | `nu = 2/a` where the prose said `2a` | right |
| S345 | a flag counting **parity-forced** zeros as failures | right |
| S351 | a line meant to divide by a polynomial that was a **no-op** | right |
| S357 | radical test depth `k ≤ 6`, and testing **inconsistency when the origin is always a solution** | right |
| S359 | a ten-minute **timeout** from not escalating incrementally | right |

Seven for seven. The argument was never the thing that broke. Three rules follow, and they are
narrower and more useful than "check your work":

1. **Never put a conclusion in a print statement.** Report measurements; write the verdict
   afterwards, against the numbers.
2. **Test the right predicate.** Ask what the degenerate locus *is* before choosing between
   inconsistency and radical membership — §1 makes that decidable.
3. **Run the test on a control whose answer you already know.** Every one of these was caught
   that way, or would have been: `Λ = u⁻¹+u` is not in the nullcone, `M = 1` has a proof, `β=0`
   vanishes by parity.

The reason this matters beyond bookkeeping: in a fleet where results are cited by other agents
within hours, a wrong *instrument* produces a confidently-worded artefact that reads exactly
like a finding. Three of the six above would have been published as results if the control had
not fired.

---

*Related: THM-1600 (the normalisation), THM-1590 (where it was unavailable), THM-1550, THM-1510,
THM-549 and THM-1430 (the repo's older instances of the same fold), and S353's Vieta note.
Script: `04-computation/` — the check in `05-knowledge/results/normalise_pattern_klein_S361.out`.
Not a HYP reservation.*
