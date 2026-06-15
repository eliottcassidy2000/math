# Two Lyapunovs, the tournament as a zero-sum game, and the E₇ 56-rep as the exceptional symplectic tournament

**Source:** kind-pasteur-2026-06-15-S5. Dispatch: consider E₇ and use the Lyapunov
ideas of arXiv:2606.11773 as inspiration. That paper is Orabona, *Last-Iterate
Convergence of Optimistic Multiplicative Weight Update (OMWU)* — its "Lyapunov" is
a **Lyapunov FUNCTION** (a monotone KL/Bregman potential certifying last-iterate
convergence to a Nash equilibrium of a two-player zero-sum game), **not** the
Lyapunov EXPONENT (random-product growth rate) the repo used for `γ_pent` (THM-488)
and `φ` (HYP-614).

## Two senses of Lyapunov, both native to tournaments

| | Lyapunov EXPONENT | Lyapunov FUNCTION |
|---|---|---|
| object | growth rate of random products / a recurrence | a monotone potential `V≥0`, `V↓` |
| repo home | `γ_pent≈0.206` (random-sign), `φ` = HYP-614 growth | **this session: OMWU on the tournament game** |
| certifies | exponential growth / entropy | convergence to equilibrium |

The dispatch's paper supplies the *function* sense, and it has a clean tournament
home that the repo had not used: **a tournament IS a symmetric zero-sum game.**

## The tournament game (verified)

For a tournament with skew-adjacency `S = A − Aᵀ` (antisymmetric), the symmetric
zero-sum game with payoff `S` has value 0 and a Nash equilibrium = the classical
**bipartisan set** (the optimal mixed strategy; unique with odd support on
odd-vertex tournaments, Laffond–Laslier–Le Breton). Running OMWU
`x_{t+1} ∝ x_t·exp(η(2Sx_t − Sx_{t-1}))` converges last-iterate to it, with the
Lyapunov potential `KL(x*, x_t) ↓ 0` (verified, `04-computation/` this session):

- **3-cycle = rock–paper–scissors:** Nash = uniform `(⅓,⅓,⅓)`; KL→0.
- **Transitive T₃:** Nash = `(1,0,0)` — a *pure* strategy on the dominator; KL starts
  at `log 3` and →0.
- **Regular T₅, Paley T₇:** Nash = uniform; KL→0.

So the tournament-game Nash equilibrium is **concentrated for transitive/dominance,
spread (uniform) for regular/balanced** — the equilibrium support is a *tournament
solution*, and it tracks the H-gradient exactly: transitive (`H=1`) ↦ pure Nash,
Paley (`H=max`) ↦ uniform Nash. This is last session's "concentrated dominance vs
balanced core" reappearing as the bipartisan equilibrium, and the OMWU Lyapunov
function certifies the convergence in both regimes.

## The dynamics' frequencies are the skew spectrum (the determinant lens)

The OMWU/game dynamics rotate at frequencies = the imaginary eigenvalues `±iμ_j` of
the antisymmetric payoff `S` (the skew spectrum). That is exactly the repo's
**determinant-lens** object: `det(I+S) = ∏_j (1+μ_j²)` (THM-468/472), and the
power sums `tr(S^{2t}) = Σμ_j^{2t}` (THM-507). So:

> The OMWU last-iterate convergence rate and oscillation frequencies of the
> tournament game are functions of the **skew spectrum** — the same `{μ_j}` the
> determinant lens and the spectral-OCF chain (THM-133) already study. The
> Lyapunov-FUNCTION analysis and the Lyapunov-EXPONENT analysis meet on one object,
> the skew spectrum.

## E₇: the exceptional symplectic tournament (the inspired direction, HYP-2530)

E₇'s **56-dimensional minuscule representation** carries an E₇-invariant
**symplectic (antisymmetric) form ω**, splitting the 56 weights into **28 antipodal
pairs** (the Freudenthal triple system; with the quartic Cartan invariant). An
antisymmetric form is the payoff of a symmetric zero-sum game — so:

> **The E₇ 56-rep is a "continuous/exceptional tournament":** the symplectic form ω
> is its skew-adjacency, the 56 weights are the strategies, the 28 antipodal pairs
> are the `λ ↔ −λ` "opponent pairs" (the complement-involution of an ordinary
> tournament). It realizes the repo's long-standing **tournament-56 bridge**
> (the LRC/unit-distance-56 threads) through E₇'s exceptional structure.

Predictions (HYP-2530), by analogy with the verified small cases:
1. The Nash equilibrium of the E₇ game is **Weyl-invariant = uniform on the 56
   weights** (they form a single Weyl orbit — the maximally-symmetric "Paley-like"
   case ⟹ uniform bipartisan set).
2. The OMWU frequencies are the eigenvalues of ω on the 56-rep, which (by the
   Coxeter/Kostant structure of minuscule reps) should be the **E₇ exponents**
   `{1,5,7,9,11,13,17}` scaled by the Coxeter number `h=18` — i.e. `μ` from
   `e^{2πi m_j/18}`. The determinant-lens `∏(1+μ²)` of this "tournament" would then
   be an E₇-invariant cousin of the **quartic Cartan invariant**.
3. The Lyapunov potential `KL(uniform, x_t)` certifies last-iterate convergence of
   OMWU on the E₇ game (the paper's theorem applied to the exceptional payoff ω).

This is the inspired EXTENSION: the OMWU Lyapunov-function machinery, the
tournament-as-game bridge, and the repo's skew-spectrum/determinant lens, lifted to
the exceptional 56-rep where E₇'s symplectic form is the "skew-adjacency." Honest
scope: the small-tournament OMWU/Lyapunov/skew-spectrum core is **verified**; the
existence of the E₇-invariant symplectic form is **sourced** (Freudenthal triple
system; arXiv:1311.0341, 1005.1275, nLab E₇); the three E₇ predictions are a
**flagged research direction**, not computed here — building ω explicitly from the
3₂₁ Gosset-polytope weights and running OMWU is the concrete next step.

Cross-links: THM-468/472/507 (the skew spectrum / determinant lens = the OMWU
frequencies), THM-133 (spectral-OCF), THM-488/HYP-614 (the Lyapunov-EXPONENT sense),
[[dominance-dodge-is-tournament-condensation-and-it-is-orthogonal-to-L-kps]]
(transitive↦pure / balanced↦uniform Nash = last session's dominance dichotomy),
the tournament-56 bridge reflections, HYP-2530 (the E₇ symplectic-tournament game).

Sources: Orabona, *Last-Iterate Convergence of OMWU*, arXiv:2606.11773;
[A Symplectic Representation of E₇](https://arxiv.org/pdf/1311.0341);
[Freudenthal triple systems by root system methods](https://arxiv.org/pdf/1005.1275);
[E₇ (nLab)](https://ncatlab.org/nlab/show/E%E2%82%87).
