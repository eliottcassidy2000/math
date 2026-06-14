# The LRC singular series — the Pollock transfer, completed

**Source:** kind-pasteur-2026-06-14-S6 (THM-501). The standing handoff from the
Pollock session (HYP-2489: "the LRC deficit is a circle-method singular series")
is now a concrete object.

## What was abstract is now a number

Two sessions ago I argued, by shape-matching, that the LRC(14) covering deficit
`D(q,S)` *looks like* a Hardy–Littlewood count: a main term plus a singular-series
fluctuation plus a finite exceptional set. This session makes the singular series
explicit. Expand the deficit in additive characters; the deviation from the
independence main term `(6/7)^{13}` is a sum over **additive resonances**
`Σ_{v∈T} t_v v ≡ 0 (mod q)` of the speeds. A resonance from a *non-zero* relation
`Σ t_v v = m` dies once `q > |m|`; only the *exact* integer relations survive to
the limit, and the Dirichlet weight collapses to the band's sinc kernel
`s(t)=sin(πt/7)/(πt)`. So

```
L(S) := lim_{q→∞} D(q,S)/q   exists
      = (6/7)^{13} + Σ_{exact integer relations Σ t_v v = 0}  (sinc-weighted)
```

is the **singular series of the LRC covering**, and it is computable from the
speeds' integer additive-relation lattice — the same Sidon/B_h object the repo's
additive-relation ladder (THM-446) is built on. The whole picture is now one
sentence: *loneliness density = the singular series = `(6/7)^{13}` minus what the
speeds' small additive relations cancel.*

## The three regimes are the three relation-densities

- **Sidon / generic speeds** (no small relations): `L ≈ (6/7)^{13} ≈ 0.135` — the
  bands behave independently, witnesses are dense, the config is trivially loose.
- **Arithmetic-progression cores** `d·{1,…,12}` (the evaders): the AP is the
  *most resonant* set — relations like `7 − 2·14 + 21 = 0` abound — so `L` is
  maximally suppressed, to `≈ 0.026–0.030`. These are exactly the hardest known
  LRC(14) configs (the evaders), and now we see *why*: an AP minimizes the singular
  series.
- **Tight** `14·{1,…,13}`: `L = 0` exactly. The singular series vanishes precisely
  at tightness — `L` is the asymptotic density of the loneliness safe set.

So the evaders are not a sporadic obstruction; they are the **extremizers of the
singular series**, the place the AP structure pushes the loneliness density lowest.
And it stays positive (`0.026 > 0`).

## The proof target, in its proper form

C'(14) — hence LRC(14) — follows from **`L(S) > 0` for every primitive
multiple-of-14 speed set**. This is the textbook circle-method endgame: *the
singular series is bounded below by a positive constant, so the main term dominates
and representations exist.* It is exactly how the 2025 Pollock proof
(Basak–Dong–Saettone–Zaharescu) and Brady's octahedral theorem close — and the
boundary I drew earlier (LRC is a *multiplicative-character* / covering problem, not
the repo's additive-tournament machinery) is respected: `L(S)` is built from
character sums over `Z/q`, with the sinc kernel of the band, and its lower bound is
a Weil/Pólya-Vinogradov-flavoured question on the relation lattice.

The honest gap is the analytic one the circle method always faces: the singular
series' defining sum is only *conditionally* convergent (the sinc weights `~1/t`,
the Dirichlet `L¹`-norm `~ log q`) — which is precisely why the naive
Pólya–Vinogradov "main term beats `√q`" bound failed for the over-correlated
(AP-core) configs two sessions ago. The resolution is the same one Linnik's method
supplies in Pollock: reduce to a ternary-form / lattice-point count and bound it
with Bombieri-along-curves, rather than a crude character-sum tail. The extremal
configs are now known (dilated APs), so the lower bound only has to be proved where
it is tightest.

## The arc

Three sessions built one bridge: *Pollock's method is the circle method with an
explicit singular series and a finite check* → *the LRC deficit has that shape* →
*the shape is the limit `L(S)`, a real number, vanishing at tightness, extremized by
APs, conjecturally `> 0`*. Efficiency-becomes-proof, one more time: computing the
deficit fast (the character expansion) **was** discovering its structure (the
singular series), and the structure **is** the proof strategy (bound `L` below).

Cross-links: THM-501 (the singular series), THM-497 (the covering reformulation /
over-correlated regime), THM-492 (band criterion), THM-446 (the additive-relation
ladder = the relation lattice controlling `L`), HYP-2489 (the transfer, now
concrete), [[pollock-as-the-bounded-arity-currency-and-the-cycle-spectrum-onset-kps3]]
(the Pollock method this completes), [[efficiency-becomes-proof-kps4]].
