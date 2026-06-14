# Efficiency becomes proof — the repo runs on speedup-theorems, and its history holds the speedups

**Source:** kind-pasteur-2026-06-13-S4. Dispatch: think about tricks to speed up
the computations, get inspired by the repo's history, and see how efficiency
becomes mathematical insight and proof.

## The principle

Almost every fast algorithm in this repo is a theorem wearing work clothes. A
speedup replaces a naive enumeration with a closed form or a recursion; the
*reason* the closed form is correct is a mathematical identity, and that identity
is the insight. The repo's engineering mandate already says this ("every rank
trick is a library AND a theorem"); stated as a slogan: **the correctness proof
of a speedup is its mathematical content.** A catalogue of our own:

| naive computation | speedup | the theorem it IS |
|---|---|---|
| enumerate Hamiltonian paths (factorial) | `H = I(Ω,2)` | the OCF (THM-002 / Grinberg–Stanley): a representation count has a closed form in the conflict graph |
| enumerate cyclic triples (`O(n³)` triples) | `c3 = C(n,3) − ΣC(s_v,2)` | Kendall–Babington-Smith (used in THM-462): `c3` is a quadratic in the scores |
| enumerate directed 5-cycles (`O(n⁵)`) | `c5 = tr(A⁵)/5`, `O(n³)` | THM-118: closed ≤5-walks can't repeat a vertex in a 2-cycle-free digraph |
| `O(m²)` skew-Hadamard transform | `O(m log m)` butterfly | THM-451: the skew-Walsh matrix is a self-similar tensor (the chirality split) |
| Burnside sum over all `n!` permutations | cycle-type recursion | the Mode-B / metagraph recursion (CLAUDE.md): `Fix(σ)` depends only on cycle type |
| optimize `M(S)` over the circle (continuum) | residue check `va mod q ∉ B_q` | the band criterion (THM-492): looseness is a finite congruence test |
| sweep all shells `q` for an LRC witness | `q ≤ 13` ⟹ pure divisibility | the band-0 lemma (THM-497 B): the cheapest clocks are divisor tests |
| Gaussian elimination over ℚ (fill-in, overflow) | small-prime modular rank | mod_rank: rank is stable mod a good prime (the engineering library) |

Each row: a slow loop on the left, a fast formula on the right, and a *named
theorem* in the middle that is the whole point. The repo doesn't optimize code as
an afterthought to the math — the optimization *is* the math, found by asking
"why could this possibly be faster?" and discovering the structure that answers.

## This session's instance — and its humbling lesson

Last session (THM-498) I found the `c5` cycle-count spectrum has a gap (`c5=10`
forbidden at `n=6`) by an exhaustive `O(n⁵)` enumeration, and left `n=7` OPEN
because `2^{21}` tournaments × the slow counter was infeasible. The fix was not a
new trick — it was **THM-118, proved in this repo on 2026-03-07**: `c5 = tr(A⁵)/5`
exactly (a closed ≤5-walk in a tournament can't revisit a vertex without a
forbidden 2-cycle; the identity fails first at `k=6`, the figure-eight of two
triangles). That is `O(n³)`. With it, the `2^{21}` sweep runs in seconds, and the
`n=7` spectrum falls out exhaustively: `[0,42] ∖ {34,37,38,39,40,41}` — the
forbidden set *migrates to the top* of the range, and `c5=10` is no longer
forbidden at `n=7`. So:

> The speedup that unblocked yesterday's open computation had been sitting in our
> own canon for six months. "Get inspired by the repo's history" is not a
> metaphor — the history literally contains the efficiency the present needs.

The lesson generalizes a recurring failure mode (MISTAKE-020's truncation, the
LRC scripts' slow counters): **before brute-forcing an invariant, ask whether the
repo already proved it equals a trace, a score polynomial, or a recursion.** The
inventory above is the first place to look.

## Efficiency becomes *proof*, not just speed

The deeper half of the dispatch: a speedup doesn't only compute faster — it
**reframes the question into a domain where it can be proved.** Two instances from
this session:

1. **The trace formula turns a combinatorial gap into a spectral exclusion.**
   `c5=10` forbidden at `n=6` looked like a brute-force accident. Through THM-118,
   `c5·5 = tr(A⁵) = Σλ_i⁵`, the 5th power-sum of the spectrum; the achievable
   `tr(A⁵)` at `n=6` is every multiple of 5 up to 45, then 55, 60, *skipping 50*.
   So the gap **is** the statement "no 6-vertex tournament has `Σλ_i⁵ = 50`" — a
   skew-spectral exclusion, living on the same `det(I+S)` / power-sum object as the
   determinant lens (THM-468/472) and the spectral-OCF chain (THM-133). A proof of
   the gap is now a power-sum-realizability question, not a search. That is
   efficiency becoming the *form* a proof must take.

2. **The "effective asymptotic + finite check" template (last session's Pollock
   method) is the ultimate speedup-as-proof:** prove a main term dominates for
   all `N` past an explicit threshold, and the infinite verification collapses to a
   finite one. The whole content is the threshold-lowering; the asymptotic *is* the
   proof. Our additive-basis DP is the finite-check half; the circle method is the
   asymptotic half. The same shape is what the LRC resource bound (t-0124) needs.

The unifying picture: **trace formulas, score polynomials, recursions, and
effective asymptotics are all the same move** — replace "enumerate the objects" by
"evaluate the invariant of the structure," where the structure (the spectrum, the
score sequence, the singular series) is lower-dimensional than the objects. The
dimension you save is exactly the structure theorem you prove. Efficiency and
insight are not two activities; the search for one is the discovery of the other.

## Coda (kps5): every reframe has a boundary, and finding it is also insight

Pushing the spectral reframe one session further produced its *edge*, which is as
valuable as the reframe. The map of tournament invariants splits cleanly:

- **spectral side** (functions of the spectrum): `c3, c5 = tr(A^k)/k`, and
  `d = det(I+S) = ∏(1+μ_j²)` — verified A-spectral (constant on every cospectral
  class at n=6).
- **non-spectral side** (need the conflict graph `Ω`): `H`, `c7`, the disjointness
  tower `α_2, α_3, …`.

The boundary is exact and mechanical: at `n ≤ 6`, `H = 1 + 2(c3+c5) + 4D` where
`D = α_2` = vertex-disjoint triangle pairs (THM-499). `α_1 = c3+c5` is the
power-sum/spectral layer; `D` is the *first* piece of conflict-graph data the
spectrum cannot see, and it appears precisely at `n=6` (at `n=5` there is no room
for two disjoint triangles, so `H` is spectral there). Cospectral tournaments —
eigenvalue twins — differ exactly in `D`, and `H` resolves them while `d` cannot.

Two payoffs of locating the boundary. (i) It explains a standing empirical fact:
the determinant lens's `d ⊥ H` (Pearson ≈ 0) is the orthogonality of a *spectral*
coordinate (`d`) and a *non-spectral* one (`H`) — they measure opposite sides of
the same boundary. (ii) It tells the forbidden-value program exactly how far the
spectral method reaches: the cycle gaps (`c5=10`) are provable on the spectral
side (a finite score-stratification: the regular class achieves `c5∈{6,8,11,12}`,
skipping 10), but the `H`-gaps `{7,21}` are not — they require the `α_2`
conflict-graph argument. Conflating the two was an overreach; the boundary is the
correction. **A reframe is a tool with a domain; knowing where the tool stops is
the same kind of knowledge as knowing it works.**

Cross-links: THM-118 (the trace identity, the speedup used here), THM-498 (the
cycle-spectrum gaps, spectrally reframed + n=7 exhaustive), THM-499 (the boundary:
H non-spectral, H=1+2(c3+c5)+4D, d spectral), THM-002 (OCF), THM-462 (Kendall–BS /
four-square), THM-451 (skew-Walsh butterfly), THM-492/497 (band criterion / band-0),
THM-133/468 (the skew spectrum / determinant lens), the Pollock reflection
[[pollock-as-the-bounded-arity-currency-and-the-cycle-spectrum-onset-kps3]].
