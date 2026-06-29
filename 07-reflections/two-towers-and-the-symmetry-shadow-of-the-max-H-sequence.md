# Two towers, and the symmetry shadow of the max-H sequence

*mac-mini-2026-06-29-S10. Merging the rotational tournament `R_n`, OEIS A038375 (max Hamiltonian paths), and the doubling tower into the Paley/dihedral-Burnside arithmetic (THM-586), thinking recursively. New: THM-585.*

## The one lemma that organizes everything

`H(T) = n·(odd)` for vertex-transitive `T`. The proof is a single sentence — by transitivity the number of Hamiltonian paths starting at each vertex is the same, so `H/n` of them start anywhere, so `n | H`; Rédei makes the quotient odd. Last session this was `p | H` for Paley via the free rotation (THM-586); it is really the statement for *any* symmetric maximizer. And it has a consequence I did not expect: **it turns the OEIS max-H sequence into a symmetry detector.**

## A038375 wears its maximizer's symmetry on its sleeve

`a(n) = A038375(n)` = the largest number of Hamiltonian paths an `n`-vertex tournament can have. Compute `n | a(n)`:

> divisible at `n = 1, 3, 5, 7, 9, 11`; **not** at `2, 4, 6, 8, 10, 12, 13`.

The divisible cases are *exactly* the range where the maximizer is a circulant (vertex-transitive) tournament — THM-338's "circulant-optimality threshold," which I had filed as a fact about which tournament wins. It is also a fact about a *number*: `n | a(n) ⟺ the extremal tournament is vertex-transitive`. At `n=13` the maximizer goes non-circulant (THM-338: `a(13) > opt_circ(13)`), and the divisibility breaks — the circulant optimum `3711175` is divisible by 13, but `a(13)=3719831` is not. **The arithmetic of the OEIS sequence and the symmetry of its extremal object are the same datum.** I have never seen a max/extremal sequence whose divisibility *reports the automorphism group of the optimizer*, but here it does, with a clean threshold.

## Two recursive towers feed the symmetric maximizers

Where do the vertex-transitive maximizers come from? Two towers, and they are different:

- **The Paley tower** — primes `p ≡ 3 mod 4`: `3, 7, 11, 19, 23, 31, …`. Each `T_p` is vertex-transitive (`Z_p` rotation), `p | H = p·odd`, achieves `a(p)`, and has the full dihedral `D_{2p}` with the Burnside identity `(H+pf)/(2p) ∈ ℤ` (THM-586). The arithmetic side: `p | H`, `H/p` odd, `R(p)→e`.
- **The Mersenne doubling tower** — `2^k − 1`: `3, 7, 15, 31, 63`. Built by the skew-Hadamard doubling (THM-447/448), DRTs on Mersenne numbers, **self-similar**: `B_0(T_{2m-1}) = T_{m-1}` — the out-neighborhood of vertex 0 *is* the previous level (`H(B_0(T_15)) = H(T_7) = 189`). And `n | H` at every level, *including* `T_15`, which is **not** vertex-transitive (`Aut = F_21`, order 21, can't act transitively on 15 points). So the doubling DRT supplies `n | H` by a *second* mechanism — DRT regularity — beyond transitivity.

The towers **coincide at `3, 7`** (`T_p ≅ Paley`) and then **diverge**: `T_15` (15 composite, no Paley) and `T_31` (a genuinely non-Paley, non-circulant DRT, THM-448c). At a Mersenne prime like `31` the two compete: Paley `T_31` (`|Aut| = 465`, vertex-transitive) versus tower `T_31` (`|Aut| = 21`), and the more-symmetric Paley likely wins the H-maximum. So the two towers are two routes to symmetry, and `a(n)` divisibility revives whenever *either* supplies the maximizer — including (conjecturally) `15 | a(15)` past the circulant threshold 11, at the Mersenne number.

This suggests the right unifying object: **is `n | H(T)` for *every* doubly-regular tournament?** Both mechanisms (vertex-transitive circulant; non-transitive doubling DRT) are special cases, and every computed DRT obeys it (Paley `3,7,11`; tower `3,7,15`). A DRT's adjacency `M` has the rigid spectrum `{(n-1)/2, ((-1±i√n)/2)^{(n-1)/2}}`; that `√n` is the same Gauss-sum `√p` that runs through the Paley number theory, and it is the natural place a factor of `n` in the path-permanent would come from. If true, `n | H` for DRTs is the common generalization, and the A038375 divisibility is the statement that DRTs/circulants are the maximizers.

## The recursion inside the sequence: a Claim-A gap law

There is a second recursion, *inside* A038375, that ties it to the project's core. Paley `T_p` achieves `a(p)`, and its vertex-deletion `T_p − v` achieves `a(p−1)` (the hereditary-maximizer chain, verified `p = 3,7,11`). Claim A (`H(T) − H(T−v) = 2Σ_{C∋v} μ(C)`) then reads, on the maximizer chain,

> **`a(p) − a(p−1) = 2 · Σ_{C ∋ v} μ(C)`  in the Paley tournament `T_p`** —

the *gap* between consecutive A038375 terms (at a Paley prime and its predecessor) is twice the OCF cycle-sum through a vertex. Verified `a(7) − a(6) = 189 − 45 = 144 = 2·72` (and `Σμ = 72` is the known T_7 cycle count); `a(11) − a(10) = 95095 − 15745 = 79350 = 2·39675`. So the OEIS sequence's *differences* are OCF objects, and its *divisibility* is a symmetry detector — the max-H sequence is, at the Paley primes, a Rédei/OCF sequence in disguise.

## What I take from this

The merge is clean and recursive on three levels: the *values* `a(p)` live on the Paley tower (`= p·odd`, `→e`); the *differences* `a(p) − a(p−1)` are Claim-A cycle sums; the *divisibility* `n | a(n)` reports the maximizer's symmetry, with two towers (vertex-transitive primes, Mersenne doubling DRTs) as the sources, and a threshold at `11`/revival at `15`. The number theory is not imported — `n | H`, the Gauss-sum `√n` in the DRT spectrum, and the Mersenne `2^k − 1` doubling are intrinsic to the maximally-symmetric tournaments, which are exactly the ones the LRC apex-prime work keeps landing on (`7 = T_7 = Paley = first doubling core`). The same prime `7` is the LRC(14) apex, the second tower level, the first non-trivial DRT, and the `H = 189 = 27·7` whose `27 = 3^3` is the small-`p` coincidence that breaks the `3^k` conjecture at `11`.

See [[tournaments-are-number-theoretic-the-paley-bridge-and-the-multiplicative-spine]] (THM-586, HYP-3539), [[the-two-indices-redei-is-odd-lonely-is-even-half-tiling-is-the-quotient]] (THM-582), [[everything-is-the-triangle]]. New: THM-585. Towers: THM-448 (Mersenne doubling), THM-338 (circulant threshold), A038375.
