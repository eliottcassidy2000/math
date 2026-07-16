# HYP-7027 — arXiv:2607.14068 (even-uniform hypergraph Moore bound) × this repo, + the T1540 cut⊕cycle probe

**Status:** CLAIMED / IN PROGRESS (death-star-2026-07-16-S21; owner directive: run the T1540
probe + deeply integrate the paper). Verify-first until results land.

## The paper (read in full, exact statements)

Bandeira–Kunisky–Nizić-Nikolac–Pesenti–Wang, "The even-uniform hypergraph Moore bound"
(2607.14068, July 15 2026): every k-uniform hypergraph (k ≥ 4 even) with
m ≥ 64n(n/ρ)^{k/2−1} hyperedges has an EVEN COVER (hyperedge set covering every vertex an
even number of times = nonzero F₂ kernel vector of the incidence system) of size
≤ 4kρ·log n — Feige's conjecture, no polylogs. Proof: Kikuchi graph K_ℓ(H) on ([n] choose ℓ),
S ~ T iff SΔT ∈ H (edge colored by SΔT), ℓ = max(k/2, ⌈ρ⌉); dense core; "unpaired girth"
(cycles with some color odd); BFS balls grow ×16 unless a short unpaired cycle exists; the
internal-edge bound |E(N_j,N_j)| ≤ ℓ|N_j| via **Lemma 2.3 (diagonal polynomial system /
one-inclusion, Frankl–Wilson/HLW-type)**: Hamming-1 adjacent labels + degree-≤d diagonal
polynomials ⟹ avg deg ≤ 2d; the labels = **color palettes** pal(S) ∈ F₂^H (odd-visited
colors on a walk from the root), f_S(X) = ∏_{v∈S}(1{v∈S₀} + Σ_{E∋v}X_E), and the
fixed-weight slice makes 1{S⊆T} = 1{S=T}. (Ack: core innovation found by GPT-5.6 Sol;
Claude Opus 4.8 + Fable 5 also used — the field's honest-AI-science norms, noted.)

## The bridge to this repo (each PRECISE, to be written into the reflection)

1. Kikuchi = a weight-ℓ hypercube slice with H-XOR moves = our waggly-layer machinery
   (d-strata of Q_m) restricted to a level set; THM-584's level parity + the Krawtchouk
   frame (OPEN-Q-040) are the same spectral geometry.
2. Even covers = the F₂ cycle space — the repo's FIRST-CLASS even graphs (E_n metagraph,
   A002854); paired/trivial vs unpaired Kikuchi cycles = our silent vs expressive mutations.
3. pal(S) = tiling coordinates anchored at a root (tournament = base path + tile XOR).
4. The Moore bound = a density ⟹ short-F₂-dependency forcing theorem = the LDPC
   rate–distance tradeoff — engineering mandate hooks: circulant LDPC (THM-125/QR_p),
   `mod_rank_library` (kernel search IS sparse F₂ rank), the [72,36,16] gauge thread.
5. LRC face: THM-890's error spectrum = small additive relations Σkᵢfᵢ = ℓe; the F₂ shadow
   of relations = even covers of the speed/divisibility hypergraph. NAMED NEW QUESTION: does
   covering-saturation (the divisibility density of covering 13-sets) FORCE short relations
   (Moore-bound style) — i.e. "covering ⟹ loose OR short-relation-coherent" — a structural
   forcing theorem shaped for the |P| ≤ 6 strata?

## The T1540 probe (this session's computation)

The miss-pattern movie = the wall-crossing event stream around the circle (states = section
vectors σ ∈ Z₇⁶, colors = walls). Sharp claim to test: **the movie's return spectrum
(repeated states at gap Δx: simultaneous approximation of the fᵢ) and THM-890's relation
spectrum (Σkᵢfᵢ = ℓe: the dual lattice) are Khintchine-TRANSFERENCE DUALS** — cut side =
returns/states, cycle side = relations/palettes; trivial returns (zero palette XOR) =
silent loops. Compute both spectra exactly on far/compact/planted cores; verify the duality
quantitatively; measure the movie's palette rank vs the relation-lattice rank.

-> THM-890/889, T1540, THM-584, OPEN-Q-040, LEM-020, engineering-synthesis; death-star-S21.
