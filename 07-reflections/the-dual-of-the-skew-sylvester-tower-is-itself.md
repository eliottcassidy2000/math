# The dual of the skew-Sylvester tower is itself (and its mirror is the other heterotic string)

**Source:** mac-mini-2026-06-14-S1. Dispatch: "what tournament object is dual to the
skew-Sylvester tower?" — the extended-Hamming slot left open in the Reed–Muller
reflection (reed-muller-on-the-tiling-cube.md). Canon: THM-480 (skew-tower row code =
d⁺ ladder), THM-447/448 (the tower), THM-451 (chirality), HYP-2409, T818.

## The question, and why it has a clean answer

The skew-Sylvester tower (THM-447) doubles a tournament by
H_{2m} = [[H, H], [−Hᵀ, Hᵀ]] (H + Hᵀ = 2I, skew type preserved), seeded at the
bordered Paley₇. Its cores are doubly regular tournaments on the Mersenne orders
2^k − 1 (THM-448). The user asks for its **dual**. The answer comes in four layers,
and the headline is the cleanest possible: **it is its own dual.**

## Layer 1 — Code level: SELF-DUAL (the tower is its own dual)

THM-480: in the tournament gauge (binary code C(H) = row span of (J−H)/2 over 𝔽₂),
the tower's row code is the **self-dual d⁺ ladder**

    ê₈ = RM(1,3)  →  d₁₆⁺  →  d₃₂⁺  →  (conj. d_{2^k}⁺, HYP-2409),

every level a doubly-even self-dual (Type II) [2^k, 2^{k-1}, 4] code. Self-dual means
C = C⊥. **So the dual code of the skew tower is the skew tower.** The orthogonal-
complement operation — the thing "dual" usually means for a code — returns the same
object. This is not a coincidence to be explained away; it is the whole point of the
dispatch slogan "codes with k = n/2 are self-dual," realized by skewness.

## Layer 2 — The RM-duality resolution: skewness fuses the two ends

Where does "extended Hamming" (the dual I asked about) go? It is hiding in the gauge
choice:

- The **Sylvester** (symmetric ±1) gauge spans RM(1,m) — the **biorthogonal LOW end**
  of the Reed–Muller filtration (dim m+1).
- Its code dual is RM(m−2,m) — the **extended-Hamming HIGH end**.
- The **skew/tournament** gauge sits at neither end. It lands on the **self-dual
  MIDDLE**: RM(1,m) ⊊ d⁺ ⊊ RM(m−2,m), with d⁺ self-dual exactly halfway
  (dim 2^{m−1}, between m+1 and 2^m − m − 1).

So the two questions — "what is the dual of the Sylvester tower?" (answer: extended
Hamming) and "what is the dual of the skew tower?" (answer: itself) — are one picture.
**Skewness is the operation that moves you off the dualizable ENDS and onto the
self-dual FIXED POINT of the RM duality involution RM(r,m) ↔ RM(m−r−1,m).** The skew
tower is the geometric mean of the Sylvester code and its extended-Hamming dual.

## Layer 3 — The genuine tournament-level dual: TRANSPOSE / CHIRALITY

The duality that acts *nontrivially* on the tower is transpose Tᵀ (= arc reversal =
the complement of the skew matrix S ↦ −S = Sᵀ). THM-451:

- Orders ≤ 8: the tower is **achiral** — skew₈ ≅ skew₈ᵀ (Hadamard-equivalent to its
  own transpose). Here self-dual and self-transpose coincide.
- Order ≥ 16: the tower is **CHIRAL** — skew₁₆ ≇ skew₁₆ᵀ. The transpose-split pair is
  exactly {had.16.3, had.16.4}, the unique transpose-distinct pair among Hall's five
  order-16 Hadamard classes.

And those two classes carry precisely the two Type II [16,8] codes:
**d₁₆⁺ (indecomposable — the tower's branch) and e₈⊕e₈ (decomposable)** — the Milnor
isospectral pair, equivalently the two even unimodular lattices in rank 16,
equivalently the two heterotic string theories **SO(32) (= D₁₆⁺) and E₈×E₈**. The
doubling selects the **indecomposable SO(32)/D₁₆⁺ branch**; its dual/mirror is the
**E₈×E₈ branch**. So:

> **The tournament object dual to the skew-Sylvester tower is its transpose (chiral
> mirror); for orders ≤ 8 the mirror IS the tower (achiral, self-dual = self-
> transpose); at order 16 the mirror splits off as the OTHER heterotic lattice
> E₈⊕E₈, which the doubling does not build.**

The code is self-dual at every level, but the *tournament* is self-transpose only up
to order 8; the chirality is where the dual becomes a genuinely different object —
and it is the most famous "isospectral but inequivalent" pair in mathematics/physics.

## Layer 4 — The even/odd sibling (the symmetric dual construction)

One level out: the skew tower is the **odd** (q ≡ 3 mod 4, χ odd, tournament, skew-
Hadamard) Paley tower. Its even sibling is the **symmetric conference / Paley-graph
tower** (q ≡ 1 mod 4, χ even, graph, symmetric-Hadamard) — the E_n even-graph side
of the project's even/odd duality. That is the dual *construction* (symmetric vs
alternating), as opposed to the dual *code* (self) or dual *tournament* (transpose).

## What this leaves open (HYP-2409, the persistence question)

THM-480 verified the d⁺ identification rigorously at orders 8, 16 and at the
weight-enumerator level at 32. The honest open core: does the tower's row code stay
**self-dual Type II** (the achievable half) and, harder, **indecomposable d⁺** (vs a
decomposable e₈⊕e₈-type) for ALL k? The self-dual Type II half is provable from the
double-and-glue structure (dim = n/2 every level via THM-451's GF(2)-rank n/2;
self-orthogonal + doubly-even from the (b,b) pair-doubling + Hᵀ glue). The
indecomposability — the tower always picking the SO(32)/connected branch — is the
deep part; the number of Type II codes explodes with length, so a structural
"the doubling glue is always connected" argument is needed (the weight-4 support
graph of the glued code stays connected). That argument, if completed, says: **the
skew tower rides the SO(32) branch forever — its dual is itself, never decomposing
into mirror halves.**

## The one-line answer

**The skew-Sylvester tower is self-dual: its row code is its own dual (the d⁺ ladder
at the k=n/2 fixed point of Reed–Muller duality, between the Sylvester biorthogonal
code and its extended-Hamming dual). The only duality that moves it is transpose-
chirality, which at order 16 reveals its mirror to be the other heterotic lattice
E₈×E₈ — the doubling itself always rides the indecomposable SO(32)/D₁₆⁺ branch.**
