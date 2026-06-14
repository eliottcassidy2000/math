# The Self-Dual Middle: where parity lives (BSD, Hodge, RH, and the tournament)

**Source:** mac-mini-2026-06-11-S2. Dispatch: merge Birch–Swinnerton-Dyer and Hodge
into the project's RH / Goldbach / Lemoine / Fermat-polygonal web. Canon: T786,
HYP-2420..2423, THM-490 (realizability semigroup), OPEN-Q-094. Builds on the
project's existing carrier: HYP-2218/2219 (Goldbach–Lemoine reflection law),
HYP-2044/2049/2051 (doubled primes), HYP-2227 (NT carrier), THM-174 (det=Pf²),
THM-442 (H²−Pf²=8Q), THM-343/075 (forbidden H), THM-477 (blue code).

This dispatch is the project's most numerology-prone yet. The discipline of this
note: **one rigorous spine, one heuristic frame, and an honest list of where the
analogy breaks.** Every claim is tagged RIGOROUS (a literal shared mechanism) or
FRAME (a structural rhyme, not a theorem).

---

## I. The rigorous spine: an alternating pairing forces a square (Pfaffian)

There is ONE theorem of linear algebra standing under all four objects:

> A non-degenerate **alternating** form on a finite/free module decomposes into
> hyperbolic planes; its determinant/order is a perfect **square**, with square
> root the **Pfaffian**.

Its incarnations are not analogies — they are the same fact in different categories:

| setting | the alternating object | the square |
|---|---|---|
| **Tournaments** (project) | S = A − Aᵀ, the skew ±1 matrix (the *odd part*) | det S = Pf(S)² (THM-174); det(I+2A)=Pf² for odd n |
| **BSD** | Cassels–Tate pairing on Ш(E) — Cassels proved it ALTERNATING for elliptic curves | \|Ш(E)\| is a perfect square (when finite) |
| **Hodge / abelian varieties** | intersection form on odd-degree H^k; Weil pairing on torsion | symplectic structure; Galois image ⊂ GSp |

The crux everywhere is **alternating** (⟨x,x⟩=0 for all x), strictly stronger than
antisymmetric. And the *defect* from alternating is itself a ±1 / doubling datum —
again the same in two settings:

- **BSD (Poonen–Stoll 1999):** for higher-genus Jacobians the polarization-induced
  Cassels–Tate pairing can be only antisymmetric, not alternating; then \|Ш\| can be
  **2 × a square**, not a square. The obstruction is a single canonical Z/2 class c
  (~13% of genus-2 hyperelliptic curves over Q have non-square Ш).
- **Tournaments (THM-442):** H² − Pf² = 8Q. The Pfaffian Pf is the square skeleton;
  Q (the #P-hard correction) measures departure from "H is itself the Pfaffian."

**RIGOROUS** that THM-174 = "alternating ⟹ square." **FRAME** (flagged, not a
theorem) that the *defect* "8Q" and Poonen–Stoll's "factor 2" are the same phenomenon
— they are different constants (8 vs 2) doing different jobs; do not identify them.
The honest statement is: *both are corrections measuring how far a clean square
captures the full invariant.*

The one place this spine must NOT be pushed: the BSD **regulator** is the determinant
of the **symmetric positive-definite** height pairing — its determinant is a generic
positive real, NOT a square. The square/Pfaffian structure of BSD lives on the Ш side
and the \|E_tors\|² denominator, never on the regulator. (Mapping THM-174 to the
regulator would be a symmetry-type error.)

## II. The heuristic frame: every duality involution has a fixed middle carrying ±1

Each object carries an involution and a distinguished **fixed locus** where a ±1 /
parity invariant lives. This is a FRAME — a way to organize, not a theorem.

| object | involution | fixed "middle" | the ±1 / parity there |
|---|---|---|---|
| RH | s ↔ 1−s (ξ(s)=ξ(1−s)) | Re(s)=½ critical line | sign of the explicit-formula error / Liouville parity |
| BSD | functional equation L ↔ Λ | center s=1 | root number w(E)=∏ₚwᵥ = (−1)^rank (parity conj.) |
| Hodge | conjugation H^{p,q}↔H^{q,p} | (p,p) diagonal | symmetric(+1, even wt) vs alternating(−1, odd wt) form |
| Goldbach/Lemoine | swap p↔q | diagonal q+q = Lemoine's 2q | E↔3E−O reflection (HYP-2218); fixed O=3E/2 |
| project | complement T↔Tᵒᵖ; grid σ | SC classes; blue code k=n/2 | Rédei parity (H odd); blue/black; SC/NS |

The user's intuition — "oddness and doubling relate to Re(s)=½" — has a **rigorous
kernel** here, exactly one: the Liouville function λ(n)=(−1)^{Ω(n)} satisfies
**λ(2n) = −λ(n)** (Ω(2n)=Ω(n)+1), and RH ⟺ Σ_{n≤x} λ(n) = O(x^{½+ε}). *Doubling
flips the parity that the ½-exponent governs.* This is the literal home of "doubling
does the work" — and it is exactly the project's HYP-2051 ("2p is the one-dyadic lift
that lets a prime cycle cross parity"). **RIGOROUS.** Everything else in the ½/doubling
story is FRAME. (Trap to avoid: Pólya's λ-sum ≤ 0 and Mertens \|M\|<√x are both FALSE;
only the weak O(x^{½+ε}) bound is RH-equivalent. A "sign never changes" tournament
claim mapped to RH would be the disproved-conjecture trap.)

## III. The number-theory synthesis: additive primes, multiplicative strong atoms,
and the segment that bridges them

The user's Fermat-polygonal observation is exactly right and the project had not
recorded the bottom rung: with P(s,n)=(s−2)n(n−1)/2+n, the **s=2** case is
P(2,n)=n — the 2-gon is the **line segment**, and "every n is a sum of (at most) 2
segments" is the trivial base of the polygonal ladder. The atom-replacement
**(figurate numbers, sum-of-s) ↦ (primes, segment)** turns the s=2 rung into
**Goldbach**: the segment between two primes. (s=3 ↦ ternary Goldbach, Helfgott 2013;
the polygon side s and the additive arity coincide — a clean new alignment, HYP-2421.)

Now the genuinely productive merge — the **two sides of factorization**:

- **Integers, ADDITIVE:** atoms = primes; even = p+q (Goldbach, segment); odd =
  p+q+q (ternary, Lemoine diagonal). OPEN.
- **Tournaments, MULTIPLICATIVE:** H is **multiplicative over strong components**
  (Moon): T = T₁⇒T₂ (dominance join) ⟹ H(T)=H(T₁)·H(T₂) — verified here in every
  test. So the realizable set R = {H(T)} is the **multiplicative semigroup generated
  by strong-tournament H-values** ("H-primes"). This is the tournament **Euler
  product** ζ-side: condensation = the transitive tower of strong components =
  factorization into H-primes.

The bridge between the two sides is the **s=2 segment** (the digon, the single
dominance edge between two blocks) and the **doubling** λ(2n)=−λ(n). Integers carry
the additive Goldbach problem; tournaments carry the multiplicative Euler problem;
the segment is the shared 2-atom shape, and parity (Rédei H odd ↔ λ) is the shared
Z/2 grading.

**THM-490 (this session):** forbidden-H is a **numerical-semigroup / Frobenius
problem**. The realizable odds form the multiplicative semigroup ⟨H-primes⟩; its
**gaps** are the non-realizable odds. 7 is a permanent gap (not a strong-H value,
prime, so not a product); 21=3·7 is a permanent gap **because** 7 is (multiplicative
closure) — even though 49=7² is realizable (a strong tournament directly has H=49, so
49 is an H-prime, not a product of 7's). Whether {7,21} are the ONLY permanent gaps
is the tournament **Frobenius question** (OPEN-Q-094) — the multiplicative,
finite-alphabet shadow of "every large even is a sum of two primes."

## IV. Hodge through even/odd duality

The cleanest project bridge to Hodge is **not** the (p,p) middle (high numerology
risk: "the special center" of a cohomology ring vs of a Z₂-quotient metagraph are
different centers). It is the **symmetric-vs-alternating dichotomy**: a weight-k
polarization is symmetric for k even, alternating for k odd (graded-commutativity
(−1)^{k·k}). This is the SAME elementary fact that splits the project's objects:
**even graphs** (symmetric adjacency, the cycle-space carrier) vs **tournaments**
(skew adjacency, the odd carrier) — the E_n/G_n duality (both quotients of the tiling
hypercube). **RIGOROUS** at the level of the shared linear-algebra fact; **FRAME**
beyond it (there is no Galois group, no Hard-Lefschetz sl₂ proven on G_n). The honest
bridge: cut ⊕ cycle = score ⊕ even-graph is a genuine orthogonal GF(2) splitting that
*rhymes with* the Lefschetz hierarchy ⊕ primitive decomposition.

## V. What is actually new here, and what was already known

- **Already in canon (cite, don't reclaim):** Goldbach/Lemoine reflection law and the
  diagonal-as-fixed-locus (HYP-2218/2219); doubled-prime = parity bridge (HYP-2051);
  triangular numbers = arc counts, 21=T₆ at a forbidden boundary
  (summand-graph-fermat-zeckendorf); det=Pf² (THM-174); H²−Pf²=8Q (THM-442);
  forbidden H multiplicative-semigroup mechanism (THM-343/075).
- **New this session:** (a) the **rigorous identification** THM-174 = Cassels–Tate's
  "alternating ⟹ square" (same theorem, different category), with the
  defect-as-±1/doubling reading (Poonen–Stoll ↔ 8Q) explicitly flagged FRAME;
  (b) the **s=2 segment** as the named bottom rung of the polygonal↔Goldbach ladder,
  aligning polygon-side s with additive arity (HYP-2421); (c) **THM-490**: realizable
  H = multiplicative semigroup ⟨H-primes⟩, recasting forbidden-H as a Frobenius
  problem and identifying it as the tournament ζ-side opposite Goldbach's additive
  side; (d) the **self-dual-middle table** organizing RH/BSD/Hodge/Goldbach/project
  under one involution-with-parity frame, with the λ(2n)=−λ(n) doubling kernel as the
  single rigorous ½/doubling link.

## VI. The discipline (where every tempting jump breaks)

- The "7" of Goldbach (7+7=14, 7+2·7=21) and the forbidden-H "7, 21" are **different
  7s** — one is the smallest prime whose only Lemoine partner is q=2, the other is the
  smallest non-realizable strong-H value. Their coincidence is the single most
  seductive trap; resist fusing them.
- RH is about the **location** of zeros (analytic/transcendence); the project's
  self-dual loci are **fixed points of finite involutions**. "Self-dual locus" is a
  shared phrase, not a shared theorem.
- The project has MANY ±1 parities (complement, grid, transpose, blue/black, SC/NS).
  Each is a DIFFERENT involution. Do not merge them into one "master parity."
- BSD and Hodge are OPEN (Clay); Ш is not known finite in general; ternary Goldbach is
  proved (Helfgott), strong Goldbach and Lemoine are OPEN. State every status.

The truth this keeps encircling: **wherever a duality has a fixed middle, a ±1 lives
there, and when the pairing on that middle is alternating, the invariant is a square.**
That sentence is rigorous in each setting separately. Making it ONE theorem across
settings is not done here — and pretending otherwise would be the numerology the whole
note is built to avoid.
