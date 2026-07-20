# Odd-valued functions, even graphs, and the two parities that never meet — until S_n makes them

**death-star-2026-07-20-S61b** (HYP-8300; owner: see the relation between odd-valued
functions and tournament-adjacent ideas, and how both relate to *even* concepts — even
graphs and even functions; chase the high-leverage question). The answer is that "odd" and
"even" split into **two independent parities** — the parity of a function's *values* and
the parity of a function under the *complement involution* — and H, the odd-cycle census,
sits at their crossing, while the even-graph dual lives on a *third* parity (vertex-degree).
S_n is what ties them together, both globally and fiberwise (THM-1445).

## 1. H is an odd-valued even function

Two facts about the same object H(T) = #{directed Hamiltonian paths of T}, verified exactly
n = 3..6:

- **Odd-valued** (Rédei): H(T) is odd for every tournament T.
- **Even function**: H(T) = H(Tᵒᵖ) — H is invariant under the complement/negation involution
  R (reverse all arcs). Equivalently (HYP-534, reconfirmed here: max |odd-degree Walsh
  coefficient| = 0) H has **only even-degree Walsh coefficients** on the orientation cube.

So H is an **odd-valued even function**: the *values* have one parity (odd), the *function*
has the opposite parity-flavor (even under R). These are orthogonal — the value-parity is
about H mod 2 (≡ 1, the trivial character), the function-parity is about how H transforms
under R (fixed). A function can be even and odd-valued with no contradiction, and H is the
canonical example.

## 2. Two "evens," and they are genuinely different

The owner named *even graphs* and *even functions* in one breath; the sharp point is that
they are **different parities**:

- **Even function** = R-invariant = Walsh support on **even-cardinality** edge sets
  (|S| even). This is the *complement* parity. H lives here.
- **Even graph** = every vertex has **even degree** = element of the cycle space
  Z ⊆ F₂^{E}. This is the *degree* parity, the dual of the cut space.

These do not coincide: at n = 3 the triangle is an even graph with an **odd** number of
edges (|S| = 3), so it is *not* in the even-function (even-cardinality) support; conversely
the even-function support at n = 3 is {∅} ∪ {the three 2-edge sets}, none of which except ∅
is an even graph. H's Fourier support is **not** the even graphs — I checked, and the
attractive guess is false. "Even" is two words here: cardinality-parity (functions) and
degree-parity (graphs), and they must not be conflated.

## 3. Where they meet: S_n, globally and fiberwise (THM-1445)

They meet numerically, through the symmetric group.

- **Global bridge** (verified n=3..6): Σ_T H(T) = **n! · 2^{C(n−1,2)}**. The odd-valued H,
  summed over all tournaments, is |S_n| times the **even-graph count** 2^{C(n−1,2)}
  (= dim of the cycle space of K_n = #tilings = #switching classes). Each directed Ham path
  (n! of them) sits in 2^{C(n−1,2)} tournaments; that is the whole identity.
- **Fiberwise bridge (THM-1445, PROVED):** Σ_{T ∈ switching class} H(T) = **n!** for *every*
  switching class. Because a switching class is a coset of the cut space, and the cut space
  restricted to any spanning tree — a Hamiltonian path is a spanning tree — is a bijection
  onto F₂^{n−1}, each of the n! directed Ham paths is realized in **exactly one** member of
  each class. Summing counts the paths: n!.

Since switching classes **are** the tilings (mac-mini THM-474 / kp THM-1430), and at **odd
n** correspond one-to-one to **even graphs / two-graphs** (opus THM-1430, E_n = the
two-graphs), THM-1445 says:

  **The even-graph (switching) projection of the odd function H is the constant n!.**

This is the high-leverage sentence. All of H's tournament-specific content — the odd-cycle
collection, the H-spectrum, the {7,21} gaps — lives **inside the switching-class fiber**
(the tiling / 2-adic direction of S60's "two arithmetics"), and is **completely invisible to
the even-graph quotient**, which sees only |S_n|. The odd function is "purely fiber" over an
even base that is trivial. That is exactly why the even-graph side (E_n, two-graphs) and the
tournament side (G_n, H) are *dual* rather than *equal*: the duality is the fibration whose
base (even graphs) washes H out to n! and whose fiber carries everything.

## 4. The Jacobian counterexample is the same odd-fiber-over-even-base

The owner asked how these relate to the JC work, and they rhyme precisely:

- The Alpoge counterexample F has an **odd fiber** — the generic fiber is 3 points (the
  triple collision), and opus THM-1350 forces the fiber size odd (a Rédei-flavored parity).
  And **3 = H(C₃)**, the smallest odd H-value above 1: the JC fiber is literally the
  Hamiltonian-path count of the 3-cycle, the A₁ weight-triple = oriented 3-cycle of HYP-8160.
- Meanwhile **det JF = −2 is even**. So the witness is an **odd fiber over an even
  determinant** — the identical odd(fiber)/even(base) shape as H (odd values, even function;
  odd fiber-content over the even-graph base).

So one picture spans all three theaters: **an odd fiber sits over an even base**, with the
odd content invisible to the even quotient. Rédei's odd H over the even-graph base (THM-1445,
projection = n!); the JC map's odd 3-fiber over its even determinant; the odd-cardinality
cycle-space generators over the even-cardinality complement parity. The recurring number is
the small odd prime **3** (= H(C₃) = JC fiber = the first nontrivial odd) sitting over the
small even prime **2** (= |det| = the complement involution = the doubling of S60).

## 5. The high-leverage question

THM-1445 shows H is flat (= n!) on the even-graph quotient, so **H's parity content is not
an even-graph invariant** — it cannot be read off E_n. The leverage question is therefore:
*within the fiber*, what determines the odd value H(T)? The fiber is a coset of the cut space,
coordinatized by a base path (S61's ∩Γ = {0} result: the base path is exactly what the iso
quotient forgets, and there is no base-path-free star refinement finer than switching). So
the odd content of H lives precisely in the **base-path / tiling coordinates within a
switching class**, and the natural next object is the one THM-1415 named: the **S_n-submodules
of the cut space over F₂** — the representation-theoretic decomposition of the fiber on which
H varies. That is where "why is H this particular odd number" must be answered, and it is the
concrete high-leverage target this synthesis points to. (The JC analogue: what determines
*which* three points the odd fiber is — the trisection/Chebyshev modulus of THM-1335 — is the
same "content lives in the fiber" question in the polynomial theater.)

## 6. Honesty and credit

§1–§3 verified exactly (n=3..6); THM-1445 has a full proof. The two-parities distinction
(§2) and the false Fourier-support guess are recorded so the conflation is not repeated.
Concurrency credited and *not* duplicated: opus THM-1430 (switching graphs = two-graphs =
E_n; even-graph bridge odd-n only), kp/mac-mini THM-1430/474 (tilings = switching classes),
HYP-534 (H even function), HYP-3808 (fiber-parity checksum), opus THM-1350 (JC odd fiber),
HYP-8160 (JC 3-fiber = oriented 3-cycle), S60 (two arithmetics), S61 (∩Γ = {0}). §4 is a
structural rhyme (odd fiber / even base), not a theorem connecting the problems; §5 is a
posed target, not a result.

## Cross-links
THM-1445 · HYP-534 · HYP-3808 · opus THM-1430 / THM-1350 · kp THM-1430 / THM-1415 · THM-474 ·
THM-466 · HYP-8160 · S60 two-arithmetics · S61 ∩Γ={0} · CONJ-001 (Claim A, the parity core).
