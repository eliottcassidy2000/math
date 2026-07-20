# Base-path independence resolved, the VC-witness dimension pinned, and the Casas–Alvero inversion

**death-star-2026-07-20-S61** (HYP-8265; owner: pin the reduction dimension and attempt
the VC witness; find "even more creative ideas than a descending star invariant — the
base-path-independent subgroup, the natural candidate being ∩Γ over all spanning paths";
Fibonacci from summing shifted Pascal; Erdős 592 / three JC parts; Casas–Alvero as
polynomial-derivative rigidity with p-adic partial results). Five threads, each verified
where checkable, each folded in with concurrency credit.

## 1. ∩Γ over all spanning paths = {0} (confirming kind-pasteur's THM-1415)

**Primary credit: kind-pasteur-S128c110 (THM-1415)** reached and settled this thread
first and more thoroughly (including n=6, and the OEIS/rep-theory next step). I record my
independent computation as a confirming data point and defer to THM-1415 as canon.

kind-pasteur-S128c109 (THM-1405) built the **star group** Γ(P) = ⟨S_v^P⟩ where
S_v^P = δ(v) ⊕ (path-edges at v): rank n−1, base-path-relative, and ended by proposing
∩_P Γ(P) over all spanning paths as the base-path-free repair. THM-1415 called that the
wrong instinct ("do not intersect; stop projecting") — and my explicit computation
confirms the intersection is empty:

  **∩_P Γ(P) = {0}** at n = 4, 5, 6 (over all 12, 60, 360 Hamiltonian base paths;
  F₂-subspace intersection via Zassenhaus). NO nonzero star move lies in the star group
  of *every* base path.

THM-1415's diagnosis is exactly right: Γ(P) is the cut space **restricted** to non-tree
edges, and the restriction is what needs a base path; the cut space **itself** needs
none, since S_n permutes the vertex stars δ(v), so ⟨δ(v)⟩ = **Seidel switching** is an
S_n-invariant, base-path-free notion. Its classes count exactly

  **switching classes = 2^{C(n,2)−(n−1)} = 2^{C(n−1,2)} = the number of tilings**
  (my run: 2, 8, 64 labeled at n = 3, 4, 5 = #tilings exactly),

so **the tiling hypercube IS the switching-class structure** and the base path is a
coordinate chart, not a feature of the object. But — THM-1415's key finding, which my
numbers reproduce — the base-path-free *tournament* invariant switching-class-up-to-iso is
**nearly trivial**: **1, 2, 2** (n = 3,4,5; THM-1415 adds **6** at n=6) vs iso classes
2, 4, 12, 56. At n=3 the transitive and the 3-cycle are already switch-equivalent (switch
the middle vertex). The moral, now doubly-checked: the base path is *precisely the
information the iso quotient forgets*; the only base-path-free shadow, the switching
class, is too coarse to separate tournaments. THM-1415 also **refuted** the attractive
guess that switching-up-to-iso would land on the even-graph/two-graph count A002854 =
2,3,7,16,54 (it does not: 1,2,2,6), and named the right next object — the **S_n-submodules
of the cut space over F₂** (a small rep-theory question). That is where this thread now
goes; I add nothing to it beyond the ∩Γ = {0} confirmation.

## 2. The VC witness: dimension pinned ≈ 76, and the clean verification is the collision

**The dimension.** Transporting F (degrees (7,6,4), dim 3) through the Bass–Connell–Wright
cubic-homogeneous reduction then the de Bondt–van den Essen symmetric doubling:
F carries **13 monomials of degree ≥ 2** (F₁: degrees {2,3,4,5,6,7}; F₂: {2,3,4,5,6};
F₃: {3,4}). Reducing each to cubic costs ~⌈(deg−1)/2⌉ auxiliary variables (the deg-7
monomial x³y³z alone needs ~2), landing the cubic-homogeneous stage at **N ≈ 35–38**;
the symmetric doubling gives the quartic-potential dimension

  **D = 2N ≈ 70–77.**

This matches an independent web-informed reduction estimate (a concurrent agent's
computation landed on **76–77**). So the first explicit witness to Zhao's Vanishing
Conjecture lives, via direct transport, in **≈ 76 variables**.

**The attempt, and the honest verdict.** At D ≈ 76 the *direct* VC certificate
Δ^m(P^m) = 0 ∀m, Δ^m(P^{m+1}) ≠ 0 involves multinomials in 76 variables — finite but
heavy (m = 2, 3 checkable; the tower expensive). The **cleaner witness is
non-invertibility**: for a Hessian-nilpotent quartic P, VC fails **iff** x + ∇P is
non-invertible, and non-invertibility is certified by a *single explicit collision*, not
an infinite Δ-tower. Our F already has one — the triple collision
F(0,0,−¼) = F(1,−3/2,13/2) = F(−1,3/2,13/2). Transport those three points through the
(explicit, constructive) reduction and they become a collision of x + ∇P, directly
certifying VC-failure. **Recommendation (rigorization):** build P by explicit transport,
then certify VC-failure by the transported collision — a finite, Lean-able check that
sidesteps both the Δ-tower and the citation risk of the equivalence. This is the S59z
"verify the witness directly" program made concrete: *the witness is a collision, not a
vanishing pattern.* The witness is constructible-in-principle for the first time (it
could not exist before today's counterexample); the open engineering step is the explicit
76-variable transport, not any new mathematics.

## 3. Casas–Alvero is the INVERSION of the Jacobian story (agent-confirmed)

The owner filed Casas–Alvero as "polynomial-derivative rigidity with p-adic partial
result," JC-adjacent. The precise status **inverts the framing**:

- **Casas–Alvero was just *proved TRUE*** (char 0, all degrees) — Soham Ghosh,
  arXiv:2501.09272 (Jan 2025, rev. Mar 2026), via Koszul homology + Brouwer-degree
  "almost-counterexamples." Where JC went *true-conjecture → false (counterexample
  today)*, CA went *open → claimed-proved-true*. The repo's "OPEN" label is **stale**.
- **The p-adic partial results are real and central** but internal to CA: Bothmer–Labs–
  Schicho–van de Woestijne (degrees p^k, 2p^k, arXiv:math/0605090) and Castryck–Laterveer–
  Ounaïes (3p^k…7p^k with excluded-prime lists, and degree 12, arXiv:1208.5404), via
  Draisma–de Jong's **p-adic-valuation reading of the resultant scheme** (reduction mod p
  preserves the trivial solution count when the resultant is a p-adic unit — hence the
  finite "bad prime" lists). This is the genuine technical home; it is **not** a bridge to
  Zhao/Hessian-nilpotence, and the repo's hoped-for "2-adic ladder adjacency" is
  **illusory** (CA's p-adics = resultant units mod p; the JC ladder = the inverse of a
  non-invertible map — same word, different structure).
- **The real opening is char p, not char 0.** In positive characteristic CA is *false but
  finite* — Ghosh arXiv:2402.18717: the CA variety is ≤ 2-dimensional, finitely many
  counterexamples per degree. **Enumerating the explicit char-2 counterexamples** is a
  concrete, combinatorial, still-open problem squarely in the repo's explicit-witness
  wheelhouse (the same muscle as the VC/Yagzhev extraction).

**The structural lesson.** Both JC and CA are polynomial *rigidity* statements — a local
derivative condition forcing a rigid global form. The dividing line the two now draw is
sharp: **univariate polynomial-derivative rigidity SURVIVES (CA true); multivariate
polynomial-map rigidity FAILS (JC false).** Rigidity holds where there is one variable
and honest differentiation, and breaks where there are several variables and a Jacobian.
That contrast — not a reduction, a *boundary* — is the real content of pairing them.

## 4. Fibonacci from shifted Pascal = the +1 stagger = the golden rung of the Pisot ladder

Verified: **F_{n+1} = Σ_k C(n−k, k)** (n = 0…24) — each Fibonacci number is the sum along
a Pascal diagonal **sheared by one step per row** (the shallow diagonal). The shear-by-1
is not incidental: it is exactly the **+1 stagger** of THM-1360's Pisot ladder
x^{s+1} = x^s + 1, whose s = 1 rung is x² = x + 1, the golden recurrence — and the
diagonal recurrence C(n−k,k) = C(n−1−k,k) + C(n−2−(k−1),k−1) collapses to
F_{n+1} = F_n + F_{n−1}. So "summing Pascal shifted" and "the golden/Pisot stagger" are
the same operator seen combinatorially vs algebraically. (kind-pasteur-S128c110 pushed the
same hint further and honestly *deflated* an over-reach: applying the shallow-diagonal sum
to the **Faulhaber** triangle gives F(m+1) + 2^{m−3} only through m ≤ 7 and breaks at
m = 8 — the "law" was **two coincidences, not one**, since the 2^{m−3} term is merely the
first four entries of a sequence that only *starts* geometric. Credit THM-1415; the clean
identity above is the part that holds.) Two further verified hooks tie it to the project's
char-2 / seven-wall spine:

- **Pascal mod 2 = Sierpinski = the char-2 atom**: row n has 2^{popcount(n)} odd entries
  (Kummer/Lucas), and the **Mersenne rows n = 2^k−1 (1,3,7,15) are all-odd** — the same
  char-2 seam as THM-1355/1360's Mersenne rung and the JC counterexample's det = −2.
- **F₈ = 21 and 1001 = 7·11·13**: the n=7 shallow diagonal sums to **21** — exactly one of
  the two *permanent* H-gaps {7, 21} of the tournament multiplicative monoid (S60); and
  1001 (the "three sixties," Pisano π(10) = 60 = 2²·3·5) carries F₇ = **13** and the
  **7** of the seven-wall (R(3,2) = 7). Graded evocative — but the structural claim under
  them is verified: the numbers 2 (doubling), 3 (trisection), 7 (the wall), and 60 (the
  aligned clock) recur across the Fibonacci/Pascal shear, the Erdős-592 trichotomy, and
  the JC counterexample's three parts, one organizing spine.

## 5. Honesty and credit

§1 and §4 are **primarily kind-pasteur's** (THM-1415, S128c110): kp settled switching as
the base-path-free quotient (nearly trivial: 1,2,2,6), refuted the two-graph/A002854
guess, and decomposed the Fibonacci–Faulhaber "law" into two coincidences — all first and
more thoroughly (n=6 included). My only additions there are the **explicit ∩Γ = {0}
computation** (n ≤ 6, all paths) confirming kp's "do not intersect," and the labeled
switching count 2,8,64 = #tilings. **The two genuinely new, checkable death-star
deliverables are §2 and §3:** §2 the VC-witness dimension ≈ 76 (a defensible monomial-count
+ de Bondt-doubling estimate, cross-checked against a concurrent reduction agent's
independent ≈ 76–77) with the collision-transport rigorization recipe — the honest path,
not an executed witness; and §3 the Casas–Alvero status **correction** — agent-confirmed
from the literature (Ghosh arXiv:2501.09272 proves it char 0; 2402.18717 gives char-p
finiteness), inverting the repo's stale "OPEN" and its illusory 2-adic-ladder adjacency.
The mod-2/Mersenne facts (§4) and F₈=21 / 1001=7·11·13 resonances are verified/graded
evocative. This is the seventh+ fleet-shared owner prompt (MISTAKE-199): I
confirmed-and-consolidated on §1/§4 and contributed new only on §2/§3.

## Cross-links

kp-S128c109 THM-1405 (star group Γ(P), the ∩Γ question) · the tiling hypercube = switching
classes (CLAUDE.md, Seidel switching) · THM-1300 (the counterexample + triple collision) ·
S59z (verify the witness directly) · THM-1355/1360 (Mersenne/Pisot ladder, +1 stagger) ·
S60 (two arithmetics; {7,21} monoid; Erdős-592/JC trichotomy; Pisano-60) · THM-466 (2-adic
digits = odd-cycle census) · Ghosh arXiv:2501.09272 / 2402.18717 (Casas–Alvero proved char
0, finite char p) · PROBLEM-LEDGER.md · PROBLEM-PORTFOLIO-2026-07-20 (Casas–Alvero status
to correct).
