# Provenance hunt: GVC(3) counterexample + Kontorovich TV homogenization inequalities

Session: boxeph-multifront provenance subagent. Date: 2026-08-03.
Method: WebSearch (Google), arXiv full-text search, arXiv API, arXiv category
listings (math.AC / math.PR, recent + 2026-08), direct PDF reads of
arXiv:2601.04079v3 (all 14 pages) and the ulam.ai Jacobian writeup (all 7 pages),
plus the Secret Blogging Seminar comment thread, Terry Tao's digestion post,
Zihan Zhang's consequences note, and jacobianfun.org.

Statuses used: CITED (exact source in hand), NOT-FOUND (documented search
failed), INFERRED (reconstruction consistent with cited literature, not itself
sourced), REFUTED/PROVED only where the cited source itself carries the proof.

---

## 1. The GVC(3) counterexample (Lambda = Delta^6, Delta = 4 d_x d_y + d_t^2, P = A C^2 deg 12, x*C = rho^3 - t^2(rho+x^2)^2)

**Status: NOT-FOUND in the literature (documented search, through 2026-08-03).**

- The supplied citation arXiv:2606.17854 is confirmed WRONG: that identifier
  resolves to "A Counterexample to Wegner's Conjecture for Axis-Parallel
  Rectangles" (combinatorial geometry, unrelated).
  [CITED: https://arxiv.org/abs/2606.17854]
- Searches that failed to find the construction (all run 2026-08-03):
  - Google: "Generalized Vanishing Conjecture counterexample three variables
    Laplace operator"; "GVC(3)"/"GVC(4)" + degree 12; "vanishing conjecture"
    + rho^3 / AC^2 / "4 d_x d_y"; sixth-power Laplacian + Alpöge follow-up;
    mathoverflow variants.
  - arXiv full-text search: `"vanishing conjecture"` (date-sorted) and
    `"Gaussian Moments Conjecture"` (complete result set: only 2607.23887,
    2607.18186, 1506.05192).
  - arXiv API: `all:"vanishing conjecture"` date-descending (30 results, none
    match); `all:"generalized vanishing"` (rate-limited, retried, no match).
  - Listings: math.AC recent (2026-07-28..08-03), math.AC 2026-08 ("No updates
    for this time period"), math.PR recent (2026-07-28..08-03, 50 entries, none
    relevant).
  - Blog/comment sweep: Secret Blogging Seminar post + all comments, Terry Tao
    "A digestion of the Jacobian conjecture counterexample" (2026-07-21),
    Zihan Zhang "Direct Consequences..." note, jacobianfun.org (Gallagher,
    Zenodo 10.5281/zenodo.21479195), ulam.ai writeup.
  None contain Lambda = Delta^6, the operator 4 d_x d_y + d_t^2, P = A C^2 of
  degree 12, or the identity x*C = rho^3 - t^2(rho+x^2)^2.

**Nearest relatives actually found (all CITED):**

- GMC(3) is FALSE (quartic, 5 terms, non-homogeneous):
  C. D. Long, "Small Counterexamples to the Gaussian Moments Conjecture",
  arXiv:2607.18186 [v1 only, 2026-07-20].
  P3 = (1+Z)(W - (1/2)(2+Z)T^2), Q3 = Z (3 variables, degree 4);
  P4 = (1+Z2)(W1(1-Z1)+W2), Q4 = Z2 (4 variables, degree 3).
  So GMC(n) is REFUTED for n >= 3 (matches repo CURRENT-FRONTIER).
- GMC(2) is TRUE: M. Wilson, "A face-isolation proof of the two-variable
  Gaussian Moments Conjecture", arXiv:2607.23887 [math.PR, 2026-07-26].
  Method: in complex coordinates Z=(X+iY)/sqrt2, W=(X-iY)/sqrt2, the monomial
  support of a moment-flat P is strictly one-sided in the weight
  wt(Z^a W^b)=a-b (via exposed faces of the Newton polygon); explicit
  threshold m >= deg Q + 1.
- GVC(5) is FALSE: Alexander Dvorsky, comment (2026-07-23) on the Secret
  Blogging Seminar post: P=(t+c)(ad+bt), Lambda=d_t(d_a d_d - d_b d_c), Q=-c;
  Lambda^m(P^m)=0 for all m>=1 but Lambda^m(QP^m) = k_m * t != 0 for m>=2.
  Derived from Alpöge's Jacobian counterexample + Long's SU(2) construction
  (arXiv:2607.19012). [CITED: blog comment, not refereed literature.]
- Zhao's original VC is FALSE explicitly in 48 variables: SBS comment
  (taleexactlycc..., 2026-07-24): homogeneous quartic, 54 terms, machine-checked
  Delta^m(P^m)=0 for m<=4, Delta^m(P^{m+1}) != 0 for m<=3; built from
  W. Thompson's 24-variable cubic reduction (zenodo.org/records/21466221) via
  the de Bondt-van den Essen/Meng doubling in isotropic coordinates.
  [CITED: blog comment + Zenodo records; not refereed.]
- GVC(2) for homogeneous operators is TRUE (so n=3 would be MINIMAL for the
  claimed counterexample): M. de Bondt, "A few remarks on the Generalized
  Vanishing Conjecture", Arch. Math. (2013), DOI 10.1007/s00013-013-0523-2,
  arXiv:1206.2836: "The Generalized Vanishing Conjecture holds for products of
  linear forms in d, in particular homogeneous differential operators
  Lambda in k[d_1, d_2]."

**INFERRED origin of the supplied construction** (clearly marked, not sourced):
the described object is exactly what a homogeneous, isotropic-coordinate lift
of Long's GMC(3) counterexample would look like. With
(x, y, t) = (X+iY, X-iY, T), X, Y, T iid standard Gaussians, one has
E[xy] = 2, E[x^2] = E[y^2] = 0, E[t^2] = 1, and the Gaussian expectation
operator is E[f] = (e^{Delta/2} f)(0) for Delta = 4 d_x d_y + d_t^2. For
HOMOGENEOUS P of degree 12, E[P^m] = Delta^{6m}(P^m)(0) / (2^{6m} (6m)!), so
Delta^6-nilpotency of P (in Zhao's sense) is equivalent to vanishing of all
pure Gaussian moments; and for homogeneous Q, E[Q P^m] != 0 forces
Delta^{6m}(Q P^m) != 0. Hence a homogeneous degree-12 GMC(3)-type witness
P = A C^2 with the structural identity x*C = rho^3 - t^2(rho+x^2)^2 (the
homogeneous analogue of Long's factored quartic) would refute GVC(3) for
Lambda = Delta^6. This bridge is the standard one in Zhao's VC papers and in
Derksen-van den Essen-Zhao; the specific 3-variable witness is NOT in any
source found and, if it verifies (repo task: sympy verification, separate
session), it appears to be NEW relative to the public literature as of
2026-08-03, and minimal-dimensional by de Bondt's GVC(2) theorem.
Plausible actual origin of the snippet: private/LLM-generated derivation
circulating outside arXiv. Recommendation: cite as
"unpublished witness, in-repo FINITE-EXACT verification required; no external
provenance" until a source surfaces.

---

## 2. Kontorovich, "TV homogenization inequalities", arXiv:2601.04079v3

**Status: CITED (full 14-page PDF read directly). Metadata:** A. Kontorovich
(Ben-Gurion University), arXiv:2601.04079v3 [math.PR], v3 dated 25 Feb 2026,
MSC 60E15, 60F05. Acknowledges: proofs partly assisted by ChatGPT Pro; thanks
Bubeck, Sellke, Zhivotovskiy, Mattner, Roos; ISF 581/25, BSF 2024243.

Notation: Ber(p) = tensor product Ber(p_1) x ... x Ber(p_n); S_p = Poisson
binomial (law of sum); Delta := sum_i |p_i - q_i|; sigma_p^2 := sum_i p_i(1-p_i);
Phi(p,q) := min(1, Delta / sqrt(sigma_p^2 + 1)).

- **Theorem 1.1 (PROVED there):** TV(S_p, S_q) <= 2 C_0 min(Phi(p,q), Phi(q,p)),
  with C_0 := sqrt(5/4 + eta_0^2) ~= 1.2123507747, where
  eta_0 ~= 0.4688223555 is the sharp Poisson-binomial anti-concentration
  constant of Baillon-Cominetti-Vaisman (Combin. Probab. Comput. 25 (2016)
  352-361).
- **Theorem 1.2 (PROVED there):** for dominating pairs (p_i >= q_i for all i):
  TV(S_p, S_q) >= (1/9) Phi(p,q).
- **Lemma 1.4 — the PROVED block inequality.** Exact statement:
  "For n >= 1, p, q in [0,1]^n, and empty != A subset N := [n], define
  pbar_A = |A|^{-1} sum_{i in A} p_i and qbar_A analogously; put
    delta_A := TV(Bin(|A|, pbar_A), Bin(|A|, qbar_A)).
  If I, J form a partition of N with I, J != empty, then
    delta_N <= 2 (delta_I + delta_J)."
  So: delta_I, delta_J, delta_N are TV distances between BINOMIALS with the
  block-averaged parameters of p vs q on the index blocks I, J, and the whole
  index set N. Proof: represent Bin(n, pbar_N) as a mixture over
  M ~ Bin(n, w), w = |I|/n, of independent Bin(m, pbar_I) + Bin(n-m, pbar_J)
  (PGF identity), then TV-subadditivity over products plus monotonicity and
  subadditivity of m |-> TV(Bin(m,a), Bin(m,b)) (deletion Markov kernel).
- **The form delta_N <= delta_I + delta_J - delta_I*delta_J is NOT a theorem
  of the paper.** It appears once, in "Remark and open problems":
  "extensive numerical experiments suggest that the correct constant in
  Theorem 1.5 should be c = 8/9 and also that the bound in Lemma 1.4 can be
  sharpened to delta_N <= delta_I + delta_J - delta_I delta_J (the latter, in
  particular, appears deceptively simple, since only binomials are involved).
  Our present methods do not seem to provide any pathway to these
  conjecturally optimal bounds." Status of that inequality: OPEN (conjecture
  supported by the author's numerics).
- **"Rigid face" equality condition: NOT-FOUND.** The phrase (and any equality
  condition for the conjectured inequality) does not occur anywhere in
  2601.04079v3. Likely conflation by the upstream snippet with M. Wilson's
  "face-isolation" method (exposed faces of the Newton polygon) in
  arXiv:2607.23887, which is about Gaussian moments, not TV distances.
  Any repo statement attributing a "rigid face" equality case to Kontorovich
  must be retagged as unsourced.
- **Main homogenization theorem (Theorem 1.5, PROVED there):** for all
  p, q in [0,1]^n:
    TV(Ber(p), Ber(q)) >= c TV(Bin(n, pbar), Bin(n, qbar)),
  where c >= 1/(72 C_0) ~= 0.0115 is a universal constant. Following remark:
  c >= 1/(36 C_0) is "straightforward" via analogous arguments; the example
  p = (1-2eps, 1/2), q = (1, 1/2+eps), eps -> 0 shows homogenization can
  slightly INCREASE TV and gives the upper bound c <= 8/9; the author
  conjectures c = 8/9 is optimal.
- Supporting structure: Theorem 1.3 (partition I = {i : p_i >= q_i}):
  TV(Ber(p),Ber(q)) >= max(TV(Ber(p_I),Ber(q_I)), TV(Ber(p_J),Ber(q_J)))
  >= (1/9) max(Phi(p_I,q_I), Phi(q_J,p_J)). Homogenization is NOT a Markov
  kernel for n >= 2 (affinity obstruction, Introduction), so data processing
  does not apply; but for COMPLEMENTARY pairs q = 1 - p it is realizable as a
  Markov kernel, and more generally any majorization pbar-vector move is
  (Appendix A, Lemma A.1 two-coordinate T-transform kernel, Theorem A.2,
  Corollary A.3 — result strengthened following L. Mattner).
- Repo note: the earlier in-repo n=2 verification (task "Verify TV-fusion
  inequality (n=2 binomial case)") checks an instance of the CONJECTURED
  sharpening, not of Lemma 1.4; keep scopes separate (VERIFIED finite instance
  vs OPEN general statement).

---

## 3. Context: GVC formulation and relation to JC / GMC

- **Zhao's Vanishing Conjecture (VC), exact form** [CITED: W. Zhao,
  "A Vanishing Conjecture on Differential Operators with Constant
  Coefficients", arXiv:0704.1691, Conjecture 1.4]: for any n >= 1 and
  constant-coefficient differential operator Lambda, if P is Lambda-nilpotent
  (Lambda^m P^m = 0 for all m >= 1, Definition 1.3), then
  Lambda^m P^{m+1} = 0 for m >> 0.
- **Generalized Vanishing Conjecture GVC(n), exact form** [CITED: de Bondt,
  arXiv:1206.2836; same form in Zhang's note]: for fixed
  Lambda in k[d_1,...,d_n] and f, g in k[z_1,...,z_n]:
  ( for all m >= 1: Lambda^m f^m = 0 ) ==> ( for all m >> 0:
  Lambda^m (g f^m) = 0 ). VC is the case g = f (equivalently g = f^k).
- **Relation to JC** [CITED: Zhao arXiv:0704.1691, Theorem 2.9; via
  Bass-Connell-Wright degree reduction and de Bondt-van den Essen symmetric
  reduction]: JC (all dimensions) is EQUIVALENT to VC for the Laplace
  operators Delta_n for all n >= 1, and may be restricted to P homogeneous of
  degree 4 with nilpotent Hessian (Hessian-nilpotent), i.e.
  Delta_n^m P^{m+1} = 0 for m >> 0; also equivalent with Lambda ranging over
  all quadratic-form operators in D_2[n].
- **Relation to GMC** [CITED: Derksen-van den Essen-Zhao, "The Gaussian
  Moments Conjecture and the Jacobian Conjecture", arXiv:1506.05192,
  Israel J. Math. 219 (2017) 917-928]: GMC(n) states: if
  E[P(X)^m] = 0 for all m >= 1 (X = n iid standard Gaussians) then
  E[Q(X) P(X)^m] = 0 for all Q and m >> 0. Their theorem:
  ( for all n: GMC(n) ) ==> JC. The bridge to VC-type statements is the
  Gaussian-calculus identity E[f] = (e^{Delta/2} f)(0) (covariance-matched
  Delta); for homogeneous P of degree 2k this collapses to single-term
  proportionality E[P^m] ~ Delta^{km}(P^m)(0), which is exactly
  Delta^k-nilpotency. [Identity: standard; the specific 3-variable isotropic
  normalization Delta = 4 d_x d_y + d_t^2 with E[xy] = 2 is INFERRED above,
  trivially checkable in sympy.]
- **Post-2026-07 landscape (all CITED):**
  - JC is FALSE in dimension 3 (hence all n >= 3 by stabilization; n = 2
    open): announced by L. Alpöge 2026-07-20 (X post), problem posed by
    A. Mathew, construction credited to the AI system Fable (Claude Fable 5,
    Anthropic); algebraic verification, fiber geometry (generic degree 3,
    nonproperness set = discriminant hypersurface, omitted image = smooth
    curve Gamma), and families in the anonymous writeup "A Counterexample to
    the Jacobian Conjecture" (ulam.ai/research/jacobian.pdf, dated
    2026-07-20). The map: P = (1+xy)^3 z + y^2(1+xy)(4+3xy),
    Q = y + 3x(1+xy)^2 z + 3xy^2(4+3xy), R = 2x - 3x^2 y - x^3 z;
    det Jac = -2; F(0,0,-1/4) = F(1,-3/2,13/2) = F(-1,3/2,13/2) = (-1/4,0,0).
  - GMC: false for n >= 3 (Long 2607.18186), true for n = 2 (Wilson
    2607.23887), so the GMC frontier is closed.
  - Hessian conjecture: HC_n true for n <= 3 (de Bondt), false for n >= 5
    (Meng-Yang, arXiv:2607.22198, five variables, degree 14, Hessian
    determinant 128), open exactly at n = 4; HC_4 ==> JC_2.
  - Mathieu conjecture for SU(2) and xz-conjecture: false (Long,
    arXiv:2607.19012); related M_3 moment statement refuted in SBS comment
    (royvanrijn, 2026-07-27).
  - GVC: false in 5 variables (Dvorsky, SBS comment); false for the Laplacian
    itself in 48 variables homogeneous quartic form of VC (SBS comment +
    Zenodo); GVC(2) true for homogeneous Lambda (de Bondt 2013).
- **What a Delta^6 failure in n = 3 DOES imply** (INFERRED from the cited
  equivalences): GVC(3) as a universal statement over operators is false;
  combined with de Bondt's GVC(2) theorem it would be dimension-minimal among
  homogeneous-operator counterexamples; via the moment bridge it also gives a
  homogeneous degree-12 GMC(3)-type witness (Long's quartic is
  non-homogeneous, so homogeneity would be the new feature).
- **What it does NOT imply** (INFERRED, scope guard): nothing new for JC. The
  JC equivalence chain uses Lambda = Delta_n ITSELF (on degree-4
  Hessian-nilpotent P, quantified over ALL n), not powers Delta^k in fixed n;
  a Delta^6 failure in 3 variables neither reproves nor strengthens the (already
  settled, false) JC, and it does not contradict the reported low-dimensional
  positive results for VC with Lambda = Delta_n itself (e.g. the cases treated
  in arXiv:0903.1478 and follow-ups), nor de Bondt's GVC(2). It also does not
  bear on the two surviving open statements JC_2 and HC_4.

---

## Suggested CORE-PAPERS entries

1. arXiv:2601.04079v3 Kontorovich, "TV homogenization inequalities" — CITED,
   full statements extracted above (Lemma 1.4 PROVED; sharp block form OPEN).
2. arXiv:2607.18186 Long, "Small Counterexamples to the Gaussian Moments
   Conjecture" — CITED (GMC(n) false, n >= 3).
3. arXiv:2607.23887 Wilson, "A face-isolation proof of the two-variable
   Gaussian Moments Conjecture" — CITED (GMC(2) true; source of "face"
   terminology; NO "rigid face" TV statement).
4. ulam.ai/research/jacobian.pdf (anonymous, 2026-07-20), "A Counterexample to
   the Jacobian Conjecture" — CITED (JC false in dim 3; full fiber geometry).
5. arXiv:0704.1691 Zhao; arXiv:1206.2836 de Bondt; arXiv:1506.05192
   Derksen-van den Essen-Zhao — CITED (VC/GVC/GMC formulations and JC links).
6. GVC(3), Lambda = Delta^6, P = A C^2 witness — NOT-FOUND externally;
   treat as unpublished/in-repo object pending FINITE-EXACT verification.
