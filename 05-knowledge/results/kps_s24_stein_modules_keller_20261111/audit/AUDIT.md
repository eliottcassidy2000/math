# Independent audit: Stein modules of Keller pairs (kind-pasteur-2026-11-11-S24; THM-5138, HYP-9668)

Audited: commit `86ab0707d6` (`README.md`, `PREREGISTRATION.md`, `scripts/`, `runs/`, `THM-5138-...md`,
`HYP-9668-...md`), and, at the coordinator's request, the uncommitted `ADDENDUM-5.2-DRAFT.md` (with
`scripts/stein_centres.py`, `runs/sn_centres.txt`). Date: 2026-11-11. Every program in this folder is the auditor's
(`tw.py` is an engine written from scratch: another basis of the ring of the slice, other primes, other filtrations).
The author's scripts were read for bugs and re-run only to check that the stored outputs reproduce (`repro/`: all
seven stored run files reproduce, timings aside). Nothing outside `audit/` was written.

## 0. Verdict

**Main packet: PASS WITH CORRECTIONS.** The proved statements (Propositions 2.1-2.3, Theorem 3.1 given the cited
comparison theorem, Theorem 4.1, Proposition 5.1, Proposition 6.1) are correct, and every number of section 6 was
reproduced by independent code (no number was found wrong). The corrections concern: one part of the hypothesis that is false as stated (M1),
missing credit and an unjustified "new" (M2, M4), failure claims on the slice that the computations do not support
(M3), heuristics stated as facts, one of them false on the toy (M5), and a remark that contradicts a canon theorem
the note itself cites (M6). Theorem 4.1 can be sharpened substantially (section 2, m1).

**Addendum (Theorem 5.2, Corollary 5.3): PASS WITH CORRECTIONS** (section 7). The proof is correct. It forces the
wording changes the coordinator anticipated, and it proves more than it states (detection at every very general
centre).

## 1. MAJOR findings

### M1. HYP-9668 (E) is false as stated; the note's own example `W = x^2 y` refutes it

*Location.* HYP-9668: title ("(E) Rank-one potentials: if H^2(C^2, d - dW ^) is one-dimensional, the kernel of its
functional is a Mathieu-Zhao subspace"), status ("(E): ... nothing is known about the Mathieu property there"),
section "(E) Rank-one potentials" ("this happens for x^p y^q", "The question is whether the Mathieu property ... is a
property ... of every one-dimensional twisted cohomology"). README section 8, bullet "Rank-one potentials"
("Whether the kernel of such a functional is Mathieu-Zhao is a natural extension of GMC(2) ...; for monomial
potentials it is a weighted-charge statement of DvdK type").

*Finding.* A Mathieu-Zhao subspace `M` of a commutative algebra that contains 1 is the whole algebra (take `a = 1`:
`a^m = 1` lies in `M` for all `m`, so `b = b a^m` lies in `M` for every `b`). For `W = x^2 y`

    1 = 2 (d/dy - W_y)(y) - (d/dx - W_x)(x)          [ (d/dy - x^2) y = 1 - x^2 y,  (d/dx - 2xy) x = 1 - 2 x^2 y ]

so `1` lies in `M = Im(d/dx - W_x) + Im(d/dy - W_y)`, while `M` has codimension one: the functional
`x^k y^l -> l! [k = 2l+1]` of the note vanishes on `M` and equals 1 on `x`. So `M` is not a Mathieu-Zhao subspace.
The same holds for every `x^p y^q` with `gcd(p, q) = 1`, `(p, q) != (1, 1)` (the class of 1 is zero by weight).
Moreover "this happens for `x^p y^q`" needs `gcd(p, q) = 1`: for `x^2 y^2` the fibre has two components and the
dimension is 2 (for `x^3 y^3`: 3). Checked: `a9_rank_one.py` (identity exact; class of 1 zero for `x^2 y`, `x^3 y^2`,
`x y^2`, nonzero for `xy` and for `x + x^2 y`), `a4_controls.py` (`x^2 y^2 -> 2`, `x^3 y^3 -> 3`, `x^3 y^2 -> 1`).

*Correction.* HYP-9668 (E): "REFUTED as stated (W = x^2 y: 1 lies in the kernel). Repaired question: if
`H^2(C^2, d - dW ^)` is one-dimensional **and the class of 1 is not zero** (e.g. Broughton's `x + x^2 y`, functional
value -1 on 1), is the kernel of the functional a Mathieu-Zhao subspace?" Replace "this happens for `x^p y^q`" by
"for `x^p y^q` with `gcd(p, q) = 1`". Delete "for monomial potentials it is a weighted-charge statement of DvdK
type" in README section 8 (for monomial potentials other than `xy` the statement is false).

*Side remark (not needed for the audit; novelty not checked beyond one search).* The same three lines contradict the
"General Image Conjecture" (M2) in characteristic zero for `A = C`, `n = 2`, `q = x^2 y`; van den Essen's survey
(arXiv:1006.5801, section 4, Comment 1) says that no counterexample in characteristic zero was known in 2010. If this
is recorded anywhere it should be cited; if not, it deserves its own sourced note, not a line in this packet.

### M2. Missing credit: K-GMC is an instance of Zhao's Image Conjecture; the "new conjectures" are not new as conjectures

*Location.* README line 9-10 ("what is new is the ladder itself ..."), section 0 item 4 ("The three weaker rungs are
new conjectures"), section 1 (no mention), section 5 ("The conjectures between GMC(2) and JC(2)"), section 8;
HYP-9668 (title, "The ladder"); THM-5138 status. No literature search is named for any novelty claim (AGENTS.md).

*Finding.* I read the source (arXiv:0902.0210; text extracted with pdftotext):

- W. Zhao, *Images of commuting differential operators of order one with constant leading coefficients*, J. Algebra
  324 (2010) 231-247. Conjecture 1.3 (Image Conjecture): for a commuting set `C` of operators `h(z) + sum c_i d_i`
  on `A[z]`, `Im C` is a Mathieu subspace. Theorem 2.5: it suffices to treat `C = {d_i - d_i(q)}`, `q` in `A[z]`,
  i.e. "the image of `d - dq ^` on `(n-1)`-forms is a Mathieu subspace". Theorem 2.10 / Corollary 2.11: for `n`
  operators in `n` variables over a field, `K[z]/Im C` is holonomic, hence finite dimensional.
- For a Keller pair with `{f, g} = 1`, `du ^ dv = dx ^ dy` and `M_F = Im(d/dx - W_x) + Im(d/dy - W_y)`,
  `W = (f - a)(g - b)`. So **K-GMC is literally the case `A = C`, `n = 2`, `q = (f - a)(g - b)` of Zhao's
  conjecture**, HYP-9668 (E) is its case "one-dimensional quotient", and "`S_F` is finite dimensional" is Zhao's
  Corollary 2.11. GMC as a statement about images of first-order operators is Zhao (J. Pure Appl. Algebra 214
  (2010), Hermite case `d/dx - 2x`) and Derksen-van den Essen-Zhao (Israel J. Math. 219 (2017): `GMC(2n) => SIC(n)`
  and the all-`n` implication to JC), both already in `CORE-PAPERS.md`.
- The nearest published equivalence of JC(2) with a Mathieu statement about a Keller pair: van den Essen, Wright,
  Zhao, *Images of locally finite derivations of polynomial algebras in two variables*, J. Pure Appl. Algebra 215
  (2011) 2130-2134 (arXiv:1004.0521), Theorem 4.3: JC(2) iff `Im D` is a Mathieu subspace for every divergence-free
  derivation `D` of `k[x,y]` with `1 in Im D` (proof through Y. Stein 1989 and Cerveau). The note's Corollary 4.2 and
  the draft Corollary 5.3 are statements of the same type for the twisted pair `(d_u - V, d_v - U)`.
- In the repo: the Fock / vacuum-Mathieu programme HYP-8350 (boxeph-2026-07-20-S166;
  `INDEX-HISTORICAL-THROUGH-2026-07-21.md:251`; `THM-1830-conditional-map-what-assuming-the-2d-nullcone-conjecture-buys.md`
  lines 70-76): "GMC(2) IS a vacuum-Mathieu statement for A_1", with the named micro-step "express ... vacuum-moment
  pairings for an A_1-endomorphism phi and identify which GMC(2)-instance controls surjectivity". The README refers
  to it only as "THM-1830 lists one as a superseded programme" (and THM-1830 is an ID shared by four files).

Not found in the literature (one sweep by a sub-agent, sources listed in section 6): the product potential
`(f - a)(g - b)`, the count `b_1 = c + sum delta_i (deg_u + deg_v)`, the bound `2d - 1`, the finite-chaos criterion.
They look new, and elementary given the cited theorems.

*Correction.* Add to README section 1 a bullet with the three references above and HYP-8350 by name; in section 0
item 4, section 5 and HYP-9668 replace "new conjectures" by "K-GMC is the instance `q = (f - a)(g - b)` of the
Image Conjecture of Zhao (2010) for two variables over C (a conjecture that fails for general `q`, see (E)); CV and
SN are consequences of it singled out here"; name the literature search in the status of THM-5138.

### M3. "SN, CV, K-GMC fail" on the slice is not supported; only finite chaos is shown to fail there

*Location.* README section 0 item 5 ("the lowest rung already fails" for both maps), section 6.2 last sentence ("SN,
CV, K-GMC and finite chaos all fail, in both directions"); THM-5138 title, part (E) ("... and one-sided Stein
nilpotence fails"); HYP-9668 status (D) ("on the slice ... it fails in both directions (classes of w V^m, w U^m
nonzero for m <= 5)") and "Why detection is plausible".

*Finding.* SN is a statement about all large `m`. The runs show nonzero normal forms for `m <= 5` (author) and, with
the audit engine, for `m <= 9` (`w V^m`, `x V^m`) and `m <= 12` (`w U^m`, `x U^m`) (`a7_slice_stability.out`). This
is evidence, not a failure of SN, CV or K-GMC; and a nonzero normal form is a certificate of non-membership only
relative to the truncation (section 4). What is proved on the slice: finite chaos fails at every centre where
`dim S_F > 3` (Theorem 4.1 and the cited dimension 11). On the toy the failure of SN and K-GMC is proved (residues),
for centres with `a != 0`.

*Correction.* Section 6.2: "finite chaos fails (Theorem 4.1: 11 > 3). No vanishing is observed: the classes of
`w V^m`, `x V^m` (`m <= 9`) and `w U^m`, `x U^m` (`m <= 12`) are nonzero modulo the generators of weight at most 70
(FINITE-EXACT relative to the truncation); the failure of SN, CV, K-GMC at the centre (2, 3) is not proved."
THM-5138 title (E): "... and no one-sided Stein nilpotence is observed (classes nonzero for m <= 5); on C^* x C its
failure is proved". HYP-9668 (D): "fails" -> "no vanishing observed for m <= 5 (audit: m <= 9, 12)". If the draft
Theorem 5.2 is adopted and stated for etale maps of surfaces (section 7), one may add the proved statement: "on the
slice, SN fails at every centre outside a countable union of curves".

### M4. The cited comparison theorem was already in the repo, as the known mechanism for Keller pairs

*Location.* README section 1: "CITED (not in the repo before): A. Dimca, M. Saito ..." and "A grep of canon and notes
for 'twisted de Rham', 'Dimca', 'bifurcation set' returned nothing".

*Finding.* `05-knowledge/results/logic_lenses_audit_A_20261012.md:32` (commit 55871a3968, 2026-10-12, the audit of
THM-4620): "Cerveau, J. Algebra 195 (1997) ..., Prop. 1.2, is the de Rham interpretation of surjectivity;
Dimca-Saito, Compositio 85 (1993) 299-309, is behind it." Also `CORE-PAPERS-PLANAR-JC.md:16` (Bonnet-Dimca). The
grep claim is true for "twisted de Rham" and "bifurcation set" and false for "Dimca". The use of Dimca-Saito for a
Jacobian pair is Gurjar-Masuda-Miyanishi 2013 (their abstract and Lemmas 3.7, 3.10: `H^(i+1)(Omega^., D_f)` is the
reduced `H^i` of the general fibre, applied to `D = f_y d_x - f_x d_y`), i.e. the case `W = f` of Theorem 3.1.

*Correction.* "CITED (already cited in the repo: `logic_lenses_audit_A_20261012.md:32`, as the theorem behind
Cerveau 1997 and Gurjar-Masuda-Miyanishi 2013, Lemma 3.7, for the potential `W = f`)"; correct the grep sentence;
credit GMM in the statement of Theorem 3.1 ("the case W = f is GMM 2013").

### M5. Heuristics stated as facts; one of them is false on the toy

*Location.* (i) README section 0 item 5: "The extra functionals of a counterexample are not Gaussian at all; they sit
on the non-properness curve, and the obstruction they create does not depend on the degree of the map." (ii)
HYP-9668, "Why detection is plausible": "at most `d` can be local Gaussians at the preimages of the centre; the
others come with the missing sheets. In the two examples on other surfaces these extra functionals are exactly what
breaks (C)". (iii) README section 7: "one residue functional ... for each puncture"; "It vanishes for all large m
only if ... which fails at the branch points"; "leaves no room for that when a is generic". (iv) README section 6.2:
"the Stein space is 3 + 8: ... eight further functionals come with the missing sheets".

*Finding.* (ii) is false for the toy. There `dim S_F = 2` equals the number of preimages of the centre, and
`lambda_+- = (1/2) L_(q_+-)`, where `L_q(h) = sum_i i! h_ii` is the local Gaussian functional at the preimage `q`
(`h = sum h_ij U^i V^j` the Taylor expansion at `q`; the sum is finite because `B` is polynomial in `V`). Checked
symbolically on `x^k y^l`, `|k| <= 5`, `l <= 4` (`a2_toy_residues.py`, item 4). So the toy has **no** extra
functional, and SN is broken by the two local Gaussians themselves: `L_q(x^(-2) V^m) = (d/du)^m (1/u)` at `u = a`,
the Taylor series of `1/u` at `a` does not terminate. For the slice the computations do not say which functionals
break SN. "At most `d` can be local Gaussians" is not proved (no local Gaussian functional on `B` is constructed for
`X = C^2`). (iii) "One residue functional for each puncture" is wrong when a puncture is ramified, which is the case
in both examples: on the toy the curve `v = v_0` has one finite puncture (`x = 0`, ramification 2) and the fibrewise
twisted `H^1` has dimension 2; on the slice two sheets escape through one puncture (`x ~ L^(-1/2)`).
(i) is the heuristic of section 7 placed in the summary without its label. (iv) `11 = 3 + 8` is a count (three
nondegenerate critical points of `W`, the preimages of the centre, against the dimension), not a proved
decomposition of the space of functionals; the same count on the toy is `2 + 0`.

*Correction.* (iv): "11 = 3 + 8, where 3 is the number of critical points of `W` (HEURISTIC reading: three local
Gaussians and eight functionals attached to the missing sheets)". (i): prefix with "Heuristically (section 7)". (ii): "On the toy the two Stein functionals are the
local Gaussians at the two preimages of the centre, and they break (C) because `1/u` is singular on the missed line;
on the slice (dimension 11, three preimages) the computations do not identify the functionals that break (C)."
(iii): "of dimension `sum_i delta_i deg_u S_i` (missing sheets; a puncture of ramification `e` carries `e` of them)";
"should fail"; "should leave no room".

### M6. Remark 3.4 contradicts THM-4122, and the restriction on the centre is empty for Keller pairs

*Location.* README Remark 3.4 ("It does not see the degrees of f and g: replacing (f, g) by (f + g^N, g) changes the
potential, not the cover"); Theorem 3.1 ("all centres outside finitely many lines"); section 0 item 5.

*Finding.* `THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction.md`, statement (4) (cited in
README section 1): every component `S_i` of the non-properness curve of a planar Keller map has normalization `A^1`,
and its two target coordinates have pole orders `(rho_i deg f / G, rho_i deg g / G)` at infinity, `rho_i >= 1`,
`G = gcd(deg f, deg g)` (Nguyen Van Chau 2004). Hence `deg_v S_i = rho_i deg f / G`, `deg_u S_i = rho_i deg g / G`,

    dim S_F = c + ((deg f + deg g)/G) * sum_i delta_i rho_i,

so the count does see the degrees, through the reduced pair; and replacing `(f, g)` by `(f + g^N, g)` changes the
bidegrees of the `S_i` and the count. Second consequence: each `S_i` goes to `(infinity, infinity)` only, so the
closure of `S` never contains `(a, infinity)` or `(infinity, b)`: **for a Keller pair on `C^2` every centre is
allowed.** (The condition matters on other surfaces: on the slice `b = 0` is excluded; verified, `a3`, `a6`.)

*Correction.* Remark 3.4: "By THM-4122 the bidegree of `S_i` is `rho_i (deg g, deg f)/gcd(deg f, deg g)`, so
`dim S_F = c + ((deg f + deg g)/gcd) sum_i delta_i rho_i`: the count depends on the degrees only through the reduced
degree pair." Theorem 3.1: add "(for a Keller pair on `C^2` this holds for every centre, by THM-4122)". Section 0
item 5: delete or mark as heuristic "does not depend on the degree of the map".

## 2. MINOR findings

- **m1. Theorem 4.1 is correct but far from sharp; the sharp form needs only `S_F != 0` and no condition on the
  centre.** *Theorem 4.1' (audit).* `X` connected smooth affine surface, `F : X -> C^2` etale, `K = Frac(B)`,
  `k = C(u, v)`. If `F` has finite chaos at `(a, b)`, then `S_F = C [1]`; and if `[1] != 0`, then
  `Tr_(K/k)(h)` lies in `C[U, V]` for every `h` in `B`, `B` is integral over `C[u, v]`, and `F` is an isomorphism
  onto `C^2`. So finite chaos at `(a, b)` holds iff `S_F = 0` or `F` is an isomorphism.
  *Proof.* `Theta : B -> S_F (x) P` is onto (the note's argument) and `P = C[U, V]` (standard module). For a linear
  form `l` on `S_F`, `tau = (l (x) 1) Theta : B -> C[U, V]` is `C[U, V]`-linear and commutes with `d_u, d_v`; extend
  it to `K -> k`. The trace form of `K/k` is nondegenerate, so `tau(h) = Tr(kappa h)` for one `kappa` in `K`; the
  trace commutes with the derivations, so `Tr((d_u kappa) h) = d_U tau(h) - tau(d_u h) = 0` for all `h`, hence
  `d_u kappa = d_v kappa = 0` and `kappa` is a constant. With `h = 1`: `tau = l([1]) Tr / d`. If `[1] = 0` then
  `Theta = 0`, so `S_F = 0`. If not, `Theta(h) = [1] (x) Tr(h)/d`, so `S_F = C [1]` and `Tr(B)` lies in `C[U, V]`;
  the power sums `Tr(h^j)` are polynomials, so `h` is integral over `C[u, v]`; `F` is finite and etale, a covering
  of `C^2`, an isomorphism. QED.
  For a Keller pair `S_F != 0` at **every** centre (`b_1(Y) = c + sum_p delta_p >= 1` by the cited theorem), so
  "finite chaos at one centre, any centre, iff automorphism", and the comparison theorem enters only through
  `S_F != 0`. Consistency: toy at `a != 0`: dimension 2, finite chaos fails (Proposition 6.1); toy at the centre
  `(0, b)`: `S_F = 0` (machine `a5_toy.out`; by hand `x^(i+1) y^j = j x^i y^(j-1)` modulo `M`), finite chaos holds
  trivially for a map that is not an isomorphism; open immersion `C^* x C -> C^2`: 1 at `(2, 3)`, 0 at `(0, 3)`.
  So the inequality `dim S_F <= d` has no counterexample, but it does not decide the double cover (2 <= 2), while
  4.1' does. Recommended: replace the rank comparison by 4.1'.
- **m2.** Theorem 4.1: the hypothesis "`S_F` finite dimensional" is not needed (a `d`-dimensional space cannot map
  onto `S_F (x) Frac(A)` otherwise); it holds anyway (Zhao, Corollary 2.11, or holonomicity).
- **m3.** Theorem 3.1: `c = 1` for generic centres (Bertini: the linear system `<1, f, g, fg>` has no base point and
  maps `C^2` onto a surface, so its general member `(f - a)(g - b) = t` is irreducible). On the slice the
  irreducibility of the fibre was checked for the centres (2, 3), (0, 0), (5, -7) (`a3`, item 7). Proposition 5.1
  carries "PROVED" (and Corollary 4.2 no status) without saying that "automorphism <=> finite chaos" rests on the
  CITED theorem: write "PROVED from the CITED comparison theorem" for that equivalence.
- **m4.** THM-5138 (E), toy: add "for centres with `a != 0`" (at `a = 0`, `S_F = 0`).
- **m5.** THM-5138 status (E) says the two-dimensional space is "PROVED by hand"; the README proves `>= 2` by hand
  and labels the dimension FINITE-EXACT. A hand proof of `<= 2` exists and should be added: modulo `M'`,
  `D_y(x^(i-1) y^(j+1)) - D_x(x^i y^j) = 2a x^(i-2) y^(j+1) - (i-j-1) x^(i-1) y^j - 2b x^(i+1) y^j` lowers the degree
  in `y` (`a != 0`), and `D_y(x^i) = a x^(i-1) - x^(i+1)` reduces `x^i` to `1, x` (`a2`, item 6).
- **m6.** README section 6.1: "`c = 2` at every truncation" -> "`c(N, N') = 2` for `N < N'`" (`c(12, 12) = 3` in
  `runs/s2_toy.txt`).
- **m7. Citations for the smooth affine case.** Sabbah, Nagoya Math. J. 141 (1996), Theorem 1.1 and the paragraph
  after it, is the statement used (relative cohomology of `(U, f^(-1)(t))`, "proved in [Dimca-Saito], ... when
  `U = C^n`, but the proof extends easily"); its Corollary 1.2 (`chi = r [chi(U) - chi(f^(-1)(t))]`) is the Euler
  characteristic form of the count and should be cited for Remark 3.3. Hien, Invent. Math. 178 (2009) is the period
  pairing for general flat connections and does not state the pair `(X, W^(-1)(t))`; cite Hien-Roucairol, Bull. SMF
  136 (2008) (period pairing for `e^g` on smooth affine varieties) or Fresan-Jossen instead. Dimca-Saito: confirmed
  (Theorem 0.3: "holds for any polynomial f"; `D_f = d - df ^`; their Example 2.11 is `x^2 y + x`). The statements
  are as the note uses them (reported by a sub-agent from the Numdam scan and from Braverman, GMM restatements; I did
  not read the proofs).
- **m8.** README section 1: THM-2633 says the inertia groups "normally generate" the monodromy group, not "generate".
  Cite colliding IDs with slugs (AGENTS.md): THM-1830 (four files), THM-1345 (two files).
- **m9. Names.** `[h V^m] = [d_u^m h]` is iterated Gaussian integration by parts (Stein's lemma); "Cameron-Martin"
  fits the exponential form `[h e^(sV)]`. "Stein kernel" has another meaning in Stein's method (Ledoux-Nourdin-
  Peccati), and "Stein" in the Jacobian literature is Yosef Stein (1989), whom the note cites two lines above its
  "Stein operators": say once that the name refers to Charles Stein's lemma.
- **m10.** Corollary 3.2: with Zoladek (Topology 47 (2008)) the generic degree of a counterexample is at least 6
  (quoted as known by Makar-Limanov, arXiv:2106.06869; the gcd results of the same paper are disputed by
  Guccione-Guccione-Valqui, arXiv:1401.1784; not in `CORE-PAPERS-PLANAR-JC.md`), which gives 11 instead of 9.
  Optional; if used, add the source to CORE-PAPERS first.
- **m11.** Section 5, "the Mathieu property is blind to UV": accurate only as "`UV` is not null". K-GMC does reach
  large radial exponents at nonzero charge (`P = V(1 + c (UV)^j)` gives `[h U^(jk) V^(m+jk)] = 0` for
  `m >= N_j(h)`, `k <= m`); what these families do not reach is charge zero. Reword.
- **m12.** PREREGISTRATION outcomes: S4 "where the truncation is exact" and Note 2 "shows that the kernel is complete
  in those weights" -> "where the truncation is stable" (section 4). S4 should read "CONFIRMED relative to the
  truncation".
- **m13.** HYP-9668, "A proof of GMC(2) that only used integration by parts in the coordinates ... would prove (A)":
  no such proof exists, since every identity of Proposition 2.2 holds on the toy, where K-GMC fails. Any proof must
  use `B = C[U, V]`. Delete or reword.
- **m14.** README section 1, Session267: the description is accurate but omits what is closest to this note:
  "Gaussian rotation invariance, or integration by parts, gives E[R P] = 0" and the one-dimensional rotational
  quotient (`planar/PROOF.md` sections 1-2), the "positive cone theorem" (a homogeneous sector of the Mathieu
  property, section 3), and the conditional implication G5. They concern the standard Gaussian, not a Keller pair;
  CONCURRENT is the right label; add one sentence.
- **m15.** HYP-9668 (A)-(D) do not say for which centres ("a centre", "a generic centre"). After Theorem 5.2 this
  quantifier decides whether the statement is open or equivalent to JC(2) (section 7).
- **m16. Observation (FINITE-EXACT, `a11_slice_centre00.out`).** At the excluded centre (0, 0) of the slice the
  potential `uv` is invariant under `(x, y) -> (-x, -y)` and the classes of the even elements vanish from some `m`
  on in the direction `d_u`: `[w V^m] = 0` for `m = 4, 6, 8`, `[x V^m] = 0` for `m = 1, 3, 5, 7, 9`,
  `[y V^m] = 0` for `m = 3, 5, 7, 9`, `[w^2 V^m] = 0` for `m = 6, 8`, while the odd elements and all `h U^m` stay
  nonzero. At (2, 0) nothing vanishes. A single-centre detection statement has to exclude or survive such centres.

## 3. What was recomputed

Engine `tw.py`: exact integer ring arithmetic; ring of the slice in the torus-weight basis `x^e q^n`, `y^e q^n`
(`q = 1 + xy`, `n` in Z), self-tested against sympy; elimination modulo 2147483647 and 1000000007 (author:
2147483629, checked prime). With the filtration (1, 2, 1) the spaces `V_N` are the author's, so `c(N, N')` and the
number of pivots must agree exactly: they do.

| quantity | author | audit |
|---|---|---|
| `xy`; `x^3 + y^3`; `x + x^2 y`; `x^2 y` | 1; 4; 1; 1 | 1; 4; 1; 1 (two primes; second weighting for the last two) |
| `fg`, `(x + y^2, y + (x+y^2)^3)`; `(f-2)(g-3)`, `(x + y^2, y + (x+y^2)^2)` | 1; 1 | 1; 1 |
| `x^2 y^2`; `x^3 y^3`; `x^3 y^2` | not computed | 2; 3; 1 |
| class of 1 for `x^2 y` | not computed | zero (exact identity) |
| toy, centre (2, 3): dimension | 2 | 2 (three filtrations, two primes; `<= 2` and `>= 2` by hand) |
| toy: `lambda_+-` kill the kernel | by hand | exact, `x^k y^l`, `-7 <= k <= 7`, `l <= 5`, symbolic `(a, b)` |
| toy: `lambda_+-(x^(-2) V^m)` | `(-1)^m m!/(2 a^(m+1))` | same, symbolic, `m <= 8` |
| toy: class of `x^(-2) V^m` | nonzero, `m <= 6` | `= ((-1)^m m!/2^(m+1)) [1] + 0 [x]` exactly, `m <= 8` (machine against formula) |
| toy: `x^(-1) V^m`; `x U^m`, `x^(-2) U^m` | nonzero; zero | nonzero; zero (`m <= 6`) |
| toy: `[U^i V^j] = i! [i=j] [1]` | 15 of 15 | 21 of 21 (`i + j <= 5`) |
| toy, centres (5, -7); (0, 3) | not computed | 2; 0 |
| THM-1300: `det J`; slice `u`, `v`; Jacobian | -2; formulas; `-2 w^3` | same (sympy, from the canon formulas) |
| master identity; `S_c`; bidegree | as stated; (3, 2) | same; irreducible, `disc_v = (u^2 - 12)^3`; boundary points `(inf, 0)`, `(inf, inf)` |
| fibre sizes of the slice map | "two sheets missing" (not tested by the scripts) | 3 at 12 rational points off `S_1` (4 of them on the conic `27v^2 - 9uv + 8 = 0`, the other factor of the discriminant); 1 at 16 rational points of `S_1`; numerically two preimages escape near `S_1` |
| `S_1` and the fibre images | 5; 2, 3, 4 | 5 (two centres); 2, 3, 4; centre (0, 0): 4 |
| slice (2, 3): dimension, pivots | 11 for `N <= 36`; 2999, 3866 at `N' = 40, 46` | 11 for `4 <= N <= 60`; 2999, 3866, 5924, 8414 at `N' = 40, 46, 58, 70`; 11 with two other filtrations |
| other centres | 11 at (5,-7), (1,1), (-3,4); 9 at (0,0) | same; 11 at (7,2), (-1,-2), (0,3); 9 at (2,0), (4,0) |
| `W = u, v, u + v` | 2, 4, 6 | 2, 4, 6 (two primes); `u - 2v`: 6; `2u + 3`: 2 |
| slice: `[U^i V^j]`, `[1]` | 15 of 15; nonzero | 15 of 15; nonzero |
| `w V^m`, `w U^m` (`m <= 5`); `x V^m`, `x U^m` (`m <= 4`) | nonzero | nonzero; also `m <= 9` (V) and `m <= 12` (U); normal forms identical at `N' = 46, 58, 70` |
| identity (d) on the slice | 21 checked, 0 failures | 27 checked, 0 failures |
| Hermite classes of `w`, `x`; of 1 | 22 nonzero each; only (0,0) | 25 nonzero each; only (0,0) |
| `runs/sn_centres.txt` (draft) | nonzero, `m <= 3`, six centres | nonzero, `m <= 4`, same centres |

No number of the packet was found wrong.

## 4. The truncation: what the runs prove

`c(N, N') = dim V_N / (M(N') cap V_N)` is exact (modulo the prime) and is an **upper bound** for the rank of `V_N`
in `S_F`. So the runs prove: the image of `V_60` in `S_F` is spanned by the 11 monomials of weight at most 4. They
do not prove `dim S_F >= 11`, nor that a nonzero normal form is a nonzero class: a nonzero normal form proves
`h not in M(N')`, and `h not in M` only if `M(N') cap V_N = M cap V_N` for `N` the weight of `h`. Given the cited
theorem (`dim S_F = 11`; the count was re-derived, including connectedness of `Y`, `delta = 2`, the 5 points,
`b_2(X) = 1` and the injectivity of `H^1(X) -> H^1(Y)`, which follows from `dw ^ dW != 0`: `a3`, items 5-7), this
condition is equivalent to `V_4 + M = B`, which the computation establishes up to weight 60 and no further. Evidence
that the truncation is saturated: the plateau 11 for all `4 <= N <= N' - 10` at `N' = 46, 58, 70`, and identical
normal forms of 35 test elements at the three values. The label FINITE-EXACT is right; "exact" (pre-registration)
is not. For the toy everything is proved by hand.

## 5. Proof-check notes

- **2.1, 2.2, 2.3.** Correct. (d): the generating-function proof is right, `e^(sV + tU - st)` has coefficients
  `H_(n,m)` with `n` the exponent of `U` and `m` that of `V` (Ito's polynomials), so `d_u^m d_v^n` goes with
  `H_(n,m)`; checked for `n, m <= 6` and on the slice. For an automorphism `M_F` is the kernel of
  `U^i V^j -> i! [i = j]`, and K-GMC is literally `GMC2Main.gmc2` (for `P, Q` in `C[X_0, X_1]`: `E(P^m) = 0` for
  `m >= 1` implies `E(Q P^m) = 0` for `m >= N`) after the substitution `X_0 = U`, `X_1 = V`.
- **3.1.** Correct given the citation. Checked: the Euler characteristic count; the intersection number of a curve
  of bidegree `(p, q)` with the (1, 1)-curve `H_t` is `p + q`, all intersections affine and transversal for generic
  `t` when the closure of `S` avoids the two base points; `THM-4620 (5)` does contain "`R = K[f]`, `L = K(f)`,
  `d = 1`" (geometrically connected generic fibre); `d = 1` gives an automorphism; `>= 2d - 1`. See M6 (every centre
  allowed) and m3 (`c = 1`). Remark 3.3: the exact sequence
  `H^1(X) -> H^1(Y) -> H^2(X, Y) -> H^2(X) -> H^2(Y) = 0` gives the formula. The note uses the injectivity of
  `H^1(X) -> H^1(Y)` in both examples without an argument; it holds: `dx/x` has a residue at the puncture `x = 0`
  of `Y`, and `dq/q` restricts to a nonzero class because `q` is an invertible function that is not constant on the
  (connected) fibres of `W` (`dw ^ dW != 0`), so it has a zero or a pole at a puncture. Add one sentence.
- **4.1.** `Theta` is well defined (values in polynomials) precisely by finite chaos, and it is the Taylor expansion
  of the class of `h` in the centre (section 7). Line by line: `Theta(d_u h) = d/ds Theta(h)`; the Leibniz step
  `d_u^m d_v^n (V h) = V d_u^m d_v^n h + n d_u^m d_v^(n-1) h` (checked on three maps) gives
  `Theta(V h) = (t + d/ds) Theta(h)`, `Theta(U h) = (s + d/dt) Theta(h)`; the six commutators of
  `s + d/dt, t + d/ds, d/ds, d/dt` are those of `U, V, d_u, d_v`; `P` is free of rank one over `C[U, V]`; the
  Submodules lemma is correct as stated (it is Kashiwara's lemma for a module on which the two derivations are
  locally nilpotent); surjectivity follows from the constant terms; `B (x)_A Frac(A) = Frac(B)` has dimension `d`
  because `B` is a domain algebraic over `A`. The conclusion needs a lower bound from the comparison theorem; as
  written it uses `2d - 1 > d` correctly; with 4.1' only `S_F != 0`.
- **5.1.** Correct: `[P_c^m] = 0` by 2.2 (b); the sets `Z_N` are Zariski closed (zeros of a polynomial with
  coefficients in `S_F`), increasing, and cover `C`; a countable union of finite sets does not; both halves; `k = 0`
  gives SN.
- **6.1.** Correct (the derivative identity under the residue, the independence, the closed formula). See M5 and m5.
- **Section 7.** No proof depends on it. See M5 (iii).
- **Scripts.** No bug found in `stein.py`, `stein_runs.py`, `stein_more.py`, `slice_pred.py`, `c3check.py`.
  `slice_pred.py` derives the curve and the intersection number but does not test the fibre sizes.

## 6. Pre-registration, credit, literature

*Pre-registration.* `git show d701fad1cd:.../PREREGISTRATION.md` (88 lines) is byte-identical to the first 88 lines
of the current file; `90e0a6855d` appends the S6 addendum (13 lines); `86ab0707d6` appends the outcomes. Nothing
above the outcomes was changed. `stein_runs.py` first appears in `90e0a6855d`, `stein_more.py` (the S6 script) in
`86ab0707d6`, after the S6 registration. Both registration commits are ancestors of `origin/main`, and Session268's
commit `e1a5a4841a` (12:43:30) sits on top of `90e0a6855d`, so the registration was on `origin/main` before the
outcome commit (12:47); that it preceded the runs cannot be read from git. The working claims W1-W5 are all in the
final note; W4 gained the rung CV. Outcome labels are right, with m12.

*Repo.* Not found before this packet: the operators `d_u - V`, `d_v - U`; a functional on `C[x,y]` defined through
a Keller pair; "JC(2) iff ... unique"; a Betti count for `(f-a)(g-b)`; the slice `F_1 = c` with its formulas (the
quartic `L` itself is THM-1310). Present and to be credited: HYP-8350 (M2); Dimca-Saito in the THM-4620 audit (M4);
THM-4122 (M6); Y. Stein 1989 in THM-4626. The credits to THM-4620 / THM-4626 are accurate.

*Literature* (sub-agent sweep; I read Zhao 2010 and van den Essen's survey myself): Dimca-Saito 1993 and Sabbah 1996
as cited (m7); Zhao 2010, van den Essen-Wright-Zhao 2011, Derksen-van den Essen-Zhao 2017 (M2); Gurjar-Masuda-
Miyanishi 2013 (M4); Bavula 2024 (JC iff the twisted module is 1-generated over `P_n`; cyclic over `A_n` always; no
de Rham or exponential twist there); Jelonek 1993 (the non-properness set is empty or a hypersurface), Nguyen Van
Chau 2004 (polynomial components); Orevkov 1987, Domrina-Orevkov 1998, Domrina 2000 (`d >= 5`), Zoladek 2008
(`d >= 6`); Ito's Hermite polynomials and their generating function as in the note.

## 7. Addendum: Theorem 5.2 (draft section 5.2)

**Verdict: PASS WITH CORRECTIONS.** The theorem and its proof are correct.

(1) *The family.* `[D_a, d_v - f + a] = [d_u, -f] + [d/da, a] = 0`, `[D_a, d_u - g + b] = 0`, symmetrically for
`D_b`, `[D_a, D_b] = 0`, and `R` is a `C[a, b]`-submodule: `G` is a module over the Weyl algebra of the plane of
centres. Right exactness gives the fibre `C[x,y]/M_c`. A simpler description, which I recommend stating: the map
`sum k_ij a^i b^j -> sum (f - d_v)^i (g - d_u)^j k_ij` is well defined (the two operators commute), kills `R`, and
induces an isomorphism

    G = B  with  a = f - d_v,  b = g - d_u,  D_a = d_u,  D_b = d_v

(relations checked on three maps, `a1`, item E). So `G` is Bavula's twisted module pulled back along the linear
automorphism `(U, V, d_U, d_V) -> (U - d_V, V - d_U, d_U, d_V)` of the Weyl algebra.
(2) *The cited facts* are used correctly. `C[a,b,x,y] e^(-(f-a)(g-b))` is holonomic (twist of the polynomial ring by
an automorphism of the Weyl algebra); `G` is the top cohomology of its direct image because, with `{f, g} = 1`, `R`
is the image of the relative twisted differential (Proposition 2.1 with parameters), and `d/da + (g - b)` is the
induced action. With the isomorphism above, holonomicity of `G` also follows from Bavula's theorem (already CITED in
THM-4620 (6)) without Bernstein. A holonomic module is a bundle with integrable connection on a dense open `Omega`,
and its fibre at `c` in `Omega` is `S_c`.
(3) *Flatness.* `D_a (1 (x) h) = 1 (x) d_u h`, and the pairing with a flat section of the dual gives
`d/da pi_h = pi_(d_u h)`, `d/db pi_h = pi_(d_v h)`. On the toy this is visible in closed form:
`lambda_+-` are flat, `d/da lambda_+(x^k y^l) = ((k+l)/2) lambda_+(x^(k-2) y^l) = lambda_+(d_u x^k y^l)`, and
`lambda_+(x^(-2)) = 1/(2a)` is not a polynomial in `a`.
(4) *Polynomiality.* Correct: a terminating Taylor series at one point of the disc makes `a -> pi_h(a, b_0)` a
polynomial; an uncountable `E_N` has an accumulation point in `Delta_2`; the identity theorem in `b` for each fixed
`a`; the same in the other variable.
(5) *Finite chaos.* Correct. The last step is simpler with 4.1' (m1): any centre of `Delta` will do, and only
`S_c != 0` is used (and by M6 every centre is allowed anyway).
(6) *Corollary 5.3.* All implications hold. `1 => 4`: for an automorphism `(f, g)` and any `(a, b)`, `(f - a, g - b)`
is an automorphism with the same dual fields and `M_c` is its Stein kernel at the centre 0; GMC(2) as proved in
canon applies after `X_0 = f - a`, `X_1 = g - b`. The rungs: the stated equivalence is correct, and it says that
`K-GMC_D` is "no Keller pair of generic degree `2 <= d <= D`" in other words; `K-GMC_4` is exactly the cited
theorems and carries no further information. Say so, and do not call the rungs intermediate statements.
(7) *The example.* `d_v = x d/dy` is locally nilpotent on `C[x, 1/x, y]`, so SN holds in that direction at every
centre. "Violates it in the direction `d_u`": for centres with `a != 0` (Proposition 6.1); at `a = 0`, `S_c = 0` and
SN holds trivially in both directions. Add the condition.

**The theorem proves more than it states.** (i) *Very general centres.* If `(f, g)` is not an automorphism, there is
one polynomial `h_0` such that the closed sets `A_N = {c in Omega : [d_u^m h_0]_c = [d_v^m h_0]_c = 0 for m >= N}`
are all proper (otherwise SN would hold on `Omega`), and SN can hold only at centres of
`(C^2 minus Omega) union A_1 union A_2 union ...`, a countable union of proper algebraic subsets. So **SN (a
fortiori CV or K-GMC) at one centre outside a countable union of curves, depending on the pair, forces an
automorphism**: HYP-9668 (D) is proved for a very general centre, and items 4-6 of Corollary 5.3 can be stated "at
one very general centre". What stays open is detection at a prescribed centre (for instance at every centre).
Observation m16 shows that special centres do behave differently on the slice. (ii) *Surfaces.* The proof works for
an etale map from a connected smooth affine surface (`B = F_+ O_X` is holonomic); the conclusion is "`S_c = 0` for
generic `c`, or `F` is an isomorphism". On the slice (`dim S_c = 11` generically, CITED) this gives the proved
replacement for M3: SN fails at every centre outside a countable union of curves. (iii) *One direction (sketch, NOT verified).* If SN
holds in the direction `d_u` only, on an open set, step (3) gives `d_u^(N+1) h in M_c` for all `c` in `Delta`, i.e.
`d_u^(N+1) h` lies in the kernel `K'` of `h -> sigma_h` at the generic point (the `C[a, b]`-torsion of `G = B`). If
`K' = 0`, `d_u` is locally nilpotent and the pair is an automorphism (THM-4620 (5)(e)). `K' = 0` holds on the toy
(the map `h -> lambda_+(h)` is injective) and for open immersions (the socle of `B` is `C[x,y]`). In general I
expect it from the period over the Lefschetz thimble through a preimage of the centre (`L(h U^i V^j) = 0` for all
`i, j` forces all Laurent coefficients of `h` on the vanishing circles to vanish), but I did not write a proof. If
it holds, "one direction over an open set" is settled as well, and the example of (7) is exactly the case where
`d_v` is locally nilpotent.

**Required changes of wording (the coordinator's question: yes).**

- README section 0 item 4 and section 5, "The three weaker rungs are new conjectures that lie between the theorem
  GMC(2) and the open problem JC(2); whether any of them implies JC(2) ... is open" -> "For a fixed centre, K-GMC, CV
  and SN are open for every Keller pair that is not an automorphism (none is known). With the centre free they are
  not intermediate: by Theorem 5.2 each of them, required at every centre of a nonempty open set (or at one very
  general centre), is equivalent to JC(2). What remains open besides JC(2): detection at a prescribed centre, and
  detection from one direction."
- README section 5, last paragraph ("Detection asks whether SN (the weakest) already forces an automorphism ...") ->
  restate with the centre quantifier; section 9 ("it does not prove detection") -> "detection is proved for very
  general centres (Theorem 5.2), not for a prescribed centre nor for one direction".
- HYP-9668: title and status. (A)-(C) must say "at a prescribed centre"; "None of (A)-(C) is known to imply JC(2)"
  -> "at every centre of an open set each of (A)-(C) is equivalent to JC(2) (THM-5138, Theorem 5.2); at a prescribed
  centre none is known to imply it"; (D) -> "PROVED at very general centres; OPEN at a prescribed centre and in one
  direction"; (E) -> M1.
- Corollary 5.3, "JC(2) <=> the Gaussian moments property holds in the coordinates of every Keller pair, at every
  centre": add that this is the Image Conjecture of Zhao for the potentials `(f - a)(g - b)` (M2), a planar
  counterpart of Zhao's Theorems 3.6-3.7 (which need all dimensions), and that the conjecture fails for general
  potentials (M1).
- THM-5138 status for the new part: "PROVED from CITED facts on holonomic modules (or from Bavula's theorem through
  `G = B`) and Theorem 4.1".

## 8. Not verified

- The proofs of Dimca-Saito, Sabbah, Hien-Roucairol, Nguyen Van Chau, Bavula, Zoladek (statements only; Dimca-Saito
  and Sabbah through a sub-agent's reading of scans).
- That the truncation is saturated in all weights (section 4); equality of the ranks over `Q` with the ranks modulo
  the two primes.
- That the pre-registration was pushed before the runs (commit order and ancestry only).
- The sketch (iii) of section 7, and whether the counterexample of M1 to the general Image Conjecture is recorded in
  the literature.
- THM-2022 and its Lean build (taken as canon, as instructed); only the statement was compared.

## 9. Files

`tw.py` (engine); `a1_hand_identities.py` (sections 2, 4, family module); `a2_toy_residues.py` (Proposition 6.1,
local Gaussians, hand bound); `a3_slice_derivation.py` (THM-1300, slice, fibre sizes, intersections,
irreducibility); `a4_controls.py`; `a5_toy.py`; `a6_slice.py` (`dims 40 46`, `classes 46 5`);
`a7_slice_stability.py`; `a8_slice_hermite.py`; `a9_rank_one.py` (M1); `a10_slice_centres_sn.py`;
`a11_slice_centre00.py`; the `.out` files next to them; `repro/` (the author's scripts re-run). Reproduce with
`python <script>` in this folder; the longest run is about 30 seconds.
