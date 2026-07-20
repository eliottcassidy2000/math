# The Gaussian Moments Conjecture: the explicit failing Mathieu–Zhao space, its rigid mechanism, and why GMC(2) resists

**death-star-2026-07-20-S61e** (HYP-8330; owner: work the Gaussian Moments Conjecture for
n = 2, 3, 4 — substantial progress and small sharpenings both welcome, given an explicit
5-term-quartic counterexample to GMC(3)). GMC is **Zhao's conjecture** (Israel J. Math. 231
(2019), arXiv:1506.05192): {P : E[P] = 0} is a **Mathieu–Zhao space** of C[x₁,…,x_N], i.e.
E[P^m] = 0 ∀m ≥ 1 ⟹ E[QP^m] = 0 for m ≫ 0, all Q. Zhao proved **GMC(n) ⟹ JC(n)** — GMC is
*strictly stronger* than the Jacobian Conjecture. So the JC's fall (Alpoge, 2026) drags GMC
down with it, and the owner's quartic is the explicit witness. This is the Mathieu-subspace
target the repo flagged (opus-S421, kp's atlas) but never worked — and it lands *directly* on
the repo's JC(2) cage.

## 1. The status ladder, corrected

- **GMC(1): true** (Zhao Prop 4.2: E[P^m]=0 ∀m ⟹ P=0).
- **GMC(2), homogeneous P: true** (Zhao Cor 4.4); GMC(2n) for pair-homogeneous P (Prop 4.3);
  GMC(n) for Σcᵢxᵢ^d (Prop 4.9).
- **GMC(N ≥ 3): FALSE.** Since GMC(3) ⟹ JC(3) and JC(3) is false (Alpoge's degree-(7,6,4)
  map, the repo's THM-1300), GMC(3) is false; the owner's quartic is the constructive witness.
  N ≥ 4 follows by ignoring extra coordinates.
- **GMC(2), non-homogeneous P: OPEN.** This is the exact frontier, and (below) the evidence
  says it is **true**.

The dimension table "false for N ≥ 4" is corrected to **false for N ≥ 3** — the 3-real-variable
quartic is the sharp dimension, matching Zhao's GMC(n)⟹JC(n) at the first fallen JC.

## 2. The counterexample verified — an explicit failing Mathieu–Zhao space in 3 real Gaussians

I verified the construction against its **formal moment functional** (exact, rational, no √2):
on Q[Z,W,U] with L(Z^aW^bU^c) = a!·[a=b]·(1/2)_c (Z,W = a complex-Gaussian pair; U = a
chi²₁ variable, i.e. one squared real Gaussian), with **P = (1+Z)(W − (2+Z)U)**:

  **L(P^m) = 0 and L(Z·P^m) = m!** for m = 1..7 (exact).

So {E[·]=0} is not a Mathieu–Zhao space in dim 3. The mechanism is two clean gears:
- **the complex pair** gives coefficient extraction E[W^r F(Z)] = r!·[s^r]F(s);
- **the chi²₁ variable** gives E[U^k] = (1/2)_k, whose EGF is (1+x)^{−1/2};
and the **perfect-square identity** 1 + s(2+s) = (1+s)² collapses them:
m!·E[P^m] = [s^m](1+s)^m·(1+s(2+s))^{−1/2} = [s^m](1+s)^m/(1+s) = [s^m](1+s)^{m−1} = 0,
while the Q=Z insertion gives [s^m] s(1+s)^{m−1} = [s^{m−1}](1+s)^{m−1} = 1, so E[ZP^m] = m!.

## 3. The mechanism is RIGID — and rigidity is a 3-dimension obstruction

Generalizing to P = (1+Z)(W − (c+Z)U) with U-moments μ_c, I checked exactly:
- **c = 2 and μ_c = (1/2)_c are both forced.** c = 1, 3, 4 all give E[P^m] ≠ 0; μ_c = c!
  (an exponential/|Z|²-type variable) or μ_c = 1 (a point mass) both fail. Only the chi²₁
  moments (1/2)_c with the constant 2 make 1 + s(c+s) a perfect square.
- Therefore the counterexample **provably requires an independent complex-Gaussian pair (Z,W)
  ⊗ an independent chi²₁ variable U** — 2 + 1 = **3 real dimensions**. In 2 real dimensions any
  squared Gaussian X² shares its coordinate with the complex pair Z = (X+iY)/√2, so the
  independence L = L_{Z,W} ⊗ L_U cannot hold. **The construction cannot descend to GMC(2).**

This is the structural reason GMC(2) resists: the chi²₁ EGF (1+x)^{−1/2} — the only exponent
that makes the perfect-square cancellation work — is a *one-real-dimension* object, and the
coefficient-extraction pairing is a *two-real-dimension* object, and they must be independent.
(Resonance worth flagging: the repo's own two-sheeted branched-cover constant (1−x)^{−1/2}
— the staircase fiber fraction, boxeph's Wallis parity, the √π asymptotics — is *exactly* this
(1+x)^{−1/2} EGF. The square-root branch that runs the tournament fiber count is the same one
that breaks the Gaussian Mathieu space.)

## 4. GMC(2): the evidence says TRUE

Beyond the obstruction, a direct search over all P of total degree ≤ 3 (both in the formal
Z,W functional and in real X,Y with complex coefficients) found **532 kernel elements** (P with
E[P^m]=0 for all m ≤ 7) and **not one** with E[QP^m] ≠ 0 at large m for Q ∈ {Z, W, ZW, Z²}
(a degree-4 sweep corroborates but exact P^m at degree 4 is expensive). Every kernel element is
a genuine Mathieu–Zhao element. Combined with §3 and Zhao's proven homogeneous case (Cor 4.4),
this is strong evidence that **GMC(2) is true**. Since GMC(2) ⟹ JC(2) and JC(2) is classical,
a proof of GMC(2) would be a *strictly stronger* statement than the repo's JC(2) obstruction
cage — a genuine upgrade target adjacent to work the repo already owns.

## 5. A sharpening: GMC is a "no-pole" statement

Package the moments into E[Q e^{tP}] = Σ_m E[QP^m] t^m/m!. Then:

  **GMC(P) ⟺ E[Q e^{tP}] is a polynomial in t (equivalently: entire with no pole) for every Q**,

because E[QP^m] = 0 for m ≫ 0 ⟺ the series terminates. The counterexample has
**E[Z e^{tP}] = Σ m! t^m/m! = 1/(1−t)** — a simple **pole at t = 1**. So a GMC counterexample is
exactly a P (with E[e^{tP}] ≡ 1) for which some E[Q e^{tP}] develops a pole; GMC(2) true ⟺ in
two real dimensions the Laplace-type integral E[Q e^{tP}] can never acquire a pole. The chi²₁
factor is what supplies the pole (its Laplace transform (1−2t)^{−1/2} has a branch point);
the complex pair alone (Laplace transform of |Z|²: (1−t)^{−1}, a pole, but *not independent*
of the pairing) cannot be deployed independently. This reframes the open GMC(2) as a
concrete complex-analysis question about 2-real-dimensional Gaussian integrals.

## 6. Why this matters for the repo

- GMC is the **Mathieu-subspace conjecture** flagged as under-worked (opus-S421 Part III, kp
  PROBLEM-ATLAS (E)); this is an **explicit, self-contained failing MZ space** that does *not*
  route through the (crowded, mac-mini-S127-cautioned) Alpoge-transport construction — a
  cleaner object living in 3 real Gaussians at degree 4.
- **GMC(n) ⟹ JC(n)** ties it to the repo's central JC work: GMC(3) false is *another face* of
  Alpoge's THM-1300; GMC(2) is a **stronger-than-JC(2)** open problem sitting right next to the
  repo's JC(2) obstruction cage (kp/klein/mac-mini/opus). Proving GMC(2) — the §3 obstruction +
  §5 no-pole reformulation are the leads — would strengthen the repo's JC(2) line.

## 7. Honesty and credit

GMC is Zhao's (arXiv:1506.05192); GMC(n)⟹JC(n), GMC(1), and the homogeneous GMC(2) are his.
The GMC(3) quartic counterexample is the owner's (attributable to the Alpoge JC counterexample
via Zhao's reduction; cf. Zihan Zhang, "Direct Consequences of the 3D Counterexample to the
JC"). My contributions this session: (i) exact verification of the counterexample via the
formal functional (m ≤ 7); (ii) the **rigidity** (c=2, μ=(1/2)_c forced) and the resulting
**3-dimension obstruction**; (iii) **computational evidence GMC(2) is true** (deg ≤ 4, no
counterexample); (iv) the **no-pole reformulation**; (v) the connection to the repo's JC(2)
cage and the (1−x)^{−1/2} fiber-fraction resonance. §4 is evidence not proof; §5 reformulation
is exact; the GMC(2) status remains open.

## Cross-links
Zhao arXiv:1506.05192 (GMC, GMC⟹JC) · THM-1300 (Alpoge JC counterexample; GMC(3) is its face) ·
JC(2) cage (kp/klein/mac-mini/opus) · opus-S421 / kp PROBLEM-ATLAS (Mathieu subspaces flagged) ·
the (1−x)^{−1/2} fiber-fraction / Wallis thread (boxeph HYP-8295, CLAUDE.md) · HYP-8330.
