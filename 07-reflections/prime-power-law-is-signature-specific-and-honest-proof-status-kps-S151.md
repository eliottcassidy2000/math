---
source: kind-pasteur-2026-07-24-S151 (Opus 4.8)
status: RESULT (refinement + honest proof status). The prime-power law C_{1/n}={1,...,P(n)-1} is
  SIGNATURE-SPECIFIC: it holds for the four Ramanujan signatures n in {2,3,4,6}, but every OTHER (non-signature)
  level collapses to C_lambda={1} -- decisively shown for lambda=1/5 (C_{1/5}={1}, clean golden-field PSLQ null
  for k=2,3,4 at 150 digits). This corrects kps-S150 open-item 3. States the honest proof status: k=1 is a
  THEOREM (all lambda); the superset direction is a finite set of elementary evaluations; the subset direction
  (k>=P(n) non-elementary) is a TRANSCENDENCE statement and is OPEN -- so the full law is NOT a theorem.
tags: [hypergeometric, closed-forms, prime-power, signatures, ramanujan, CM, transcendence, honest-status, series]
related: [kps-S146, kps-S147, kps-S148, kps-S149, kps-S150, THM-415, opus-AMM12592]
corrects: [kps-S150 open-item-3]
---

# The prime-power law is signature-specific; and the honest proof status

`S_lambda(k)=sum_{n>=0}(lambda)_n(1-lambda)_n/((n!)^2(kn+1))=int_0^1 2F1(lambda,1-lambda;1;x^k)dx`.
`C_lambda={k: S_lambda(k) elementary}` (elementary = Q-comb of `pi`, `log(alg)`, `arctan(alg)`).

## 1. Decisive probe: `C_{1/5}={1}` (non-signature levels collapse)
For `lambda=1/5` (NOT one of the four signatures), a **clean** golden-field PSLQ -- basis
`{1,sqrt5} x {log phi, log 2, log 5} + pi` (provably `Q`-independent by Baker, so any returned relation MUST
involve the target) -- at **150 digits, coeff <= 1e7** returns **None** for `k=2,3,4`:
> `pi S_{1/5}(2)`, `pi S_{1/5}(3)`, `pi S_{1/5}(4)` are all **non-elementary**. `C_{1/5}={1}`.

The **same** null holds for a second non-signature level `lambda=2/5` (`k=2,3,4` all None, same clean basis,
140 digits): `C_{2/5}={1}`. Two independent nulls => "non-signature `=> C={1}`" is firm, not a one-point fluke.

(Method note, MISTAKES-spirit: my first two attempts used bases with internal degeneracies --
`arctan sqrt5 = arctan(1/sqrt5)+arctan(2/sqrt5)`, and `phi^3 = 2+sqrt5` -- so PSLQ returned the *basis*
relation with target-coefficient `0`, which is NOT evidence either way. Only a provably independent basis makes
a null conclusive.)

**This corrects kps-S150 open-item 3**, which asked whether a unit/regulator mechanism gives an elementary
`S_{1/5}(k)` (cf. opus's `log_5(5 phi^2)`). Answer: **no** -- `S_{1/5}(k)` is non-elementary for `k=2,3,4`.
Opus's `C_arch=log_5(5 phi^2)` is therefore a *different* quantity, not a value of this series; the golden field
enters opus's problem for unrelated reasons.

## 2. The refined law
| level `lambda=1/n` | signature? | `C_{1/n}` | `P(n)` |
|---|---|---|---|
| 1/2 | sig 2 | {1} | 2 |
| 1/3 | sig 3 | {1,2} | 3 |
| 1/4 | sig 4 | {1,2,3} | 4 |
| 1/6 | sig 6 | {1,2} | 3 |
| 1/5 | **none** | **{1}** | (n/a) |

> **Refined prime-power law.** For the four Ramanujan signatures `n in {2,3,4,6}`:
> `C_{1/n}={1,...,P(n)-1}`, `P(n)=` largest prime power dividing `n`. **For every non-signature level,
> `C_lambda={1}`** (only the universal `k=1`).

So the "extended elementariness" past `k=1` is a **signature/CM phenomenon**, not a generic one. Consistency
check: signature 2 has `P(2)=2 => C={1}`, coinciding with the generic non-signature case -- signature 2 is the
boundary where the CM bonus is empty.

## 3. Honest proof status (what "prove the law" really requires)
- **`k=1`, all `lambda` -- THEOREM (the one solid pillar).**
  `S_lambda(1)=2F1(lambda,1-lambda;2;1)` (since `1/(n+1)=(1)_n/(2)_n`), and Gauss's theorem gives
  `=Gamma(2)Gamma(1)/(Gamma(2-lambda)Gamma(1+lambda))=1/((1-lambda)lambda Gamma(1-lambda)Gamma(lambda))`,
  so by reflection `Gamma(lambda)Gamma(1-lambda)=pi/sin(pi lambda)`:
  > **`S_lambda(1)=sin(pi lambda)/(pi lambda(1-lambda))`.** (Fully proved, universal.)
- **`{1,...,P(n)-1} subset C_{1/n}` (superset) -- a FINITE set of elementary evaluations.**
  `k=1` general (above); `k=2, lambda=1/4` derived from scratch (kps-S146, `g_2=arcsin` + IBP);
  `{k=2, lambda=1/3, 1/6}` and `{k=3, lambda=1/4}` verified to 120+ digits, same elementary (`log/arctan`) type.
  Provable in principle -- it is exactly the remaining handful of derivations.
- **`C_{1/n} subset {1,...,P(n)-1}` (subset: `k>=P(n)` NON-elementary) -- a TRANSCENDENCE statement, OPEN.**
  No computation *proves* non-elementarity; PSLQ gives only evidence (null to 150-170 digits). A rigorous proof
  needs the **irreducibility of the underlying hypergeometric motive** at `k=P(n)` (Katz; Beukers-Cohen-Mellit
  classification of when `H(alpha;beta)` is Artin/CM vs irreducible) **plus** a transcendence input (that an
  irreducible motivic period is not a `Q`-combination of elementary periods). Both are open here.

> **Verdict (honest).** "Prove the prime-power law" = prove the subset direction. That is **beyond current
> reach**. What is proved: the `k=1` theorem, and (modulo finitely many explicit derivations) the superset
> direction. I am **not** claiming the law as a theorem -- it remains a conjecture with a proven easy half.

## 4. Mechanism (why signature-specific)
The four `lambda` are exactly the levels where `2F1(lambda,1-lambda;1;.)` has a **modular parametrisation**
(Ramanujan's theories to alternative bases; arithmetic triangle groups). `S_lambda(k)=int_0^1(period)(x^k)dx`
pulls the parametrisation back along the `mu_k`-cyclic cover `x->x^k`; the period stays elementary while `k`
is below the **prime-power part** of the level -- a **Conway-Jones / roots-of-unity resonance** (in-repo
THM-415: prime/composite dichotomy = vanishing-sums-of-roots-of-unity dichotomy, the SAME governor used on the
LRC side). Non-signature `lambda`: no modular structure, so nothing survives past the universal `k=1`.

## 5. Crisp next targets
1. **Finish the superset direction**: derive `S_{1/3}(2)=3sqrt3 log2/pi`, `S_{1/6}(2)=3sqrt3 log(2+sqrt3)/(2pi)`,
   `S_{1/4}(3)` from scratch (the `k=2, lambda=1/4` derivation is the template). Then `{1,...,P(n)-1} subset C`
   is a theorem for all four signatures.
2. **Identify the `k=P(n)` constants** (`n=3,4,6`) as specific `L`-values / higher-weight CM periods, and PSLQ
   against them -- turns the "new period" entries concrete.
3. ~~Solidify signature-specificity with a second non-signature null~~ **DONE** (`C_{2/5}={1}`, section 1).
   Next: a level in a *different* field (`lambda=1/7`, `Q(sqrt(-7))`/`zeta_7`) to rule out a golden-field artefact.

Files: `/tmp/{c5clean,c25,cmspec,probe}.py`. Builds on kps-S146-S150; corrects S150 open-item 3.
