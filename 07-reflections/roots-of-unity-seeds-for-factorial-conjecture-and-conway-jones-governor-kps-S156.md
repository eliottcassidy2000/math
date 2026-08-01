---
source: kind-pasteur-2026-07-24-S156 (Opus 4.8)
status: RESULT (clean theorem + reformulation + capacity control) on the Factorial-Conjecture radial search.
  PROVES the roots-of-unity seed L((sum zeta^j X_j)^k)=k!*[m|k] -- so the FC low-moment structure IS governed by
  vanishing sums of roots of unity (Conway-Jones), the SAME governor as the series thread. Reformulates an FC(m)
  counterexample as a polynomial whose pushforward of Exp^m has all holomorphic moments zero (the seed gives only
  m-fold symmetry; the leak is at m|k). Exact "capacity control" for FC(2): the degree-2 antisymmetric family
  kills L(f^2),L(f^4) but is FORCED to L(f^6)!=0 -- consistent with FC(2) true. FC(3) leak-closing left open.
  Engages (and corrects) three inspiration threads honestly.
tags: [factorial-conjecture, roots-of-unity, conway-jones, periods, pushforward, capacity, series, honest-status]
related: [kps-S153, kps-S154, kps-S155, THM-415, THM-2801]
---

# Roots-of-unity seeds for FC, and the Conway-Jones governor

Search target (kps-S155): a radial/diagonal `f in C[X_1..X_m]`, `f!=0`, with all factorial moments
`L(f^k)=0` (`L(X^a)=prod a_i! = E[X^a]`, `X_i` iid `Exp(1)`). FC(m) says none exists.

## 1. The roots-of-unity seed -- THEOREM
> For `zeta=e^{2 pi i/m}` and `f = sum_{j=0}^{m-1} zeta^j X_{j+1}`:
> **`L(f^k) = k! * [m | k]`** (`= k!` if `m|k`, else `0`).

**Proof.** `L(f^k)=E[f^k]=k! * h_k(1,zeta,...,zeta^{m-1})` (multinomial expansion; each `E[X^a]=prod a_i!`
cancels the `1/prod a_i!`). The complete homogeneous symmetric polynomial at the `m`-th roots of unity has
generating function `sum_k h_k t^k = prod_j 1/(1-zeta^j t) = 1/(1-t^m)`, so `h_k(1,...,zeta^{m-1}) = [m|k]`. QED.
(Verified numerically `m=2,3,4`.) `m=2` is `f=x-y`: `L((x-y)^k)=k![2|k]` (`0,2,0,24,...`).

**Reading.** The seed kills every moment *except multiples of `m`*. That "vanish unless `m|k`" is a **vanishing
sum of roots of unity** -- **Conway-Jones / Lam-Leung** territory, the exact governor the series thread runs on
(kps-S149-S153) and the LRC thread cites (THM-415). So FC's low-moment structure and the series' closed-form
locus share one arithmetic engine. The pushforward measure `mu_f = f_*(Exp^m)` is **`zeta_m`-rotation invariant**
(that is what `E[W^k]=0` unless `m|k` means).

## 2. Reformulation (what a counterexample must be)
`L(f^k)=int w^k d mu_f`, so:
> **FC(m) counterexample `<=>` a polynomial `f` whose pushforward `mu_f=f_*(Exp^m)` has ALL holomorphic moments
> `int w^k d mu_f = 0` (`k>=1`)** -- i.e. `mu_f` is moment-balanced (a "fully rotationally symmetric" pushforward).
The seed reaches only `m`-fold (`zeta_m`) symmetry; **closing the `m|k` leak = upgrading `zeta_m`-symmetry to full
rotational symmetry through a polynomial map.** Whether a polynomial can do this is exactly the KZ/rigidity
question (kps-S154): a linear `f` gives `mu_f` with generating function `1/(1-t^m)` (never balanced), so a
counterexample must be **nonlinear**.

## 3. Capacity control (FC(2), exact)
Antisymmetric `f=(x-y)+a(x^2-y^2)+b(x^2 y-x y^2)` (odd moments auto-zero). Exact moment polynomials `L(f^2)`,
`L(f^4)`, `L(f^6)` in `(a,b)` computed. The variety `L(f^2)=L(f^4)=0` has **8 nonzero solutions**, and at every
one `|L(f^6)| ~ 10^4-10^5 != 0`.
> **Two parameters kill two even moments but are FORCED to leave the third nonzero -- the tower TERMINATES.**
This is a concrete instance of the friend's "capacity" picture (how many cancellations a plane curve `f=const`
supports), and it is consistent with **FC(2) true**. The `L(f^6)!=0` values are exactly the *reliable*
certificate direction: **certifying NON-vanishing**, never vanishing (kps-S151 discipline).

## 4. FC(3): open
The seed leaks at `3|k`; corrections must keep the cyclic `omega`-weight (so non-multiples of 3 stay auto-zero)
and kill `L(f^3),L(f^6),...`. Low-degree least-squares/exact searches were run but are computationally heavy and
**did not resolve it** -- honestly, no FC(3) counterexample found and none excluded here. The sharp target stands:
a **nonlinear, cyclic-weight** `f` in 3 vars whose `mu_f` is fully rotationally symmetric. FC(3) "may be false"
(friend) would be this object; FC(3) true would be a capacity obstruction one dimension up from Sec 3.

## 5. The three inspiration threads -- engaged, with honest corrections
1. **"S(k) `<->` FC; I proved `S(k>=4)` irreducible."** Correction (kps-S153): I did **not** prove that -- those
   non-elementarity claims were false-negative-prone and are *disciplined evidence, not proof*. The transferable,
   *reliable* half is real and used here: PSLQ/exact certificates of **NON-vanishing** (find one `m` with
   `L(f^m)!=0`), exactly the `L(f^6)!=0` control in Sec 3. The "edge-value / no-pole / monodromy" view is a shared
   lens, but the irreducibility/vanishing *direction* is the unreliable one -- keep it as evidence.
2. **AMM golden `<->` Paley `<->` periods; "edge vs exponential resummation."** The FC generating function
   `Phi_f(t)=int e^{tf-|x|}` is the exponential evaluation of the FC period family -- genuine framing. Caveat:
   `Phi_f(t)` **diverges for real `t>1`** (integrand `~e^{(t-1)x}`), so `"Phi_f == 1"` is *formal*; the moments
   `L(f^{mj})` (asymptotic series `sum ~ 1/(1-t^m)` for the seed) are the real objects. Use the edge/exponential
   duality as heuristic, not identity.
3. **The `n=2 | n>=3` wall is universal.** Real pattern: JC false `n>=3` (THM-1300), GMC/SIC false `n>=2/3`
   (THM-2801), FC(2) true-ish / FC(3) maybe-false (Sec 3-4). But rigor differs across members: JC/SIC walls are
   **proved**; the FC wall and the series' `"C_{1/4}={1,2,3}, irreducible at 4"` are **evidence/conjecture**
   (the latter retracted-to-evidence, kps-S153). The unifying "capacity holds until a threshold" shape is genuine
   (Sec 3 makes it explicit for FC(2)); do not upgrade the conjectural members to theorems.

## 6. Status
Delivered: the roots-of-unity seed theorem (FC `<->` Conway-Jones), the pushforward reformulation, and an exact
FC(2) capacity control. Open: FC(3) leak-closing (the nonlinear cyclic-weight rotational-symmetry search). This
is the concrete, honest instantiation of the kps-S154 JC<->period edge on the FC side.

Files: `/tmp/{tower,fast2,fc2solve}.py`. Builds on kps-S153/S154/S155.
