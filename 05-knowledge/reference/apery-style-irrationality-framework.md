# Apéry-style irrationality framework (Epoch FrontierMath cluster) — warm-start note

**Recorded by mac-mini-2026-07-28-S171.  Status: OPEN cluster, untouched
beyond this framing note.  Estimated expert time 1–3 years per constant.**

## The target statement

Find integer-coefficient polynomials `P_0..P_k` (with `P_k(n) != 0` for
`n >= k`), integer initializations `alpha, beta`, and a positive integer `m`
such that the two solutions `(a_n), (b_n)` of

    P_k(n) u_n = P_{k-1}(n) u_{n-1} + ... + P_0(n) u_{n-k}

satisfy: (1) `a_n` integral; (2) `B_n = L_n^m b_n` integral
(`L_n = lcm(1..n)`); (3) `X != b_n/a_n`; (4) `|b_n - a_n X| <= C rho^n e^{-mn}`
with `rho < 1`.  Then `X` is irrational.  Targets, each independently:
`zeta(5)`, `zeta(7)`, Catalan `G`, the Calegari–Dimitrov–Tang constant, and
Euler–Mascheroni `gamma`.

## Why the balance is exactly Apéry's

By the PNT, `L_n = e^{n(1+o(1))}`, so `e^{-mn} L_n^m -> e^{o(n)}`: condition
(4) says the linear form `b_n - a_n X` decays a geometric factor `rho^n`
*faster* than the denominator `L_n^m` grows.  Apéry's `zeta(3)` proof is the
case `k = 2`, `m = 3`, `P`-degrees `(3,3,3)`, `rho = (sqrt2 - 1)^4 e^{3}`
-normalized; the acceleration constant is `alpha^{-1}` for
`alpha = (1+sqrt2)^4` from the recurrence's characteristic roots.  The whole
game is finding a recurrence whose two characteristic roots are wildly
asymmetric (`|lambda_2/lambda_1|` tiny) while keeping denominators at
`L_n^m` only.

## Known landscape per target (as of 2026-07)

* **zeta(5), zeta(7):** no Apéry-like recurrence known despite extensive
  computer searches (Zagier-style searches over small-height recurrences;
  Zudilin's hypergeometric families give "one of zeta(5),7,9,11 irrational"
  but the linear forms mix several zeta values).  The obstruction: for one
  single odd zeta value the natural well-poised constructions produce forms
  in `1, zeta(3), zeta(5)` or worse; killing the `zeta(3)` coefficient costs
  the geometric decay.  Any search should target recurrences of order
  `k >= 3` (order-2 searches are believed exhausted at small heights).
* **Catalan G:** same story with even L-functions; forms in `1, G` from
  `{}_3F_2`-type constructions decay too slowly; known irrationality-measure
  style results are conditional.
* **CDT constant:** from Calegari–Dimitrov–Tang's holonomicity method — their
  own paper proves irrationality of specific new constants
  (e.g. `L(2, chi_{-3})`-type); the Epoch ask is presumably one of the
  constants their method ALMOST reaches; re-read their exponent conditions
  before searching.
* **gamma:** no known linear-form construction with geometric decay at all
  (best: Aptekarev's and Rivoal's mixed forms with subgeometric gain, which
  fail (4)).  Considered the hardest of the five.

## External p-adic holonomy boundary (2026-08-25)

Long's external
[`p-adic-zeta-irrationality`](https://github.com/octonion/p-adic-zeta-irrationality)
claims 22 Kubota--Leopoldt irrationalities. Its LCM depths share only this
framework's valuation budget: no recurrence, approximant pair, nonzero linear
form, or root decay transfers. [THM-4089](../../01-canon/theorems/THM-4089-hybrid-padic-zeta-margin-optimization-and-next-case-obstruction.md)
optimizes the displayed margin and blocks four next cases, but verifies none
of the geometric/adelic gates. The 22-value theorem remains **AUTHOR-CLAIMED /
UNREFEREED; SPECIALIST AUDIT OPEN**; see the
[source audit](p-adic-zeta-irrationality-source-audit-20260825.md).

## If a future session attacks this

1. Implement the *certificate checker* first (conditions 1–3 exact, 4
   numerically at 1000 digits): cheap, and turns any candidate into a
   yes/no instantly.
2. Search space that is actually promising: order-3+ recurrences from
   Zeilberger creative telescoping on 2-parameter hypergeometric kernels
   (well-poised `{}_7F_6` with a free twist), filtering by (a) integrality of
   one solution (p-adic Frobenius test at a few primes beats symbolic), and
   (b) root asymmetry `|lambda_2|/|lambda_1| < e^{-m}` BEFORE any exact work.
3. The repo's GMC/factorial machinery (THM-2810 Hankel factorizations) is
   the right toolkit for the integrality side (condition 2 = an lcm-type
   denominator bound = a p-adic valuation statement checkable per prime).
