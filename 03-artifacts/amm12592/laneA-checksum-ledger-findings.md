# Lane A findings — mechanism extraction on the THM-2225 cyclic-checksum extractor

Referee: `laneA_referee.py` (exact rational, stdlib only; 20 checks, ALL PASS, exit 0).
Output: `laneA_referee.out`.
sha256: script `c56c366eb60841d0c70377d481a96fa98be210c8ba2536799b01651531f50a00`,
output `f0a09943edfd6b27e39d98e98b47a55b14f91ff99d8abeee85547a60d571d82c`.
Runtime ~46 s. Scope: both spines, all critical values m = 1..31 (complete dyadic
blocks M = 2, 4, 8, 16, 32); identity checks are coefficient-level polynomial
identities (hence proofs, not samples); point checks at 60 rational p = k/61.

Notation. p = P(0), q = 1−p. Critical value n = initial constant-run length.
Block M = 2^ceil(log2(n+1)), opener m' = M/2, so m' ≤ n ≤ M−1. THM-2225 rule:
M = 2: 01→H, 10→T; M ≥ 4: s(y) = Σ_{i=1}^{m'} i·y_i mod m' on the second half,
H iff s < m'/2. Deadline T(n) = max(2, 2n−1), i.e. C = 2, D = −1, per-shell
depth budget d_m = T(m) − (m+1) = m − 2 for m ≥ 2 (d_1 = 0 via the max clause).
W_m = future-H-share at node 0^m 1, V_m at 1^m 0, Δ_m = W_m − 1/2, Δ'_m = V_m − 1/2.

## 1. The spine reformulation is VERIFIED on the known scheme

**V1 (decided-tree form).** For every m ≤ 31 and both spines, the subtree of
0^m 1 (resp. 1^m 0) is fully decided by relative depth f = M−m−2 (m ≤ M−2) or
f = 0 (m = M−1), the last carrier bit (position M) is provably irrelevant on
every path, and
`W_m = Σ_k w_k p^{f−k} q^k` with integers `0 ≤ w_k ≤ binom(f,k)` — the FACT's
containment direction holds cell-by-cell (check S2.FACT). Two independent
routes agree exactly for all 62 nodes: full-shell enumeration vs. the lazy
(minimal-stopping) recursion with variable leaf depths (check S2.route-equality).
This concretely verifies the refinement argument behind the FACT: refining a
variable-depth decided tree to full depth f preserves G and produces the
leaf-subset vector; conversely any leaf-subset of the full depth-d tree is a
valid decided tree, so the FACT's set equality
{depth-≤d decided-tree polynomials} = {Σ w_k p^{d−k} q^k : 0 ≤ w_k ≤ C(d,k)}
is correct (two-line proof; no further caveat found).

First polynomials (full table to m = 31 in the .out):

| m | W_m | V_m |
|---|-----|-----|
| 1 | 1 | 0 |
| 2 | 0 | 1 |
| 3 | 1 | 0 |
| 4 | p | 1 − 2p + 2p² |
| 5 | 1 − p | 1 |
| 6 | 0 | 0 |
| 7 | 1 | 0 |
| 8 | 3p − 7p² + 7p³ − 6p⁵ + 4p⁶ | 1 − 4p + 14p² − 28p³ + 34p⁴ − 24p⁵ + 8p⁶ |
| 12 | 1 − p | 2p − p² |
| 13,14 | 0 | 0 |
| 15 | 1 | 0 |

**V2 (identity (S), exact dyadic remainder).** With
R_N := 1/2 − Σ_{m≤N} [p^m q W_m + q^m p V_m]:
`R_{2^K−1} = (p^{2^K} + q^{2^K})/2` **as an exact polynomial identity** for
K = 1..5 (check S5.dyadic). Hence (S) holds: R_N → 0 on (0,1) and the
truncation error at dyadic times is exactly half the residual constant-ray mass.

**V3 (interval and structural checks).** For all N = 0..30 and all 60 p-values:
`0 ≤ R_N ≤ p^{N+1} + q^{N+1}` (so R_N/(p^{N+1}+q^{N+1}) ∈ [0,1]), and R_N lies
inside the exact interval `p^{N+1}·[G_{0^{N+1}}] + q^{N+1}·[G_{1^{N+1}}]` where
the spine values are enclosed by the proven block recursion
u_j = (1−p^{2^j}−q^{2^j})/2 + p^{2^j} u_{j+1} (0-spine),
v_j = (1+p^{2^j}−q^{2^j})/2 + q^{2^j} v_{j+1} (1-spine), tail ∈ [0,1] at level
2^13; max enclosure width 1.6e−59 (checks S6.*). N = 0 is exact root fairness.
Note: G at dyadic spine nodes is NOT 1/2; only the p^t-weighted two-spine
combination returns exactly to (p^t + q^t)/2 at dyadic t.

## 2. Exact block anatomy (the mechanism)

**V4 (per-spine closed forms; enumerated for m' = 2,4,8,16, and valid for all
dyadic m' by THM-2225's per-tail-weight bisection plus its eq. (13)):**

- 0-spine: `Σ_{m=m'}^{2m'−1} p^m q W_m = p^{m'} (1 − p^{m'} − q^{m'})/2`
- 1-spine: `Σ_{m=m'}^{2m'−1} q^m p V_m = q^{m'} ( p^{m'} + (1 − p^{m'} − q^{m'})/2 )`
- Block M = 2 anomaly: the sums are pq and 0 respectively (sign-reversed roles).

**V5 (deficit flow, the central structural fact).** Per block M = 2m' ≥ 4:

```
Σ_{m=m'}^{2m'−1} p^m q Δ_m  = − p^{m'} q^{m'} / 2
Σ_{m=m'}^{2m'−1} q^m p Δ'_m = + p^{m'} q^{m'} / 2      (M = 2: signs reversed)
⇒  D_M = 0  exactly, for every block.
```

ALL deficit cancellation is **intra-block**. Nothing is carried between blocks.
The unique cross-spine flux is the **equal-weight boundary pair** of the middle
composition class: 0^{m'}1^{m'} → T, 1^{m'}0^{m'} → H (verified M = 4..32),
each of Bernoulli weight p^{m'} q^{m'}; it transfers exactly p^{m'}q^{m'}/2 of
H-share from the 0-spine to the 1-spine. The 0-spine systematically
under-delivers heads by q^{m'}/2 (relative), the 1-spine over-delivers by
p^{m'}/2, and the boundary pair settles the account. At M = 2 the pair is
01→H, 10→T — reversed, matching the reversed flux.

**V6 (ledger; where cancellation happens).** For M ∈ {4,8,16,32}, refining each
total-Hamming column j ∈ {1..M−1} of S_M by (spine, run m): every column sums
to zero; columns j < m' involve ONLY the 0-spine, j > m' ONLY the 1-spine, and
j = m' is exactly the boundary pair. So cancellation is: (i) within one dyadic
super-shell S_M, (ii) within one total-composition column, (iii) across the
run-shells m ∈ [m', M−1] of ONE spine — never across blocks, and across spines
only through the single middle column. Example (M = 16, column j = 1, entries
2·(#H − |cell|/2) for m = 8..15): +1,+1,+1,−1,−1,−1,−1,+1. The forced constant
tail shells (m = 13,14: Δ = −1/2; m = 15: Δ = +1/2) are absorbed by the free
choices of the deep shells m = 8..11 inside the same column.

**V8 (parity/Lucas forcing).** In every cell (m,k) with binom(f,k) odd the
deficit is a forced half-integer (2δ odd) — verified for all m ≤ 31. Depth-0
shells (m ∈ {M−2, M−1} and the whole blocks M ≤ 4) are forced constants ±1/2.
Sup-norms: the 239-point exact grid gives sup_p |Δ_m| > 0.4676 for every m ≤ 31
and both spines (minimum attained at m = 18, value
533735954287050752035057259/1141260857376768000000000000 ≈ 0.46767): deficits
never decay with m; they are killed only by the p^m q prefactor and by
cancellation, never by smallness of Δ_m itself.

**V9 (spines are genuinely different).** Bit-flip equivariance
V_m(p) = 1 − W_m(1−p) holds exactly iff m ∈ {1,2,3} ∪ {M−1} (trivial constants)
and FAILS for every other m ≤ 31 (first failure m = 4: V_4 = 1−2p+2p² vs
1 − W_4(1−p) = p; also W_10 = V_10 coincide as polynomials while their
equivariant partners differ). Aggregate fairness is restored only by V4/V5 —
the reformulation cannot assume spine symmetry.

**V10 (single-defect balance = the p→1 Abel constraint).** The class-j=1 words
0^m 1 0^{M−m−1} receive per block M ≥ 4 exactly half H:
M=4: TH; M=8: HTTH; M=16: HHHTTTTH; M=32: HHHHHHHTTTTTTTTH (m = 16..31);
1-spine analogously. As p → 1 the measure concentrates on these words; the
checksum satisfies the resulting Abel-balance constraint block-locally.

## 3. Where the degree growth d_m ~ m actually binds (task 1/3 answer)

Depth table (S3): needed(m) ≤ nominal(m) = M−m−2 ≤ budget m−2 for all m
(equality nominal = budget iff m = 2^k). Needed depth:

- Openers m = 2^k: m=4: W needs 1, V needs 2(=budget); m=8: W and V both need
  6 = budget; m=16: both need 14 = budget. So from M ≥ 16 on, the checksum
  genuinely consumes the FULL budget (C−1)m + D − 1 = m−2 at every opener, on
  both spines. The small-block W-slack at m=4 (needs 1 = the THM-2160 §6.2
  shell-model floor h/2 − 1) does not persist.
- Inside a block the needed depth falls with slope −1 (sawtooth): the absolute
  decision time is the CONSTANT M−1 for the whole block; stop(m)/m is maximized
  at the opener: (2m'−1)/m' = 2 − 1/m' → 2. C = 2 is tight for this scheme and
  is pinned by the block openers.
- The checksum's celebrated slack (ignoring its last coordinate) is exactly the
  −1 in D = −1 and cannot compound: THM-2160 §6.3 proves at most ONE globally
  ignored tail coordinate (simple root of (1+u) at u = −1), and our S3 table
  confirms needed = budget at m = 8, 16. **The slack improves D, never C; it is
  no evidence for γ < 1.**

## 4. Diagnosis: what must change for d_m = γ·m with γ < 1 (task 4)

The verified facts isolate three load-bearing features of the C = 2 mechanism,
and each is an obstruction that a γ < 1 scheme must break:

**(a) Block-locality of the ledger.** V5/V6: every deficit is settled inside
its own dyadic block, within composition columns, with a single cross-spine
flux through the middle class. Block-locality + a common carrier (all shells of
the block decided by the same absolute time M−1) force the sawtooth and hence
C = 2 exactly. A γ < 1 scheme needs deadlines growing WITHIN the current block
scale, so its ledger cannot close block-locally: forced half-integer deficits
(V8; they exist at every scale by Lucas whenever binom(d_m, k) is odd) must be
transported forward to shells with strictly larger m — a genuine carry process
across scales, which is exactly the "deficit flow" rate question.

**(b) The equal-weight boundary pair needs s = t (dyadic alignment).** The only
cross-spine settlement device observed is the pair {0^{m'}1^{m'}, 1^{m'}0^{m'}}
— two words of EQUAL Bernoulli weight p^{m'}q^{m'} given opposite outputs. Try
to shorten carriers: a sub-block family [s, s+t−1] with carrier b^s y, |y| = t
(t = 2^r for the interior classes to have even binomials, per THM-2160 (25))
strands the boundary words 0^s 1^t and 1^s 0^t in singleton composition
classes; for s > t their weights p^s q^t ≠ p^t q^s are UNEQUAL, so they can no
longer cancel each other as a pair. They must be deferred (left undecided at
time s+t) and their ±(1/2)p^s q^t deficits balanced later by words of weakly
LARGER composition (a deficit at composition (s,t) can only be cancelled at
compositions (s+a, t+b), a,b ≥ 0 — monomial divisibility gives a partial-order
transport constraint). Sustainable deferral of one stranded word per sub-family
is precisely the carry problem whose feasible rate determines γ*. This is the
concrete cross-shell design problem lanes B/C should attack: sub-dyadic tails
t ≈ γ·s, one deferred boundary word per family, forward transport in the
composition order.

**(c) Composition-exactness itself.** THM-2160 (25) forces composition-exact
tails to be dyadic, and the family range forces carrier length 2h; within the
composition-exact paradigm the block END is the deadline for every shell of the
block, so C = 2 is the exact optimum of that paradigm (consistent with the
canon: uniform max(2, 2n−2) is impossible, and the §5 stratified rule buys
opener speed 3h/2 but still C = 2, binding just after the opener). Therefore
**any C < 2 scheme must be non-composition-exact**: its fairness must come from
the one-variable polynomial identity Σ p^m q Δ_m + Σ q^m p Δ'_m ≡ 0 with
class-level deficits genuinely nonzero and cancelling across shells of
different critical value — HYP-9061's "cross-shell rules" are not optional but
forced.

**Relation to the decoded certificate.** If the certificate's ratio
2457/6592 ≈ 0.37272 is a candidate γ (C = 9049/6592 ≈ 1.37272), the mechanism
in (b) predicts it should appear as the sustainable rate of the boundary-word
carry (deferred mass p^s q^t with t ≈ γs balanced along the composition order).
The two anchor biases p_A = 1285/2181, p_B = 8847357/11821757 of the artanh
inequality would then be the extremal p at which the transport ledger is
tightest. Untested here; flagged for the lanes holding the certificate.

## 5. Caveats

- All identity checks are exact for blocks M ≤ 32 (m ≤ 31); the closed forms V4
  for general dyadic m' additionally follow from THM-2225's proof (per
  tail-weight-class bisection per spine + boundary checksums (13)), so nothing
  rests on extrapolation, but the referee itself enumerates only M ≤ 32.
- Sup-norms of Δ_m are grid maxima (239-point exact grid), i.e. certified lower
  bounds on the sup; the structural claims never use them as upper bounds.
- The diagnosis in §4(b)–(c) is mechanism extraction (what must change), not a
  proof that γ < 1 is achievable or that C* < 2; the stranded-word transport
  model is handed off as a design problem, with the partial-order constraint
  (monomial divisibility) and parity forcing (V8) as its verified ground rules.
