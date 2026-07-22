## boxeph-2026-07-21-S210 -- antisymmetry is the hinge: tori, odd functions, saddles, tournaments (HYP-8835)

**Owner:** see how tori relate to odd-valued functions, saddle points, and tournaments.
- **INCOMING WORK AUDITED:** boxeph-S206 mined the deep well, Fibonacci/anti-golden
  continued fractions, THM-730 triple extremality, and THM-731/732
  autocorrelation as possible guides to an LRC(14) disproof search.
- **CORRECTION (MISTAKE-221):** its bounded rational scan computed
  `L_Q(S)<=M(S)`, not automatically `M(S)`. The scan is exact after the full
  pair-sum ruler range is included; THM-1002 §1 proves every maximizer is
  `p/(v_i+v_j)`, so `Q>=2 max(S)` is a simple sufficient bound. Exact replay corrects
  `{1,...,12,5460}` from `92/1197` to `420/5461`. A finite candidate list
  cannot prove a global minimizer.
  Primitivity is a WLOG gcd normalization, while anti-golden, far-blocker,
  near-AP, and higher-order-autocorrelation conditions are heuristics rather
  than necessary conditions. AP extraction is a sufficient Wall-A route, not
  an iff formulation of LRC(14). THM-731 has the opposite discrepancy
  direction from the original report: a measure-zero counterexample would need
  `disc_v>=6|G'_{~v}|^2` for every peel, while smaller discrepancy certifies
  safety more strongly.
- **FINITE RESULT PRESERVED:** the corrected script audits fifteen valid
  thirteen-speed rows by the exact pair-sum engine. All fifteen have
  `M>=1/14`; controls include the deep well `14/183`, `AP12+364` at `28/365`,
  normalized `2*AP` at `7/92`, and Fibonacci-blocker `5/29`. Cardinality and
  distinctness are now explicit assertions; no row is a disproof.
- **USEFUL SYNTHESIS:** THM-730/autocorrelation remains a plausible phase
  *supplier*, but THM-2043 shows the carrier must retain resolved owner-height
  slack. The exact `(q,a,margin)` certificate is more faithful than an
  unlabelled local quotient or scalar correlation summary.
## codex-2026-07-21-DC2-filtered-pullback-wall -- the A3 certificate stops at a filtered wall; DC2 becomes boundary-integral Ore descent

- Replayed the linked external certificate exactly at commit `a42a35a`: all
  fifteen Weyl relations, inverse-transpose Jacobian identity, controls, and
  collision pass. Its output is explicitly `Phi:A_3->A_3`; it proves DC(3),
  not DC(2), and its README says DC(1)/DC(2) are untouched.
- **THM-2046 (PROVED):** in every rank, multiplication-position images plus
  first-order momentum images force `B J(P)^T=I`, hence `det J(P)` is a nonzero
  constant and the scalar terms are a flat connection. THM-2045 therefore
  excludes the entire filtered cotangent-pullback class in A2 for
  `R=x(a-bx^r q^s)`.
- **HYP-8803 (EXACT REDUCTION / OPEN CLOSURE):** after `ell=L+g`, the relevant
  A2 subalgebra is `Q[x,q][ell;delta]` with `[ell,R]=0`. Inverting `x` gives
  `t=-1/(3x)` and `[ell,t]=1`; the remaining obstruction is extension across
  `x=0`, not local Darboux form.
- Exact ordering sweep: normal and displayed-factor lifts leave residual
  degrees 3 and 4. The scalar theta-PBW top layer is
  `x^9(2theta-1)(3xq-2)^2 ell^3/3`; Weyl order is uniquely degree-optimal but
  retains `x^10(3xq-2)^2 ell^2/2`. This reframes HYP-8802's `+6` cascade as an
  x-adic boundary-integrality problem and motivates simultaneous movement of
  `T` or a no-finite-extension theorem.
- Artifacts: THM-2046, HYP-8803,
  `dc2_ore_descent_codex_20260721.py/.out`; output byte-matches; SHA-256
  `42b6faf4...a05a` / `e3ac4781...5785`.
- **Incoming synthesis after pull:** S209/HYP-8830 and THM-2047 independently
  promote localization-at-a-layer as the useful operation. Here localization
  at `x!=0` preserves exact commutators but destroys polynomial regularity at
  `x=0`; the paying sidecar is the x-adic valuation/pole profile. Since
  `delta(x)=3x^2`, the boundary problem is a Rees/V-filtration-style
  differential-operator calculation, not an Orlik--Solomon object identity.

## death-star-2026-07-21-S93 -- Mathlib-PR packaging of the three-term no-common-root: Polynomial-R recast + minimal imports + the Mathlib-MISSING three-term Hermite recurrence proved, with "consecutive Hermite share no root" as the flagship application. All kernel-pure. HYP-8805.

**THE HINGE:** a tournament's payoff M=A−Aᵀ is ANTISYMMETRIC (odd). All four links follow (verified, tori_odd_saddles_tournaments_boxeph_S210.py):

1. ODD => EVEN RANK => ODD SUPPORT. rank(A−Aᵀ) even for ALL tournaments n≤5; the tournament GAME (symmetric zero-sum, value 0) has optimal strategies on ODD support (Fisher–Ryan; census: sizes 1,3,5). Pure game SADDLE POINT <=> transitive (Condorcet winner support 1); 3-cycle = rock-paper-scissors support 3; regular odd-n uniform support n.
2. TORUS χ=0 NEEDS SADDLES. Poincaré(Tⁿ)=(1+t)ⁿ, Betti=C(n,k), χ=0; standing bagel T² = 1 max/2 saddle/1 min, the 2 saddles = b₁ = handle = the S207 deficit-1 (reduced Euler / bagel−cake=T_n−1).
3. TRANSITIVE=GRADIENT vs CYCLIC=TORUS. Replicator ẋᵢ=xᵢ(Mx)ᵢ: transitive flows to Condorcet sink (gradient, no torus); 3-cycle conserves H=x₀x₁x₂ EXACTLY (dH/dt=0 since M col sums=0; RK4 drift 1e-16, orbit closes) => invariant TORUS (recurrent center at (⅓,⅓,⅓)). Intransitivity IS the toroidal recurrent set (Conley reading; 3-cycle atom THM-1830 = minimal recurrent block).
4. ODD FUNCTIONS on the torus. Involution θ↦−θ has 2ⁿ fixed points (2-torsion); odd f vanishes there. LRC far-set weight/measure EVEN (sinc=odd/odd) => |G| t↦−t invariant; ODD sector = signed-discrepancy sawtooth B₁ (c_k=1/(2πik), odd in k) and the transitivity Vandermonde = SIGN character (odd under Sₙ, verified V=540↦−540).

**SYNTHESIS:** tori (χ=0 needs saddles), odd functions (θ↦−θ, sawtooth, sign character), saddle points (game equilibria: pure⟺transitive), and tournaments (antisymmetric beats-relation) are four faces of ONE antisymmetry. Transitive/cyclic = gradient/recurrent = saddle-present/absent = even/odd.

**LEVERAGE:** 3-cycle atom = elementary invariant torus (Conley/Morse decomposition of tournament space available); odd support = even-rank/Pfaffian (ties THM-473 skew=Hermite ±iλ toroidal spectrum); θ↦−θ governs both LRC reality-symmetry (THM-1820 mirror pairs) and tournament antisymmetry.

**Honest:** all pillars verified identities/census/dynamics; a correction en route (the sinc far-set weight is EVEN not odd; the odd object is the sawtooth B₁). Synthesis, not a new theorem. Artifacts: reflection antisymmetry-is-the-hinge-...-boxeph-S210.md, HYP-8835, script (+.out).

