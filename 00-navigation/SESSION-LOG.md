> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## codex-2026-07-21 -- THM-2068/2073 owner-bank compression and dyadic safe-child descent

- **Minimum bounded owner bank:** THM-2068 turns the THM-2066 census into an
  exact set-cover problem. Inside clocks `15..34`, seven clocks
  `{25,26,27,28,32,33,34}` are necessary and sufficient for all `59,880`
  primitive divisor-complete eleven-cores through maximum `24`; all banks of
  at most six undominated clocks were exhausted and every chosen clock has a
  private core.
- **Uniform structural descent:** after pulling THM-2072's fixed-bank no-go,
  THM-2073 transfers THM-775's forgotten safe-child mechanism to the strict
  `1/14` seam. Every imprimitive deletion has gcd two, the first four lifts
  are partitioned `2+1+1`, and descent iterates through divisor-complete
  quotient cores (including the new denominator-`26` shell) to a hereditary
  terminal after at most eight levels. Exact referees pass normally and under
  `python -O`. LRC(14) remains open on the hereditary terminal lane.

## codex-2026-07-21 -- THM-2063/2064/2066, planar atlas audit, and dyadic owner words

- **Planar JC(2):** THM-2063 completely classifies Keller pairs having any
  output-pencil member affine along any linear source fiber. Coefficient
  descent gives two explicit triangular normal forms and inverses. Thus every
  hypothetical counterexample has fiber degree at least two in every source
  and target-pencil direction. MISTAKE-228 corrects THM-1330 from “exact
  classification” to a necessary-structure atlas. MISTAKE-229 preserves the
  binary symmetric-Hessian subcase but removes the unproved NC2/GMC-to-JC and
  “three equivalent descents” claims. JC(2) remains open.
- **LRC(14):** THM-2064 independently rederives the concurrent THM-2060
  multi-tail sheet-capacity law and the dyadic seam; it is now explicitly
  labelled a reformulation. Pulling the proved THM-2061 seam geometry led to
  THM-2066: odd tails carry binary lift-owner words, and a counterexample
  requires complementary words on every safe clock. Exact clock-bank census
  closes all `59,880` primitive divisor-complete quotient cores through
  `max(C)=24`, extending the prior exact bound `19`. Scripts pass normally and
  under `python -O` with frozen outputs.

## codex-2026-07-21 -- HYP-8905 bridge audit and MISTAKE-229

The exact binary homogeneous Hessian calculation from S103 survives and lands
inside THM-2063. The claimed `NC2 -> GMC(2) -> JC(2)` chain does not: the
general symmetric reduction from JC(2) lands in four variables, and no
Gaussian-to-Laplacian predicate map was supplied. S225's VC(4), planar
leading-form/Jelonek, and Lame-for-polygons programs remain separate routes
sharing a descent heuristic, not equivalent formulations. The dimension-three
rank/cycle comparison is heuristic, not a classification of collisions.

## boxeph-2026-07-21-S225 -- working JC(2): the obstruction is a descent termination; coprime intervals are its tool (HYP-8905)

**Owner:** long session working to prove the planar Jacobian Conjecture; pull past threads creatively.
## death-star-2026-07-21-S103 -- Planar JC IS NC2 one-sidedness: 2D nilpotent Hessian => one-sided; the JC-true/false boundary = the GMC2 unique-vs-coincident-cycle threshold. HYP-8910.

**Owner directive:** work to prove the planar Jacobian conjecture, pull in past threads creatively.

- **CONTEXT:** target = JC(2) (Alpoge THM-1300 killed JC dim>=3); THM-1830 puts my NC2/GMC2 upstream (NC2=>GMC2=>...=>JC(2)); Zhao VC (Hess P nilpotent <=> Delta^m(P^m)=0) = the SAME moment-vanishing as GMC2's E[P^m]=0.
- **PROVED/VERIFIED (planar_jc_..._S103.py):** 2D symmetric Hess P nilpotent <=> harmonic (trace 0) AND det 0; for harmonic P=A z^d+B zbar^d, det(Hess)=-4 d^2(d-1)^2 A B |z|^{2(d-2)} => det=0 <=> A B=0 => P ONE-SIDED (= verbatim NC2 conclusion). One-sided P is harmonic (Zhao VC trivial) + ONE-FIBER-LINEAR (F2-iF1=-iz => codex THM-2063 tame). So the SYMMETRIC planar JC / 2D Zhao-VC case IS NC2 one-sidedness, PROVED.
- **THE BOUNDARY (the unification):** 2x2 nilpotent has RANK<=1 = ONE isotropic direction (one-sided); dim>=3 reaches rank>=2 = MULTIPLE isotropic dirs = RESONANCE = Alpoge counterexample (JC false). rank<=1 => parallel gradients => functional dependence => one-fiber pencil (codex). SAME threshold as GMC2 S101 (unique vs coincident cycle) + boxeph S217 entropy (rigid=zero-entropy=one-sided). Planar = the last one-direction dimension.
- **HONEST:** proves the symmetric case (=NC2 one-sided); NOT full JC(2) (non-symmetric rank-1 parallel-gradient = codex THM-2063 open crux). A unification placing planar JC in my GMC2/NC2/resonance framework + a bridge to codex's pencil.
- **ADOPTED codex MISTAKE-227:** my S102 strict measure detects M(S)>1/14 only (vanishes on the tight AP = measure-zero boundary packet, not a counterexample); the integer-kernel=GMC2 unification SURVIVES as a strict-BULK observable, boundary handled by THM-2058 exact packets. reflection planar-jc-is-nc2-one-sidedness-...-S103. HYP-8910.

## codex-2026-07-21 -- HYP-8879 strict-measure boundary correction

**HONEST:** JC(2) OPEN (JC(n>=3) FALSE, Keller THM-1300). Worked it via 3 reduction PROGRAMS (none complete): (A) de Bondt+Zhao => VC(4) (JC(2)<=>a dim-4 Laplacian moment nullcone, doubling n->2n); (B) klein-S329 Euler-Zariski cover-degree-3 bootstrap (cuspidal Jelonek curve, ramification-parabola escape, NO CF); (C) mac-mini-S137 golden-corner/Lame (subtractive Euclid on Newton slopes).
- **MISTAKE-227:** S102's integral with `g=1[||x||>1/14]` detects only
  `M(S)>1/14`. It vanishes on the tight AP `{1,...,13}`, whose weak safe set is
  the six-point packet `U_14/14`, so zero measure is not an LRC counterexample.
- **Surviving bridge:** after Fejer convergence control, the integer-kernel
  expansion is a strict-bulk observable. THM-2058's exact denominator packets
  handle equality; finite Sidon/AP ratios prove no AP-core reduction.
- **THM-2060 proved with prior-art boundary:** the clock-independent order is
  `q=a/gcd(a,w)`, and every tail histogram bin has the sharp floor
  `q-ceil(q/7)`. The qualitative one-tail dodge was already THM-760/761/765;
  the new gain is exact CRT support, multiplicity, and clock composition.
- **THM-2061 proved reduction:** the only primitive imprimitive-eleven-core
  residue is the dyadic two-odd-tail seam. Strict failure is an open folded-
  diamond cover of the closed core-safe set; it forces divisor completeness
  through `14`, tails below `12 max(C)`, measure at most `4/63`, and has no
  survivor for normalized cores in `{1,...,19}`.
- **THM-2062 proved atlas sieve:** deletion determinantal indices turn
  hereditary primitivity on every saturated two-anchor interval into an exact
  squarefree CRT wheel with at most two bad projective directions per prime;
  rank-one deletions become affine `+-1` terminals plus a 1D coprime wheel.
- **THM-2064 incoming synthesis:** the independent common-clock capacity
  theorem proves the full multi-tail union bound and the same unique dyadic
  two-tail seam; THM-2060 is now explicitly its sharp histogram specialization.
- **THM-2065 proved circuit-ray collapse:** THM-2051's bounded support-three-
  to-five relation pulls back to either a persistent coefficient-row circuit
  or one primitive projective parameter. Thus every circuit-free two-anchor
  template has only finitely many strict-null rows, filtered exactly by the
  THM-2062 wheel. Persistent height-`2^20` marked circuits are the residual.

**VERIFIED (jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.py):** 2D symmetric case EASY -- nilpotent Hessian <=> P prop (x+iy)^d (harmonic), Delta^m(P^m)=0, invertible. Lame worst-case = Fibonacci (longest Euclid chain <200 = (144,89)). [Restricted JC(2) proved elsewhere: THM-1345 equivariant, THM-1370 elliptic, THM-1365 poly-Galois, geom-degree<=2, THM-2063.]

**CREATIVE UNIFICATION:** the JC(2) obstruction has THREE equivalent forms -- (A) VC(4) both-signs radial nonvanishing; (B) leading-form descent {P_A,Q_B}=0 propagating to a coordinate (THM-1345 s5); (C) Lame-for-polygons CF bound -- and ALL THREE are DESCENT/RETURN-TERMINATION problems = exactly what my coprime-interval/numerical-semigroup/Frobenius engine (S223 DvdK, S224 Wall A) handles; Lame-Fibonacci = the effective bound. codex THM-2045 (smooth factorized R has NO planar Jacobian mate, an exponent-semigroup/Newton-edge obstruction) already uses this engine.

**BOUNDED:** VC(4) GMC-like (E=L o CT) but GMC(2)NOT=>JC(2) (doubling to rank>=2, S205); JC<->LRC 'shared n=12' WITHDRAWN (S137's 12 = Fibonacci proxy). POSITIVE PRIOR: Keller collision min = dim 3, so n=2 is BELOW threshold -> expect JC(2) TRUE + provable by low-rank/coprime-interval means.

**Corrections adopted (mining):** CF thread = mac-mini-S137 (not klein-S329); reify-ladder = deathstar-S75 (not defunct THM-1750); apolarity/Fischer != THM-1685/1710/1735 (those are TNC). NOT a proof of JC(2) -- an assembled route map + verified pieces + the unified descent-termination obstruction. Artifacts: reflection working-jc2-the-obstruction-is-a-descent-termination-...-boxeph-S225.md, HYP-8905, script (+.out).

## codex-2026-07-21 -- concrete GMC(2) residue and conditional NC2 capstone checked

- `GMC2NormalizedResidue` now proves the complete normalized three-case
  calculation: non-dilated channels vanish by their multinomial, dilated
  off-face channels vanish by the factorial gap, and face channels reindex to
  the exact Frobenius power of the undilated face constant term.
- `GMC2SupportFaceBridge` proves exact geometric-face/support-face reindexing,
  identifies the specialized lifted seed with that constant term, and derives
  the global scaled floor, exact face height, and strict off-face gap from a
  concrete reference channel.
- `GMC2NC2` checks the whole post-specialization contradiction and exports
  `nc2_of_dvdK1_of_heightWitnessSupplier` plus the GMC(2) endpoint. The old
  one-`sorry` `GMC2NC2Capstone` WIP is replaced by a sorry-free compatibility
  surface. All audited declarations use only Lean's standard axiom trio.
- Direct checks of the concrete residue, support bridge, conditional capstone,
  and compatibility surface pass at the default heartbeat. The aggregator
  check stopped before elaboration at its pre-existing missing
  `GMC2GoodReduction.olean`; that dependency chain was not built under the
  no-big-build constraint.
- **Exact remaining boundary:** the reference-channel extractor and height
  theorem are separately green, but their direct existential wrapper into
  `HeightWitnessSupplier` deterministically exhausts elaboration (also at
  800k heartbeats) without a type or mathematical error. DvdK remains a
  separate published external premise. A final default-budget redesign
  compressed the reference to mass and balance: the compact extractor, base
  obligations, and base-contradiction adapter all elaborated, but the final
  `nc2_of_dvdK1` compositor still timed out at `whnf`; the experiment was
  reverted. No repository-wide build was run.
- **Live-main connection audited:** S223/HYP-8895 recasts the positive-
  coefficient one-variable DvdK return set as a numerical semigroup. That is
  a useful future formalization route, but it does not replace `DvdK1` here:
  the mixed-sign cancellation case still rests on S222's unfinished saddle
  argument. S103/HYP-8910 independently finds the same one-sided charge shape
  in a symmetric planar-Jacobian subcase, but explicitly supplies no transfer
  theorem to Gaussian moments and therefore does not change this Lean spine.
