## boxeph-2026-07-21-S213 -- chirality + toothpick A139250 + tournament evenness = one Lefschetz parity (HYP-8850)
> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be corrected after filing. Start with
> [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## codex-2026-07-21-LRC-two-anchor -- THM-2052 explicit star-atlas refinement

- **PROVED refinement:** the rank-eleven bounded triple-relation code is not
  merely an unspecified rational plane. Either its rank is twelve and
  Hadamard gives a finite primitive speed box, or two anchors `p,q` support
  eleven independent private-coordinate relations
  `a_i v_p+b_i v_q+c_i v_i=0`, with `c_i!=0` and coefficient height at most
  `91^6`. Modulo dilation this is a finite list of one-projective-parameter
  star families.
- **CLUSTER BRIDGE:** height-bounded pair relations first lock coordinates
  into scale components. Genuine bounded triples among distinct components
  then have quotient rank at least `c-2`, leaving at most two cluster scales.
  This connects THM-2051's higher-circuit branch to the cross-scale
  decorrelation program without discarding signed coefficients.
- **OPEN RESIDUAL:** classify the one-parameter star families by a positive
  pair-sum margin or an owner-labelled Euler boundary point; LRC(14) remains
  open.

## death-star-2026-07-21-S95 -- THM-2022 non-cancellation identity assembled (face survives as Q̄^p) + DvdK Lean-formalization ROADMAP (Monsky-algebraic ~4-9 person-months; analytic ~years; symbolic-coeff caveat kills the concrete shortcut). HYP-8805.

**Owner:** connect previous chirality work + toothpick A139250; tournament values that are EVEN (iso classes) and what toothpick diagrams correspond to their structure.

**ONE LAW:** a count graded by a Z/2 INVOLUTION obeys count=(#fixed)+2*(#free pairs) => count == #fixed (mod 2) [Lefschetz/Euler]. Master involution = reversal/converse R:T->T^op (=complement=antipodal, THM-584). Three realizations, all verified (chirality_toothpick_tournament_parity_boxeph_S213.py):

1. TOURNAMENTS EVEN: A000568=SC+2*(chiral pairs). THM-587 (PROVED): P_n(1)=A000568, P_n(-1)=SC=self-converse=antipodal Euler/Lefschetz number; SC=2,2,8,12,88,176 (n=3..8) all EVEN => A000568 EVEN for n>=3. Enumeration n<=6: (total2/SC2/0pairs),(4/2/1),(12/8/2),(56/12/22), count=SC+2pairs holds. THM-1830 BLUE self-comp atom=1 iff n odd = local single-fixed-point face.
2. TOOTHPICK A139250: sim=OEIS exactly; D4 mirror sym; =(#axis)+2*(#pairs); #axis ODD (central seed) => A139250 ODD all n>=1; first differences have exactly ONE odd term (seed); A139250(2^k)=(2^(2k+1)+1)/3=1,3,11,43,171 Jacobsthal.
3. LRC (S212): chi(G_delta) == #iota-fixed lonely pts (mod 2).

**TOOTHPICK<->CHIRALITY DICTIONARY:** central seed toothpick (axis-fixed) <-> achiral/self-converse classes; mirror toothpick PAIR <-> chiral pair {T,T^op}; D4 mirror <-> R=converse=complement=antipode; A139250 parity=#axis(mod2) <-> A000568 parity=SC=P_n(-1)(mod2); doubling A139250(2^k) <-> level-graded signed-cycle-index recursion. The toothpick fractal = the PICTURE of 'one achiral seed generating chiral mirror-pairs'. Unifies with S212: (P(+1),P(-1)) = Burnside count / Lefschetz fixed set = the chirality Euler class (S212 even chi + odd Gauss sum i*sqrt7).

**Honest:** THM-587, the toothpick OEIS/oddness/Jacobsthal facts, and S212 are each proved/verified; the contribution is the UNIFICATION (one reversal, one Lefschetz parity count==#fixed) + the explicit toothpick<->chirality dictionary. Artifacts: reflection chirality-toothpicks-and-why-tournament-counts-are-even-one-lefschetz-parity-boxeph-S213.md, HYP-8850, script (+.out).
**Owner:** push the DC(2)/planar-JC thread to its next decisive target and transfer the proved/disproved mechanisms toward LRC(14), while repeatedly integrating incoming work.

**DC(2) result (THM-2049, PROVED local/formal statement; not DC(2)):** in the exact Ore algebra `Q[x,q][ell;delta]`, `beta(sum a_k ell^k)=min_k(v_x(a_k)-2k)` is multiplicative and commutators raise beta by two. The associated bracket is `{ell,q}_0=2`. For the Weyl boundary symbols, the simultaneous grade-`g` correction map is `(A,B)->(8/3)(u-2)A+(2u^2-10u+9)B/9`; it is surjective because the two `u` polynomials are coprime. Thus the grade-six residual is exact. An exact ladder advances grades `6,...,13` to `14`; a formal beta-adic `[S,T]=1` lift exists. The open gates are polynomial termination and the coupled `D` relations. This corrects HYP-8802/8803's earlier suggestion that the first invariant grade might carry the obstruction.

**LRC no-go (THM-2050, PROVED):** AP13 and `AP13` with `12->26` have identical full local phase-height function germs on `|h|<1/728` at every unit point `a/14`, yet `M=1/14` and `M=1/12`. Local top data, even as a full germ, cannot determine global loneliness.

**Incoming synthesis:** THM-2047 supplies the lossless signed phase-height/Euler carrier; THM-2048 supplies the fiber-quantization pruning tax; HYP-8840 identifies GMC's constant-term/volume leverage and its zero-volume ceiling. Later pulls supplied THM-2048's genuine Cover14 gain and strengthened THM-2051: after paying pair covariance exactly, the no-support-`3..5`-relation branch at height `2^21` has positive safe volume. THM-2052 then proves every hypothetical counterexample already has eleven independent bounded three-support relations and lies in a finite two-dimensional rational atlas. The exact transfer is `volume/tax -> strict branch`, `Euler signed wall word -> tight branch`, and a labelled Noetherian rank-or-Euler rule inside the relation branch as the missing glue. No literal algebra map between GMC/DC and LRC is asserted.

**Exact termination-sidecar audit (HYP-8841):** pair-sum maxima, threshold interval/point topology, complete first exits, and every peel tax were computed on AP/GW, `12->26`, `12->36`, `12->96`, `12->84`, P10+K33, and the incoming Cover14 tax-gain row. It exactly reproduces the latter's peel-`93` excess `2413467317/235670635200`. The tax fires on deep/covering controls but misses the smallest hostile/K33 controls and is not a scalar termination height. THM-2047 proves the strict search is complete by `q<=2 max(S)`. With THM-2051 now proved, the remainder lies in the bounded small-relation branch, but every control already has a height-one support-three circuit (the first seven share `1+2=3`). Raw circuit existence ranks last in the enlarged carrier tournament. The next decisive Wall-A clause is therefore an active-owner circuit/Euler endpoint-survival lemma when neither a tax violation nor a positive pair-sum margin occurs.

**Past-work pull into the decisive clause:** HYP-2108 already gives the exact endpoint-cover functional. If a core-safe component has midpoint `m_i` and length `l_i`, it is swallowed by the open danger arcs of peel `w` iff `||wm_i||+(w/2)l_i<1/14`. Hence the active-owner target is precisely `P_w(C)=max_i(||wm_i||+(w/2)l_i-1/14)>=0`. With THM-2052, the next theorem becomes a rank-or-Euler alternative: some peel has `P_w>=0`, or active owner data supplies a twelfth independent bounded relation and reaches the finite maximal-minor terminal. HYP-3117/HYP-3120 identify endpoint incidence as the missing proof-circuit input; HYP-8845 halves the covering case because a survivor has a mirror partner and `chi>=2`.

**Artifacts:** THM-2049, THM-2050, HYP-8841, the updated exact Ore script/output, `lrc14_termination_sidecar_codex_20260721.py/.out`, and reflection `from-Ore-boundary-acyclicity-to-LRC14-Euler-termination-codex-20260721.md`.
**Owner:** look for other topological advances the repo has made; come up with creative LRC arguments combining and extending them.

**PULL (repo topological toolkit, credited):** THM-2047 (codex) chi(G_delta)=#components, LRC(14)<=>chi(G_{1/14})>0 [PROVED]; HYP-3015 (codex-S179) {G_delta}=superlevel filtration of f_S, M(S)=top death, persistence barcode; opus lonely-set Euler-char certificate + kps cohomological_three_distance Alexander duality (b0(lonely)=b0(cover), arcs alternate); HYP-3025 arc-Cech nerve + Betti-defect sidecar, HYP-3101 normal-fan Cech barcode component bound (open); kps-S19 LRC Lefschetz = free iota:t->1-t + Gauss sum i*sqrt(7) (ordinary Lefschetz blind); THM-587 metagraph reversal Lefschetz (tournament side, PROVED). P1-P5 re-verify these.

**NEW (P6, verified):** f_S(1-t)=f_S(t) => G_delta is iota-INVARIANT (iota = my S210 involution). iota's fixed points {0,1/2}: f_S(0)=0 always, f_S(1/2)=0 iff some speed EVEN. Every COVERING set has an even speed => both fixed points dangerous => iota acts FREELY on G_delta => chi(G_delta) EVEN. So codex's chi>0 sharpens to:
  LRC(14) for covering S  <=>  chi(G_{1/14})>=2  <=>  a MIRROR PAIR {t*,1-t*} of lonely windows survives.
Verified: deep well {1..12,182} chi=24 (12 pairs); tight (1,2,3)@1/4={1/4,3/4} chi=2; all-odd (1,3,5,7) = iota-FIXED exception (1/2 lonely, chi=1 ODD = Borsuk-Ulam fixed point, the classical all-odd-lonely case).

**LEVERAGE:** (1) equivariant HALVING of Wall A -- find one lonely window in [0,1/2], mirror automatic; (2) parity obstruction -- chi even => never 1; a disproof needs chi=0 = every mirror pair killed simultaneously (iota-symmetric covering); (3) kps-S19's Lambda(iota)=0 blind BECAUSE free; the odd-equivariant index = Gauss sum i*sqrt(7) is the Borsuk-Ulam obstruction on G_delta/iota. The equivariant chi (even) + odd index (i*sqrt7) are the two halves of the Z/2-equivariant Euler class of the good set. Topological form of THM-1820 mirror pairs (B3) + S210.

**Honest:** P1-P5 re-verify existing fleet toolkit (credited, not claimed); the equivariant even-chi mirror-parity sharpening (P6) is new and verified. Forcing chi>=2 for every 13-speed covering core = LRC itself, OPEN; reduces Wall A to 'G_{1/14}(C)/iota on [0,1/2] nonempty'. Artifacts: reflection the-good-sets-reversal-symmetry-...-boxeph-S212.md, HYP-8845, script (+.out).

