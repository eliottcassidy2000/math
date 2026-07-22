## boxeph-2026-07-21-S212 -- the good set's reversal symmetry: an equivariant (mirror-parity) sharpening of codex's chi(G_delta) criterion (HYP-8845)

- **PROVED (THM-2051, strengthened):** THM-965 gives the sharp pair covariance
  floor `-6/637`, attained exactly at reduced ratio `1:13`. Paying all `78`
  pairs exactly and applying the whole-product argument only above support two
  proves: if a thirteen-speed row has no exact support-three-through-five
  relation of coefficient height at most `2^21`, then its continuous quintic
  Bonferroni functional is positive, hence it has a positive-measure strict
  lonely set.
- **MECHANISM:** approximate the centered danger indicator by a degree-`H`
  Fejer polynomial. Relation dissociation kills every finite product constant
  term. A BV translation estimate controls each whole signed centered-product
  tail in `L1`. After the pair payment, the reserve is `1188/16807`; the
  higher-support error is below `25791480/559473817`, leaving the exact margin
  `96283836/3916316719>0`. The original all-support `2^20` version also holds.
- **SYNTHESIS:** this bypasses THM-946's still-open absolute `T4/T5`
  strip/slab estimates for the coarse alternative and converts HYP-8841's
  ambient termination problem into descent on a finite union of bounded-
  height genuine higher-relation hyperplanes. Twelve independent relation rows
  would bound primitive speeds by `5^6 H^12` via maximal minors and Hadamard;
  the missing supplier must force new independent rows. This does not classify
  the structured branch;
  the tight AP is deliberately relation-rich. THM-940 is only the discrete
  analogue and is not silently identified with the continuous proof.
## codex-2026-07-21-DC2-LRC14-termination -- local acyclicity versus finite Euler termination (THM-2049/2050, HYP-8841)

**Owner:** push the DC(2)/planar-JC thread to its next decisive target and transfer the proved/disproved mechanisms toward LRC(14), while repeatedly integrating incoming work.

**DC(2) result (THM-2049, PROVED local/formal statement; not DC(2)):** in the exact Ore algebra `Q[x,q][ell;delta]`, `beta(sum a_k ell^k)=min_k(v_x(a_k)-2k)` is multiplicative and commutators raise beta by two. The associated bracket is `{ell,q}_0=2`. For the Weyl boundary symbols, the simultaneous grade-`g` correction map is `(A,B)->(8/3)(u-2)A+(2u^2-10u+9)B/9`; it is surjective because the two `u` polynomials are coprime. Thus the grade-six residual is exact. An exact ladder advances grades `6,...,13` to `14`; a formal beta-adic `[S,T]=1` lift exists. The open gates are polynomial termination and the coupled `D` relations. This corrects HYP-8802/8803's earlier suggestion that the first invariant grade might carry the obstruction.

**LRC no-go (THM-2050, PROVED):** AP13 and `AP13` with `12->26` have identical full local phase-height function germs on `|h|<1/728` at every unit point `a/14`, yet `M=1/14` and `M=1/12`. Local top data, even as a full germ, cannot determine global loneliness.

**Incoming synthesis:** THM-2047 supplies the lossless signed phase-height/Euler carrier; THM-2048 supplies the fiber-quantization pruning tax; HYP-8840 identifies GMC's constant-term/volume leverage and its zero-volume ceiling. Later pulls supplied THM-2048's genuine Cover14 gain and promoted THM-2051: the no-small-relation branch has positive safe volume, so every hard row lies on a support-`2..5`, coefficient-height-`<=2^20` circuit hyperplane. The exact transfer is `volume/tax -> strict branch`, `Euler signed wall word -> tight branch`, and a labelled Noetherian deletion rule inside the circuit branch as the missing glue. No literal algebra map between GMC/DC and LRC is asserted.

**Exact termination-sidecar audit (HYP-8841):** pair-sum maxima, threshold interval/point topology, complete first exits, and every peel tax were computed on AP/GW, `12->26`, `12->36`, `12->96`, `12->84`, P10+K33, and the incoming Cover14 tax-gain row. It exactly reproduces the latter's peel-`93` excess `2413467317/235670635200`. The tax fires on deep/covering controls but misses the smallest hostile/K33 controls and is not a scalar termination height. THM-2047 proves the strict search is complete by `q<=2 max(S)`. With THM-2051 now proved, the next decisive Wall-A clause is owner-labelled endpoint survival inside the bounded small-relation branch when neither a volume-tax violation nor a positive pair-sum margin occurs. The proof-carrier tournament is transitive; signed threshold topology wins and raw unit germs lose.

**Artifacts:** THM-2049, THM-2050, HYP-8841, the updated exact Ore script/output, `lrc14_termination_sidecar_codex_20260721.py/.out`, and reflection `from-Ore-boundary-acyclicity-to-LRC14-Euler-termination-codex-20260721.md`.
**Owner:** look for other topological advances the repo has made; come up with creative LRC arguments combining and extending them.

**PULL (repo topological toolkit, credited):** THM-2047 (codex) chi(G_delta)=#components, LRC(14)<=>chi(G_{1/14})>0 [PROVED]; HYP-3015 (codex-S179) {G_delta}=superlevel filtration of f_S, M(S)=top death, persistence barcode; opus lonely-set Euler-char certificate + kps cohomological_three_distance Alexander duality (b0(lonely)=b0(cover), arcs alternate); HYP-3025 arc-Cech nerve + Betti-defect sidecar, HYP-3101 normal-fan Cech barcode component bound (open); kps-S19 LRC Lefschetz = free iota:t->1-t + Gauss sum i*sqrt(7) (ordinary Lefschetz blind); THM-587 metagraph reversal Lefschetz (tournament side, PROVED). P1-P5 re-verify these.

**NEW (P6, verified):** f_S(1-t)=f_S(t) => G_delta is iota-INVARIANT (iota = my S210 involution). iota's fixed points {0,1/2}: f_S(0)=0 always, f_S(1/2)=0 iff some speed EVEN. Every COVERING set has an even speed => both fixed points dangerous => iota acts FREELY on G_delta => chi(G_delta) EVEN. So codex's chi>0 sharpens to:
  LRC(14) for covering S  <=>  chi(G_{1/14})>=2  <=>  a MIRROR PAIR {t*,1-t*} of lonely windows survives.
Verified: deep well {1..12,182} chi=24 (12 pairs); tight (1,2,3)@1/4={1/4,3/4} chi=2; all-odd (1,3,5,7) = iota-FIXED exception (1/2 lonely, chi=1 ODD = Borsuk-Ulam fixed point, the classical all-odd-lonely case).

**LEVERAGE:** (1) equivariant HALVING of Wall A -- find one lonely window in [0,1/2], mirror automatic; (2) parity obstruction -- chi even => never 1; a disproof needs chi=0 = every mirror pair killed simultaneously (iota-symmetric covering); (3) kps-S19's Lambda(iota)=0 blind BECAUSE free; the odd-equivariant index = Gauss sum i*sqrt(7) is the Borsuk-Ulam obstruction on G_delta/iota. The equivariant chi (even) + odd index (i*sqrt7) are the two halves of the Z/2-equivariant Euler class of the good set. Topological form of THM-1820 mirror pairs (B3) + S210.

**Honest:** P1-P5 re-verify existing fleet toolkit (credited, not claimed); the equivariant even-chi mirror-parity sharpening (P6) is new and verified. Forcing chi>=2 for every 13-speed covering core = LRC itself, OPEN; reduces Wall A to 'G_{1/14}(C)/iota on [0,1/2] nonempty'. Artifacts: reflection the-good-sets-reversal-symmetry-...-boxeph-S212.md, HYP-8845, script (+.out).

