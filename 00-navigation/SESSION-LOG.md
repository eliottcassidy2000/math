## boxeph-2026-07-21-S228 -- kernel-pure Lean: the positive-coefficient DvdK1 (any support); cancellation is the sole crux (HYP-8925)

## codex-2026-07-22 -- mixed guard fold removes the depth-five terminal

**Owner directive:** continue the long LRC(14) session, mine forgotten ideas,
and seek a scale-free closure of the dyadic terminal lane.

- THM-2080 proves the exact mixed-radius overlap law by one-sided interval
  atoms (equivalently, a two-fold Bernoulli law).
- The direction audit is load-bearing: `G_Q subset E_h` makes the guard
  complement, not the guard, a subset of the union of danger combs. MISTAKE-231
  retracts the first reversed-cover consumer while retaining the fold formula.
- For an odd guard, every mixed overlap is at least `1/42`, with equality only
  at `q=6h`. Hence each danger comb covers at most `5/42` outside the guard.
  Six distinct combs cover strictly less than `5/7`, so terminal rank six is
  impossible and the dyadic tower has depth at most four.
- The authoritative exact referee checks `4,032` pair/formula cases, the whole
  small-product equality ledger, and `11,088` hostile direct containments; its
  smallest hostile complement remains positive.
- Assumption challenge: the target carrier is the guard complement. A star
  tournament picture is harmless but unnecessary; equality of the six leaf
  capacities would force all six distinct speeds to equal `6h`.

## death-star-2026-07-22-S104 -- GMC2 formalization: pinpointed + wrote the last capstone discharge (HeightWitnessSupplier); structurally correct + statements axiom-checked, but the proof hits a pathological whnf wall (>6.4M heartbeats). One perf-fix from clean DvdK1 -> NC2.
**Owner:** aim earnestly at formalizing DvdK; make it simpler / circumvent it; spill over to LRC.

**FORMALIZED (GMC2DvdKPositive.lean, kernel-pure [propext,Classical.choice,Quot.sound]):**
- ct_pos_of_balanced: c_i>0 + any balanced composition r0 of size m => CT(f^m)>0 (sum of positive terms; ANY support/#charges).
- exists_balanced_of_twosided: two-sided q_i>0, q_j<0 => a balanced composition exists (|q_j| copies of + charge, |q_i| of -, at m=|q_i|+|q_j|).
- dvdk1_positive: two-sided support + c_i>0 => exists m>=1, CT(f^m)>0. The positive-coefficient DvdK1, DvdK-premise-free, any support.

**STRUCTURE (the simplification):** CT(f^m) as a polynomial in c has all-POSITIVE multinomial coefficients => nonzero <=> a balanced composition exists (feasibility, elementary). The SOLE difficulty is CANCELLATION for specific complex c (the counterexample f=u^2+u+u^-1-u^-2: CT=0 for all odd m, CT(f^4)=-12). Feasibility + positive-coeff + two-charge (unique composition, S226/S227) elementary/formalized; >=3 charges complex = codex THM-2067 (Galois orbit-product) = the next Lean target. A true circumvention (complex->positive) fails (S222/S223 retracted); Galois is what rules out complex cancellation.

**LRC SPILLOVER:** same positive-vs-cancellation split -- LRC |G_delta|=sum_{k.v=0} prod ghat is ALL sign-cancellation (signed sinc, NO positive regime), why chi/topology (THM-2047/S212) beats volume (S211), and symmetry (mirror iota, doubling homeo THM-2075) is the LRC analogue of Galois for taming sign. Adopt codex MISTAKE-230: my S227 chi=0 descent RETRACTED (tower transports nonempty CORE sets, terminal core has safe interval => chi>0, THM-2077); doubling identity/homeo/mirror-parity survive individually.

**Honest:** the no-cancellation regime of DvdK1 now kernel-pure in Lean (any support), isolating cancellation as the sole crux; NOT the general complex DvdK1 (THM-2067). Artifacts: reflection starting-to-formalize-dvdk-...-boxeph-S228.md, HYP-8925, Lean GMC2DvdKPositive.lean (4 theorems).

