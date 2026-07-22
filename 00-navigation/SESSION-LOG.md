## boxeph-2026-07-21-S227 -- doubling homeomorphism + mirror-parity (LRC reduction); full two-charge DvdK in Lean (HYP-8920)

**Owner:** complete GMC(2) formalization; LRC math -- doubling as a continuous bijection + the unique safe-child condition.

**LRC (mapping the owner's prompts to codex's frontier):** 'making doubling a continuous bijection' = THM-2075 (doubling D:t->2t is a HOMEOMORPHISM between consecutive safe sets, chi/#comp/endpoint INVARIANT per address-sheet); 'the unique hild-safe condition' = THM-2073 (unique safe-child law, dyadic tower). VERIFIED (doubling_homeomorphism_meets_mirror_parity_boxeph_S227.py): phi_{2C}(t)=phi_C(2t) => G_{2C}=D^{-1}G_C; D 2-to-1 on S^1 but homeomorphism on each binary-address half (deck tau:t->t+1/2, distinct from my mirror iota); chi even for dyadic-seam covering sets.
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
  terminal. THM-2076's Haar-capacity lemma forces terminal size at least six,
  sharpening depth to at most five. Exact referees pass normally and under
  `python -O`. THM-2075 then proves that doubling is a homeomorphism along the
  tower: component/Euler counts persist, lengths and measure halve exactly,
  each component carries one constant binary address, and every endpoint has
  an inherited terminal-core owner. THM-2078 uses the terminal guard-height
  bound and an exact denominator-`8192` bitset audit to close every nontrivial
  tower with terminal maximum at most `24`: all `4,484,931` cores of ranks
  `6,...,10` were filtered, `30,594` were hereditary/divisor-complete, and no
  allowed guard survived even the necessary rational grid. LRC(14) remains
  open on the unbounded hereditary terminal lane and its address assignment.
- **Hostile parity correction:** MISTAKE-230 retracts HYP-8920's claimed
  `chi(G_S)=0 -> chi(G_terminal)=0` descent. THM-2075 starts at `G_C`, not at
  the empty `G_S`; the original tails kill both outer children. Mirror
  evenness of the nonempty terminal survives, but it does not discharge tail
  coverage. The S227 two-charge DvdK Lean work is independent and unaffected.

**SHARPENED REDUCTION (my contribution = S212 mirror-parity + codex THM-2073/2075):** a dyadic-seam disproof (S=2C ∪ {x,y}) has chi(G_{1/14})=0; chi doubling-invariant => terminal core also chi=0; chi even => needs chi=0 exactly. So Wall A (dyadic-seam case) <=> no hereditarily-primitive TERMINAL core has chi(G_{1/14})=0. Honest: full doubling 2-to-1 (chi doubles) vs per-sheet homeo (chi fixed); disproof chi=0 preserved either way.

**GMC(2) LEAN:** extended GMC2DvdKTwoCharge.lean with exists_nonzero_ct_pair' (SWAP orientation, index 0 = -n charge). Both pair orientations + dvdk1_pair KERNEL-PURE ([propext, Classical.choice, Quot.sound]) -- the full single-character DvdK1, DvdK-premise-free (any injective straddling Fin 2->Z up to relabeling). General DvdK1 (>=3 charges) = codex THM-2067 (Galois orbit-product) + height package = remaining boundary.

**Honest:** an LRC reduction (composing frontier theorems, verified in parts, not a Wall A proof) + a completed kernel-pure single-character DvdK1 leaf. Artifacts: reflection doubling-homeomorphism-plus-mirror-parity-...-boxeph-S227.md, HYP-8920, script (+.out), Lean GMC2DvdKTwoCharge.lean.

