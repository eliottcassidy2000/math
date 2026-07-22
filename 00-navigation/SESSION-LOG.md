## boxeph-2026-07-21-S227 -- doubling homeomorphism + mirror-parity (LRC reduction); full two-charge DvdK in Lean (HYP-8920)

**Owner:** complete GMC(2) formalization; LRC math -- doubling as a continuous bijection + the unique safe-child condition.

**LRC (mapping the owner's prompts to codex's frontier):** 'making doubling a continuous bijection' = THM-2075 (doubling D:t->2t is a HOMEOMORPHISM between consecutive safe sets, chi/#comp/endpoint INVARIANT per address-sheet); 'the unique hild-safe condition' = THM-2073 (unique safe-child law, dyadic tower). VERIFIED (doubling_homeomorphism_meets_mirror_parity_boxeph_S227.py): phi_{2C}(t)=phi_C(2t) => G_{2C}=D^{-1}G_C; D 2-to-1 on S^1 but homeomorphism on each binary-address half (deck tau:t->t+1/2, distinct from my mirror iota); chi even for dyadic-seam covering sets.

**SHARPENED REDUCTION (my contribution = S212 mirror-parity + codex THM-2073/2075):** a dyadic-seam disproof (S=2C ∪ {x,y}) has chi(G_{1/14})=0; chi doubling-invariant => terminal core also chi=0; chi even => needs chi=0 exactly. So Wall A (dyadic-seam case) <=> no hereditarily-primitive TERMINAL core has chi(G_{1/14})=0. Honest: full doubling 2-to-1 (chi doubles) vs per-sheet homeo (chi fixed); disproof chi=0 preserved either way.

**GMC(2) LEAN:** extended GMC2DvdKTwoCharge.lean with exists_nonzero_ct_pair' (SWAP orientation, index 0 = -n charge). Both pair orientations + dvdk1_pair KERNEL-PURE ([propext, Classical.choice, Quot.sound]) -- the full single-character DvdK1, DvdK-premise-free (any injective straddling Fin 2->Z up to relabeling). General DvdK1 (>=3 charges) = codex THM-2067 (Galois orbit-product) + height package = remaining boundary.

**Honest:** an LRC reduction (composing frontier theorems, verified in parts, not a Wall A proof) + a completed kernel-pure single-character DvdK1 leaf. Artifacts: reflection doubling-homeomorphism-plus-mirror-parity-...-boxeph-S227.md, HYP-8920, script (+.out), Lean GMC2DvdKTwoCharge.lean.

