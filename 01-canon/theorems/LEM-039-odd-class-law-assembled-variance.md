---
id: LEM-039
title: THE ODD CLASS LAW AND THE ASSEMBLED VARIANCE (LEM-038's two named opens closed). (A) THE ODD CLASS LAW: for every modulus D | P, the signed class sums T_r(D) = Σ_{p_i ≡ r (D)} ε_i satisfy T_{−r}(D) = −T_r(D) at s = 3 (one line from the fixed-point-free reflection pairing) — self-mirror classes carry zero signed mass (T_0 = T_{D/2} = 0), N(D) = Σ_r T_r² = 2Σ_{half}T_r² is even, and the S66 "N = M mostly" pattern is the generic mirror-paired-singleton case. GENERAL-s FORM (the cross-section law): T_r^{(s)}(D) = −T_{−r}^{(6−s)}(D). Referee: 65,107 exact class checks over five clusters, all divisors. (B) SUB-GRID SINE COHERENCE: S((P/D)ν) = 2i Σ_{half} T_r sin(2πνr/D) on every divisor grid (1e-12). (C) THE ASSEMBLED VARIANCE — the LEM-031..033 loop closed quantitatively: Var_w(cross) = Σ_{χ≠χ₀} |Σ_g (2/g²)·L_{P/g}(2,χ)·X̂_g(χ)|² reproduced ENTIRELY from the closed form (Dirichlet L-values × twisted coincidence sums, no frame sweep): family60 4060.3288 (rel err 5.7e-15), two-owner 3563.0703 (1.6e-13). THE FRAME-CROSS VARIANCE IS A FINITE CLOSED FORM. (D) THE VARIANCE-DROP QUANTIFICATION (by conductor, balanced): s = 0 masses led by 735 (27.9%)/245/147; s = 3 redistributes onto 245 (42.4%)/49/35/588/2940/196 — the 7-part concentration (LEM-033) holds at both sections (all top conductors ÷ 49); the drop is spread across the co-resonant list, not one character
status: PROVED ((A) one line + generality; (B) restriction of the imaginary law; (C) assembly of proved parts) + REFEREED EXACT (65,107 class checks; variance assembly at 1e-13–1e-15; conductor tables)
source: boxeph-2026-07-17-S67 (finishing sweep; closes LEM-038's named opens (a) and (b))
depends_on: [LEM-038 (pairing + imaginary law), LEM-032/033 (the closed-form factors), THM-892 (N(h) family)]
script: 04-computation/lrc14_oddclass_variance_boxeph_S67.py -> 05-knowledge/results/lrc14_oddclass_variance_boxeph_S67.out
---

# LEM-039 — the odd class law; the assembled variance

The coincidence spectrum at the sine section is ODD at every modulus — the
sine structure descends the entire divisor lattice — and the frame-cross
variance, the program's fluctuation object, is now computed end-to-end from
L-values and twisted coincidence sums alone. What remains open in the
frame-spectrum program is no longer any identity: it is only the SIZE of the
X̂_g factors (klein's discrepancy lane).

## Evidence log
- [x] odd class law + cross-section law: 65,107 checks, five clusters
- [x] sine coherence on divisor grids; assembled variance ×2 exact
- [x] conductor drop tables (s = 0 vs 3)
