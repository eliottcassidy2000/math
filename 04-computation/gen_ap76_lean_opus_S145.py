"""
gen_ap76_lean_opus_S145.py -- generates TournamentH7/LRCAP76Certificate.lean:
the tight AP76 certificate (muGood (1/7) AP76 >= 2314528732/40290957525) as a
24-interval instantiation of the AP44 pattern (boxeph-S2's handoff), discharging
the AP76Certificate Prop that LRCTailDiameter's k=13 diam<=75 assembly consumes.
"""
from fractions import Fraction as F

CELLS = [
    (0, 1, 1, 76, F(0), F(2, 175)),
    (12, 73, 1, 6, F(78, 469), F(1, 6)),
    (1, 6, 12, 71, F(1, 6), F(76, 455)),
    (15, 76, 1, 5, F(99, 497), F(1, 5)),
    (1, 5, 15, 74, F(1, 5), F(97, 483)),
    (18, 73, 1, 4, F(40, 161), F(1, 4)),
    (1, 4, 19, 75, F(1, 4), F(125, 497)),
    (25, 76, 1, 3, F(169, 511), F(1, 3)),
    (1, 3, 25, 74, F(1, 3), F(167, 497)),
    (29, 73, 2, 5, F(95, 238), F(2, 5)),
    (2, 5, 29, 72, F(2, 5), F(188, 469)),
    (37, 75, 1, 2, F(253, 511), F(1, 2)),
    (1, 2, 38, 75, F(1, 2), F(258, 511)),
    (43, 72, 3, 5, F(281, 469), F(3, 5)),
    (3, 5, 44, 73, F(3, 5), F(143, 238)),
    (49, 74, 2, 3, F(330, 497), F(2, 3)),
    (2, 3, 51, 76, F(2, 3), F(342, 511)),
    (56, 75, 3, 4, F(372, 497), F(3, 4)),
    (3, 4, 55, 73, F(3, 4), F(121, 161)),
    (59, 74, 4, 5, F(386, 483), F(4, 5)),
    (4, 5, 61, 76, F(4, 5), F(398, 497)),
    (59, 71, 5, 6, F(379, 455), F(5, 6)),
    (5, 6, 61, 73, F(5, 6), F(391, 469)),
    (75, 76, 1, 1, F(173, 175), F(1)),
]

def lit(fr):
    if fr.denominator == 1:
        return f"{fr.numerator}"
    return f"{fr.numerator} / {fr.denominator}"

def gen():
    tot = sum(d - c for (_, _, _, _, c, d) in CELLS)
    assert tot == F(2314528732, 40290957525)
    L = []
    A = L.append
    A("/-")
    A("  TournamentH7.LRCAP76Certificate — THE TIGHT AP₇₆ DENSITY-FLOOR CERTIFICATE.")
    A("  opus-2026-07-07-S145 (HYP-5197), executing boxeph-S2's handoff with their bridge:")
    A("  the 24 roof-superlevel intervals of the q ≤ 6 Farey-76 cells, each in Good(1/7)(AP₇₆)")
    A("  by one `good_of_roof_gt` call, pairwise disjoint (sorted), summing EXACTLY to")
    A("  2314528732/40290957525 — the published exact value of μ_{1/7}(AP₇₆).")
    A("  This discharges `TailDiameter.AP76Certificate`, making the k=13 diameter ≤ 75 leg")
    A("  (`hlarge_floor_of_diam_le`) unconditional.")
    A("  Cell data verified exactly: 04-computation/gen_ap76_lean_opus_S145.py.")
    A("  Kernel-pure: no sorry, no native_decide.")
    A("-/")
    A("import Mathlib")
    A("import TournamentH7.LRCFareyRoofBridge")
    A("import TournamentH7.LRCTailDiameter")
    A("")
    A("namespace LonelyRunner")
    A("namespace AP76Certificate")
    A("")
    A("open TailDiameter FareyRoofBridge MeasureTheory")
    A("open scoped ENNReal")
    A("")
    # interval lemmas
    for i, (p, q, pp, qq, c, d) in enumerate(CELLS, start=1):
        side_R = (c == F(p, q))  # node at left end -> crossing at d -> roof uses hxu
        A(f"/-- Interval {i}: Farey-76 cell `({p}/{q}, {pp}/{qq})`, roof-superlevel"
          f" `({lit(c)}, {lit(d)})`. -/")
        A(f"theorem I{i}_good :")
        A(f"    Set.Ioo ({lit(c)} : ℝ) ({lit(d)}) ⊆"
          f" Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by")
        A("  intro x hx; obtain ⟨hxl, hxu⟩ := hx")
        roof_h = "hxu" if side_R else "hxl"
        A(f"  refine ⟨good_of_roof_gt (p := ({p} : ℤ)) (q := ({q} : ℤ))"
          f" (p' := ({pp} : ℤ)) (q' := ({qq} : ℤ))")
        A("    (k := (76 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))")
        A("    (by norm_num) (by norm_num) (by norm_num) (by norm_num)")
        A(f"    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu])"
          f" (by push_cast; nlinarith [{roof_h}]),")
        lo = "le_of_lt hxl" if c == 0 else "by nlinarith [hxl]"
        hi = "le_of_lt hxu" if d == 1 else "by nlinarith [hxu]"
        A(f"    {lo}, {hi}⟩")
        A("")
    # disjointness helper
    A("/-- disjointness of two `Ioo` when the first ends at or before the second starts. -/")
    A("private theorem ioo_disj {a b c e : ℝ} (h : b ≤ c) :")
    A("    Disjoint (Set.Ioo a b) (Set.Ioo c e) := by")
    A("  apply Set.disjoint_left.mpr")
    A("  intro x hx1 hx2")
    A("  exact absurd (lt_of_lt_of_le hx1.2 (le_trans h (le_of_lt hx2.1))) (lt_irrefl x)")
    A("")
    # main certificate on Icc 1 76
    A("/-- **THE AP₇₆ CERTIFICATE on `{1,…,76}`.**  The 24 disjoint roof intervals sum to")
    A("the exact value. -/")
    A("theorem ap76_sum :")
    A("    ENNReal.ofReal ((2314528732 : ℝ) / 40290957525) ≤")
    A("      muGood (1 / 7) (Finset.Icc (1 : ℤ) 76) := by")
    for i, (_, _, _, _, c, d) in enumerate(CELLS, start=1):
        A(f"  set I{i} := Set.Ioo ({lit(c)} : ℝ) ({lit(d)}) with hI{i}")
    # pairwise disjointness haves
    A("  -- pairwise disjointness (intervals are sorted)")
    for j in range(2, 25):
        for i in range(1, j):
            A(f"  have d{i}_{j} : Disjoint I{i} I{j} := ioo_disj (by norm_num)")
    # union
    ustr = "I1"
    for i in range(2, 25):
        ustr += f" ∪ I{i}"
    A(f"  -- union ⊆ Good ∩ [0,1]")
    A(f"  have hu : {ustr} ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 76) ∩ Set.Icc (0 : ℝ) 1 := by")
    subchain = "I1_good"
    for i in range(2, 25):
        subchain = f"(Set.union_subset {subchain} I{i}_good)"
    A(f"    exact {subchain[1:-1] if subchain.startswith('(') else subchain}")
    # measure chain
    A("  -- measure of the disjoint union = sum of lengths (peel from the right)")
    vols = " + ".join(f"volume I{i}" for i in range(1, 25))
    A(f"  have hvol : volume ({ustr}) = {vols} := by")
    rws = []
    for k in range(24, 1, -1):
        # Disjoint (I1 ∪ ... ∪ I_{k-1}) I_k
        if k == 2:
            dis = "d1_2"
        else:
            dis = f"d1_{k}"
            for i in range(2, k):
                dis = f"({dis}.union_left d{i}_{k})"
        rws.append(f"measure_union {dis} measurableSet_Ioo")
    A("    rw [" + ",\n        ".join(rws) + "]")
    A(f"  have hbound : {vols}")
    A("      ≤ muGood (1 / 7) (Finset.Icc (1 : ℤ) 76) := by")
    A("    rw [← hvol]; exact measure_mono hu")
    A("  refine le_trans ?_ hbound")
    A("  simp only [" + ", ".join(f"hI{i}" for i in range(1, 25)) + ", Real.volume_Ioo]")
    A("  rw [" + ",\n      ".join("← ENNReal.ofReal_add (by norm_num) (by norm_num)"
                                  for _ in range(23)) + "]")
    A("  exact ENNReal.ofReal_le_ofReal (by norm_num)")
    A("")
    # translate to Icc 0 75 and discharge the Prop
    A("/-- **`AP76Certificate` HOLDS** — the Prop `LRCTailDiameter` consumes, on `{0,…,75}`")
    A("via translation.  With `hlarge_floor_of_diam_le`, the k=13 diameter ≤ 75 leg is now")
    A("unconditional. -/")
    A("theorem ap76_certificate : TailDiameter.AP76Certificate := by")
    A("  unfold TailDiameter.AP76Certificate")
    A("  have hE : (Finset.Icc (1 : ℤ) 76).image (fun e => e - 1) = Finset.Icc (0 : ℤ) 75 := by")
    A("    ext n; simp only [Finset.mem_image, Finset.mem_Icc]")
    A("    constructor")
    A("    · rintro ⟨a, ⟨h1, h2⟩, rfl⟩; omega")
    A("    · intro h; exact ⟨n + 1, ⟨by omega, by omega⟩, by omega⟩")
    A("  have htr : muGood (1 / 7) (Finset.Icc (0 : ℤ) 75)")
    A("      = muGood (1 / 7) (Finset.Icc (1 : ℤ) 76) := by")
    A("    have := muGood_translate (1 / 7) (Finset.Icc (1 : ℤ) 76) 1")
    A("    rw [hE] at this; exact this")
    A("  rw [htr]; exact ap76_sum")
    A("")
    A("/-- The k=13 `hlarge` floor, UNCONDITIONAL: every integer family inside a translated")
    A("window of diameter ≤ 75 has `μ_{1/7} ≥ m_P`. -/")
    A("theorem hlarge_floor_diam75_unconditional")
    A("    (E : Finset ℤ) (m D : ℤ) (hD : D ≤ 75)")
    A("    (hE : ∀ e ∈ E, e - m ∈ Finset.Icc (0 : ℤ) D) :")
    A("    ENNReal.ofReal ((14249 : ℝ) / 252252) ≤ muGood (1 / 7) E :=")
    A("  TailDiameter.hlarge_floor_of_diam_le ap76_certificate E m D hD hE")
    A("")
    A("end AP76Certificate")
    A("end LonelyRunner")
    return "\n".join(L) + "\n"

if __name__ == "__main__":
    src = gen()
    out = "lean/TournamentH7/TournamentH7/LRCAP76Certificate.lean"
    with open(out, "w", encoding="utf-8") as f:
        f.write(src)
    print(f"wrote {out}: {len(src.splitlines())} lines")
