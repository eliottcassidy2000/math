/-
  TournamentH7.LRCEChannelCert247 -- GENERATED (HYP-8010).
  Kernel-exact loneliness value  sSup (margin '' [0,1]) = 4/247  for the
  family {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,61,240} (61 speeds), witness t0 = 70/247.
  F_4(61), THM-1286's second D=4 gate member (predicted-then-found, S59c).
-/
import TournamentH7.LRCEChannelCert247Checks

namespace LonelyRunner
namespace EChannelCert

open GridAttainment TournamentH7.LRCWitness

set_option maxRecDepth 8000 in
set_option maxHeartbeats 4000000 in
/-- All moduli pass (flat Bool assembly, no case dispatch). -/
theorem moduli247_ok : moduli247.all (certCheckS l247 4 247) = true := by
  simp only [moduli247, List.all_cons, List.all_nil, Bool.and_eq_true]
  exact ⟨chk247_2, chk247_3, chk247_4, chk247_5, chk247_6, chk247_7, chk247_8, chk247_9, chk247_10, chk247_11, chk247_12, chk247_13, chk247_14, chk247_15, chk247_16, chk247_17, chk247_18, chk247_19, chk247_20, chk247_21, chk247_22, chk247_23, chk247_24, chk247_25, chk247_26, chk247_27, chk247_28, chk247_29, chk247_30, chk247_31, chk247_32, chk247_33, chk247_34, chk247_35, chk247_36, chk247_37, chk247_38, chk247_39, chk247_40, chk247_41, chk247_42, chk247_43, chk247_44, chk247_45, chk247_46, chk247_47, chk247_48, chk247_49, chk247_50, chk247_51, chk247_52, chk247_53, chk247_54, chk247_55, chk247_56, chk247_57, chk247_58, chk247_59, chk247_60, chk247_61, chk247_62, chk247_63, chk247_64, chk247_65, chk247_66, chk247_67, chk247_68, chk247_69, chk247_70, chk247_71, chk247_72, chk247_73, chk247_74, chk247_75, chk247_76, chk247_77, chk247_78, chk247_79, chk247_80, chk247_81, chk247_82, chk247_83, chk247_84, chk247_85, chk247_86, chk247_87, chk247_88, chk247_89, chk247_90, chk247_91, chk247_92, chk247_93, chk247_94, chk247_95, chk247_96, chk247_97, chk247_98, chk247_99, chk247_100, chk247_101, chk247_102, chk247_103, chk247_104, chk247_105, chk247_106, chk247_107, chk247_108, chk247_109, chk247_110, chk247_111, chk247_112, chk247_113, chk247_114, chk247_115, chk247_116, chk247_117, chk247_118, chk247_119, chk247_120, chk247_122, chk247_241, chk247_242, chk247_243, chk247_244, chk247_245, chk247_246, chk247_247, chk247_248, chk247_249, chk247_250, chk247_251, chk247_252, chk247_253, chk247_254, chk247_255, chk247_256, chk247_257, chk247_258, chk247_259, chk247_260, chk247_261, chk247_262, chk247_263, chk247_264, chk247_265, chk247_266, chk247_267, chk247_268, chk247_269, chk247_270, chk247_271, chk247_272, chk247_273, chk247_274, chk247_275, chk247_276, chk247_277, chk247_278, chk247_279, chk247_280, chk247_281, chk247_282, chk247_283, chk247_284, chk247_285, chk247_286, chk247_287, chk247_288, chk247_289, chk247_290, chk247_291, chk247_292, chk247_293, chk247_294, chk247_295, chk247_296, chk247_297, chk247_298, chk247_299, chk247_301, chk247_480, trivial⟩

theorem chk247_mem : ∀ S ∈ moduli247, certCheckS l247 4 247 S = true := by
  have h := moduli247_ok
  rw [List.all_eq_true] at h
  exact h

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2694004 in
/-- Every pair sum is a listed modulus (Bool contains sweep). -/
theorem sums_covered247 :
    ∀ vi ∈ l247, ∀ vj ∈ l247, (vi + vj) ∈ moduli247 := by
  have h : (l247.all fun vi => l247.all fun vj =>
      moduli247.contains (vi+vj)) = true := by decide
  simp only [List.all_eq_true] at h
  intro vi hvi vj hvj
  exact List.contains_iff_mem.mp (h vi hvi vj hvj)

theorem check_247 : certCheck l247 4 247 = true := by
  simp only [certCheck, List.all_eq_true]
  intro vi hvi vj hvj
  exact chk247_mem _ (sums_covered247 vi hvi vj hvj)

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2000000 in
/-- The e-channel certificate. -/
theorem cert_247 : Cert v247 4 247 := by
  have h := cert_of_check l247 4 247 check_247 v247
    (by decide) (by decide) (by decide)
  exact_mod_cast h

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2000000 in
/-- **Kernel-exact loneliness value**: `sSup (margin '' [0,1]) = 4/247`. -/
theorem member_247_exact :
    sSup (margin v247 '' Set.Icc (0:ℝ) 1) = (4 : ℝ) / (247 : ℝ) := by
  have hD : ((4:ℤ) : ℝ) = (4:ℝ) := by norm_num
  have hQ : ((247:ℤ) : ℝ) = (247:ℝ) := by norm_num
  rw [← hD, ← hQ]
  apply margin_sSup_eq_of_cert v247 (by decide) 4 247 (by norm_num) (by norm_num)
    ((70 : ℝ) / 247) (by constructor <;> norm_num)
  · intro i m
    have hband : 4 ≤ (v247 i * 70) % 247 ∧
        (v247 i * 70) % 247 ≤ 247 - 4 := by
      fin_cases i <;> decide
    have h := RungFloor.rung_floor_single 247 4 (v247 i) 70 (by norm_num) hband m
    have hc : ((247 : ℤ) : ℝ) = (247 : ℝ) := by norm_num
    rw [hc] at h
    convert h using 2 <;> norm_num
  · exact cert_247

end EChannelCert
end LonelyRunner

#print axioms LonelyRunner.EChannelCert.check_247
#print axioms LonelyRunner.EChannelCert.cert_247
#print axioms LonelyRunner.EChannelCert.member_247_exact
