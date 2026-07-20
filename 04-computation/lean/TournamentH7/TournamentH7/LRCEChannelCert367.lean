/-
  TournamentH7.LRCEChannelCert367 -- GENERATED (HYP-8010).
  Kernel-exact loneliness value  sSup (margin '' [0,1]) = 4/367  for the
  family {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,91,360} (91 speeds), witness t0 = 53/367.
  F_4(91), THM-1286's third D=4 gate member (discovered S59c).
-/
import TournamentH7.LRCEChannelCert367Checks

namespace LonelyRunner
namespace EChannelCert

open GridAttainment TournamentH7.LRCWitness

set_option maxRecDepth 8000 in
set_option maxHeartbeats 4000000 in
/-- All moduli pass (flat Bool assembly, no case dispatch). -/
theorem moduli367_ok : moduli367.all (certCheckS l367 4 367) = true := by
  simp only [moduli367, List.all_cons, List.all_nil, Bool.and_eq_true]
  exact ⟨chk367_2, chk367_3, chk367_4, chk367_5, chk367_6, chk367_7, chk367_8, chk367_9, chk367_10, chk367_11, chk367_12, chk367_13, chk367_14, chk367_15, chk367_16, chk367_17, chk367_18, chk367_19, chk367_20, chk367_21, chk367_22, chk367_23, chk367_24, chk367_25, chk367_26, chk367_27, chk367_28, chk367_29, chk367_30, chk367_31, chk367_32, chk367_33, chk367_34, chk367_35, chk367_36, chk367_37, chk367_38, chk367_39, chk367_40, chk367_41, chk367_42, chk367_43, chk367_44, chk367_45, chk367_46, chk367_47, chk367_48, chk367_49, chk367_50, chk367_51, chk367_52, chk367_53, chk367_54, chk367_55, chk367_56, chk367_57, chk367_58, chk367_59, chk367_60, chk367_61, chk367_62, chk367_63, chk367_64, chk367_65, chk367_66, chk367_67, chk367_68, chk367_69, chk367_70, chk367_71, chk367_72, chk367_73, chk367_74, chk367_75, chk367_76, chk367_77, chk367_78, chk367_79, chk367_80, chk367_81, chk367_82, chk367_83, chk367_84, chk367_85, chk367_86, chk367_87, chk367_88, chk367_89, chk367_90, chk367_91, chk367_92, chk367_93, chk367_94, chk367_95, chk367_96, chk367_97, chk367_98, chk367_99, chk367_100, chk367_101, chk367_102, chk367_103, chk367_104, chk367_105, chk367_106, chk367_107, chk367_108, chk367_109, chk367_110, chk367_111, chk367_112, chk367_113, chk367_114, chk367_115, chk367_116, chk367_117, chk367_118, chk367_119, chk367_120, chk367_121, chk367_122, chk367_123, chk367_124, chk367_125, chk367_126, chk367_127, chk367_128, chk367_129, chk367_130, chk367_131, chk367_132, chk367_133, chk367_134, chk367_135, chk367_136, chk367_137, chk367_138, chk367_139, chk367_140, chk367_141, chk367_142, chk367_143, chk367_144, chk367_145, chk367_146, chk367_147, chk367_148, chk367_149, chk367_150, chk367_151, chk367_152, chk367_153, chk367_154, chk367_155, chk367_156, chk367_157, chk367_158, chk367_159, chk367_160, chk367_161, chk367_162, chk367_163, chk367_164, chk367_165, chk367_166, chk367_167, chk367_168, chk367_169, chk367_170, chk367_171, chk367_172, chk367_173, chk367_174, chk367_175, chk367_176, chk367_177, chk367_178, chk367_179, chk367_180, chk367_182, chk367_361, chk367_362, chk367_363, chk367_364, chk367_365, chk367_366, chk367_367, chk367_368, chk367_369, chk367_370, chk367_371, chk367_372, chk367_373, chk367_374, chk367_375, chk367_376, chk367_377, chk367_378, chk367_379, chk367_380, chk367_381, chk367_382, chk367_383, chk367_384, chk367_385, chk367_386, chk367_387, chk367_388, chk367_389, chk367_390, chk367_391, chk367_392, chk367_393, chk367_394, chk367_395, chk367_396, chk367_397, chk367_398, chk367_399, chk367_400, chk367_401, chk367_402, chk367_403, chk367_404, chk367_405, chk367_406, chk367_407, chk367_408, chk367_409, chk367_410, chk367_411, chk367_412, chk367_413, chk367_414, chk367_415, chk367_416, chk367_417, chk367_418, chk367_419, chk367_420, chk367_421, chk367_422, chk367_423, chk367_424, chk367_425, chk367_426, chk367_427, chk367_428, chk367_429, chk367_430, chk367_431, chk367_432, chk367_433, chk367_434, chk367_435, chk367_436, chk367_437, chk367_438, chk367_439, chk367_440, chk367_441, chk367_442, chk367_443, chk367_444, chk367_445, chk367_446, chk367_447, chk367_448, chk367_449, chk367_451, chk367_720, trivial⟩

theorem chk367_mem : ∀ S ∈ moduli367, certCheckS l367 4 367 S = true := by
  have h := moduli367_ok
  rw [List.all_eq_true] at h
  exact h

set_option maxRecDepth 8000 in
set_option maxHeartbeats 8976604 in
/-- Every pair sum is a listed modulus (Bool contains sweep). -/
theorem sums_covered367 :
    ∀ vi ∈ l367, ∀ vj ∈ l367, (vi + vj) ∈ moduli367 := by
  have h : (l367.all fun vi => l367.all fun vj =>
      moduli367.contains (vi+vj)) = true := by decide
  simp only [List.all_eq_true] at h
  intro vi hvi vj hvj
  exact List.contains_iff_mem.mp (h vi hvi vj hvj)

theorem check_367 : certCheck l367 4 367 = true := by
  simp only [certCheck, List.all_eq_true]
  intro vi hvi vj hvj
  exact chk367_mem _ (sums_covered367 vi hvi vj hvj)

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2000000 in
/-- The e-channel certificate. -/
theorem cert_367 : Cert v367 4 367 := by
  have h := cert_of_check l367 4 367 check_367 v367
    (by decide) (by decide) (by decide)
  exact_mod_cast h

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2000000 in
/-- **Kernel-exact loneliness value**: `sSup (margin '' [0,1]) = 4/367`. -/
theorem member_367_exact :
    sSup (margin v367 '' Set.Icc (0:ℝ) 1) = (4 : ℝ) / (367 : ℝ) := by
  have hD : ((4:ℤ) : ℝ) = (4:ℝ) := by norm_num
  have hQ : ((367:ℤ) : ℝ) = (367:ℝ) := by norm_num
  rw [← hD, ← hQ]
  apply margin_sSup_eq_of_cert v367 (by decide) 4 367 (by norm_num) (by norm_num)
    ((53 : ℝ) / 367) (by constructor <;> norm_num)
  · intro i m
    have hband : 4 ≤ (v367 i * 53) % 367 ∧
        (v367 i * 53) % 367 ≤ 367 - 4 := by
      fin_cases i <;> decide
    have h := RungFloor.rung_floor_single 367 4 (v367 i) 53 (by norm_num) hband m
    have hc : ((367 : ℤ) : ℝ) = (367 : ℝ) := by norm_num
    rw [hc] at h
    convert h using 2 <;> norm_num
  · exact cert_367

end EChannelCert
end LonelyRunner

#print axioms LonelyRunner.EChannelCert.check_367
#print axioms LonelyRunner.EChannelCert.cert_367
#print axioms LonelyRunner.EChannelCert.member_367_exact
