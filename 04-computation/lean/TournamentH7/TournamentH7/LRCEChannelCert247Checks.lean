/-
  TournamentH7.LRCEChannelCert247Checks -- GENERATED (HYP-8010).
  Defs + 181 per-modulus certificate decides for the 4/247 member
  (family {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,61,240}, worst S = 480).  Stable: build once.
-/
import TournamentH7.LRCEChannelCert

namespace LonelyRunner
namespace EChannelCert

/-- The family as a `Fin 61` function. -/
def v247 : Fin 61 → ℤ := fun i =>
  if (i : ℕ) ≤ 58 then (i : ℕ) + 1 else if (i : ℕ) = 59 then 61 else 240

/-- The value list. -/
def l247 : List ℕ := [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,61,240]

/-- The 181 distinct pair-sum moduli. -/
def moduli247 : List ℕ := [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 122, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 301, 480]

section PerModulusChecks247
set_option maxRecDepth 8000
theorem chk247_2 : certCheckS l247 4 247 2 = true := by decide
theorem chk247_3 : certCheckS l247 4 247 3 = true := by decide
theorem chk247_4 : certCheckS l247 4 247 4 = true := by decide
theorem chk247_5 : certCheckS l247 4 247 5 = true := by decide
theorem chk247_6 : certCheckS l247 4 247 6 = true := by decide
theorem chk247_7 : certCheckS l247 4 247 7 = true := by decide
theorem chk247_8 : certCheckS l247 4 247 8 = true := by decide
theorem chk247_9 : certCheckS l247 4 247 9 = true := by decide
theorem chk247_10 : certCheckS l247 4 247 10 = true := by decide
theorem chk247_11 : certCheckS l247 4 247 11 = true := by decide
theorem chk247_12 : certCheckS l247 4 247 12 = true := by decide
theorem chk247_13 : certCheckS l247 4 247 13 = true := by decide
theorem chk247_14 : certCheckS l247 4 247 14 = true := by decide
theorem chk247_15 : certCheckS l247 4 247 15 = true := by decide
theorem chk247_16 : certCheckS l247 4 247 16 = true := by decide
theorem chk247_17 : certCheckS l247 4 247 17 = true := by decide
theorem chk247_18 : certCheckS l247 4 247 18 = true := by decide
theorem chk247_19 : certCheckS l247 4 247 19 = true := by decide
theorem chk247_20 : certCheckS l247 4 247 20 = true := by decide
theorem chk247_21 : certCheckS l247 4 247 21 = true := by decide
theorem chk247_22 : certCheckS l247 4 247 22 = true := by decide
theorem chk247_23 : certCheckS l247 4 247 23 = true := by decide
theorem chk247_24 : certCheckS l247 4 247 24 = true := by decide
theorem chk247_25 : certCheckS l247 4 247 25 = true := by decide
theorem chk247_26 : certCheckS l247 4 247 26 = true := by decide
theorem chk247_27 : certCheckS l247 4 247 27 = true := by decide
theorem chk247_28 : certCheckS l247 4 247 28 = true := by decide
theorem chk247_29 : certCheckS l247 4 247 29 = true := by decide
theorem chk247_30 : certCheckS l247 4 247 30 = true := by decide
theorem chk247_31 : certCheckS l247 4 247 31 = true := by decide
theorem chk247_32 : certCheckS l247 4 247 32 = true := by decide
theorem chk247_33 : certCheckS l247 4 247 33 = true := by decide
theorem chk247_34 : certCheckS l247 4 247 34 = true := by decide
theorem chk247_35 : certCheckS l247 4 247 35 = true := by decide
theorem chk247_36 : certCheckS l247 4 247 36 = true := by decide
theorem chk247_37 : certCheckS l247 4 247 37 = true := by decide
theorem chk247_38 : certCheckS l247 4 247 38 = true := by decide
theorem chk247_39 : certCheckS l247 4 247 39 = true := by decide
theorem chk247_40 : certCheckS l247 4 247 40 = true := by decide
theorem chk247_41 : certCheckS l247 4 247 41 = true := by decide
theorem chk247_42 : certCheckS l247 4 247 42 = true := by decide
theorem chk247_43 : certCheckS l247 4 247 43 = true := by decide
theorem chk247_44 : certCheckS l247 4 247 44 = true := by decide
theorem chk247_45 : certCheckS l247 4 247 45 = true := by decide
theorem chk247_46 : certCheckS l247 4 247 46 = true := by decide
theorem chk247_47 : certCheckS l247 4 247 47 = true := by decide
theorem chk247_48 : certCheckS l247 4 247 48 = true := by decide
theorem chk247_49 : certCheckS l247 4 247 49 = true := by decide
theorem chk247_50 : certCheckS l247 4 247 50 = true := by decide
theorem chk247_51 : certCheckS l247 4 247 51 = true := by decide
theorem chk247_52 : certCheckS l247 4 247 52 = true := by decide
theorem chk247_53 : certCheckS l247 4 247 53 = true := by decide
theorem chk247_54 : certCheckS l247 4 247 54 = true := by decide
theorem chk247_55 : certCheckS l247 4 247 55 = true := by decide
theorem chk247_56 : certCheckS l247 4 247 56 = true := by decide
theorem chk247_57 : certCheckS l247 4 247 57 = true := by decide
theorem chk247_58 : certCheckS l247 4 247 58 = true := by decide
theorem chk247_59 : certCheckS l247 4 247 59 = true := by decide
theorem chk247_60 : certCheckS l247 4 247 60 = true := by decide
theorem chk247_61 : certCheckS l247 4 247 61 = true := by decide
theorem chk247_62 : certCheckS l247 4 247 62 = true := by decide
theorem chk247_63 : certCheckS l247 4 247 63 = true := by decide
theorem chk247_64 : certCheckS l247 4 247 64 = true := by decide
theorem chk247_65 : certCheckS l247 4 247 65 = true := by decide
theorem chk247_66 : certCheckS l247 4 247 66 = true := by decide
theorem chk247_67 : certCheckS l247 4 247 67 = true := by decide
theorem chk247_68 : certCheckS l247 4 247 68 = true := by decide
theorem chk247_69 : certCheckS l247 4 247 69 = true := by decide
theorem chk247_70 : certCheckS l247 4 247 70 = true := by decide
theorem chk247_71 : certCheckS l247 4 247 71 = true := by decide
theorem chk247_72 : certCheckS l247 4 247 72 = true := by decide
theorem chk247_73 : certCheckS l247 4 247 73 = true := by decide
theorem chk247_74 : certCheckS l247 4 247 74 = true := by decide
theorem chk247_75 : certCheckS l247 4 247 75 = true := by decide
theorem chk247_76 : certCheckS l247 4 247 76 = true := by decide
theorem chk247_77 : certCheckS l247 4 247 77 = true := by decide
theorem chk247_78 : certCheckS l247 4 247 78 = true := by decide
theorem chk247_79 : certCheckS l247 4 247 79 = true := by decide
theorem chk247_80 : certCheckS l247 4 247 80 = true := by decide
theorem chk247_81 : certCheckS l247 4 247 81 = true := by decide
theorem chk247_82 : certCheckS l247 4 247 82 = true := by decide
theorem chk247_83 : certCheckS l247 4 247 83 = true := by decide
theorem chk247_84 : certCheckS l247 4 247 84 = true := by decide
theorem chk247_85 : certCheckS l247 4 247 85 = true := by decide
theorem chk247_86 : certCheckS l247 4 247 86 = true := by decide
theorem chk247_87 : certCheckS l247 4 247 87 = true := by decide
theorem chk247_88 : certCheckS l247 4 247 88 = true := by decide
theorem chk247_89 : certCheckS l247 4 247 89 = true := by decide
theorem chk247_90 : certCheckS l247 4 247 90 = true := by decide
theorem chk247_91 : certCheckS l247 4 247 91 = true := by decide
theorem chk247_92 : certCheckS l247 4 247 92 = true := by decide
theorem chk247_93 : certCheckS l247 4 247 93 = true := by decide
theorem chk247_94 : certCheckS l247 4 247 94 = true := by decide
theorem chk247_95 : certCheckS l247 4 247 95 = true := by decide
theorem chk247_96 : certCheckS l247 4 247 96 = true := by decide
theorem chk247_97 : certCheckS l247 4 247 97 = true := by decide
theorem chk247_98 : certCheckS l247 4 247 98 = true := by decide
theorem chk247_99 : certCheckS l247 4 247 99 = true := by decide
theorem chk247_100 : certCheckS l247 4 247 100 = true := by decide
theorem chk247_101 : certCheckS l247 4 247 101 = true := by decide
theorem chk247_102 : certCheckS l247 4 247 102 = true := by decide
theorem chk247_103 : certCheckS l247 4 247 103 = true := by decide
theorem chk247_104 : certCheckS l247 4 247 104 = true := by decide
theorem chk247_105 : certCheckS l247 4 247 105 = true := by decide
theorem chk247_106 : certCheckS l247 4 247 106 = true := by decide
theorem chk247_107 : certCheckS l247 4 247 107 = true := by decide
theorem chk247_108 : certCheckS l247 4 247 108 = true := by decide
theorem chk247_109 : certCheckS l247 4 247 109 = true := by decide
theorem chk247_110 : certCheckS l247 4 247 110 = true := by decide
theorem chk247_111 : certCheckS l247 4 247 111 = true := by decide
theorem chk247_112 : certCheckS l247 4 247 112 = true := by decide
theorem chk247_113 : certCheckS l247 4 247 113 = true := by decide
theorem chk247_114 : certCheckS l247 4 247 114 = true := by decide
theorem chk247_115 : certCheckS l247 4 247 115 = true := by decide
theorem chk247_116 : certCheckS l247 4 247 116 = true := by decide
theorem chk247_117 : certCheckS l247 4 247 117 = true := by decide
theorem chk247_118 : certCheckS l247 4 247 118 = true := by decide
theorem chk247_119 : certCheckS l247 4 247 119 = true := by decide
theorem chk247_120 : certCheckS l247 4 247 120 = true := by decide
theorem chk247_122 : certCheckS l247 4 247 122 = true := by decide
theorem chk247_241 : certCheckS l247 4 247 241 = true := by decide
theorem chk247_242 : certCheckS l247 4 247 242 = true := by decide
theorem chk247_243 : certCheckS l247 4 247 243 = true := by decide
theorem chk247_244 : certCheckS l247 4 247 244 = true := by decide
theorem chk247_245 : certCheckS l247 4 247 245 = true := by decide
theorem chk247_246 : certCheckS l247 4 247 246 = true := by decide
theorem chk247_247 : certCheckS l247 4 247 247 = true := by decide
theorem chk247_248 : certCheckS l247 4 247 248 = true := by decide
theorem chk247_249 : certCheckS l247 4 247 249 = true := by decide
theorem chk247_250 : certCheckS l247 4 247 250 = true := by decide
theorem chk247_251 : certCheckS l247 4 247 251 = true := by decide
theorem chk247_252 : certCheckS l247 4 247 252 = true := by decide
theorem chk247_253 : certCheckS l247 4 247 253 = true := by decide
theorem chk247_254 : certCheckS l247 4 247 254 = true := by decide
theorem chk247_255 : certCheckS l247 4 247 255 = true := by decide
theorem chk247_256 : certCheckS l247 4 247 256 = true := by decide
theorem chk247_257 : certCheckS l247 4 247 257 = true := by decide
theorem chk247_258 : certCheckS l247 4 247 258 = true := by decide
theorem chk247_259 : certCheckS l247 4 247 259 = true := by decide
theorem chk247_260 : certCheckS l247 4 247 260 = true := by decide
theorem chk247_261 : certCheckS l247 4 247 261 = true := by decide
theorem chk247_262 : certCheckS l247 4 247 262 = true := by decide
theorem chk247_263 : certCheckS l247 4 247 263 = true := by decide
theorem chk247_264 : certCheckS l247 4 247 264 = true := by decide
theorem chk247_265 : certCheckS l247 4 247 265 = true := by decide
theorem chk247_266 : certCheckS l247 4 247 266 = true := by decide
theorem chk247_267 : certCheckS l247 4 247 267 = true := by decide
theorem chk247_268 : certCheckS l247 4 247 268 = true := by decide
theorem chk247_269 : certCheckS l247 4 247 269 = true := by decide
theorem chk247_270 : certCheckS l247 4 247 270 = true := by decide
theorem chk247_271 : certCheckS l247 4 247 271 = true := by decide
theorem chk247_272 : certCheckS l247 4 247 272 = true := by decide
theorem chk247_273 : certCheckS l247 4 247 273 = true := by decide
theorem chk247_274 : certCheckS l247 4 247 274 = true := by decide
theorem chk247_275 : certCheckS l247 4 247 275 = true := by decide
theorem chk247_276 : certCheckS l247 4 247 276 = true := by decide
theorem chk247_277 : certCheckS l247 4 247 277 = true := by decide
theorem chk247_278 : certCheckS l247 4 247 278 = true := by decide
theorem chk247_279 : certCheckS l247 4 247 279 = true := by decide
theorem chk247_280 : certCheckS l247 4 247 280 = true := by decide
theorem chk247_281 : certCheckS l247 4 247 281 = true := by decide
theorem chk247_282 : certCheckS l247 4 247 282 = true := by decide
theorem chk247_283 : certCheckS l247 4 247 283 = true := by decide
theorem chk247_284 : certCheckS l247 4 247 284 = true := by decide
theorem chk247_285 : certCheckS l247 4 247 285 = true := by decide
theorem chk247_286 : certCheckS l247 4 247 286 = true := by decide
theorem chk247_287 : certCheckS l247 4 247 287 = true := by decide
theorem chk247_288 : certCheckS l247 4 247 288 = true := by decide
theorem chk247_289 : certCheckS l247 4 247 289 = true := by decide
theorem chk247_290 : certCheckS l247 4 247 290 = true := by decide
theorem chk247_291 : certCheckS l247 4 247 291 = true := by decide
theorem chk247_292 : certCheckS l247 4 247 292 = true := by decide
theorem chk247_293 : certCheckS l247 4 247 293 = true := by decide
theorem chk247_294 : certCheckS l247 4 247 294 = true := by decide
theorem chk247_295 : certCheckS l247 4 247 295 = true := by decide
theorem chk247_296 : certCheckS l247 4 247 296 = true := by decide
theorem chk247_297 : certCheckS l247 4 247 297 = true := by decide
theorem chk247_298 : certCheckS l247 4 247 298 = true := by decide
theorem chk247_299 : certCheckS l247 4 247 299 = true := by decide
theorem chk247_301 : certCheckS l247 4 247 301 = true := by decide
theorem chk247_480 : certCheckS l247 4 247 480 = true := by decide
end PerModulusChecks247

end EChannelCert
end LonelyRunner
