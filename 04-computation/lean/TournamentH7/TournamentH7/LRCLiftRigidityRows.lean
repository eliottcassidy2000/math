/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S77)
-/
import TournamentH7.LRCKernelGate13

/-!
# The hdich lift-rigidity finite leg: all 144 single lifts, kernel-checked strictly loose

Every single lift `r -> r + 13k` (r, k in 1..12) of the full-residue base {1..12} is
STRICTLY loose (M > 1/13), each by kernel `decide` at its exact max-min witness --
standard axioms only.  Floor: 14/169 at liftRow_143 = {1..11,168}, the n=13 deep well
(killer 13^2 - 1, witness 14/169 = 14/13^2; MISTAKE-104).  With residue pinning (S75),
the sieve, and the 13-band window (S76), this closes hdich's finite leg.
-/

namespace LonelyRunner
namespace LiftRigidity

open KernelGate13

theorem liftRow_0 : ∃ t : ℝ, StrictLonely13 (![14,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 16)

theorem liftRow_1 : ∃ t : ℝ, StrictLonely13 (![27,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 2) (den := 29)

theorem liftRow_2 : ∃ t : ℝ, StrictLonely13 (![40,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_3 : ∃ t : ℝ, StrictLonely13 (![53,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_4 : ∃ t : ℝ, StrictLonely13 (![66,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_5 : ∃ t : ℝ, StrictLonely13 (![79,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_6 : ∃ t : ℝ, StrictLonely13 (![92,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_7 : ∃ t : ℝ, StrictLonely13 (![105,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_8 : ∃ t : ℝ, StrictLonely13 (![118,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_9 : ∃ t : ℝ, StrictLonely13 (![131,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_10 : ∃ t : ℝ, StrictLonely13 (![144,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_11 : ∃ t : ℝ, StrictLonely13 (![157,2,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 14)

theorem liftRow_12 : ∃ t : ℝ, StrictLonely13 (![1,15,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 9) (den := 19)

theorem liftRow_13 : ∃ t : ℝ, StrictLonely13 (![1,28,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 15) (den := 32)

theorem liftRow_14 : ∃ t : ℝ, StrictLonely13 (![1,41,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_15 : ∃ t : ℝ, StrictLonely13 (![1,54,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_16 : ∃ t : ℝ, StrictLonely13 (![1,67,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_17 : ∃ t : ℝ, StrictLonely13 (![1,80,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_18 : ∃ t : ℝ, StrictLonely13 (![1,93,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_19 : ∃ t : ℝ, StrictLonely13 (![1,106,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_20 : ∃ t : ℝ, StrictLonely13 (![1,119,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_21 : ∃ t : ℝ, StrictLonely13 (![1,132,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_22 : ∃ t : ℝ, StrictLonely13 (![1,145,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_23 : ∃ t : ℝ, StrictLonely13 (![1,158,3,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 15)

theorem liftRow_24 : ∃ t : ℝ, StrictLonely13 (![1,2,16,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 6) (den := 17)

theorem liftRow_25 : ∃ t : ℝ, StrictLonely13 (![1,2,29,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 6) (den := 17)

theorem liftRow_26 : ∃ t : ℝ, StrictLonely13 (![1,2,42,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_27 : ∃ t : ℝ, StrictLonely13 (![1,2,55,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_28 : ∃ t : ℝ, StrictLonely13 (![1,2,68,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_29 : ∃ t : ℝ, StrictLonely13 (![1,2,81,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_30 : ∃ t : ℝ, StrictLonely13 (![1,2,94,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_31 : ∃ t : ℝ, StrictLonely13 (![1,2,107,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_32 : ∃ t : ℝ, StrictLonely13 (![1,2,120,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_33 : ∃ t : ℝ, StrictLonely13 (![1,2,133,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_34 : ∃ t : ℝ, StrictLonely13 (![1,2,146,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_35 : ∃ t : ℝ, StrictLonely13 (![1,2,159,4,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 16)

theorem liftRow_36 : ∃ t : ℝ, StrictLonely13 (![1,2,3,17,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 19)

theorem liftRow_37 : ∃ t : ℝ, StrictLonely13 (![1,2,3,30,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 5) (den := 19)

theorem liftRow_38 : ∃ t : ℝ, StrictLonely13 (![1,2,3,43,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_39 : ∃ t : ℝ, StrictLonely13 (![1,2,3,56,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_40 : ∃ t : ℝ, StrictLonely13 (![1,2,3,69,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_41 : ∃ t : ℝ, StrictLonely13 (![1,2,3,82,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_42 : ∃ t : ℝ, StrictLonely13 (![1,2,3,95,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_43 : ∃ t : ℝ, StrictLonely13 (![1,2,3,108,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_44 : ∃ t : ℝ, StrictLonely13 (![1,2,3,121,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_45 : ∃ t : ℝ, StrictLonely13 (![1,2,3,134,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_46 : ∃ t : ℝ, StrictLonely13 (![1,2,3,147,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_47 : ∃ t : ℝ, StrictLonely13 (![1,2,3,160,5,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 17)

theorem liftRow_48 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,18,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 19)

theorem liftRow_49 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,31,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 19)

theorem liftRow_50 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,44,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_51 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,57,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_52 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,70,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_53 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,83,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_54 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,96,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_55 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,109,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_56 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,122,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_57 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,135,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_58 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,148,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_59 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,161,6,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 18)

theorem liftRow_60 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,19,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 4) (den := 23)

theorem liftRow_61 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,32,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 7) (den := 44)

theorem liftRow_62 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,45,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_63 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,58,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_64 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,71,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_65 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,84,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_66 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,97,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_67 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,110,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_68 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,123,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_69 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,136,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_70 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,149,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_71 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,162,7,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 3) (den := 19)

theorem liftRow_72 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,20,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_73 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,33,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_74 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,46,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_75 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,59,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_76 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,72,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_77 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,85,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_78 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,98,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 15) (den := 104)

theorem liftRow_79 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,111,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_80 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,124,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_81 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,137,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_82 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,150,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_83 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,163,8,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 7)

theorem liftRow_84 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,21,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_85 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,34,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_86 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,47,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_87 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,60,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_88 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,73,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_89 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,86,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_90 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,99,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_91 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,112,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 44) (den := 117)

theorem liftRow_92 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,125,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_93 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,138,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_94 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,151,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_95 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,164,9,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 8)

theorem liftRow_96 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,22,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_97 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,35,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_98 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,48,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_99 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,61,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_100 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,74,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_101 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,87,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_102 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,100,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_103 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,113,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_104 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,126,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 29) (den := 130)

theorem liftRow_105 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,139,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_106 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,152,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_107 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,165,10,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 9)

theorem liftRow_108 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,23,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_109 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,36,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_110 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,49,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_111 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,62,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_112 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,75,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_113 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,88,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_114 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,101,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_115 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,114,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_116 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,127,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_117 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,140,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 43) (den := 143)

theorem liftRow_118 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,153,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_119 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,166,11,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 10)

theorem liftRow_120 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,24,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_121 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,37,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_122 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,50,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_123 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,63,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_124 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,76,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_125 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,89,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_126 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,102,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_127 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,115,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_128 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,128,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_129 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,141,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_130 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,154,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 71) (den := 156)

theorem liftRow_131 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,167,12] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 11)

theorem liftRow_132 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,25] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_133 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,38] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_134 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,51] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_135 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,64] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_136 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,77] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_137 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,90] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_138 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,103] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_139 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,116] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_140 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,129] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_141 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,142] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_142 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,155] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 1) (den := 12)

theorem liftRow_143 : ∃ t : ℝ, StrictLonely13 (![1,2,3,4,5,6,7,8,9,10,11,168] : Fin 12 → ℤ) t :=
  strictLonely13_of_kernelWitness (by norm_num) (by decide)
    (num := 14) (den := 169)

end LiftRigidity
end LonelyRunner
