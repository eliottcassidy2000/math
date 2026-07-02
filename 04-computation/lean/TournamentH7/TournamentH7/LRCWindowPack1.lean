/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S47 kernel edition)
-/
import TournamentH7.LRCKernelGate

/-!
# Bounded-window census band (0, 18] -- kernel-gate edition [WIP: DO NOT root-register]

All 966 primitive covering 13-tuples with max <= 18, witnesses exact, checks TRUE
(python-mirrored) and equivalent to the generator's rational verification.

**CRASH FORENSICS (S47, unresolved -- fresh-context debug needed):** lean.exe dies with
0xC0000005 (access violation) elaborating this file -- reproduced at 3 rows with purged
artifacts -- while the IDENTICAL row form (kernelRow_AP/kernelRow_2) compiles green inside
LRCKernelGate.lean.  Facts: crash also occurred for the native_decide edition at 966 rows
(originally attributed to parallel load); a corrupt-olean cascade from the first crash
poisoned the ROOT build (fixed by purging artifacts -- purge before debugging!).  Next
probes: (1) copy winRow_0 verbatim INTO LRCKernelGate.lean (isolates file-vs-form);
(2) binary-diff this file's bytes vs the gate file around the theorem syntax;
(3) try `decide +kernel` / set_option maxRecDepth; (4) check Windows stack-size lore.
-/

namespace LonelyRunner
namespace WindowPack

theorem winRow_0 : Lonely 14 (![1,2,3,4,5,6,8,9,10,11,12,13,14] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_1 : Lonely 14 (![1,2,3,4,5,6,8,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_2 : Lonely 14 (![1,2,3,4,5,6,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_3 : Lonely 14 (![1,2,3,4,5,6,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_4 : Lonely 14 (![1,2,3,4,5,7,8,9,10,11,12,13,14] : Fin 13 → ℤ)
    ((((4 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_5 : Lonely 14 (![1,2,3,4,5,7,8,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_6 : Lonely 14 (![1,2,3,4,5,7,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((3 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_7 : Lonely 14 (![1,2,3,4,5,7,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_8 : Lonely 14 (![1,2,3,4,5,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((4 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_9 : Lonely 14 (![1,2,3,4,5,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_10 : Lonely 14 (![1,2,3,4,5,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_11 : Lonely 14 (![1,2,3,4,5,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_12 : Lonely 14 (![1,2,3,4,5,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_13 : Lonely 14 (![1,2,3,4,5,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_14 : Lonely 14 (![1,2,3,4,5,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_15 : Lonely 14 (![1,2,3,4,5,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((4 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_16 : Lonely 14 (![1,2,3,4,5,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_17 : Lonely 14 (![1,2,3,4,5,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_18 : Lonely 14 (![1,2,3,4,5,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_19 : Lonely 14 (![1,2,3,4,5,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_20 : Lonely 14 (![1,2,3,4,6,7,8,9,10,11,12,13,14] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_21 : Lonely 14 (![1,2,3,4,6,7,8,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_22 : Lonely 14 (![1,2,3,4,6,7,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_23 : Lonely 14 (![1,2,3,4,6,7,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_24 : Lonely 14 (![1,2,3,4,6,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_25 : Lonely 14 (![1,2,3,4,6,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_26 : Lonely 14 (![1,2,3,4,6,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_27 : Lonely 14 (![1,2,3,4,6,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_28 : Lonely 14 (![1,2,3,4,6,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_29 : Lonely 14 (![1,2,3,4,6,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_30 : Lonely 14 (![1,2,3,4,6,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_31 : Lonely 14 (![1,2,3,4,6,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_32 : Lonely 14 (![1,2,3,4,6,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_33 : Lonely 14 (![1,2,3,4,6,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_34 : Lonely 14 (![1,2,3,4,6,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_35 : Lonely 14 (![1,2,3,4,6,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_36 : Lonely 14 (![1,2,3,4,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_37 : Lonely 14 (![1,2,3,4,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_38 : Lonely 14 (![1,2,3,4,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_39 : Lonely 14 (![1,2,3,4,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_40 : Lonely 14 (![1,2,3,4,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_41 : Lonely 14 (![1,2,3,4,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_42 : Lonely 14 (![1,2,3,4,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_43 : Lonely 14 (![1,2,3,4,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_44 : Lonely 14 (![1,2,3,4,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((3 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_45 : Lonely 14 (![1,2,3,4,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_46 : Lonely 14 (![1,2,3,4,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_47 : Lonely 14 (![1,2,3,4,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((3 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_48 : Lonely 14 (![1,2,3,4,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_49 : Lonely 14 (![1,2,3,4,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_50 : Lonely 14 (![1,2,3,4,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_51 : Lonely 14 (![1,2,3,4,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_52 : Lonely 14 (![1,2,3,4,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_53 : Lonely 14 (![1,2,3,4,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_54 : Lonely 14 (![1,2,3,4,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_55 : Lonely 14 (![1,2,3,4,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_56 : Lonely 14 (![1,2,3,4,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_57 : Lonely 14 (![1,2,3,4,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((9 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_58 : Lonely 14 (![1,2,3,4,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_59 : Lonely 14 (![1,2,3,4,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_60 : Lonely 14 (![1,2,3,4,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (24 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_61 : Lonely 14 (![1,2,3,5,6,7,8,9,10,11,12,13,14] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_62 : Lonely 14 (![1,2,3,5,6,7,8,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_63 : Lonely 14 (![1,2,3,5,6,7,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_64 : Lonely 14 (![1,2,3,5,6,7,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_65 : Lonely 14 (![1,2,3,5,6,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_66 : Lonely 14 (![1,2,3,5,6,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_67 : Lonely 14 (![1,2,3,5,6,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_68 : Lonely 14 (![1,2,3,5,6,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_69 : Lonely 14 (![1,2,3,5,6,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_70 : Lonely 14 (![1,2,3,5,6,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_71 : Lonely 14 (![1,2,3,5,6,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_72 : Lonely 14 (![1,2,3,5,6,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_73 : Lonely 14 (![1,2,3,5,6,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_74 : Lonely 14 (![1,2,3,5,6,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_75 : Lonely 14 (![1,2,3,5,6,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_76 : Lonely 14 (![1,2,3,5,6,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_77 : Lonely 14 (![1,2,3,5,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_78 : Lonely 14 (![1,2,3,5,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_79 : Lonely 14 (![1,2,3,5,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_80 : Lonely 14 (![1,2,3,5,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_81 : Lonely 14 (![1,2,3,5,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_82 : Lonely 14 (![1,2,3,5,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_83 : Lonely 14 (![1,2,3,5,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_84 : Lonely 14 (![1,2,3,5,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_85 : Lonely 14 (![1,2,3,5,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_86 : Lonely 14 (![1,2,3,5,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_87 : Lonely 14 (![1,2,3,5,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_88 : Lonely 14 (![1,2,3,5,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_89 : Lonely 14 (![1,2,3,5,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_90 : Lonely 14 (![1,2,3,5,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_91 : Lonely 14 (![1,2,3,5,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_92 : Lonely 14 (![1,2,3,5,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_93 : Lonely 14 (![1,2,3,5,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_94 : Lonely 14 (![1,2,3,5,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_95 : Lonely 14 (![1,2,3,5,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_96 : Lonely 14 (![1,2,3,5,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_97 : Lonely 14 (![1,2,3,5,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_98 : Lonely 14 (![1,2,3,5,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_99 : Lonely 14 (![1,2,3,5,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_100 : Lonely 14 (![1,2,3,5,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_101 : Lonely 14 (![1,2,3,5,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_102 : Lonely 14 (![1,2,3,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_103 : Lonely 14 (![1,2,3,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_104 : Lonely 14 (![1,2,3,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_105 : Lonely 14 (![1,2,3,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_106 : Lonely 14 (![1,2,3,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_107 : Lonely 14 (![1,2,3,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_108 : Lonely 14 (![1,2,3,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_109 : Lonely 14 (![1,2,3,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_110 : Lonely 14 (![1,2,3,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_111 : Lonely 14 (![1,2,3,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_112 : Lonely 14 (![1,2,3,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_113 : Lonely 14 (![1,2,3,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_114 : Lonely 14 (![1,2,3,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_115 : Lonely 14 (![1,2,3,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_116 : Lonely 14 (![1,2,3,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_117 : Lonely 14 (![1,2,3,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_118 : Lonely 14 (![1,2,3,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_119 : Lonely 14 (![1,2,3,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_120 : Lonely 14 (![1,2,3,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_121 : Lonely 14 (![1,2,3,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_122 : Lonely 14 (![1,2,3,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_123 : Lonely 14 (![1,2,3,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_124 : Lonely 14 (![1,2,3,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_125 : Lonely 14 (![1,2,3,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_126 : Lonely 14 (![1,2,3,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_127 : Lonely 14 (![1,2,3,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_128 : Lonely 14 (![1,2,3,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_129 : Lonely 14 (![1,2,3,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_130 : Lonely 14 (![1,2,3,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_131 : Lonely 14 (![1,2,3,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_132 : Lonely 14 (![1,2,3,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_133 : Lonely 14 (![1,2,3,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_134 : Lonely 14 (![1,2,3,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_135 : Lonely 14 (![1,2,3,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_136 : Lonely 14 (![1,2,3,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_137 : Lonely 14 (![1,2,3,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_138 : Lonely 14 (![1,2,3,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_139 : Lonely 14 (![1,2,3,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_140 : Lonely 14 (![1,2,3,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_141 : Lonely 14 (![1,2,3,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_142 : Lonely 14 (![1,2,3,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_143 : Lonely 14 (![1,2,3,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_144 : Lonely 14 (![1,2,3,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_145 : Lonely 14 (![1,2,3,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_146 : Lonely 14 (![1,2,4,5,6,7,8,9,10,11,12,13,14] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_147 : Lonely 14 (![1,2,4,5,6,7,8,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_148 : Lonely 14 (![1,2,4,5,6,7,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_149 : Lonely 14 (![1,2,4,5,6,7,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_150 : Lonely 14 (![1,2,4,5,6,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_151 : Lonely 14 (![1,2,4,5,6,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_152 : Lonely 14 (![1,2,4,5,6,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_153 : Lonely 14 (![1,2,4,5,6,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_154 : Lonely 14 (![1,2,4,5,6,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_155 : Lonely 14 (![1,2,4,5,6,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_156 : Lonely 14 (![1,2,4,5,6,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_157 : Lonely 14 (![1,2,4,5,6,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_158 : Lonely 14 (![1,2,4,5,6,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_159 : Lonely 14 (![1,2,4,5,6,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_160 : Lonely 14 (![1,2,4,5,6,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_161 : Lonely 14 (![1,2,4,5,6,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_162 : Lonely 14 (![1,2,4,5,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_163 : Lonely 14 (![1,2,4,5,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_164 : Lonely 14 (![1,2,4,5,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_165 : Lonely 14 (![1,2,4,5,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_166 : Lonely 14 (![1,2,4,5,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_167 : Lonely 14 (![1,2,4,5,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (26 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_168 : Lonely 14 (![1,2,4,5,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_169 : Lonely 14 (![1,2,4,5,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_170 : Lonely 14 (![1,2,4,5,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_171 : Lonely 14 (![1,2,4,5,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_172 : Lonely 14 (![1,2,4,5,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (26 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_173 : Lonely 14 (![1,2,4,5,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_174 : Lonely 14 (![1,2,4,5,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_175 : Lonely 14 (![1,2,4,5,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_176 : Lonely 14 (![1,2,4,5,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_177 : Lonely 14 (![1,2,4,5,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_178 : Lonely 14 (![1,2,4,5,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_179 : Lonely 14 (![1,2,4,5,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_180 : Lonely 14 (![1,2,4,5,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (26 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_181 : Lonely 14 (![1,2,4,5,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_182 : Lonely 14 (![1,2,4,5,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_183 : Lonely 14 (![1,2,4,5,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_184 : Lonely 14 (![1,2,4,5,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_185 : Lonely 14 (![1,2,4,5,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_186 : Lonely 14 (![1,2,4,5,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_187 : Lonely 14 (![1,2,4,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_188 : Lonely 14 (![1,2,4,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_189 : Lonely 14 (![1,2,4,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_190 : Lonely 14 (![1,2,4,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_191 : Lonely 14 (![1,2,4,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_192 : Lonely 14 (![1,2,4,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_193 : Lonely 14 (![1,2,4,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_194 : Lonely 14 (![1,2,4,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_195 : Lonely 14 (![1,2,4,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_196 : Lonely 14 (![1,2,4,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_197 : Lonely 14 (![1,2,4,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_198 : Lonely 14 (![1,2,4,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_199 : Lonely 14 (![1,2,4,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_200 : Lonely 14 (![1,2,4,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_201 : Lonely 14 (![1,2,4,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_202 : Lonely 14 (![1,2,4,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_203 : Lonely 14 (![1,2,4,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_204 : Lonely 14 (![1,2,4,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_205 : Lonely 14 (![1,2,4,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_206 : Lonely 14 (![1,2,4,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_207 : Lonely 14 (![1,2,4,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((3 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_208 : Lonely 14 (![1,2,4,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_209 : Lonely 14 (![1,2,4,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_210 : Lonely 14 (![1,2,4,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_211 : Lonely 14 (![1,2,4,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_212 : Lonely 14 (![1,2,4,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_213 : Lonely 14 (![1,2,4,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_214 : Lonely 14 (![1,2,4,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_215 : Lonely 14 (![1,2,4,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_216 : Lonely 14 (![1,2,4,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_217 : Lonely 14 (![1,2,4,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_218 : Lonely 14 (![1,2,4,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (26 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_219 : Lonely 14 (![1,2,4,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_220 : Lonely 14 (![1,2,4,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_221 : Lonely 14 (![1,2,4,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_222 : Lonely 14 (![1,2,4,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_223 : Lonely 14 (![1,2,4,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_224 : Lonely 14 (![1,2,4,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_225 : Lonely 14 (![1,2,4,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_226 : Lonely 14 (![1,2,4,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_227 : Lonely 14 (![1,2,4,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_228 : Lonely 14 (![1,2,4,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_229 : Lonely 14 (![1,2,4,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_230 : Lonely 14 (![1,2,4,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_231 : Lonely 14 (![1,2,5,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_232 : Lonely 14 (![1,2,5,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_233 : Lonely 14 (![1,2,5,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_234 : Lonely 14 (![1,2,5,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_235 : Lonely 14 (![1,2,5,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_236 : Lonely 14 (![1,2,5,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_237 : Lonely 14 (![1,2,5,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_238 : Lonely 14 (![1,2,5,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_239 : Lonely 14 (![1,2,5,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_240 : Lonely 14 (![1,2,5,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_241 : Lonely 14 (![1,2,5,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_242 : Lonely 14 (![1,2,5,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_243 : Lonely 14 (![1,2,5,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_244 : Lonely 14 (![1,2,5,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_245 : Lonely 14 (![1,2,5,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_246 : Lonely 14 (![1,2,5,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_247 : Lonely 14 (![1,2,5,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_248 : Lonely 14 (![1,2,5,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_249 : Lonely 14 (![1,2,5,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_250 : Lonely 14 (![1,2,5,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_251 : Lonely 14 (![1,2,5,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_252 : Lonely 14 (![1,2,5,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_253 : Lonely 14 (![1,2,5,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_254 : Lonely 14 (![1,2,5,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_255 : Lonely 14 (![1,2,5,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_256 : Lonely 14 (![1,2,5,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_257 : Lonely 14 (![1,2,5,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_258 : Lonely 14 (![1,2,5,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_259 : Lonely 14 (![1,2,5,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_260 : Lonely 14 (![1,2,5,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_261 : Lonely 14 (![1,2,5,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_262 : Lonely 14 (![1,2,5,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (26 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_263 : Lonely 14 (![1,2,5,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_264 : Lonely 14 (![1,2,5,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_265 : Lonely 14 (![1,2,5,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_266 : Lonely 14 (![1,2,5,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_267 : Lonely 14 (![1,2,5,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_268 : Lonely 14 (![1,2,5,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_269 : Lonely 14 (![1,2,5,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_270 : Lonely 14 (![1,2,5,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_271 : Lonely 14 (![1,2,5,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_272 : Lonely 14 (![1,2,5,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_273 : Lonely 14 (![1,2,5,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_274 : Lonely 14 (![1,2,5,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_275 : Lonely 14 (![1,2,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_276 : Lonely 14 (![1,2,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_277 : Lonely 14 (![1,2,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_278 : Lonely 14 (![1,2,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_279 : Lonely 14 (![1,2,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_280 : Lonely 14 (![1,2,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_281 : Lonely 14 (![1,2,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_282 : Lonely 14 (![1,2,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_283 : Lonely 14 (![1,2,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_284 : Lonely 14 (![1,2,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_285 : Lonely 14 (![1,2,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_286 : Lonely 14 (![1,2,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_287 : Lonely 14 (![1,2,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_288 : Lonely 14 (![1,2,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_289 : Lonely 14 (![1,2,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_290 : Lonely 14 (![1,2,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_291 : Lonely 14 (![1,2,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_292 : Lonely 14 (![1,2,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_293 : Lonely 14 (![1,2,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_294 : Lonely 14 (![1,2,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_295 : Lonely 14 (![1,2,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_296 : Lonely 14 (![1,2,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_297 : Lonely 14 (![1,2,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_298 : Lonely 14 (![1,2,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_299 : Lonely 14 (![1,2,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_300 : Lonely 14 (![1,2,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_301 : Lonely 14 (![1,3,4,5,6,7,8,9,10,11,12,13,14] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_302 : Lonely 14 (![1,3,4,5,6,7,8,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_303 : Lonely 14 (![1,3,4,5,6,7,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_304 : Lonely 14 (![1,3,4,5,6,7,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_305 : Lonely 14 (![1,3,4,5,6,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_306 : Lonely 14 (![1,3,4,5,6,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_307 : Lonely 14 (![1,3,4,5,6,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_308 : Lonely 14 (![1,3,4,5,6,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_309 : Lonely 14 (![1,3,4,5,6,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_310 : Lonely 14 (![1,3,4,5,6,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_311 : Lonely 14 (![1,3,4,5,6,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_312 : Lonely 14 (![1,3,4,5,6,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_313 : Lonely 14 (![1,3,4,5,6,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_314 : Lonely 14 (![1,3,4,5,6,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_315 : Lonely 14 (![1,3,4,5,6,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_316 : Lonely 14 (![1,3,4,5,6,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_317 : Lonely 14 (![1,3,4,5,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_318 : Lonely 14 (![1,3,4,5,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_319 : Lonely 14 (![1,3,4,5,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_320 : Lonely 14 (![1,3,4,5,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_321 : Lonely 14 (![1,3,4,5,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_322 : Lonely 14 (![1,3,4,5,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_323 : Lonely 14 (![1,3,4,5,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_324 : Lonely 14 (![1,3,4,5,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_325 : Lonely 14 (![1,3,4,5,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_326 : Lonely 14 (![1,3,4,5,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_327 : Lonely 14 (![1,3,4,5,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_328 : Lonely 14 (![1,3,4,5,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_329 : Lonely 14 (![1,3,4,5,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_330 : Lonely 14 (![1,3,4,5,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_331 : Lonely 14 (![1,3,4,5,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_332 : Lonely 14 (![1,3,4,5,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_333 : Lonely 14 (![1,3,4,5,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_334 : Lonely 14 (![1,3,4,5,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_335 : Lonely 14 (![1,3,4,5,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_336 : Lonely 14 (![1,3,4,5,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_337 : Lonely 14 (![1,3,4,5,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_338 : Lonely 14 (![1,3,4,5,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_339 : Lonely 14 (![1,3,4,5,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_340 : Lonely 14 (![1,3,4,5,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_341 : Lonely 14 (![1,3,4,5,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_342 : Lonely 14 (![1,3,4,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_343 : Lonely 14 (![1,3,4,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_344 : Lonely 14 (![1,3,4,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_345 : Lonely 14 (![1,3,4,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_346 : Lonely 14 (![1,3,4,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_347 : Lonely 14 (![1,3,4,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_348 : Lonely 14 (![1,3,4,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_349 : Lonely 14 (![1,3,4,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_350 : Lonely 14 (![1,3,4,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_351 : Lonely 14 (![1,3,4,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_352 : Lonely 14 (![1,3,4,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_353 : Lonely 14 (![1,3,4,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_354 : Lonely 14 (![1,3,4,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_355 : Lonely 14 (![1,3,4,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_356 : Lonely 14 (![1,3,4,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_357 : Lonely 14 (![1,3,4,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_358 : Lonely 14 (![1,3,4,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_359 : Lonely 14 (![1,3,4,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_360 : Lonely 14 (![1,3,4,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_361 : Lonely 14 (![1,3,4,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_362 : Lonely 14 (![1,3,4,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_363 : Lonely 14 (![1,3,4,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_364 : Lonely 14 (![1,3,4,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_365 : Lonely 14 (![1,3,4,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_366 : Lonely 14 (![1,3,4,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_367 : Lonely 14 (![1,3,4,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_368 : Lonely 14 (![1,3,4,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_369 : Lonely 14 (![1,3,4,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_370 : Lonely 14 (![1,3,4,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_371 : Lonely 14 (![1,3,4,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_372 : Lonely 14 (![1,3,4,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_373 : Lonely 14 (![1,3,4,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_374 : Lonely 14 (![1,3,4,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_375 : Lonely 14 (![1,3,4,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_376 : Lonely 14 (![1,3,4,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_377 : Lonely 14 (![1,3,4,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_378 : Lonely 14 (![1,3,4,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_379 : Lonely 14 (![1,3,4,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_380 : Lonely 14 (![1,3,4,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_381 : Lonely 14 (![1,3,4,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_382 : Lonely 14 (![1,3,4,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((4 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_383 : Lonely 14 (![1,3,4,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_384 : Lonely 14 (![1,3,4,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_385 : Lonely 14 (![1,3,4,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_386 : Lonely 14 (![1,3,5,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_387 : Lonely 14 (![1,3,5,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_388 : Lonely 14 (![1,3,5,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_389 : Lonely 14 (![1,3,5,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_390 : Lonely 14 (![1,3,5,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_391 : Lonely 14 (![1,3,5,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_392 : Lonely 14 (![1,3,5,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_393 : Lonely 14 (![1,3,5,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_394 : Lonely 14 (![1,3,5,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_395 : Lonely 14 (![1,3,5,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_396 : Lonely 14 (![1,3,5,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_397 : Lonely 14 (![1,3,5,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_398 : Lonely 14 (![1,3,5,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_399 : Lonely 14 (![1,3,5,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_400 : Lonely 14 (![1,3,5,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_401 : Lonely 14 (![1,3,5,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_402 : Lonely 14 (![1,3,5,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_403 : Lonely 14 (![1,3,5,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_404 : Lonely 14 (![1,3,5,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_405 : Lonely 14 (![1,3,5,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_406 : Lonely 14 (![1,3,5,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_407 : Lonely 14 (![1,3,5,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_408 : Lonely 14 (![1,3,5,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_409 : Lonely 14 (![1,3,5,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_410 : Lonely 14 (![1,3,5,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_411 : Lonely 14 (![1,3,5,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_412 : Lonely 14 (![1,3,5,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_413 : Lonely 14 (![1,3,5,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_414 : Lonely 14 (![1,3,5,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_415 : Lonely 14 (![1,3,5,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_416 : Lonely 14 (![1,3,5,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_417 : Lonely 14 (![1,3,5,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_418 : Lonely 14 (![1,3,5,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_419 : Lonely 14 (![1,3,5,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_420 : Lonely 14 (![1,3,5,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((13 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_421 : Lonely 14 (![1,3,5,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_422 : Lonely 14 (![1,3,5,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((13 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_423 : Lonely 14 (![1,3,5,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((13 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_424 : Lonely 14 (![1,3,5,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_425 : Lonely 14 (![1,3,5,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_426 : Lonely 14 (![1,3,5,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_427 : Lonely 14 (![1,3,5,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_428 : Lonely 14 (![1,3,5,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_429 : Lonely 14 (![1,3,5,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((13 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_430 : Lonely 14 (![1,3,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_431 : Lonely 14 (![1,3,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_432 : Lonely 14 (![1,3,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_433 : Lonely 14 (![1,3,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_434 : Lonely 14 (![1,3,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_435 : Lonely 14 (![1,3,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_436 : Lonely 14 (![1,3,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_437 : Lonely 14 (![1,3,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_438 : Lonely 14 (![1,3,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_439 : Lonely 14 (![1,3,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_440 : Lonely 14 (![1,3,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_441 : Lonely 14 (![1,3,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_442 : Lonely 14 (![1,3,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_443 : Lonely 14 (![1,3,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_444 : Lonely 14 (![1,3,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_445 : Lonely 14 (![1,3,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_446 : Lonely 14 (![1,3,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_447 : Lonely 14 (![1,3,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_448 : Lonely 14 (![1,3,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_449 : Lonely 14 (![1,3,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_450 : Lonely 14 (![1,3,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_451 : Lonely 14 (![1,3,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_452 : Lonely 14 (![1,3,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_453 : Lonely 14 (![1,3,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_454 : Lonely 14 (![1,3,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((13 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_455 : Lonely 14 (![1,3,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_456 : Lonely 14 (![1,4,5,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_457 : Lonely 14 (![1,4,5,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_458 : Lonely 14 (![1,4,5,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_459 : Lonely 14 (![1,4,5,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_460 : Lonely 14 (![1,4,5,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_461 : Lonely 14 (![1,4,5,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_462 : Lonely 14 (![1,4,5,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_463 : Lonely 14 (![1,4,5,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_464 : Lonely 14 (![1,4,5,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_465 : Lonely 14 (![1,4,5,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_466 : Lonely 14 (![1,4,5,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_467 : Lonely 14 (![1,4,5,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_468 : Lonely 14 (![1,4,5,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_469 : Lonely 14 (![1,4,5,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_470 : Lonely 14 (![1,4,5,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_471 : Lonely 14 (![1,4,5,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_472 : Lonely 14 (![1,4,5,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_473 : Lonely 14 (![1,4,5,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_474 : Lonely 14 (![1,4,5,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_475 : Lonely 14 (![1,4,5,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_476 : Lonely 14 (![1,4,5,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_477 : Lonely 14 (![1,4,5,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_478 : Lonely 14 (![1,4,5,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_479 : Lonely 14 (![1,4,5,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_480 : Lonely 14 (![1,4,5,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_481 : Lonely 14 (![1,4,5,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_482 : Lonely 14 (![1,4,5,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_483 : Lonely 14 (![1,4,5,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_484 : Lonely 14 (![1,4,5,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_485 : Lonely 14 (![1,4,5,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_486 : Lonely 14 (![1,4,5,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_487 : Lonely 14 (![1,4,5,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (26 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_488 : Lonely 14 (![1,4,5,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_489 : Lonely 14 (![1,4,5,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_490 : Lonely 14 (![1,4,5,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_491 : Lonely 14 (![1,4,5,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_492 : Lonely 14 (![1,4,5,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_493 : Lonely 14 (![1,4,5,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_494 : Lonely 14 (![1,4,5,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_495 : Lonely 14 (![1,4,5,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_496 : Lonely 14 (![1,4,5,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_497 : Lonely 14 (![1,4,5,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_498 : Lonely 14 (![1,4,5,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_499 : Lonely 14 (![1,4,5,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_500 : Lonely 14 (![1,4,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_501 : Lonely 14 (![1,4,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_502 : Lonely 14 (![1,4,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_503 : Lonely 14 (![1,4,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_504 : Lonely 14 (![1,4,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_505 : Lonely 14 (![1,4,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_506 : Lonely 14 (![1,4,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_507 : Lonely 14 (![1,4,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_508 : Lonely 14 (![1,4,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_509 : Lonely 14 (![1,4,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_510 : Lonely 14 (![1,4,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_511 : Lonely 14 (![1,4,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_512 : Lonely 14 (![1,4,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_513 : Lonely 14 (![1,4,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_514 : Lonely 14 (![1,4,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_515 : Lonely 14 (![1,4,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_516 : Lonely 14 (![1,4,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_517 : Lonely 14 (![1,4,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_518 : Lonely 14 (![1,4,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_519 : Lonely 14 (![1,4,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_520 : Lonely 14 (![1,4,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_521 : Lonely 14 (![1,4,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_522 : Lonely 14 (![1,4,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_523 : Lonely 14 (![1,4,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_524 : Lonely 14 (![1,4,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_525 : Lonely 14 (![1,4,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_526 : Lonely 14 (![1,5,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_527 : Lonely 14 (![1,5,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_528 : Lonely 14 (![1,5,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_529 : Lonely 14 (![1,5,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_530 : Lonely 14 (![1,5,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_531 : Lonely 14 (![1,5,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_532 : Lonely 14 (![1,5,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_533 : Lonely 14 (![1,5,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_534 : Lonely 14 (![1,5,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_535 : Lonely 14 (![1,5,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_536 : Lonely 14 (![1,5,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_537 : Lonely 14 (![1,5,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_538 : Lonely 14 (![1,5,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_539 : Lonely 14 (![1,5,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_540 : Lonely 14 (![1,5,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_541 : Lonely 14 (![1,5,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_542 : Lonely 14 (![1,5,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_543 : Lonely 14 (![1,5,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_544 : Lonely 14 (![1,5,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_545 : Lonely 14 (![1,5,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_546 : Lonely 14 (![1,5,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_547 : Lonely 14 (![1,5,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_548 : Lonely 14 (![1,5,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_549 : Lonely 14 (![1,5,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_550 : Lonely 14 (![1,5,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((13 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_551 : Lonely 14 (![1,5,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_552 : Lonely 14 (![1,6,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_553 : Lonely 14 (![1,6,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((10 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_554 : Lonely 14 (![1,6,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_555 : Lonely 14 (![1,6,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_556 : Lonely 14 (![1,6,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_557 : Lonely 14 (![1,6,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_558 : Lonely 14 (![1,6,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_559 : Lonely 14 (![1,7,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_560 : Lonely 14 (![2,3,4,5,6,7,8,9,10,11,12,13,14] : Fin 13 → ℤ)
    ((((1 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_561 : Lonely 14 (![2,3,4,5,6,7,8,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_562 : Lonely 14 (![2,3,4,5,6,7,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_563 : Lonely 14 (![2,3,4,5,6,7,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_564 : Lonely 14 (![2,3,4,5,6,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((1 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_565 : Lonely 14 (![2,3,4,5,6,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_566 : Lonely 14 (![2,3,4,5,6,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_567 : Lonely 14 (![2,3,4,5,6,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_568 : Lonely 14 (![2,3,4,5,6,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_569 : Lonely 14 (![2,3,4,5,6,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_570 : Lonely 14 (![2,3,4,5,6,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_571 : Lonely 14 (![2,3,4,5,6,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_572 : Lonely 14 (![2,3,4,5,6,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_573 : Lonely 14 (![2,3,4,5,6,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_574 : Lonely 14 (![2,3,4,5,6,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_575 : Lonely 14 (![2,3,4,5,6,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_576 : Lonely 14 (![2,3,4,5,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((1 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_577 : Lonely 14 (![2,3,4,5,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_578 : Lonely 14 (![2,3,4,5,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_579 : Lonely 14 (![2,3,4,5,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_580 : Lonely 14 (![2,3,4,5,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_581 : Lonely 14 (![2,3,4,5,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_582 : Lonely 14 (![2,3,4,5,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_583 : Lonely 14 (![2,3,4,5,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_584 : Lonely 14 (![2,3,4,5,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_585 : Lonely 14 (![2,3,4,5,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_586 : Lonely 14 (![2,3,4,5,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_587 : Lonely 14 (![2,3,4,5,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_588 : Lonely 14 (![2,3,4,5,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_589 : Lonely 14 (![2,3,4,5,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_590 : Lonely 14 (![2,3,4,5,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_591 : Lonely 14 (![2,3,4,5,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_592 : Lonely 14 (![2,3,4,5,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_593 : Lonely 14 (![2,3,4,5,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_594 : Lonely 14 (![2,3,4,5,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_595 : Lonely 14 (![2,3,4,5,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_596 : Lonely 14 (![2,3,4,5,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_597 : Lonely 14 (![2,3,4,5,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_598 : Lonely 14 (![2,3,4,5,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_599 : Lonely 14 (![2,3,4,5,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_600 : Lonely 14 (![2,3,4,5,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_601 : Lonely 14 (![2,3,4,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((1 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_602 : Lonely 14 (![2,3,4,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_603 : Lonely 14 (![2,3,4,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_604 : Lonely 14 (![2,3,4,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_605 : Lonely 14 (![2,3,4,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_606 : Lonely 14 (![2,3,4,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_607 : Lonely 14 (![2,3,4,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_608 : Lonely 14 (![2,3,4,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_609 : Lonely 14 (![2,3,4,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_610 : Lonely 14 (![2,3,4,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_611 : Lonely 14 (![2,3,4,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_612 : Lonely 14 (![2,3,4,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_613 : Lonely 14 (![2,3,4,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_614 : Lonely 14 (![2,3,4,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_615 : Lonely 14 (![2,3,4,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_616 : Lonely 14 (![2,3,4,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_617 : Lonely 14 (![2,3,4,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_618 : Lonely 14 (![2,3,4,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_619 : Lonely 14 (![2,3,4,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_620 : Lonely 14 (![2,3,4,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_621 : Lonely 14 (![2,3,4,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_622 : Lonely 14 (![2,3,4,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_623 : Lonely 14 (![2,3,4,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_624 : Lonely 14 (![2,3,4,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_625 : Lonely 14 (![2,3,4,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_626 : Lonely 14 (![2,3,4,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_627 : Lonely 14 (![2,3,4,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_628 : Lonely 14 (![2,3,4,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_629 : Lonely 14 (![2,3,4,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_630 : Lonely 14 (![2,3,4,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_631 : Lonely 14 (![2,3,4,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_632 : Lonely 14 (![2,3,4,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_633 : Lonely 14 (![2,3,4,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_634 : Lonely 14 (![2,3,4,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_635 : Lonely 14 (![2,3,4,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_636 : Lonely 14 (![2,3,4,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_637 : Lonely 14 (![2,3,4,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_638 : Lonely 14 (![2,3,4,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_639 : Lonely 14 (![2,3,4,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_640 : Lonely 14 (![2,3,4,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_641 : Lonely 14 (![2,3,4,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_642 : Lonely 14 (![2,3,4,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_643 : Lonely 14 (![2,3,4,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_644 : Lonely 14 (![2,3,4,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_645 : Lonely 14 (![2,3,5,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((1 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_646 : Lonely 14 (![2,3,5,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_647 : Lonely 14 (![2,3,5,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_648 : Lonely 14 (![2,3,5,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_649 : Lonely 14 (![2,3,5,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_650 : Lonely 14 (![2,3,5,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_651 : Lonely 14 (![2,3,5,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_652 : Lonely 14 (![2,3,5,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_653 : Lonely 14 (![2,3,5,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_654 : Lonely 14 (![2,3,5,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_655 : Lonely 14 (![2,3,5,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_656 : Lonely 14 (![2,3,5,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_657 : Lonely 14 (![2,3,5,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_658 : Lonely 14 (![2,3,5,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_659 : Lonely 14 (![2,3,5,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_660 : Lonely 14 (![2,3,5,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_661 : Lonely 14 (![2,3,5,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_662 : Lonely 14 (![2,3,5,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_663 : Lonely 14 (![2,3,5,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_664 : Lonely 14 (![2,3,5,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_665 : Lonely 14 (![2,3,5,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_666 : Lonely 14 (![2,3,5,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_667 : Lonely 14 (![2,3,5,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_668 : Lonely 14 (![2,3,5,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_669 : Lonely 14 (![2,3,5,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_670 : Lonely 14 (![2,3,5,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_671 : Lonely 14 (![2,3,5,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_672 : Lonely 14 (![2,3,5,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_673 : Lonely 14 (![2,3,5,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_674 : Lonely 14 (![2,3,5,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_675 : Lonely 14 (![2,3,5,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_676 : Lonely 14 (![2,3,5,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_677 : Lonely 14 (![2,3,5,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_678 : Lonely 14 (![2,3,5,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_679 : Lonely 14 (![2,3,5,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_680 : Lonely 14 (![2,3,5,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_681 : Lonely 14 (![2,3,5,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_682 : Lonely 14 (![2,3,5,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_683 : Lonely 14 (![2,3,5,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_684 : Lonely 14 (![2,3,5,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_685 : Lonely 14 (![2,3,5,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_686 : Lonely 14 (![2,3,5,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_687 : Lonely 14 (![2,3,5,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_688 : Lonely 14 (![2,3,5,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_689 : Lonely 14 (![2,3,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_690 : Lonely 14 (![2,3,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_691 : Lonely 14 (![2,3,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_692 : Lonely 14 (![2,3,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_693 : Lonely 14 (![2,3,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_694 : Lonely 14 (![2,3,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_695 : Lonely 14 (![2,3,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_696 : Lonely 14 (![2,3,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_697 : Lonely 14 (![2,3,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_698 : Lonely 14 (![2,3,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_699 : Lonely 14 (![2,3,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_700 : Lonely 14 (![2,3,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_701 : Lonely 14 (![2,3,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_702 : Lonely 14 (![2,3,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_703 : Lonely 14 (![2,3,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_704 : Lonely 14 (![2,3,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_705 : Lonely 14 (![2,3,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_706 : Lonely 14 (![2,3,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_707 : Lonely 14 (![2,3,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_708 : Lonely 14 (![2,3,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_709 : Lonely 14 (![2,3,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_710 : Lonely 14 (![2,3,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_711 : Lonely 14 (![2,3,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_712 : Lonely 14 (![2,3,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_713 : Lonely 14 (![2,3,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_714 : Lonely 14 (![2,3,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_715 : Lonely 14 (![2,4,5,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((1 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_716 : Lonely 14 (![2,4,5,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_717 : Lonely 14 (![2,4,5,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_718 : Lonely 14 (![2,4,5,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_719 : Lonely 14 (![2,4,5,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_720 : Lonely 14 (![2,4,5,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_721 : Lonely 14 (![2,4,5,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_722 : Lonely 14 (![2,4,5,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_723 : Lonely 14 (![2,4,5,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_724 : Lonely 14 (![2,4,5,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_725 : Lonely 14 (![2,4,5,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (17 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_726 : Lonely 14 (![2,4,5,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_727 : Lonely 14 (![2,4,5,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_728 : Lonely 14 (![2,4,5,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_729 : Lonely 14 (![2,4,5,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_730 : Lonely 14 (![2,4,5,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_731 : Lonely 14 (![2,4,5,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_732 : Lonely 14 (![2,4,5,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_733 : Lonely 14 (![2,4,5,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_734 : Lonely 14 (![2,4,5,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (16 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_735 : Lonely 14 (![2,4,5,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_736 : Lonely 14 (![2,4,5,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_737 : Lonely 14 (![2,4,5,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_738 : Lonely 14 (![2,4,5,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((2 : ℚ) / (15 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_739 : Lonely 14 (![2,4,5,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_740 : Lonely 14 (![2,4,5,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_741 : Lonely 14 (![2,4,5,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_742 : Lonely 14 (![2,4,5,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_743 : Lonely 14 (![2,4,5,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_744 : Lonely 14 (![2,4,5,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_745 : Lonely 14 (![2,4,5,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_746 : Lonely 14 (![2,4,5,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (26 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_747 : Lonely 14 (![2,4,5,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_748 : Lonely 14 (![2,4,5,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_749 : Lonely 14 (![2,4,5,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_750 : Lonely 14 (![2,4,5,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_751 : Lonely 14 (![2,4,5,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_752 : Lonely 14 (![2,4,5,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_753 : Lonely 14 (![2,4,5,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_754 : Lonely 14 (![2,4,5,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_755 : Lonely 14 (![2,4,5,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_756 : Lonely 14 (![2,4,5,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_757 : Lonely 14 (![2,4,5,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_758 : Lonely 14 (![2,4,5,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_759 : Lonely 14 (![2,4,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_760 : Lonely 14 (![2,4,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_761 : Lonely 14 (![2,4,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_762 : Lonely 14 (![2,4,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_763 : Lonely 14 (![2,4,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_764 : Lonely 14 (![2,4,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_765 : Lonely 14 (![2,4,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_766 : Lonely 14 (![2,4,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_767 : Lonely 14 (![2,4,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_768 : Lonely 14 (![2,4,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_769 : Lonely 14 (![2,4,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_770 : Lonely 14 (![2,4,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_771 : Lonely 14 (![2,4,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_772 : Lonely 14 (![2,4,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_773 : Lonely 14 (![2,4,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_774 : Lonely 14 (![2,4,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_775 : Lonely 14 (![2,4,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_776 : Lonely 14 (![2,4,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_777 : Lonely 14 (![2,4,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_778 : Lonely 14 (![2,4,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_779 : Lonely 14 (![2,4,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_780 : Lonely 14 (![2,4,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_781 : Lonely 14 (![2,4,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_782 : Lonely 14 (![2,4,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_783 : Lonely 14 (![2,4,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_784 : Lonely 14 (![2,4,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_785 : Lonely 14 (![2,5,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_786 : Lonely 14 (![2,5,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_787 : Lonely 14 (![2,5,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_788 : Lonely 14 (![2,5,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_789 : Lonely 14 (![2,5,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_790 : Lonely 14 (![2,5,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_791 : Lonely 14 (![2,5,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_792 : Lonely 14 (![2,5,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_793 : Lonely 14 (![2,5,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_794 : Lonely 14 (![2,5,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_795 : Lonely 14 (![2,5,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_796 : Lonely 14 (![2,5,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_797 : Lonely 14 (![2,5,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_798 : Lonely 14 (![2,5,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_799 : Lonely 14 (![2,5,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_800 : Lonely 14 (![2,5,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_801 : Lonely 14 (![2,5,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_802 : Lonely 14 (![2,5,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_803 : Lonely 14 (![2,5,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_804 : Lonely 14 (![2,5,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_805 : Lonely 14 (![2,5,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_806 : Lonely 14 (![2,5,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_807 : Lonely 14 (![2,5,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_808 : Lonely 14 (![2,5,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (28 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_809 : Lonely 14 (![2,5,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_810 : Lonely 14 (![2,5,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_811 : Lonely 14 (![2,6,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_812 : Lonely 14 (![2,6,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_813 : Lonely 14 (![2,6,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((6 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_814 : Lonely 14 (![2,6,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((5 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_815 : Lonely 14 (![2,6,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_816 : Lonely 14 (![2,6,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((7 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_817 : Lonely 14 (![2,6,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_818 : Lonely 14 (![2,7,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((8 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_819 : Lonely 14 (![3,4,5,6,7,8,9,10,11,12,13,14,15] : Fin 13 → ℤ)
    ((((1 : ℚ) / (18 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_820 : Lonely 14 (![3,4,5,6,7,8,9,10,11,12,13,14,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_821 : Lonely 14 (![3,4,5,6,7,8,9,10,11,12,13,14,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_822 : Lonely 14 (![3,4,5,6,7,8,9,10,11,12,13,14,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_823 : Lonely 14 (![3,4,5,6,7,8,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_824 : Lonely 14 (![3,4,5,6,7,8,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_825 : Lonely 14 (![3,4,5,6,7,8,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_826 : Lonely 14 (![3,4,5,6,7,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_827 : Lonely 14 (![3,4,5,6,7,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_828 : Lonely 14 (![3,4,5,6,7,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_829 : Lonely 14 (![3,4,5,6,7,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_830 : Lonely 14 (![3,4,5,6,7,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_831 : Lonely 14 (![3,4,5,6,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_832 : Lonely 14 (![3,4,5,6,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_833 : Lonely 14 (![3,4,5,6,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_834 : Lonely 14 (![3,4,5,6,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_835 : Lonely 14 (![3,4,5,6,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_836 : Lonely 14 (![3,4,5,6,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_837 : Lonely 14 (![3,4,5,6,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_838 : Lonely 14 (![3,4,5,6,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_839 : Lonely 14 (![3,4,5,6,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_840 : Lonely 14 (![3,4,5,6,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_841 : Lonely 14 (![3,4,5,6,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_842 : Lonely 14 (![3,4,5,6,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_843 : Lonely 14 (![3,4,5,6,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_844 : Lonely 14 (![3,4,5,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_845 : Lonely 14 (![3,4,5,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_846 : Lonely 14 (![3,4,5,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_847 : Lonely 14 (![3,4,5,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_848 : Lonely 14 (![3,4,5,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_849 : Lonely 14 (![3,4,5,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_850 : Lonely 14 (![3,4,5,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_851 : Lonely 14 (![3,4,5,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_852 : Lonely 14 (![3,4,5,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_853 : Lonely 14 (![3,4,5,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_854 : Lonely 14 (![3,4,5,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_855 : Lonely 14 (![3,4,5,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_856 : Lonely 14 (![3,4,5,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_857 : Lonely 14 (![3,4,5,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_858 : Lonely 14 (![3,4,5,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_859 : Lonely 14 (![3,4,5,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_860 : Lonely 14 (![3,4,5,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_861 : Lonely 14 (![3,4,5,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_862 : Lonely 14 (![3,4,5,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_863 : Lonely 14 (![3,4,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_864 : Lonely 14 (![3,4,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_865 : Lonely 14 (![3,4,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_866 : Lonely 14 (![3,4,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_867 : Lonely 14 (![3,4,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_868 : Lonely 14 (![3,4,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_869 : Lonely 14 (![3,4,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_870 : Lonely 14 (![3,4,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_871 : Lonely 14 (![3,4,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_872 : Lonely 14 (![3,4,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_873 : Lonely 14 (![3,4,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_874 : Lonely 14 (![3,4,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_875 : Lonely 14 (![3,4,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_876 : Lonely 14 (![3,4,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_877 : Lonely 14 (![3,4,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_878 : Lonely 14 (![3,4,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_879 : Lonely 14 (![3,4,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_880 : Lonely 14 (![3,4,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_881 : Lonely 14 (![3,4,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_882 : Lonely 14 (![3,4,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_883 : Lonely 14 (![3,4,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_884 : Lonely 14 (![3,4,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_885 : Lonely 14 (![3,4,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_886 : Lonely 14 (![3,4,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_887 : Lonely 14 (![3,4,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_888 : Lonely 14 (![3,4,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_889 : Lonely 14 (![3,5,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_890 : Lonely 14 (![3,5,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_891 : Lonely 14 (![3,5,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_892 : Lonely 14 (![3,5,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_893 : Lonely 14 (![3,5,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((9 : ℚ) / (19 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_894 : Lonely 14 (![3,5,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_895 : Lonely 14 (![3,5,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_896 : Lonely 14 (![3,5,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_897 : Lonely 14 (![3,5,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_898 : Lonely 14 (![3,5,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_899 : Lonely 14 (![3,5,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_900 : Lonely 14 (![3,5,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_901 : Lonely 14 (![3,5,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_902 : Lonely 14 (![3,5,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_903 : Lonely 14 (![3,5,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_904 : Lonely 14 (![3,5,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_905 : Lonely 14 (![3,5,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_906 : Lonely 14 (![3,5,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_907 : Lonely 14 (![3,5,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_908 : Lonely 14 (![3,5,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_909 : Lonely 14 (![3,5,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((11 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_910 : Lonely 14 (![3,5,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_911 : Lonely 14 (![3,5,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_912 : Lonely 14 (![3,5,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_913 : Lonely 14 (![3,5,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((13 : ℚ) / (27 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_914 : Lonely 14 (![3,5,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_915 : Lonely 14 (![3,6,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_916 : Lonely 14 (![3,6,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_917 : Lonely 14 (![3,6,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_918 : Lonely 14 (![3,6,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_919 : Lonely 14 (![3,6,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_920 : Lonely 14 (![3,6,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_921 : Lonely 14 (![3,6,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_922 : Lonely 14 (![3,7,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((12 : ℚ) / (25 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_923 : Lonely 14 (![4,5,6,7,8,9,10,11,12,13,14,15,16] : Fin 13 → ℤ)
    ((((1 : ℚ) / (20 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_924 : Lonely 14 (![4,5,6,7,8,9,10,11,12,13,14,15,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_925 : Lonely 14 (![4,5,6,7,8,9,10,11,12,13,14,15,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_926 : Lonely 14 (![4,5,6,7,8,9,10,11,12,13,14,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_927 : Lonely 14 (![4,5,6,7,8,9,10,11,12,13,14,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_928 : Lonely 14 (![4,5,6,7,8,9,10,11,12,13,14,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_929 : Lonely 14 (![4,5,6,7,8,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_930 : Lonely 14 (![4,5,6,7,8,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_931 : Lonely 14 (![4,5,6,7,8,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_932 : Lonely 14 (![4,5,6,7,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_933 : Lonely 14 (![4,5,6,7,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_934 : Lonely 14 (![4,5,6,7,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_935 : Lonely 14 (![4,5,6,7,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_936 : Lonely 14 (![4,5,6,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_937 : Lonely 14 (![4,5,6,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_938 : Lonely 14 (![4,5,6,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_939 : Lonely 14 (![4,5,6,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_940 : Lonely 14 (![4,5,6,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_941 : Lonely 14 (![4,5,6,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_942 : Lonely 14 (![4,5,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_943 : Lonely 14 (![4,5,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_944 : Lonely 14 (![4,5,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_945 : Lonely 14 (![4,5,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_946 : Lonely 14 (![4,5,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_947 : Lonely 14 (![4,5,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_948 : Lonely 14 (![4,5,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_949 : Lonely 14 (![4,6,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (21 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_950 : Lonely 14 (![4,6,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_951 : Lonely 14 (![4,6,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_952 : Lonely 14 (![4,6,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_953 : Lonely 14 (![4,6,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_954 : Lonely 14 (![4,6,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_955 : Lonely 14 (![4,6,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_956 : Lonely 14 (![4,7,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_957 : Lonely 14 (![5,6,7,8,9,10,11,12,13,14,15,16,17] : Fin 13 → ℤ)
    ((((1 : ℚ) / (22 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_958 : Lonely 14 (![5,6,7,8,9,10,11,12,13,14,15,16,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_959 : Lonely 14 (![5,6,7,8,9,10,11,12,13,14,15,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_960 : Lonely 14 (![5,6,7,8,9,10,11,12,13,14,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_961 : Lonely 14 (![5,6,7,8,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_962 : Lonely 14 (![5,6,7,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_963 : Lonely 14 (![5,6,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_964 : Lonely 14 (![5,7,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (23 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

theorem winRow_965 : Lonely 14 (![6,7,8,9,10,11,12,13,14,15,16,17,18] : Fin 13 → ℤ)
    ((((1 : ℚ) / (24 : ℚ)) : ℚ) : ℝ) :=
  KernelGate.lonely_of_kernelWitness (by norm_num) (by decide)

end WindowPack
end LonelyRunner
