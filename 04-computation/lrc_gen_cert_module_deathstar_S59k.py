#!/usr/bin/env python3
"""
death-star-2026-07-19-S59k -- HYP-8010: the certificate-module GENERATOR.

gen_member_module(family, D, Q, witness_a, tag, module_name) emits a complete
Lean module proving  sSup (margin v '' Icc 0 1) = D/Q  for the family, in the
S59j working stack:
  per-modulus Bool decides (certCheckS, one declaration per pair-sum modulus,
  memory-GC-friendly) + flat-constructor assembly (moduli_ok) +
  contains-sweep coverage (Bool decide + List.contains_iff_mem) +
  lazy reflection (cert_of_check) + element-wise floor (rung_floor_single).

It VALIDATES everything Python-side first and refuses to emit on any failure:
  V1 the witness band: every element's residue at witness_a lies in [D, Q-D];
  V2 Q is a pair sum and D/Q is the exact M (via the pair-sum evaluator if
     affordable, else the certificate itself is the upper-bound proof);
  V3 the per-modulus certificate: every (S, kk) has a killing element
     (exactly Lean's certCheckS semantics -- if V3 passes, the decides pass).

Usage (this session): python3 lrc_gen_cert_module_deathstar_S59k.py
  -> emits TournamentH7/LRCEChannelCert247.lean for F_4(61) = {1..59,61,240}.
"""
import sys, time
from math import gcd

def validate_and_data(fam, D, Q, a):
    fam = sorted(fam)
    n = len(fam)
    # V1: witness band
    for v in fam:
        r = (v * a) % Q
        if not (D <= r <= Q - D):
            raise SystemExit(f"V1 FAIL: element {v} residue {r} outside [{D},{Q-D}]")
    # V2: Q is a pair sum
    sums = sorted(set(x + y for x in fam for y in fam))
    if Q not in sums:
        raise SystemExit(f"V2 FAIL: Q={Q} is not a pair sum")
    # V3: the certificate, exactly certCheckS semantics per modulus
    t0 = time.time()
    for S in sums:
        for kk in range(S):
            ok = False
            for w in fam:
                r = (w * kk) % S
                if r * Q <= D * S or (S - r) * Q <= D * S:
                    ok = True
                    break
            if not ok:
                raise SystemExit(f"V3 FAIL: no kill at (S={S}, kk={kk})")
    print(f"  validation passed: {len(sums)} moduli, "
          f"worst S={max(sums)}, {time.time()-t0:.1f}s")
    return fam, sums

def gen_member_module(fam, D, Q, a, tag, path, extra_note=""):
    fam, sums = validate_and_data(fam, D, Q, a)
    n = len(fam)
    fam_str = ",".join(map(str, fam))
    mod_str = ", ".join(map(str, sums))
    chks = "\n".join(
        f"theorem chk{tag}_{S} : certCheckS l{tag} {D} {Q} {S} = true := by decide"
        for S in sums)
    covs = "\n".join(
        f"theorem cov{tag}_{i} : covOK{tag} {v} = true := by decide"
        for i, v in enumerate(fam))
    cov_facts = ", ".join(f"cov{tag}_{i}" for i in range(n))
    facts = ", ".join(f"chk{tag}_{S}" for S in sums)
    # the Fin-family definition: consecutive prefix then exceptions
    # (general form: nested ite over the index; here we emit a linear ite chain
    #  compressed for prefix runs)
    runs = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and fam[j+1] == fam[j] + 1:
            j += 1
        runs.append((i, j))
        i = j + 1
    # build ite chain: for a run [i..j] with fam[k] = fam[i] + (k - i):
    conds = []
    for (i, j) in runs:
        off = fam[i] - i
        if i == j:
            conds.append((f"(i : ℕ) = {i}", f"{fam[i]}"))
        else:
            conds.append((f"(i : ℕ) ≤ {j}", f"(i : ℕ) + {off}"))
    body = ""
    for idx, (c, val) in enumerate(conds):
        if idx == len(conds) - 1:
            body += f"{val}"
        else:
            body += f"if {c} then {val} else "
    hb_cov = max(2000000, 4 * n * n * len(sums))
    checks_module = f'''/-
  TournamentH7.LRCEChannelCert{tag}Checks -- GENERATED (HYP-8010).
  Defs + {len(sums)} per-modulus certificate decides for the {D}/{Q} member
  (family {{{fam_str}}}, worst S = {max(sums)}).  Stable: build once.
-/
import TournamentH7.LRCEChannelCert

namespace LonelyRunner
namespace EChannelCert

/-- The family as a `Fin {n}` function. -/
def v{tag} : Fin {n} → ℤ := fun i =>
  {body}

/-- The value list. -/
def l{tag} : List ℕ := [{fam_str}]

/-- The {len(sums)} distinct pair-sum moduli. -/
def moduli{tag} : List ℕ := [{mod_str}]

/-- vi is covered: every pair sum vi + vj is a listed modulus (opaque for assembly). -/
def covOK{tag} (vi : ℕ) : Bool := l{tag}.all fun vj => moduli{tag}.contains (vi + vj)

section PerModulusChecks{tag}
set_option maxRecDepth 8000
{chks}
{covs}
end PerModulusChecks{tag}

end EChannelCert
end LonelyRunner
'''
    module = f'''/-
  TournamentH7.LRCEChannelCert{tag} -- GENERATED (HYP-8010).
  Kernel-exact loneliness value  sSup (margin '' [0,1]) = {D}/{Q}  for the
  family {{{fam_str}}} ({n} speeds), witness t0 = {a}/{Q}.{extra_note}
-/
import TournamentH7.LRCEChannelCert{tag}Checks

namespace LonelyRunner
namespace EChannelCert

open GridAttainment TournamentH7.LRCWitness

set_option maxRecDepth 8000 in
set_option maxHeartbeats 4000000 in
/-- All moduli pass (flat Bool assembly, no case dispatch). -/
theorem moduli{tag}_ok : moduli{tag}.all (certCheckS l{tag} {D} {Q}) = true := by
  simp only [moduli{tag}, List.all_cons, List.all_nil, Bool.and_eq_true]
  exact ⟨{facts}, trivial⟩

theorem chk{tag}_mem : ∀ S ∈ moduli{tag}, certCheckS l{tag} {D} {Q} S = true := by
  have h := moduli{tag}_ok
  rw [List.all_eq_true] at h
  exact h

/-- Every pair sum is a listed modulus (flat assembly of cached per-element decides;
    S59m lesson: the monolithic 91x91x271 kernel sweep blows the build window). -/
theorem sums_covered{tag} :
    ∀ vi ∈ l{tag}, ∀ vj ∈ l{tag}, (vi + vj) ∈ moduli{tag} := by
  have h : (l{tag}.all covOK{tag}) = true := by
    simp only [l{tag}, List.all_cons, List.all_nil, Bool.and_eq_true]
    exact ⟨{cov_facts}, trivial⟩
  simp only [List.all_eq_true] at h
  intro vi hvi vj hvj
  have h2 := h vi hvi
  simp only [covOK{tag}, List.all_eq_true] at h2
  exact List.contains_iff_mem.mp (h2 vj hvj)

theorem check_{tag} : certCheck l{tag} {D} {Q} = true := by
  simp only [certCheck, List.all_eq_true]
  intro vi hvi vj hvj
  exact chk{tag}_mem _ (sums_covered{tag} vi hvi vj hvj)

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2000000 in
/-- The e-channel certificate. -/
theorem cert_{tag} : Cert v{tag} {D} {Q} := by
  have h := cert_of_check l{tag} {D} {Q} check_{tag} v{tag}
    (by decide) (by decide) (by decide)
  exact_mod_cast h

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2000000 in
/-- **Kernel-exact loneliness value**: `sSup (margin '' [0,1]) = {D}/{Q}`. -/
theorem member_{tag}_exact :
    sSup (margin v{tag} '' Set.Icc (0:ℝ) 1) = ({D} : ℝ) / ({Q} : ℝ) := by
  have hD : (({D}:ℤ) : ℝ) = ({D}:ℝ) := by norm_num
  have hQ : (({Q}:ℤ) : ℝ) = ({Q}:ℝ) := by norm_num
  rw [← hD, ← hQ]
  apply margin_sSup_eq_of_cert v{tag} (by decide) {D} {Q} (by norm_num) (by norm_num)
    (({a} : ℝ) / {Q}) (by constructor <;> norm_num)
  · intro i m
    have hband : {D} ≤ (v{tag} i * {a}) % {Q} ∧
        (v{tag} i * {a}) % {Q} ≤ {Q} - {D} := by
      fin_cases i <;> decide
    have h := RungFloor.rung_floor_single {Q} {D} (v{tag} i) {a} (by norm_num) hband m
    have hc : (({Q} : ℤ) : ℝ) = ({Q} : ℝ) := by norm_num
    rw [hc] at h
    convert h using 2 <;> norm_num
  · exact cert_{tag}

end EChannelCert
end LonelyRunner

#print axioms LonelyRunner.EChannelCert.check_{tag}
#print axioms LonelyRunner.EChannelCert.cert_{tag}
#print axioms LonelyRunner.EChannelCert.member_{tag}_exact
'''
    checks_path = path.replace(".lean", "Checks.lean")
    with open(checks_path, 'w') as f:
        f.write(checks_module)
    with open(path, 'w') as f:
        f.write(module)
    print(f"  emitted {checks_path} + {path}: {n} speeds, {len(sums)} moduli")

if __name__ == "__main__":
    fam = list(range(1, 60)) + [61, 240]
    gen_member_module(fam, 4, 247, 70, "247",
        "04-computation/lean/TournamentH7/TournamentH7/LRCEChannelCert247.lean",
        "\n  F_4(61), THM-1286's second D=4 gate member (predicted-then-found, S59c).")
