#!/usr/bin/env python3
"""death-star-S59m v2: per-element coverage decides behind an OPAQUE named
predicate covOK367 (the moduli_ok pattern -- simp unfolds only the outer list,
never the inner all), so the assembly is a flat Bool constructor."""
fam = sorted(list(range(1, 90)) + [91, 360])
covs = "\n".join(
    f"theorem cov367_{i} : covOK367 {v} = true := by decide"
    for i, v in enumerate(fam))
facts = ", ".join(f"cov367_{i}" for i in range(len(fam)))
cov_module = f'''/-
  TournamentH7.LRCEChannelCert367Cov -- GENERATED (S59m v2).
  Per-element coverage decides for the 4/367 member, behind the opaque
  predicate covOK367 so assembly-side simp cannot unfold the inner sweep.
-/
import TournamentH7.LRCEChannelCert367Checks

namespace LonelyRunner
namespace EChannelCert

/-- vi is covered: every pair sum vi + vj is a listed modulus. -/
def covOK367 (vi : ℕ) : Bool := l367.all fun vj => moduli367.contains (vi + vj)

section Cov367
set_option maxRecDepth 8000
{covs}
end Cov367

end EChannelCert
end LonelyRunner
'''
open("04-computation/lean/TournamentH7/TournamentH7/LRCEChannelCert367Cov.lean", "w").write(cov_module)

p = "04-computation/lean/TournamentH7/TournamentH7/LRCEChannelCert367.lean"
s = open(p).read()
marker = "/-- Every pair sum is a listed modulus (assembled from cached per-element decides). -/"
i0 = s.find(marker)
assert i0 >= 0, "patched marker not found"
i1 = s.find("theorem check_367", i0)
assert i1 > i0
new_cov = f'''/-- Every pair sum is a listed modulus (assembled from cached per-element decides). -/
theorem sums_covered367 :
    ∀ vi ∈ l367, ∀ vj ∈ l367, (vi + vj) ∈ moduli367 := by
  have h : (l367.all covOK367) = true := by
    simp only [l367, List.all_cons, List.all_nil, Bool.and_eq_true]
    exact ⟨{facts}, trivial⟩
  simp only [List.all_eq_true] at h
  intro vi hvi vj hvj
  have h2 := h vi hvi
  simp only [covOK367, List.all_eq_true] at h2
  exact List.contains_iff_mem.mp (h2 vj hvj)

'''
s = s[:i0] + new_cov + s[i1:]
open(p, "w").write(s)
print("v2 emitted: opaque covOK367 + 91 decides + flat assembly")
