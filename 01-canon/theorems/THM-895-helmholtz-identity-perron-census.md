---
id: THM-895
title: THE HELMHOLTZ IDENTITY AND THE PERRON CENSUS — HodgeRank-decomposing the unit arc flow of a tournament on K_n: gradient energy = x/(4n) EXACTLY (harmonic part 0), cycle energy = C(n,2) − x/(4n); so the axis IS the Helmholtz gradient energy and THM-866's +8 walk is gradient ascent. The Kendall–Wei/Perron invariant λ (spectral radius; the all-frames-at-once strength) satisfies λ = 0 ⟺ transitive (acyclic ⟺ nilpotent; count n! verified), corr(λ, −x) = corr(λ, cycle energy) = 0.9313/0.9218 (n = 5/6 labeled censuses), and λ splits every interior x-level at n = 6 (7/7) — the feedback invariant strictly refines the frame-dependent first moment
status: Helmholtz identity PROVED (two lines: Σ_{u<v}(d_u−d_v)² = nΣd² − (Σd)² = nx with potential φ = d/(2n); harmonic = 0 on complete graphs) + verified (400 random n ≤ 9); Perron facts machine-exact (full labeled censuses n = 5, 6); λ = 0 ⟺ transitive PROVED (one line)
source: mac-mini-2026-07-16-S115 (owner: Helmholtz + arborescences + Kendall–Wei self-compositions)
related:
  - THM-866 (the walk = gradient-energy ascent, now literal), THM-855 F3/F5 (per-flip and drift laws = gradient-side mechanics)
  - klein cont.9 HYP-7023/7024 (arborescences = determinant shadow of H — the curl side's tree census; complementary)
  - opus THM-894 (Kendall–Wei on the LRC resonance side)
script: 04-computation/helmholtz_perron_toothpick_macmini_S115.py -> 05-knowledge/results/helmholtz_perron_toothpick_macmini_S115.out
---

# THM-895 — the Helmholtz identity and the Perron census

**(i) Helmholtz/HodgeRank.** Give each arc of tournament T unit flow. On K_n the Hodge
decomposition has no harmonic part, and the least-squares potential is φ_v = d_v/(2n).
Gradient energy: Σ_{u<v}(φ_u − φ_v)² = [nΣd² − (Σd)²]/(4n²) = **x/(4n)** (Σd = 0).
Cycle (curl) energy: C(n,2) − x/(4n). ∎

So: **x is the gradient energy of the tournament flow** (up to 4n), the transitive
tournament is the max-gradient/min-curl point, the regular floor is max-curl, and the
tie-splitting +8 walk of THM-866 is literally gradient-energy ascent in steps of 2/n.

**(ii) Perron census (Kendall–Wei).** λ(T) = spectral radius of the adjacency = the limit
of score-of-scores-of-… (R∘R∘…; the strength defined through all others simultaneously).
- **λ = 0 ⟺ T transitive** (λ = 0 ⟺ A nilpotent ⟺ acyclic ⟺ transitive); labeled count
  at λ = 0 equals n! exactly (n = 5, 6).
- corr(λ, −x) = **0.9313** (n = 5), **0.9218** (n = 6) — and identically equal to
  corr(λ, cycle energy), which is the mechanism: λ is a curl-side invariant.
- **λ splits x-levels**: at n = 6 all 7 interior levels carry positive λ-spread; at n = 5
  only levels x = 8, 16 split (the other levels are single-orbit-thin). The
  all-frames-at-once invariant strictly refines the frame-dependent first moment from
  n = 6 on — velocity-relative/acceleration-absolute, tournament-exact.

**Reading.** The Helmholtz split organizes the invariant zoo: gradient side = scores, x,
the walk, majorization (THM-866/868-869 territory — fully solvable/graded); curl side =
λ, odd cycles, H's higher digits, arborescence discord (klein cont.9) — where the
monodromy and the lawlessness live (THM-865/466). The 0.93 correlation with its exact
complement-form is the quantitative version of "the two sides see the same tournament
from opposite ends".
