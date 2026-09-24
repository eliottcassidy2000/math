---
id: THM-4469
title: "Mahler bridge: Collatz no-divergence on an adjacent supercritical block pair is equivalent to a generalized Mahler Z-number statement"
status: >
  PROVED + INDEPENDENTLY AUDITED. Let B, B' be parity blocks of length L and
  weight a, with p = 3^a > 2q, q = 2^L, and adjacent carries R_B' = R_B + 1.
  Then some positive integer has a 3x+1 parity vector that is eventually a
  concatenation of B and B' if and only if some real xi > 0 has
  xi (p/q)^j in Z + [R_B, R_B + 1]/(p - q) for every j >= 0 (Theorem M').
  For any set A of at least two weight-a blocks with p > q, the Mahler-type
  statement Z(p/q, I_A) implies no-divergence on A^N (Theorem M). The smallest
  instance is (L, a) = (10, 7), with B = 0111101110 (R = 4726),
  B' = 1101100111 (R = 4727), p/q = 2187/1024 and interval
  [4726, 4727]/1163. Consequence: Collatz, and already no-divergence (T1) or
  the Periodicity Conjecture, imply Mahler-type exclusions whose interval
  length (alpha/(alpha-1))/p exceeds the 1/p reached by
  Flatto-Lagarias-Pollington, Dubickas and Bugeaud.
source: collatz-procgen-20260922 session, implication-atlas lane (Theorems M and M', 2026-09-24); audited and promoted by the session orchestrator 2026-09-24
depends_on: []
related:
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
  - 05-knowledge/results/collatz_procgen_20260922_hard_class.md (Proposition M, the coupled Z-number form of C2)
  - 05-knowledge/hypotheses/HYP-9123-supercritical-strips-periodicity.md
  - 05-knowledge/hypotheses/HYP-9134-mahler-bridge-smallest-instance.md
script: 04-computation/experiments/procgen_atlas_20260924_mahler_bridge.py
script_audit: 04-computation/experiments/procgen_atlas_20260924_orchestrator_check.py
output: 05-knowledge/results/procgen_atlas_20260924.out
output_audit: 05-knowledge/results/procgen_atlas_20260924_orchestrator_check.out
script_sha256: cb5ba7ee0d4e2a8cfd2bbf0fb19d9b79def81b1e82bb74fbc3fa5572953c9840
script_audit_sha256: d27b2053a9381612c8c7d96e15daa814c7e7a3f8fc79b17f6ebd17ef814f5b7b
output_sha256: 24929d0dc895b58aec4b60cceb04792ddf8ccd3456a653e32f3cc01b315ccd5e
output_audit_sha256: e93ce493b39663688f1934ea27feff4176916194dde69e1e611ad6b7950ac54c
hash_basis: raw LF bytes
audit: >
  The orchestrator re-derived both proofs line by line: the tail identity
  alpha t_j = R_(c_j)/q + t_(j+1); the digit-forcing window
  [R_B - q/(p-q), R_B + 1 + q/(p-q)], whose only integers are R_B and R_B + 1
  when p > 2q; and the cylinder fact. Independent code, written without
  reading the lane's scripts, confirms:
  * the carries 4726 and 4727;
  * the cylinder classes 990 and 187 mod 1024 (the divisibility classes,
    stable on 49 lifts);
  * the complete adjacent-pair census for L <= 20:
    (10,7):4, (16,11):5, (19,13):1, (20,14):8;
  * the T-side counts 8192, 16, 0 for 1, 2, 3 blocks on [1, 2^22];
  * exact-rational agreement of the xi-side and the T-side at depth 2 on
    56 sampled integer parts, and failure at depth 1 for 200 random
    non-cylinder starts.
  The lane's full atlas pipeline was re-run and its output reproduced byte
  for byte (peak memory 60 MB).
---

# THM-4469 -- the Mahler bridge for adjacent block pairs

**PROVED + INDEPENDENTLY AUDITED.** Atlas: [procgen_atlas_20260924_collatz_implication_atlas](../../05-knowledge/results/procgen_atlas_20260924_collatz_implication_atlas.md) §3.1.

## 1. Setting

`T` is the shortcut map `T(x) = x/2` (x even) and `(3x+1)/2` (x odd) on `Z`.
A *block* is a parity word `B` of length `L` with `a` ones, read first letter
first. On the set of `x` whose next `L` parities spell `B` (its *cylinder*),

```text
T^L(x) = (3^a x + R_B)/2^L,    R_() = 0,   R_(z e) = 3^e R_z + e 2^|z|.
```

**Cylinder fact.** The cylinder of `B` is exactly the class
`{x : 2^L | 3^a x + R_B}`. Both sets are single classes mod `2^L` (Terras's
bijection; `3^a` is a unit mod `2^L`), and the first lies in the second.

Put `p = 3^a`, `q = 2^L`, `alpha = p/q`, `rho = q/p`. For an interval `I`, let
`Z(alpha, I)` be the statement *no `xi > 0` has `xi alpha^j in Z + I` for
every `j >= 0`*.

## 2. Statements

**Theorem M.** Let `p > q` and let `A` be a set of at least two weight-`a`
blocks of length `L`. Put `I_A = [min_A R_B, max_A R_B]/(p - q)`. If some
positive integer's parity vector is `u c` with `c in A^N`, then some `xi > 0`
satisfies `xi alpha^j in Z + I_A` for all `j`. Hence `Z(alpha, I_A)` implies
that no positive integer has a parity vector eventually in `A^N`.

**Theorem M′ (exact equivalence).** Let `A = {B, B'}` with `R_B' = R_B + 1`
and `p > 2q`, and put `I = [R_B, R_B + 1]/(p - q)`. Then some positive integer
has a parity vector eventually in `{B, B'}^N` **iff** some `xi > 0` has
`xi alpha^j in Z + I` for every `j >= 0`.

## 3. Proofs

*Theorem M.*
1. Let `x = T^|u|(x_0)` and `x_j = T^(jL)(x)`. Then `x_(j+1) = (p x_j + R_(c_j))/q`.
2. Set `t_j = (1/p) sum_(i>=0) R_(c_(j+i)) rho^i` and `xi = x + t_0`.
3. Since `rho/q = 1/p`, `alpha t_j = R_(c_j)/q + t_(j+1)`. So `alpha (x_j + t_j) = x_(j+1) + t_(j+1)`, and by induction `xi alpha^j = x_j + t_j`.
4. Finally `t_j in [R_min, R_max] / (p (1 - rho)) = I_A`. ∎

*Theorem M′, direction ⟹.* This is Theorem M with `A = {B, B'}`.

*Theorem M′, direction ⟸.*
1. `|I| = 1/(p-q) < 1`, so `xi alpha^j = x_j + t_j` with unique `x_j in Z` and `t_j in I`.
2. `e_j := q x_(j+1) - p x_j = p t_j - q t_(j+1)` is an integer in `[R_B - q/(p-q), R_B + 1 + q/(p-q)]`.
3. Since `p > 2q` gives `q/(p-q) < 1`, the only integers there are `R_B` and `R_B + 1 = R_B'`.
4. Hence `x_(j+1) = (p x_j + R_(c_j))/q` with `c_j in {B, B'}`. By the cylinder fact, `x_j` lies in the cylinder of `c_j` and `T^L(x_j) = x_(j+1)`.
5. `xi alpha^j -> infinity`, so some `x_(j0)` is a positive integer, and its parity vector lies in `{B, B'}^N`. ∎

## 4. Instances and reach

* **The complete census of adjacent pairs with `p > 2q` for `L <= 20`:**

  | `(L, a)` | `alpha` | pairs |
  |---|---|---|
  | `(10, 7)` | `2187/1024` | 4 |
  | `(16, 11)` | `177147/65536` | 5 |
  | `(19, 13)` | `1594323/524288` | 1 |
  | `(20, 14)` | `4782969/1048576` | 8 |

  The smallest is `0111101110 (4726)`, `1101100111 (4727)`, with `I = [4726, 4727]/1163`.
* **Theorem M's non-adjacent instances** (`|I_A| < 1`, `alpha |I_A| < 1`) exist at `(4,3)`, `(5,4)`, `(7,5)` and many more (atlas §3.1).
* **Where proved Mahler-type results stop.**
  * Flatto–Lagarias–Pollington: every `xi > 0` has spread of `{xi (p/q)^n}` at least `1/p`.
  * Dubickas: excludes `[s, s + 1/p]` for `q < p < q^2`.
  * Bugeaud: excludes `[s, s + 1/p]` for almost every `s`.

  All of these exclude intervals of length at most `1/p`. Theorem M′'s intervals have length `(alpha/(alpha - 1))/p`: `1.88/p`, `1.59/p`, `1.49/p` and `1.28/p` at `L = 10, 16, 19, 20`. At that length the real confinement is *equivalent* to a free binary choice at each step together with 2-adic divisibility, so any exclusion there proves a HARD slice of no-divergence. Such a slice has density `a/L > log_3 2`, entropy `1/L` bits per step and bounded discrepancy, and no capacity or repetition argument reaches it.

## 5. Consequences and scope

* **One proves the other.**
  * Collatz ⟹ T1 ⟹ `Z(2187/1024, [4726, 4727]/1163)`.
  * The Periodicity Conjecture ⟹ T1 (Bernstein–Lagarias) ⟹ the same.
  * Conversely, a proof of this Mahler-type statement proves no-divergence on the class `{B, B'}^N`.
  * This is the only arrow from the Collatz family to a famous-type problem found by the 2026-09-24 atlas.
* **Mahler's own `3/2` problem does not transfer.** It is the boundary `a = L`: one block, no pairs. Its interval `[0, 1/2)` is a proper part of its digits' tail hull, so the real condition there is an extra constraint rather than an equivalent one (atlas §3.3: two-place orthogonality).
* **Heuristics.** `alpha |I| = alpha/(q(alpha - 1)) < 0.002` for these instances, so no `xi` is expected. The adjacent-pair class is the Mahler-side twin of the cube-swap number (HYP-9127): explicit, supercritical, and beyond both proved mechanisms.
* **Priority.** No claim of priority is made. The equivalence was not found in the sources read: FLP, Bugeaud 2004, Akiyama–Frougny–Sakarovitch 2008, Mahler 1968, Lagarias's surveys and Dubickas (the last through secondary accounts).

The open instance is registered as HYP-9134.
