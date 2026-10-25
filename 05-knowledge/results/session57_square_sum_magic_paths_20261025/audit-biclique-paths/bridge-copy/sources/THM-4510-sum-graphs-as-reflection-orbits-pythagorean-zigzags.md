---
id: THM-4510
title: "Sum graphs are unions of reflection matchings x -> t - x, and two reflections compose to the translation by the gap between their targets. Two targets give a Hamiltonian path only for gaps 1 and 2. Three targets give one iff the union is tight and gcd(b-a, c-b) = 1 (or 2 with a odd); the path is then one orbit of a rotation. Every primitive Pythagorean triple s^2 + t^2 = u^2 gives a square-sum chain of 1..t^2-1 and of 1..t^2 using only s^2, t^2, u^2; the owner's '15, with 8 and 9 at the ends' is the triple (3,4,5), and the first square-sum window {15,16,17} is the single tight triple (9,16,25). In the Collatz alphabet C_n the top layer is a forced zigzag translating by |3^a - 2P|, and three failing families exist at every level with endpoints affine in 2^p and 3^(a-1), governed by rho_a = 2^p/3^(a-1)"
status: >
  PROVED + INDEPENDENTLY AUDITED (Lemma R, Theorem A2, Theorem A3 with
  Corollaries A3', A3'', A3''', Lemma T, Lemma K, the Choke Lemma,
  Propositions W1-W3 for every level a); FINITE-EXACT (C_n Hamiltonian
  set to 2200, reproduced by a second code path; W_8 = [4374, 6561]
  Hamiltonian at every n; propagation alone decides every failure in
  range); CITED (Gerbicz, A090461; Weyl equidistribution); CONJECTURE B
  (propagation completeness, connection form of the inner endpoints,
  rho-scaling, right ends); ANALOGY (interval-exchange/Rauzy reading;
  |2^q - 3^a| as Collatz cycle denominators). Collatz is not addressed.
  Setting: G_S(n) on 1..n, x ~ y iff x + y is in S; M_t = {x, t - x}
  is a matching, and G_S(n) is the union of the M_t.
  (R) r_t2 o r_t1 is the translation by t2 - t1; components lie in
  orbits (x + gZ) u (t_0 - x + gZ).
  (A2) Two targets s < t never make a cycle. They give a Hamiltonian path
  of [n] iff {s,t} is {n,n+1} or {n+1,n+2}, or the gap is 2 with s odd and
  s in {n-1, n, n+1}. For C_n (powers of 2 and 3) this happens only at
  n = 3, 7, 8 (the Gersonides pairs (3,4), (8,9)).
  (A3) Three targets a < b < c give a Hamiltonian path iff (i) max degree
  <= 2, (ii) |M_a| + |M_b| + |M_c| = n - 1, and (iii)
  gcd(b-a, c-b) = 1, or 2 with a odd. The alternate vertices then follow
  the rotation x -> x + (c-b) mod (c-a).
  (A3''') For every primitive Pythagorean triple s^2 + t^2 = u^2 (s < t),
  the squares s^2, t^2, u^2 give a Hamiltonian path of [n] exactly for
  n in {t^2 - 1, t^2} (and n = 17 for (3,4,5)), following the rotation by
  s^2 on Z/t^2. At n = t^2 - 1 its ends are s^2 and t^2/2 (t even), or
  s^2/2 and s^2 (t odd). (3,4,5) gives Q_15 with ends 9 and 8. The next
  triples give three-square chains of Q_143 (ends 25, 72), Q_224 (32, 64)
  and Q_575 (49, 288).
  (A3'') Three consecutive squares form a Hamiltonian union only for
  (9,16,25) at n = 15, 16, 17. So the first square-sum window
  {15,16,17} is this single tight triple. It closes at 18 because 7 gets
  a third partner.
  (T) In C_n every target above max(P, Q) is 2P or 3Q, so the top layer
  is a forced two-target zigzag translating by |3^a - 2P|. For W_3..W_8
  these are 5, 17/47, 13, 217/295, 139/1909, 1631.
  (W1-W3) With B = 3^(a-1), P the power of 2 in (B, 2B) and rho = P/B:
  - W1: C_n has no Hamiltonian path for
    max(5P/4, 3B - 3P/4) <= n < min(3P/2, 3B - P/2) (two forced ends and
    a choke at P/4; nonempty iff 4/3 < rho < 12/7);
  - W2: similarly with a double choke at P/4 (nonempty iff
    4/3 < rho < 8/5);
  - W3: chokes at P/4 and P/8 with disjoint resolving sets (nonempty iff
    4/3 < rho < 3/2).
  They occur at levels of density 0.363, 0.263, 0.170 (Weyl). Examples:
  W_7's 1419-1535, 1714-1791, 1803-1919; new failing ranges
  [40960, 42664] at a = 10, and [334833, 393215] and [419830, 465904] at
  a = 12.
  (B, FINITE-EXACT) W_8 = [4374, 6561] is Hamiltonian at every n, the
  first fully Hamiltonian window since W_2. Every failure up to W_8 is
  refuted by unit propagation alone. C_243 has exactly one Hamiltonian
  path.
  Square sums 18-24: 18 is local (three leaves), 19-22 are small forcing
  conflicts (chokes at 3, 5 and the defect 2), and 24 is the only global
  failure.
  OPEN:
  - Conjecture B;
  - closed families for the deeper conflicts (W_5: 195-206;
    W_7: 1617-1619, 1664-1702, 2115-2119);
  - any positive Hamiltonicity statement for all levels.
source: collatz-procgen-20260922 session, sumgraph lane (2026-09-26), answering the owner's prompts on small graphs that encode arithmetic and the square-sum problem at 15 with 8 and 9 at the ends; the reflection-matching idea was the orchestrator's; audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md (Q_n, C_n, the Window Theorem, the zigzag law)
related:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (the gaps |2^p - 3^a|, as numbers only)
  - 05-knowledge/results/crossroads_family_20260926_squares.md (the n = 15 structure, first switch at 46)
note: 05-knowledge/results/procgen_sumgraph_20260926_reflection_orbits.md
scripts: 04-computation/experiments/procgen_sumgraph_20260926_{run,solver,theory}.py
script_audit: 04-computation/experiments/procgen_sumgraph_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_sumgraph_20260926.out
output_sha256: 4e1465d2beaf6bc974b451dc79732030fff6e3be632d8cff72ffcce0049b624b
output_audit: 05-knowledge/results/procgen_sumgraph_20260926_orchestrator_check.out
output_audit_sha256: 174c03427524840d0647ff0ba7ace8bb9ee465ab3b95ad7410cd8246df90919b
hash_basis: raw bytes
audit: >
  The orchestrator read the proofs of:
  - Lemma R, Theorems A2 and A3 (acyclicity via translations; tightness;
    the gcd condition from the orbit decomposition);
  - Corollaries A3', A3'', A3''';
  - Lemma T and the Choke Lemma;
  - Propositions W1-W3,
  and found them sound.
  Independent code (procgen_sumgraph_20260926_orchestrator_check.py,
  written from the note's statements without reading the lane's solver
  or theory code) confirms:
  - Theorem A2 by brute force over all target pairs, n = 3..60;
  - Theorem A3 over all 85261 target triples with n = 3..24 (758 paths);
  - Corollary A3''' for all 11 primitive triples with t <= 60 (paths
    exactly at t^2-1, t^2, plus 17 for (3,4,5), with the stated ends);
  - Corollary A3'' (j <= 25, n <= 700);
  - Lemma T for n <= 5000;
  - the W1-W3 certificates, re-derived from neighbour lists at every n of
    the ranges for a = 5, 7, 10 and at 202 n of each range for a = 12,
    and the rho-criteria for all a <= 200;
  - OR-tools CP-SAT (independent of both lanes' solvers; paths verified
    edge by edge) on W_5: PATH at 179, 180, 194, 224, 243 and NONE at
    195, 200, 206, 207, 223.
  The lane's full pipeline was re-run with --full (1400 s, 47 checks,
  peak RSS 66 MB). It is identical to the committed .out except for one
  timing field.
  Scope notes:
  - W_8's full Hamiltonicity and the propagation-completeness statements
    rest on the lane's solver. The C_n list to 2200 now has two
    independent code paths (THM-4505's lanes and this lane).
  - The IET/Rauzy reading and the Collatz-denominator coincidence are
    ANALOGY.
---

# THM-4510 — sum graphs as reflection orbits, and why 8 and 9 are (3, 4, 5)

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_sumgraph_20260926_reflection_orbits](../../05-knowledge/results/procgen_sumgraph_20260926_reflection_orbits.md).

**UPDATE 2026-09-27 (cnpos lane, [positive-Hamiltonicity note](../../05-knowledge/results/procgen_cnpos_20260926_positive_hamiltonicity.md); orchestrator-audited): an exact reduction at every level and verified constructions far beyond `W_8`.**

**Theorem F (PROVED; first return).** For every `n` in every window `W_a` (`n != 3^a`), let `T1 < T2` be the two top targets and `m = T2 - T1 = |3^a - 2^k|`.
- The forced top zigzag returns to the residual `[1, T2 - 1 - n]` as the reflection `phi(x) = T1 - x (mod m)`.
- So every solution of the *residual problem* (the `phi`-pairs plus a matching by C-edges) lifts to a Hamiltonian path of `C_n`.
- The reduction is an equivalence in the lower part of every two-power window. Elsewhere Theorem F' allows one binary choice per vertex of `[3^a - P, P]`.

Also PROVED:
- **Proposition U.** The Gersonides relations `2-1 = 4-3 = 3-2 = 9-8 = 1` make the first zone steps rotations by units mod `m`, hence cycle-free.
- **Proposition R.** A single rotation solves the residual problem only when `m + 1` or `m + 2` is a target, or two targets differ by `m`. For `a <= 300` (FINITE-EXACT) this happens only at `a = 2, 3, 5`, the Pillai coincidences `4-3 = 9-8`, `9-4 = 32-27`, `16-3 = 256-243`.
- **Proposition S.** Exchanging two path edges needs a target relation `t1 + t2 = t3 + t4`. Below `2^200` there are exactly four, all through the edges `{1,2}` or `{1,3}`.

FINITE-EXACT, with every path re-verified edge by edge:
- `C_(T1 - 1)` is Hamiltonian at every level `2 <= a <= 15`, up to `n = 14,348,906`.
- The window right ends are Hamiltonian for `4 <= a <= 14`: Conjecture B4 through `a = 14`.
- Every Hamiltonian `n` in `W_3..W_8` has a verified path. Of these, 3 (`n = 5116, 5117, 5248`) come only from the sumgraph solver's search.

OPEN:
- the residual problem at every level (Conjecture RP), and with it a positive theorem for infinitely many levels;
- `RP(59, 22)` and `RP(383, 1)` have no solution, so not every residual instance is solvable.

Collatz typing: `m = |3^a - 2^k|` coincides with the cycle denominators as numbers only (ANALOGY).

Audit: the orchestrator's own checker (procgen_cnpos_20260926_orchestrator_check.py, written from the definition of `C_n`; the lane's construction is used only to produce certificates) verified:
- the constructed paths for the clean `n` at `a = 2..11`, up to `n = 131071`;
- the window right ends at `a = 6, 7, 8, 9, 10`.

Rerun: the lane's runner in default mode (every 8th `n` of `W_8`; 376 s, peak 352 MB, 13 checks, ALL CHECKS PASSED) is identical to the committed full run (`--full`, 4784 s) up to timing and the sampled `W_8` section.

Resource note: one development run exceeded the lane's 500 MB cap (1.2 GB for about 26 s). The final runner peaks at 316-352 MB.

## 1. Every sum graph is a union of reflections

For each target `t`, the pairs `{x, t - x}` form a matching, the reflection about `t/2`. Two reflections compose to a translation by the difference of their targets. So along any chain the numbers in even positions, and those in odd positions, move by target gaps.

This forces the structure:
- **Two targets** never close a cycle. They cover `1..n` in one chain only when the targets are consecutive (Gersonides pairs such as `8, 9`) or two apart.
- **Three targets** give a chain exactly when the union is tight and the gaps are coprime. The chain is then a single orbit of a rotation: a discrete rotation whose continuous analogue would need an irrational angle.

## 2. The square-sum problem at 15

The chain `8,1,15,10,6,3,13,12,4,5,11,14,2,7,9` uses the three squares `9, 16, 25`, and `9 + 16 = 25`. It is the Pythagorean triple `(3,4,5)`: alternate entries rotate by `9` modulo `16`, and the ends are `3^2 = 9` and `4^2/2 = 8`.

Every primitive triple does the same:
- `(5,12,13)` chains `1..143` with ends `25, 72`;
- `(8,15,17)` chains `1..224` with ends `32, 64`.

What is special about 15 is only that it is the first such chain: below 15 there is none, and the whole first window `15, 16, 17` is this one triple. THM-4505's Pell identity is the statement that among these triples only `(3,4,5)` has the zigzag shape `(2k+1, 4k, 6k+1)`.

## 3. The Collatz alphabet at every scale

In `C_n` (sums that are powers of 2 or 3), the biggest vertices can only pair through `2P` and `3^a`. So they form a forced zigzag translating by the gap `|3^a - 2P|`.

Just below, "chokes" appear. A choke is a small vertex such as `P/4` whose partners are all taken by forced vertices. They switch Hamiltonicity off on explicit ranges of `n` whose ends are affine in `2^p` and `3^(a-1)`. Whether a level has these ranges depends only on `rho_a = 2^p/3^(a-1)`, the fractional part of `(a-1) log_2 3`. This is the exact form of the "fractal recursion across levels" that the Beatty word only coarsely records.
