---
id: THM-4051
title: "Order-seven exact frontier for the tournament Hamiltonian/Pfaffian strong base"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Across all 456
  isomorphism classes of tournaments of order seven, exactly 353 are strongly
  connected, and every one satisfies the open THM-1950/THM-3729 base
  H(C)>=max(E_even(C),E_odd(C;1)). There is no equality. The minimum additive
  slack is 21, at representative 38 with (H,E_even,E_odd)=(25,2,4). The
  minimum ratio H/max(E_even,E_odd) is 27/8, attained by representative 85298,
  the Paley tournament, with (189,8,56). This corrects THM-1950's sampled
  order-seven minimum 4.22 and its unsupported monotone-growth language. The
  all-order strong base and hence H>=disc remain OPEN.
source: long-precise-frontiers / 2026-08-24
audit: >
  PASS. The primary path reads the complete 456-class universe certified by
  THM-1370-h-spectrum-omits-7-21-all-n, detects strong connectivity by forward
  and converse reachability, counts Hamiltonian paths by subset DP, and
  computes the even and all-one-root odd Pfaffian energies through Bareiss
  determinants and exact Fraction inversion. It also checks all 456 converse
  symmetries. An independent path uses literal enumeration of all 7!
  permutations, dominating-cut tests, and recursive Pfaffian square sums.
  Both paths return the same 353 rows and the same semantic digest; normal and
  optimized streams byte-match the frozen outputs.
depends_on:
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-1950-h-ge-disc-reduced-to-strongly-connected
  - THM-3729-tournament-rooted-response-deletion-average
related:
  - THM-468-tournament-determinant-floor
script: 04-computation/tournament_strong_base_n7_exact_frontier.py
output: 05-knowledge/results/tournament_strong_base_n7_exact_frontier.out
script_sha256: 4d6efb8a7db92dba74429fd472697046bcff7a57b66ed18118a908238f2f2a71
output_sha256: 30e68764683f94a19e10b604e21c4bacfec9aedd1ebefe58b8f38f4944577bc6
independent_audit_script: 04-computation/tournament_strong_base_n7_exact_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_strong_base_n7_exact_independent_audit.out
independent_audit_script_sha256: af0ab20c3cf31c9065c17dc5eb9b7d13be745d362423233748cf20f77ef22e82
independent_audit_output_sha256: 4ad351682da0ae0b2a71e54698530337a23f95684dfc1308ff23944835d62e91
universe: 05-knowledge/results/n7_class_reps_boxeph_S152.txt
universe_sha256: 151070ab6234a9fd9f449f05f93feab64a04b17ad8dae71ba7b7a19679ccd262
semantic_sha256: ef45816e64458457a1b5045c599dd08475e3304a177482de21c1c1f5b61ed350
hash_basis: raw LF bytes
---

# THM-4051 -- exact order-seven frontier of the strong base

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** This theorem
closes one finite order. It does not prove an asymptotic estimate or the
all-order conjecture.

## 1. Target and universe

For a tournament `C` with skew adjacency matrix `K`, let `H(C)` be its number
of directed Hamiltonian paths. Following THM-3729, put

```text
E_even(C) = 2^(-(n-1)) sum_(S even) Pf(K[S])^2 = disc(C),
E_odd(C;1) = 2^(-(n-1)) sum_(S odd) p_1(S)^2
           = disc(C) 1^T(I+K)^(-1)1.                 (1)
```

THM-1950 reduces `H>=disc` to the open strongly-connected base

```text
H(C) >= max(E_even(C),E_odd(C;1)).                    (2)
```

The maintained file `n7_class_reps_boxeph_S152.txt` contains 456 canonical
order-seven representatives. Its completeness is inherited from the
augmentation-cover proof in
THM-1370-h-spectrum-omits-7-21-all-n (not the colliding Jacobian theorem with
the same legacy number). Exactly 353 representatives are strongly connected.

## 2. Exact result

Every one of the 353 strong representatives satisfies `(2)`, and none gives
equality. More sharply,

```text
min_C [H(C)-max(E_even(C),E_odd(C;1))] = 21,          (3)
min_C H(C)/max(E_even(C),E_odd(C;1))   = 27/8.        (4)
```

The extrema are different:

```text
quantity            canonical rep    H     E_even    E_odd
minimum slack       38                25    2         4
minimum ratio       85298             189   8         56
```

Representative `85298` is the Paley tournament on `F_7`. It is regular, so
`1^T(I+K)^(-1)1=7` and `E_odd=7 E_even=56`, exactly as THM-3729 predicts.
Across all strong order-seven classes the inverse response ranges from `2/3`
to `7`.

## 3. Two independent exact paths

The primary certificate uses:

- forward and converse reachability for strong connectivity;
- subset dynamic programming for `H`;
- Bareiss elimination for `det(I+K)`;
- exact rational Gaussian elimination for the all-one response; and
- a converse check on every one of the 456 representatives.

The independent audit shares only the frozen representative universe. It uses
literal enumeration of the `7!` vertex orders for `H`, dominating cuts for
strong connectivity, and recursive Pfaffian expansion for both energies.
Both paths yield the semantic digest

```text
ef45816e64458457a1b5045c599dd08475e3304a177482de21c1c1f5b61ed350.
```

## 4. Correction boundary and survivor

THM-1950 previously reported a sampled order-seven minimum ratio near `4.22`
and inferred that the available room was increasing. The complete quotient
universe finds the Paley value `27/8=3.375`; thus both the numerical extremum
and the monotonicity inference were sampling artifacts. The strongest
survivor is exact: `(2)` holds strictly through order seven. What remains open
is precisely the same all-order strong base, with the high-symmetry Paley
family now promoted to a mandatory hostile control for any proposed margin
theorem.

## 5. Reproduction

From the repository root:

```bash
python3 04-computation/tournament_strong_base_n7_exact_frontier.py
python3 -O 04-computation/tournament_strong_base_n7_exact_frontier.py
python3 04-computation/tournament_strong_base_n7_exact_independent_audit.py
python3 -O 04-computation/tournament_strong_base_n7_exact_independent_audit.py
```

Each normal/optimized pair must byte-match its maintained output, and the two
semantic digests must agree.
