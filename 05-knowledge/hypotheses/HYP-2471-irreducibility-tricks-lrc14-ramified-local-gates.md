# HYP-2471 - Irreducibility tricks transfer to LRC14 as ramified local-gate proof tactics

**Status:** OPEN proof-transfer program; exact exception diagnostic verified.

**Source:** codex-2026-06-13. Extends HYP-2451, HYP-2452, HYP-2469, and HYP-2470.

**Computation:** `04-computation/irreducibility_tricks_lrc_transfer_codex.py`; stored output `05-knowledge/results/irreducibility_tricks_lrc_transfer_codex.out`.

## Statement

The useful part of polynomial irreducibility theory for LRC14 is not the word
"irreducible"; it is the local-gate architecture:

```text
primitive normal form
+ local residue gate that empties all split survivors
+ valuation/Newton gate for ramified local exceptions
+ factor-capture budget for hidden allocations
+ dominance/recombination controls outside the finite window.
```

Translated to the LRC14 proof search:

```text
Gauss primitive part        -> gcd/core/carry-window normalization
mod-p irreducibility        -> Q27 / small-shell set-cover infeasibility
Eisenstein/Newton polygon   -> prime-ideal carry and ramified shell openings
Singh/Cohn factor capture   -> obligation-token budget and large-speed dominance
convolution lift            -> hidden denominator-obligation cover lift
```

Thus a normalized LRC14 counterexample should be studied by which local
obligation splits survive, not by the scalar fact "q is blocked."

## Exact HYP-2470 Diagnostic

The two four-core deletion packets that are Q27-feasible in HYP-2470 have the
same ramified shape.

```text
E1:
  D=(28,42,56,84)
  A=(91,322,350,504,936)
  row=(7,14,21,35,49,63,70,77,91,322,350,504,936)
  Q27 witness: none
  first plain witness: q=33, a=13
  Bprime(any)=(True, 322, 83/175812)
  p0=148352/2169475
  base Q27 obligations=630
  added Q27 cover counts=(91:266, 322:182, 350:182, 504:140, 936:90)
```

```text
E2:
  D=(42,56,70,84)
  A=(91,119,700,1008,1066)
  row=(7,14,21,28,35,49,63,77,91,119,700,1008,1066)
  Q27 witness: none
  first plain witness: q=31, a=14
  Bprime(any)=(True, 700, 105799/188042400)
  p0=6870869/125074950
  base Q27 obligations=504
  added Q27 cover counts=(91:266, 1008:140, 700:126, 1066:74, 119:70)
```

In both rows:

```text
7-ideal occupancy = 12/13 speeds
primitive escape count = 1
the unique non-7 added speed is divisible by 13
```

More explicitly:

```text
E1 primitive escape: 936 = 2^3 * 3^2 * 13
E2 primitive escape: 1066 = 2 * 13 * 41
```

The valuation histograms match:

```text
v_7(row):  {0:1, 1:11, 2:1}
v_13(row): {0:11, 1:2}
```

This is the LRC analogue of an Eisenstein/Newton situation: the residue face
looks stubborn, but the valuation/carry face has almost all mass in one prime
ideal plus one primitive leading token. The missing plain shells `31` and `33`
detect exactly that token.

## Proof Consequence

HYP-2470 already proves the residue-face theorem:

```text
carry-window row retaining >=8 core speeds
=> Q27 witness or plain witness q <= 41.
```

HYP-2471 sharpens the next human proof target:

```text
Q27-feasible four-deletion packet
=> ramified 7-ideal / 13-clock portal
=> q in {31,33,41} opens, or Bprime/positive measure opens.
```

The finite scan found only two such packets, but the proposed proof should not
depend on their names. It should prove that the 7-adic occupancy plus a single
13-clock primitive escape cannot maintain Q27 blocking across the missing
plain-shell layer.

## Tournament Analysis

Vertices considered before settling the quotient: runners, gaps, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, matroid circuits, denominator obligations, added-speed cover
masks, valuation primes, and proof tricks.

Selected vertices in the printed tournament are proof carriers/tricks. The
pairwise observable is

```text
(exactness, local_to_global, exception_power,
 LRC_leverage, computability, descent_power).
```

The switch/gauge is primitive normal form first; scalar q-blocking is only the
residue face. The tie Hamiltonian path is the listed trick order in the script.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
edge_flips_vs_tie_path=15/45
leader=integral_convolution_lift_ilp
```

The transitivity is itself useful: at the current LRC14 frontier, the proof
carriers are not fighting. They rank as a pipeline: integral lift/set-cover
certificates first, residue blockers next, valuation exceptions next, then
normalization and factor-capture budgets.

## Next Steps

1. Extract dual/Farkas-style certificates from the HYP-2465/HYP-2470 set-cover
   infeasibility runs, so Q27 becomes a human-readable mod-p-style certificate.
2. Prove the ramified portal lemma: a four-deletion Q27 packet with 12/13 speeds
   in the 7-ideal and a single 13-clock primitive escape opens at `q=31`, `33`,
   `41`, or Bprime.
3. Build the below-eight-core survivor ledger with fields for shell-27 class,
   13-clock debt, divisor fiber, owner/Bprime, and support-load budget.
4. Add a Cohn/Perron outside-window normalizer: very large speeds either open a
   Bprime interval or fold to an equivalent carry-window resource address.
5. Store first-witness shell sequences as a Galois-cycle-type analogue for
   below-eight-core rows.
