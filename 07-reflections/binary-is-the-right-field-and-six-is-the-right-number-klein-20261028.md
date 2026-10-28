> **CORRECTION 2026-10-28.** This reflection predates the repairs to THM-4840: regularity does not imply signed-graph realizability or TU-signability of every chosen matrix; the Hamming extension bit is edge-count parity, not common degree parity. Integer lifts can change dependencies under reduction. Use [the repaired statements](../05-knowledge/results/free_bits_seven_fiftysix_20261028/INCOMING_REPAIR.md) over the historical wording below.

# Binary is the right field, and six is the right number

*klein-2026-10-28. A REFLECTION on the synthesis session behind
[THM-4840](../01-canon/theorems/THM-4840-tutte-camion-dichotomy-for-arithmetic-square-class-banks-hamming-15-is-the-equal-parity-degree-space-of-k6-rank-6-prime-banks-unit-cayley-tournaments-c3-t5-p7-and-fano-chromatic-numbers.md),
[HYP-9438](../05-knowledge/hypotheses/HYP-9438-every-binary-matroid-is-a-prime-square-class-matroid-arithmetic-r10-and-kuratowski-banks-sum-free-chromatic-numbers-square-sum-descent-and-the-unbounded-dual-word-signature-of-cryptids.md)
and the [packet](../05-knowledge/results/unifying_binary_dependence_free_bits_7_15_56_20261028/README.md). It records why
things worked; every claim carries its file's status.*

## 1. Ask which field the dependence lives in

"Ordinary graph cycle structures cannot represent every arithmetic dependence" sounded like a limitation to be worked
around. It is a theorem with a decision procedure. Square classes are exponent parities, so the dependence is binary;
graph cycles, signed and unsigned incidence, gain graphs and unimodular matrices are all characteristic-zero objects;
Tutte's two excluded-minor theorems and Camion's signing theorem say exactly when the binary object has a
characteristic-zero shadow. The seven-prime Fano bank has none. The lesson is not "graphs fail"; it is "the field was
chosen by the arithmetic, and the graph picture is the regular special case". The same question — which field — is the
one the dipole-law reflection asked of Collatz (2-adic tests, Archimedean certificate) and the one the K₆ identity answers
for the fifteen-prime bank (binary, but one projection away from a graph).

## 2. The hidden number is six

7, 15 and 5 were offered as a mystery. They are `PG(2,2)`, `PG(3,2)` and `dim RM(1,4)`, and the object behind `PG(3,2)` is
a six-set: its fifteen points are the duads, its lines the triangles and perfect matchings of `K_6`, and the Hamming code
is "all degrees of one parity". The affine lift that session85 found (`15 → 16`, five checks) is `K_6` plus its global
degree parity, and the same lift at `7 → 8` gives `AG(3,2)` whose 56 triples are the other 56 in the directive. Once the
six is visible, the counts on six vertices that the repository already had — 15 arcs, 56 tournament classes, the strong
floor `f(6) = 15` — can be typed honestly: same vertices, no mechanism, coincidence. Typing coincidences is what keeps a
genuine identity (the Eisenstein unit reading of `C_3, T_5, P_7`) from drowning.

## 3. Compute the annihilator, not the list of guards

Three notes from three sessions said the same thing in three languages: the coarsest lawful quotient (Myhill–Nerode with
legality observable), the row-span criterion with its exact lost-bit budget, and the Boolean-ideal certificate. The common
form is one sentence: retain the annihilator of the future test algebra; a bit is free iff its test is dependent; a test
costs its new rank; nonlinear legality is ideal membership. The paradigm shift the owner asked for is a change of object
from guards to the dual space. Its payoff is a sorting of the repository's one-bit phenomena into mechanisms, and one
crisp bridge: the halting certificates of the cryptids are dual words of unbounded support, which is why bounded-scope
deciders cannot see them.

## 4. Concurrency, again

Session85 pushed a Fano–Paley note an hour before this one was finished; the 2026-09-17 and 09-22 notes already had the
difference set. The rule from the dipole nights holds: read the newest commits before typing anything as new, inherit by
name, and keep only the genuinely new part — here the ring reading and the match with the pure-blue alphabet.
