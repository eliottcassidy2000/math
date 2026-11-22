# Physical two-port branches and exact inverse reads

**2026-11-22 · PROVED in the stated families. Chan–Pak OPEN.**

The first result here is an all-word construction, not a census:
**every finite binary word compiles to an actual two-port tournament
whose ordinary Hamiltonian-path count uniquely recovers the word.**
The compiler uses linearly many vertices in the word length. Its image
has power-law lower and upper bounds but has density zero.

- [One count, two legal return moves](ONE-COUNT-INVERSE.md): legal states
  land on opposite sides of an exact linear test. An irreducible return
  matrix and an algebraic norm give a uniform scalar observable. A binary
  delimiter and width bands remove the need to supply the depth. A four-variable integer
  feasibility problem gives a polynomial-time inverse on the restricted image.
- [Two independent flanks and an exact +6 switch](TWO-FLANKS.md): both
  ports can run their own rooted history. The two middle tournaments
  give outputs in disjoint classes modulo 32; factoring one explicitly
  transformed target lists every compatible pair of actual histories.
- [Odd valuations and recurrent tags](ODD-TAGS.md): a physical common reset
  becomes rank one modulo 15 without losing its exact integer inverse.
  Two strong ordinary count readouts recover its history, while the first
  read itself separates the final digit and stays coprime to 15. A separate
  contact-word construction steers both marked coordinates jointly to
  **every vector modulo every odd modulus**, with infinitely many
  distinct strong counts at each vector.
- [Strong self-converse branch readout](SHELL-BRANCHES.md): one common outer
  wrapper recovers the full slot polynomial of any fixed-order host from
  a single ordinary count. Binary outer choices can vary the reflected-pair
  count while leaving a prescribed reversal-fixed count unchanged.
- Exact replay of the new constant identities and short literal controls:
  run `python3 check_twoport_inverse.py` from this directory. The proof,
  not the replay, covers arbitrary words and the large common suffix.

## Inheritance and transfer

The closest prior mechanism is the actual four-coordinate
[root-child compiler](../two_port_return_induction_20261103/ROOTED_RESETS.md),
joined here with the [marked shell limits](../chanpak_matching_folds_branching_20261122/SHELL-INFORMATION.md).
The hostile is a selector that only distinguishes reductions or residue
classes: it can certify the final digit without identifying an absolute
predecessor. The retained sidecar is the integral profile and the legal
controller word, including the port renaming. The closest near miss is
equal total count on two marked hosts with different future responses.

Two ideas from the requested readings are operational here.
[Pasten–Vidaux](https://arxiv.org/pdf/2611.24748), §§4–7, distinguish
local approximation from exact selection using valuation propagation and
a height inequality. The one-count inverse replaces the arithmetic
input by a nonzero algebraic norm, while retaining the local-to-global
gap as a separate proof obligation.
[The linked RH/GRH remark](https://alpo.ge/remark.pdf), pp. 1, 3–4,
suggests iterating a lift with known gain and charging its cost. Here
we can compute literal transfer amplification, bound all other modes,
and account for every added vertex. No RH/GRH assertion is assumed or
evaluated.

## Boundary

A decoder of the constructed image is not a proof that an arbitrary
odd target belongs to it. The missing step is physical branch capacity:
enough source states or lower-cost return moves to make the union of
inverse domains cover the admissible targets. The reusable design rule
is to distinguish the tests for a local tag, a legitimate marked
predecessor, and an exact absolute image.
