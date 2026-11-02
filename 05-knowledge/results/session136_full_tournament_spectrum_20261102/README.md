# Session136: structure of the full tournament spectrum

**Checkpoint 1. The full claim that H takes every positive odd value except7 and21 remains OPEN.**

## Local universality, exact holes

The strong one-apex tournament family realizes every parity-compatible residue modulo every modulus, simultaneously for compatible CRT requirements. Every such residue contains arbitrarily large counts. Its profinite closure, and therefore that of the strong and unrestricted spectra, is the odd coset `1+2 Zhat`.

The construction uses actual positive binary append words, a modular commutator, and a dyadic reset that preserves the odd part. No inverse graph operation is assumed. The explicit word bound is polynomial in the numeric modulus, not its bit length.

[Proof and scope](local_closure/PROOF.md) · [Standalone constructor and checks](local_closure/check.py).

This cannot imply exact completeness: the same family misses13 at every order. The unrestricted exclusions7/21 are inherited. No finite collection of congruence tests on H alone can certify an odd target missing; magnitude bounds or additional graph data remain essential.

## A verified warning about matrix lenses

The inherited GF(2)-similarity lens already fails from order6. Our [small literal certificate](LENS_BOUNDARY.normal.json) exhibits two adjacency matrices with rank powers6,5,4,3,2,1,0, hence the same single nilpotent Jordan block, but Hamiltonian-path counts1 and41. This rechecks a known information-loss boundary rather than claiming a new matrix-classification theorem. [Checker](lens_boundary.py).

## Provenance and scope

The conjecture is recorded by Chan–Pak in [Computational complexity of counting coincidences, Section6.5](https://arxiv.org/html/2308.10214v2#S6.SS5), crediting bof and Royle's2016 question. This session makes no literature-priority claim. Existing one-apex transfer laws and exact inverse are credited inside the proof.

[Independent proof and execution audit](independent/LOCAL_CLOSURE.md) · [Concept board](BOARD.md). Other structural lanes are being audited separately and are not dependencies of this checkpoint.
