# Circulant Hadamard {1, 4} against {7, 21}; circulant tournaments modulo p; the packing-polynomial lattices at every level

Session klein-2026-11-16-circ (2026-11-16), on the owner's directive: how the circulant Hadamard classification (orders
1 and 4 only) relates abstractly to the holes {7, 21}, and the packing-polynomial lattices. Concurrent (other machine,
same words): `packing_hadamard_flatness_20261117` — the identity `P(q) = G(2q, q)`, the Hadamard defect inside `G`, the
orbit-phase reduction mod a prime-order automorphism, the `(x − 2y)` invisibility, and the typing "allowed orders versus
excluded counts"; read first, cited, not repeated. Canon: THM-5270 (new).

## Status

- **ANALOGY, typed** ([HADAMARD.md](HADAMARD.md)). Three finite exceptional sets from two kinds of constraint: circulant
  Hadamard survivors `{1, 4}` (symmetry: flat autocorrelation), prime-count survivors among vertex-transitive tournaments
  `{3}` (symmetry: `n | H`, floors), and the excluded counts `{7, 21}` (growth: floors and a finite kernel). The honest
  tournament analogue of "circulant Hadamard only at 1, 4" is "a cyclic tournament has a prime count only at order 3";
  the holes are the single-orbit counts that orders 7 and 21 cannot have.
- **PROVED** (HADAMARD.md §2, THM-5270 (B)). For a circulant tournament of prime order `p`:
  `G_T ≡ y^p + ((p−1)/2)·x (mod p)` coefficientwise, hence `p | P_q` for every `q`; checked on every connection set of
  orders 3–13 (908 checks). Counts: 3; 15; 175, 189; 92411 … 95095; 3,669,497 … 3,711,175; `a_{p−2}` and `a_{p−4}` are
  constant across each order (the concurrent flatness), the defect first shows at `a_{p−6}`.
- **PROVED** ([LATTICES.md](LATTICES.md) §1, THM-5270 (A)). The delayed-choice lattice at every level:
  `L_n^{(k+1)} = (k+1)^{m_0} Π_{c=1}^{k} ((k+1)² − c²) ℤ` with `m_0 = 1` (odd `n`) or `2` (even), by the band-vector ℤ-basis
  of the Vandermonde kernel: `120/360`, `5040/20160`, `362880/1814400` for levels 3, 4, 5. Checked to `n = 14`.
- **FINITE-EXACT** (LATTICES.md §2). The diagonal forgets most of `G`: at order 8, 1728 diagonals for 3007 bivariate
  polynomials, and 4503 of 6880 classes share their diagonal with a class of different cycle structure (none at order
  5, 16 at order 6, 196 at order 7).
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md) — 6 predictions, 4 confirmed, 2 refuted (the diagonal is
  coarser than guessed).

## Files

- `scripts/`: `circulant.py`, `lattices.py`, `lattices2.py`, `fibres.py`, `fibres2.py`, `nearjoin_check.py` (helper block).
- `receipts/`: `circulant.out`, `lattices.out`, `lattices2.out`, `fibres.out`, `fibres2.out`, each with a byte-identical
  `.optimized.out`.
- `prediction_circ.txt`, `SHA256SUMS`. Input: `profiles8.npz` from the seams packet.

Reproduce: `python3 circulant.py` (6 min: order 13 needs five subset convolutions per connection set), `python3 lattices2.py`
(sympy), `python3 fibres2.py` (1 min), `python3 lattices.py`, `python3 fibres.py`.
