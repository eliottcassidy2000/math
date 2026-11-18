# Self-converse cores as exact computational probes (2026-11-18)

* [Exact blue oracle reduction](encoder/PROOF.md): recover H and the
  two-slot path-cover count of any tournament from two strong
  self-converse Hamilton counts. A same-order version changes one
  central chord in their common Hamilton frame; it is not a native
  line in the tiling metagraph.
* [Folded DP](fold/PROOF.md): compute fixed Hamilton paths using one
  signed endpoint per reversal pair. Compute the entire fixed
  odd-cycle-packing polynomial, and identify the nonfixed orbit
  overlap that carries the even-order mod-8 digit. For the full
  exact H use a reflected middle split that retains nonfixed terms.
* [A physical 7-adic selector](build/PROOF.md): a one-center cyclic
  mirror is stuck at even 7-adic valuation. A two-center mirror
  escapes that obstruction, hits any prescribed multiple of seven
  modulo 7^k from an explicit admissible class of source profiles,
  and sits on one mod-8 sheet. [Four-center CRT](build/CRT.md) chooses all odd residues modulo 8·7^k using a fixed transitive source bank. Since these hosts are strong and have order growing linearly in the width, Busch gives only O(log X) distinct counts <= X for this fixed bank. Includes a literal order-16 7p
  witness. This is not absolute target selection.
* [Absolute misses and a controlled exit](encoder/ABSOLUTE.md):
  the 63 obstruction and the mirrored SCC factorization are inherited;
  the small blue chart independently reproduces the complete list of
  self-converse holes below 75. The previous two-contact power family
  actually starts from an identified blue class and reaches a chiral
  class, with 65 to 63 as a concrete control.
* [A structural one-apex limit](encoder/ONE_APEX.md): strong
  self-converse sources with a transitive deletion are completely
  described and have only on the order of sqrt(X) distinct counts
  through X. [A closed growth question](encoder/RAY.md): a blue
  center of count 567 attains the previously only bounded
  self-converse growth degree along the even-exponent 7 times 3 ray.

The distinction is sharp: all Hamilton counts reduce to the blue
sector as exact *computational reads*, but not as identical witness
values. A constructive Chan–Pak approach must choose a center and
either decode a physical output or cross out of the blue sector.
The all-target proof is open.
