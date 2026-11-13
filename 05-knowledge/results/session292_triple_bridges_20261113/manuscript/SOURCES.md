# Source pins and scope

## Supplied manuscript

- Read `unified.tex` at main commit `8a79f4a20e21afe8cb27de1283d087b79539537b`.
- [Pinned source](https://github.com/startigerjln/Nielsen-and-Semita-Riemann-Proof-in-Lean4/blob/8a79f4a20e21afe8cb27de1283d087b79539537b/unified.tex).
- Retrieved UTF-8 SHA256: `91a7917a5681b1a4c54ec280261d428043dcf2e16d6f27bfc9cbae347016f9ad`.
- Source was retrieved twice, once pinned, and byte equality checked. We do not
  republish its full text here. The paper calls `appendix_body.tex` at its end;
  no audit of the imported Lean appendix, axioms, or root import was performed.
- Inspired mechanisms: reconstruction consistency; local/global limits;
  complexity of continuation; phase observables. Claimed RH/independence results
  are not dependencies of this packet.

## User's totient inequality

**CITED, not newly proved:** J.-L. Nicolas, *Small values of the Euler function
and the Riemann hypothesis*, Acta Arithmetica155(3),2012,311–321,
[doi:10.4064/aa155-3-7](https://doi.org/10.4064/aa155-3-7),
[primary full text](https://www.impan.pl/shop/en/publication/transaction/download/product/83792).
Theorem1.1 equation(1.5) and Corollary1.1 give exactly the supplied threshold
`N_120569`, whose largest prime is1591883, and constant
`exp(gamma)(4+gamma-log(4pi))`. The threshold is a primorial, not the integer120569
or the product of primes at most120569. The equivalence itself is not RH's proof.
Our fixed-support envelope is elementary calculus and uses no RH hypothesis.

## Repo lineage read

- THM-2961, `01-canon/theorems/THM-2961-ordered-continuation-poset-and-finite-group-context-tax.md`:
  full-context equality is the stable quotient; scalar labels need not suffice.
- `05-knowledge/results/collatz_procgen_20260923_trunk_rh.md`: genuine Iwasawa
  chart and precise loss of arithmetic/dynamical targets.
- `05-knowledge/results/denominator_fibre_recursion_20261004.md`: denominator
  ideals versus numerator phase; retained phase enables recursive lifts.
- `05-knowledge/results/difference_visibility_norms_20261004.md`: exact totient
  fibers and sampling-specific Euler products.
- `05-knowledge/results/definability_lenses_spectrum_diophantine_carmichael_20261105.md`:
  finite-rank path series and missing constructive numeration step; finite
  multiplicity claims retain the later scope correction.
- Session291 cluster, ports, Witt notes: exact endpoint algebra and eight-mode
  profile dictionary. Formal sums are not physical addition gadgets.

## Incoming proof-language boundary

Read `PROOFS.md`, E1–E2, in commit
`1ece157001637120f4260478719aef48c2e6d9e9` of
`05-knowledge/results/core_descent_endpoints_20261112`. Its complete-observer
collision and selected-tail separation are inherited exact results, not
recomputed in this lane. The immutable link is given beside the application.
