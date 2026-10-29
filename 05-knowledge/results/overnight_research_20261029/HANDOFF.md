# Overnight paired arithmetic research handoff

Status: active research on 2026-10-29. This file is the continuation point
for the eight hourly follow-ups requested by the user. Preserve proof
status and resume from completed checkpoints rather than restarting.

Current worktree: `codex/session-sha_pairing_f1_20261029`, based on the
clean pushed prior packet `sha_marked_join_f1_20261029`.

## Completed checkpoints: inherit, do not restart

1. [TWO_PORTS](TWO_PORTS.md) proves the complete weighted classification
   H_s+H_(s n), s|gcd(d,e), and the bounded-weight Mobius/convolution
   recurrence. It exhausts all65,536 fixed native cuts and the separate
   simultaneous-C3 subcensus. At weight50 the packed square counts
   104,060,401 cuts. Minimum weights for H25,H27,H57 are4,4,6.
2. [MINIMAL_RATIONAL_PORTS](MINIMAL_RATIONAL_PORTS.md) proves that
   (det M,P^T M^-1 P) is the exact response quotient for retained-span
   determinant queries. It preserves re-marking and global sign reversal;
   integral cokernels and clock-restricted minimality remain separate.
3. [TERNARY_TERMINATION](TERNARY_TERMINATION.md) proves universal
   positive-integer termination with n(8/5)^Q. Incoming work gives stronger
   bounds, five increasing cylinders, periodic exclusions and an exact
   Cantor-cell recurrence. No positive aperiodic exception is excluded.
4. [FAMILY_HEIGHT](FAMILY_HEIGHT.md) proves degree24 for the record4239
   j-map and Theta(X^(1/12)) distinct curves at invariant height X. Its
   named section is exact. Generic rank17 is cited, not proved here.
5. [ARITHMETIC_INTEGRATION](ARITHMETIC_INTEGRATION.md) gives the explicit
   quadratic selector from15 matchings to6 K3,3 matchings, credits the
   incoming Frobenius/Hamming/fruit results and retains the12/13 clocks.
   [ALGORITHMS](ALGORITHMS.md) reads pinned primary techniques; the
   external manuscripts' full proof/build graphs have not been audited.

Six adjacent checkers pass normally and with Python -O, with identical
receipts. See [AUDIT](AUDIT.md) and mode_receipts.json. The earlier M4
Pfaffian sign was repaired in place and recorded in MISTAKES.md. The
August rank30 certificate and weighted81/27/A000568 controls were replayed.
After integrating incoming commit4d9c3ad1ce, the cyclic-kernel, exceptional
cell and all66 high-rank coverage checkers were also replayed successfully.

## Incoming results to inherit with their own scopes

Read the current Session97 README and the current ternary-integer packet,
then their targeted proofs. Recent Session97 additions include:

- controllers/EXCEPTIONAL_CANTOR_RECURRENCE.md: occupied-cell recurrence
  N_L=N_(L-1)+N_(L-6)+3N_(L-11), a dimension about0.441171 and an integer
  sparsity bound. This is not an emptiness or global Collatz theorem.
- supersingular/CYCLIC_KERNEL_RECURRENCE.md: all odd cycle lengths,
  reciprocal factors and odd ord_p(2), including Ogg intersection
  {7,23,31,47,71}. Retain the characteristic separately from cycle length.
- supersingular/ALL_M_AFFINE_PORT_RECONSTRUCTION.md: four-port affine
  reconstruction with its m=3 hostile; this is a different query language.
- ranks/phase2/README.md: all66 listed rank-at-least30 point banks have
  local certificates; two require F3 instead of F2. Three named sections
  and a common A3(2) lattice are proved for the newer families. The full
  classical rank17 certificate belongs to determinant948, not1152.
- fractions/FULL_IMAGE.md and DENSITY.md: an arithmetic rational-cover
  guard and one-third lift density, conditional on the explicitly cited
  source-generator theorem. Real positivity alone remains insufficient.

## Next bounded obligations, in priority order

1. The finite-clock unrestricted-integer subproblem is now closed in
   [FINITE_CLOCK_PORTS](FINITE_CLOCK_PORTS.md): orbit sums of pair probes
   are legal and recover every invariant Gram entry. A unipotent infinite
   clock is an all-word hostile. Next restrict the operation set to actual
   cross-block or native cuts; the orbit sums need not belong to it.
   A promising exact next case uses A=[[0,-1],[1,-1]] on each two-port
   block. Solving A^T B A=B gives B=[[a,b],[-a-b,a]], whose determinant is
   the Eisenstein norm a^2+ab+b^2. Retain phase(a,b), source product de
   and primitive content; norm alone changes under a marked translation.
2. Replace abstract paired-group realization by an actual source-fixed
   arithmetic map. For Sha finiteness one still needs finite prime support
   and bounded primary exponent (or another sufficient argument), with
   all Selmer levels compatible. A square-order count supplies neither.
3. Develop adaptive local-character certificates: select cheap primes
   and observation moduli while keeping root labels, generators, curve
   transport and torsion guards. Compare the two F2 failures repaired
   over F3 before assuming missing information requires full reconstruction.
4. For the cheap controller, retain the ORIGINAL source, full carry and
   both clocks through resets. Try a decisive arithmetic exclusion of
   aperiodic positive points in the exact failure attractor; if it fails,
   record the missing coordinate. Do not substitute Haar probability.
5. Recover section banks for the actual newer K3 fibrations. Combine
   their exact family-height exponent with proved section/intersection
   data; no unproved rank-tail law may be promoted to an upper bound.

Each continuation should produce a bounded proof, exact experiment or
documented obstruction, then commit and push. Check incoming corrections
before extending these results. Do not repeat a completed census merely
to keep the session busy.

Anchor: marked paired groups and arithmetic realization. Niche: the
positive-integer parser certificate. Wildcard: high-rank families and
the two fifteen-element structures. Sha finiteness, Collatz and LRC(14)
remain open.

Do not submit anything to external leaderboards or message other people.
Commit and push coherent verified checkpoints. A final morning synthesis
must say exactly which obligations remain open.
