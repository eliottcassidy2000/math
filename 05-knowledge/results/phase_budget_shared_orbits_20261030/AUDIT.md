# Audit and reproduction

**Status discipline:** proofs are elementary scoped arguments in this packet,
or explicitly use the preceding audited finite certificate. Computational
checks have the universes below. Review was by the same agent; independent
algorithms are not an independent human or subagent audit. No formal Lean
verification, universal Collatz proof, A000568 complexity improvement, plane
chromatic determination or planar Jacobian solution is claimed.

Run from the repository root, with D denoting
`05-knowledge/results/phase_budget_shared_orbits_20261030`:

```text
python3 D/collatz_audit.py
python3 D/audit_potential.py
python3 D/tournament_audit.py
python3 D/palette_audit.py
python3 D/jacobian_audit.py
python3 D/ledger_audit.py
```

Replace D by the displayed directory. C++17, Python3, NumPy and SymPy are
the runtime dependencies. The five direct audit scripts also pass under
`python3 -O`; their new checks use explicit exceptions rather than assertions.
The original full leader-ledger battery was rerun normally after its repair
and completed, including40 primitive covering controls. Its older internal
assertions remain as inherited tests.

## What was checked

* **Collatz word:**255 valuation letters, integer clock271. Literal big-integer
  gcds give charge154. Exact rational arithmetic gives old deficit769/5000 and
  first-arrival deficit-60951/5000. The two calculations use separate descriptions
  of the return events. The compatible address lift agrees with every epoch.
* **Word controls:** all32767 nonempty positive compositions of clock<=15,
  plus1250 U 1^L U controls for0<=L<1250. The bounded-controller proof covers
  every L, using thirteen explicit phase divisors. The full old budget remains
  refuted, not merely unproved.
* **Potentials:** eight exported integer vectors,2*5^E states each. A different
  NumPy implementation checks38,085,920 inequalities in total, every seed and
  the unbounded-letter tail. Byte hashes and extrema are in
  [potential_audit.json](potential_audit.json). Binaries are generated in a
  temporary directory and are not needed in Git. The producer run at E=9 exits
  with code2 after printing a positive path; that expected diagnostic is not a
  successful certificate. The independently verified printed path is(2) in
  [COLLATZ.md](COLLATZ.md).
* **Tournament controls:** eight composite families through57 vertices. Signed
  DSU, independent binary elimination and explicit common-bank enumeration agree.
  Each supplied selector is checked as an actual Hamilton path in a tournament
  fixed by both parents, with endpoints retained. Six weighted intersections at
  orders9 and15 compare factorized signed sums with literal source summation.
  Inconsistent signs and arc-reversing permutations are negative controls.
* **Palette controls:** all45 blocker types acting on90 options; exact-cover
  enumeration produces precisely nine tight covers. An independent sweep of
  all2^15 simple graphs checks the six-star sublemma. Fifteen independent blockers
  explicitly show why method obstruction is weaker than six-colour forcing.
* **Jacobian controls:**140 torus and1897 finite-cyclic projections of exact
  polynomial automorphisms. Symbolic differentiation verifies the success/failure
  criterion and explicit determinant defects. Characteristic zero, nonzero
  coefficients, r,s>=2 and the finite-cyclic no-alias guard are retained.
* **722 ledger repair:** eight valid primitive/nonprimitive controls and five
  normalized all-odd hostiles. {2,6} was admitted by the original guard and failed
  its later maximum check; the new guard rejects it both normally and under -O.

## Proof review boundaries

The potential induction has the direction V(target)>=V(source)+debt; seed
upper bounds and V<=0 therefore bound every finite accumulated debt. The tail
uses min V and max V<=0, not an empirical absence of large-letter violations.
The active quotient is sufficient for the *old truncated budget*, not future
return recognition. E<=8 cannot be extrapolated past the explicit E=9 hostile.

The controller proof uses C(1,VU), not C(1,UV), and includes the empty and full
suffixes. The family consequence uses the already proved sharp E=5 result;
no weaker finite scan is substituted for it.

The tournament count is affine consistency plus free components. Homogeneous
orbit counts alone can miss inconsistent offsets. Substitution injectivity
requires disjoint nonempty named blocks; port identifications require their
own guards. A selector can exclude a desired source, so it is not a universal
cover. Pairwise orbit information also does not restore Session113's collective
transporter character.

The palette proof establishes exact coverage of the ordered-pair method, then
uses a separate circle-intersection condition. It does not infer a boundary
relation from graph minors, global colour count or the number of cross-edges.
The Jacobian formula retains every resonant monomial; distinct exponent pairs
justify the no-cancellation step. All original two-shear maps remain invertible.

Canonical corrections are recorded in `01-canon/MISTAKES.md`. The predecessor
budget status and results index are updated; earlier finite receipts remain
historical evidence with their original, bounded universes.
