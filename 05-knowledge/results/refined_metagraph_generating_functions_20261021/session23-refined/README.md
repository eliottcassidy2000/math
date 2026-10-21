# Marked fixed-degree series and exact recurrence orders

2026-10-21, Session23. PROVED from inherited finite-atom and SCC laws;
independent mathematical and computational review supplied separately.
No complete strong-center census or metagraph connectivity theorem is claimed.

New extensions:

- Every fixed fundamental black degree has a rational multivariate series
  marking order, SCC count, finitely many prime valuations of automorphism
  order, presentation count h, and blue count b. Converse half-pairing
  squares every additive marker, not just the order variable.
- Every nonempty fixed-h ordinary stratum and every nonempty fixed-deficit
  self-converse stratum has an exact reduced denominator, determined by the
  maximal feasible number of nonrigid SCCs. The argument uses positivity at
  a real pole and irreducibility of its polynomial carrier.
- At degree 3^r the exact reduced denominator is
  (1-z-z^3)^(r+1)(1-z^2-z^6)^e, with e=floor(r/2)+1 except e=2 at r=1.
  Thus the minimal eventual recurrence orders begin
  9,18,21,24,33,36,45,48,... . The cases 9 and 18 were already known;
  the all-r noncancellation law is the extension.
- Marking raw automorphism order as an exponent is not D-finite in that
  marker over Q(z), already on h=1, whereas prime-valuation marking is
  rational. The proof is a coefficient-gap obstruction.

Files:

- PROOF.md: complete definitions, proofs, inherited attribution, hostile
  examples, and the source-to-target connection contract.
- check.py: supporting exact controls, using Python and SymPy.
- normal.json and optimized.json: identical normal/-O receipts.
- source_pins.json: repository-relative inherited sources with LF SHA256
  pins and the source commit. No live web claim or external novelty claim.
- manifest.json: this packet's file hashes and reproduction information.

Reproduce from this directory:

    python3 check.py > replay.json
    python3 -O check.py > replay-optimized.json
    cmp replay.json normal.json
    cmp replay-optimized.json optimized.json

The main checker uses a declared synthetic involutive atom alphabet, with
actual K1,C3,R5,S4 signatures plus artificial adversarial symbols. It compares
16,601 literal words through order13 with independent marked convolution,
retaining six statistics, then checks denominators symbolically. Its output
is not a tournament census and does not promote invented atom signatures.

Novelty audit found the older rigid bivariate series and exact recurrence
orders at degrees1 and3; they are explicitly inherited. The general
multivariate and prime-power denominator statements are new repository
extensions within the bounded search performed, not claims of global
literature priority.
