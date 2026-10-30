# Rectangle-only counting: remove redundant choices before modular observation

2026-10-30. **PROVED / VERIFIED Lean overlay; FINITE-EXACT native controls;
CITED external developments; compiler and complexity improvement OPEN.**
Same-agent mathematical review, with Lean kernel checking of the specified
overlay and separate native finite controls. No claim of literature priority.

The pinned Rudich--Rousseau factoring engine enumerates choices on all
`G^4` grid quadruples, although its child construction reads only the
`R=binom(G,2)^2` actual rectangles. We prove that the entire child clique
table is constant on each discarded-coordinate fibre. Enumerating one
canonical representative per fibre gives a new all-depth semantic adapter:

    branches per step: q^(G^4)  ->  q^R,       q=ceil(lambda/3)+1;
    reduced value:     V_new(T)=2^E value(T);
    original value:    V_old(T)=P V_new(T);
    original leaves:   H_old=P H_new,
    P=product_levels q^(G^4-R).

The external rectangle identity is inherited. The contribution here is
the formally checked adapter to its actual child constructor, the
all-depth normalization and leaf recurrences, and the actual scheduled
parameter audit. The compiler is still the original one: no new TM2
running-time theorem or practical factoring benchmark is asserted.

This removes a specific obstruction found in the earlier
[scaled-count audit](../scaled_count_fft_carriers_20261030/SOURCE_AUDIT.md).
That audit correctly warned that the **original** engine has odd factors.
We do not divide them out of existing modular readings; we change the
producer to stop generating them. Session112's already-proved
[effective-modulus law](../session112_counting_interfaces_20261030/arithmetic/MODULAR_OBSERVERS.md)
then applies to the reduced, purely dyadic multiplier.

For any fixed odd modulus bank of lcm L, the old readings retain a count
modulo `L/gcd(L,P)`, whereas reduced readings retain it modulo L. The exact
gain is `gcd(L,P)`. At the first two positive uniform depth thresholds
of the actual `K0=48` schedule, the bank `3,5,...,39` gives:

| recursion depth | bits needed to write input bit length n | old effective modulus | new effective modulus | gain |
|---|---:|---:|---:|---:|
| 1 | 3,702 | 363,761,673,275 | 166,966,608,033,225 | 459 |
| 2 | 59,346 | 325,470,970,825 | 166,966,608,033,225 | 513 |

Here `n=(Kb(48,depth))^40` is itself the input bit length. These are
astronomical theoretical thresholds, not tested input integers or
practical speed measurements. The small bank illustrates exact
information loss; it is nowhere near sufficient to recover the enormous
counts at these thresholds. The all-size odd-bank theorem still needs a
bound on the unscaled count.

The core synthesis is precise: remove a coordinate only after proving
that every declared future consumer is independent of it. Equality of
child tables establishes that here. Equality of scalar summaries would
not. In tournaments the inherited 81/27 pair still has equal labelled
cyclic triples and unequal attachment behaviour; it passes the hostile
replay and cannot be merged by this rule.

The external circuit development points to a different operation:
Mishra--Santhanam's recent primary preprint uses independent parameters
and reparameterization in gate elimination. Incoming Session117 proves
that the provisional three-weight template cannot exceed coefficient4
just by retuning weights. These findings motivate context-aware
reparameterization, but neither a small semantic carrier nor its large
minimal state count establishes a general circuit lower bound. See the
typed transfers and concrete next obligations in [SYNTHESIS.md](SYNTHESIS.md).

Files:

- [PROOFS.md](PROOFS.md): mechanisms, quantifiers, recurrences and failures.
- [ReducedEngine.lean](ReducedEngine.lean): eight printed theorem axiom
  audits, all only `propext`, `Classical.choice`, `Quot.sound`.
- [SOURCES.md](SOURCES.md): attached PDF, primary sources, version and
  evidence-level distinctions; current Big O detail page was refreshed.
- [AUDIT.md](AUDIT.md), [certificate.json](certificate.json): 1,106,254
  local native checks, identical normal/optimized outputs, inherited81/27
  control separately reported, and the Lean receipt.
- [HANDOFF.md](HANDOFF.md): the next decisive compiler test.

This checkpoint does not solve factoring complexity beyond the imported
result, ETH, P versus NP, LRC(14), universal Collatz, or unrestricted
A000568 enumeration. It supplies a checked transformation at a specific
algorithmic interface and records what an end-to-end improvement requires.
