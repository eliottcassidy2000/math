# Reproduction, independent checks, and scope

All commands below run from this packet's directory. Python3 and its standard
library suffice for the proof controls; the large square-chart search also
uses g++ with C++17. No network access is needed for the saved numerical-source
comparisons. An explicit --refresh option is required to replace those sources.

## Main verification

    python -B replay.py
    python -B seal.py --check

The replay runs the newly integrated proof controls normally and with Python
optimization enabled. It requires exact equality of their normalized receipt
files, and stops on any failed explicit check. It does not rely on assertions
that disappear under optimization. REPLAY.json records the selected commands
and matching receipt hashes. The manifest records every packet artifact and
selected inherited proof/code dependencies, with LF-normalized text hashes.

The checks are finite controls of stated mathematical proofs. Passing a finite
prefix does not establish an unproved all-order assertion. The complete small
automata, complete cut universes and factor-pair certificates have their own
finite exhaustiveness arguments in the corresponding notes.

## Source-return and operation dictionaries

    python -B reset/uniform_return.py
    python -B reset/composition.py
    python -B reset/ready_rebase.py --inverse
    python -B reset/iterable_return.py
    python -B reset/atlas_returns.py
    python -B reset/atlas_two_birth.py
    python -B reset/classify_modules.py

The last command independently rechecks the complete ordinary census through
eight vertices before classifying all fixed transitive-module source templates.
The ready inverse uses the inherited complete one-/two-start grammar and its
proved target-dependent order bound. The atlas universes are exactly9216
prescribed-sink cuts and270336 arbitrary two-birth cut pairs; their receipts
retain every base scalar candidate and its full module response. Any *_probe
file is exploratory unless a note explicitly promotes its scoped result.

Independent root-lane graph constructions and counting implementations:

    python -B audit/retained_reset.py
    python -B audit/ready_rebase.py
    python -B audit/atlas_instructions.py

These use separate vertex names and endpoint-subset recurrence, preserve the
literal old induced graph, and verify full response coordinates. The smallest
graphs are also checked by literal Hamilton-word enumeration.

## Classification and growth

    python -B defect/two_inport.py
    python -B defect/multi_inport.py
    python -B defect/cycle_residual.py
    python -B defect/pentagon_inflation.py
    python -B arithmetic/dyadic_defect.py
    python -B metagraph/minimal_defect.py
    python -B metagraph/response_families.py
    python -B metagraph/classify_cap9.py --output metagraph/cap9.json
    python -B metagraph/audit_operations.py

The strong defect2–6 enumeration is complete at every order, with automaton
termination and exact productive languages recorded. The later canonical
prefix theorem supplies a direct finite grammar for arbitrary fixed branching
defect; it does not pretend the huge tail bank has been enumerated for every
budget. Its universal proof is reviewed in DYADIC_DEFECT.md and the independent
review records. The weak/program theorem uses positivity, an explicit logarithmic
grid and root-SCC composition, with all intermediate defect hypotheses retained.

The root controls audit/template_controls.py and audit/weak_factor.py separately
check2008 strong branching roots and3208 clone insertions, and5405 rooted SCC
factorizations. Arithmetic and metagraph lanes independently review the same
proof mechanisms. The classified strong/weak order sequences and automorphism
refinements are not unrooted A000568 terms.

## Arithmetic ranges and exact intersections

From arithmetic/:

    python -B maxima.py --output maxima_receipt.json
    python -B audit_maxima.py --output maxima_audit.json
    python -B quadratic_primes.py --output quadratic_receipt.json
    python -B compare_sources.py
    python -B chart_inverse.py
    python -B tail_squares.py
    python -B minimum_squares.py
    python -B carry_ray_probe.py
    python -B prime_source_filter.py
    python -B near_square_charts.py

The maxima producer and independent audit both cover every target2..1,000,000.
The quadratic sieve covers every root through1,000,000, with a proper-divisor
witness for each rejected candidate and full-factorization Lucas certificates
for the54110 survivors. Its complete indexed table is committed. Add
--table output.txt to either producer for a fresh complete export. The factor
proofs and the A135244/A002496 extremal equivalence are not inferred from those
finite comparisons.

TAIL_SQUARES.md proves its infinite intersection by a complete source inventory,
periodic congruence obstructions and exhaustive odd factor pairs, keeping the
Fermat primality predicate explicit. The chart inverse is complete in its
normalized grammar; a supplied max_short option is explicitly only a restricted
probe. No OEIS submission, new world record, unrestricted spectrum coverage or
prime-infinitude claim is made.

## Repository checks and preservation

The repository documentation checker and git diff --check pass at close-out.
The session used an isolated worktree and committed only its packet and explicit
navigation entries. Concurrent incoming birth-family and H13 carry work was
read, credited and integrated through clean rebases. No historical log was used
as current truth, no RESERVED statement was promoted, and no theorem number was
claimed. The final commit and remote verification are reported with delivery.
