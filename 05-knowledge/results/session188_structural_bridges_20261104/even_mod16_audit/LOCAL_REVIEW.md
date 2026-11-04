# Independent formal run-table review: cyclic even-child modulo sixteen

PASS for the local proof and complete formal run certificate. This review pins source proof d759e62510d0469dd9e864a46117f2bccc02e8256892162664a021821f07a78c and the other final hashes in local.normal.json. The separate physical/profile audit reviews the producer's large named bank and original-module controls; this artifact deliberately owns the independent algebraic route.

## Completeness and proof review

The modulo-eight normalization gives q2/H in {6,14} and q3/H in {4,12}. Forgetting ordered path-slot labels makes q_k divisible by k!, while odd induced Hamilton counts identify the parity of q4/4! with S(N,4). The Stirling recurrence makes that parity one for every even N>=4. Thus normalized q4=8, normalized q5 is 0 or 8, and every coordinate k>=6 vanishes modulo16. The displayed eight vectors are an upper domain; their realizability is neither needed nor assumed. The normalization retains the odd H scalar, and multilinearity contributes one such scalar per child.

Every Hamilton word with uniform child-port incidences decomposes uniquely into maximal child runs and the two distinct port positions. Each child contributes one complete ordered-cover coordinate indexed by its number of runs. Consequently any word using at least six runs of a child contributes zero modulo16, and the retained finite language is complete. This validates the all-order reduction, not just a search at bounded physical size.

All512 formal triples have nonzero residue banks. More strongly the same physical block query a=A,c=A union B has contrast8 for every triple. Multiplication by the three odd normalization scalars preserves8. Combined with the separately audited five-type modulo-eight theorem, this closes local C3 quotient-converse separation for all nonempty child orders and arbitrary unequal profiles. The two-query root-family consequence correctly uses exact converse-complement duality and monotonicity to select the baseline-negative member. This is source-relative detection, not unrestricted reconstruction or injectivity of all internal orbit outputs.

No arbitrary ambient lift is inferred from the local residue certificate. A separate all-cover certificate is being reviewed for that purpose.

## Independent algorithm

`automaton.py` does not import producer code or use its rational formula. Its state is the triple of run multiplicities (each between zero and five), a two-bit used-port mask, and the last label. There are five labels: the three old children and the actual distinct ports. Old labels cannot repeat consecutively; transitions follow the actual quotient/port directed arcs. Each port can occur at most once. Starting at every possible label, the dynamic program counts each legal run word exactly once. A completed word must contain both ports. At most17 labels can occur, so the finite layered program is exhaustive. Missing old children contribute q0=0 when profiles are applied.

For all64 physical mask pairs and both quotient directions, this independently constructs the complete integer run coefficients. Multiplying the contrast series by (1-XYZ)^3 recovers every displayed numerator coefficient and checks all additional coordinates through multiplicity five are zero. The rational numerator's global degree bound is justified by its explicit finite matrix expression, not inferred from this bounded diagnostic.

`audit_local.py` evaluates both entire integer source/target functionals for every512 formal triple. It compares all32,768 exact contrast entries against the retained source table, with no hash-only grouping or numerical tolerance. It retains all65,536 exact source/target values and full numerator vectors. The residue histogram is exactly64 banks with12 nonzeros,256 with24, and192 with28. The six universal queries are (1,3),(1,5),(2,3),(2,6),(4,5),(4,6). There are eight distinct full residue vectors. As an additional diagnostic, only the three q2 bits affect those C3 contrast vectors.

Normal and optimized executions agree byte-for-byte. Explicit exceptions enforce every check under `-O`.

## Replay and exact allowlist

Run `python audit_local.py > local.normal.json` and `python -O audit_local.py > local.optimized.json`, then compare the outputs.

Local audit allowlist: LOCAL_REVIEW.md, automaton.py, audit_local.py, local.normal.json, local.optimized.json, LOCAL_SHA256SUMS.

The neighboring cover-automaton and physical-cover files are for the separate all-cover ambient review and are not part of this local audit checkpoint.
