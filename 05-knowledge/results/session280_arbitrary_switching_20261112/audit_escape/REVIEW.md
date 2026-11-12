# Independent review of receiver rearming and the six-count carry

**Session280, 2026-11-12. PASS WITH MINOR CORRECTIONS APPLIED.**
No Chan–Pak completeness or eventually gapless tail is established.

## Review scope

Read the full current `escape/PROOF.md`, inspected the Python and C++
producers, independently derived the receiving updates and carry formula,
and implemented a separate backward-subset recursion for the finite seed
profiles. The independent implementation does not import the lane evaluator
or its fusion matrices. It validates the actual tournament arcs before
counting.

The inherited rooted-fusion multiplication law is treated as the stated
prior dependency; this review verifies its new applications rather than
re-auditing its entire earlier proof.

## Passed mathematical claims

1. **Endpoint-cover injection.** Splitting a Hamilton path of the deletion
   immediately after vertex i injects its paths into the E_i covers. Thus
   E_i>=h. Summing over terminal vertices counts all nonempty first slots,
   giving sum_i E_i=Q-h. The decomposition of V by the incoming neighbors of
   the root is exact. Consequently the source/sink equality characterizations
   and positive gap bounds include the singleton boundary correctly.
2. **Arbitrary same-root sandwiches.** The displayed defect
   `2(Q_L-V_L)h_R-16Q_L(V_R-h_R)` follows directly by successive left and
   right receiving updates. If the right root is not a source, the negative
   term strictly dominates. Equality is exactly left sink/right source.
   Those contexts form an ordinal sum and multiply the +2 difference by
   h_L h_R. The generalized `0<a<=b` criterion is also valid. No claim about
   differing endpoint-specific contexts or internal edge surgery is made.
3. **Physical rearming.** Independent literal profiles agree for the original
   transported three inputs and the four-edge repair. The repaired bank has
   common `(h,V,e)=(15,441,15)`, although Q,U,d differ. These are precisely the
   coordinates consumed by a common left unit context in the count update.
   Both receiving updates preserve their equality. Thus arbitrary finite
   sequences of those contexts preserve the whole three-count bank.
4. **All-depth recurrence.** Repeated C1 gives `V'=4V-15`, `H'=H+2V`,
   hence `V_k=436*4^k+5` and
   `H_i(k)=269+2i+872*(4^k-1)/3+10k`. The vertex count is9+2k.
5. **Paid carry.** For even k, `V_k-21` is divisible by30, so
   `z=(V_k-21)/15` is even. The inherited recursion gives
   `g(2z+1)-g(2z)=-z`. Therefore adjacent physical contexts at P=2z,P+1
   shift the rearmed triplet by amounts differing exactly6. The six resulting
   scalar counts are consecutive odd integers. No baseline is subtracted
   from a physical graph or confused with a realizable count.
6. **Limit of the widened bank.** The contexts differ in one final incidence
   bit, changing their own V by `-2^(n-1)`. Thus the two output triplets have
   different receivers, by `-15*2^(n-1)`. Corresponding inputs also
   have equal h,Q. Therefore ANY common physical fusion sandwich multiplies
   this nonzero defect by Q_L h_R, not only left unit contexts. Independent neutral
   global padding only multiplies V by powers of2, while the two receivers'
   ratio is strictly between1 and2. These observations correctly prevent
   claiming a renewable six-way bank from the current result.

7. **General width-doubling criterion.** For an actual r-bank with common
   h,V,e=h, the conditions V>h+2r and V=h+2r modulo2h give a positive
   even z. The same carry gives exactly2r and hence2r consecutive odd
   outputs. The stated receiver defect and equal h,Q between corresponding
   outputs are exact. This is conditional on physical bank supply, not an
   existence theorem for such a bank at each r.

## Independent finite controls

`independent.py` and its optimized run produce identical `results.json`.
All gates use explicit exceptions, not removable Python assertions.

- Nine-vertex repaired profile:
  `(H,h,Q,U,V,d,e)=(273,15,768,441,441,9,15)`.
- The two lower bank members have counts269/271 and common h15,V441,e15.
- Verified the four repair arcs and the five direct corridor arcs.
- Verified the common Hamilton frame for the five-arc move. The separate
  zero-count repair frame also passes.
- Independently counted the six literal fifteen-vertex output graphs:
  **24875,24877,24879,24881,24883,24885**.
- Independently checked the exact carry at every even k from0 through100;
  this supports, but does not replace, the symbolic all-depth argument.
- Replayed the lane Python verifier in normal and optimized modes; outputs
  agree (49 literal sandwich controls,31 unit contexts, and further controls
  listed in the lane receipt).
- Compiled/replayed the C++ four-arc search at t2: exactly58,905 candidates,
  three matches. Independently recomputed all three matching profiles: each
  has the same full displayed profile as the selected repair. This search is
  finite, and no minimality or uniqueness claim is inferred.

## Corrections applied before publication

- The draft explanation incorrectly wrote `4^k=1 mod30` for even k. For
  positive even k the residue is16. The conclusion is unchanged because
  `436*4^k+5=21 mod30`, also at k0. Current prose states that correct identity.
- The first bank edge was called a single root arc. It reverses TWO
  root-incidence arcs, `{0,5}` and `{0,6}` in the transported labels. Current
  prose records both; the second direct corridor move has five arcs.
- The optional C++ search now has argument/dimension guards rather than
  silently allowing out-of-range fixed-array use. Its declared t2 replay
  and mathematical result were unchanged.

These are local exposition and defensive-interface fixes, not a failure of
the rearming or six-count theorem. Exact reviewed hashes are in MANIFEST.json.

## Reproduction

From this directory: `python3 independent.py` (also `python3 -O independent.py`).
From the sibling escape directory: `python3 verify.py`; and
`g++ -O3 -std=c++17 search_repair.cpp -o /tmp/repair280audit && /tmp/repair280audit 2 4`.

No global odd-target coverage, growing interval width, common receiver for
all six outputs, or inter-scale low-jump path is proved.
