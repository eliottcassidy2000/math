# Independent review of the size-six extinction refutation

2026-10-19. Status: PASS for the scoped claims.

Reviewed `GILBREATH_ALPHABET_THREE_AND_SIZE_SIX.md` in full. The simultaneous
alphabet-three A/B/C induction is valid; the independent geometric shield
proof supplies a second mechanism. The failure at maximum four is explicit.

The zero-run formula uses the first one in the original spatial right sea.
Before that one, the temporal boundary trace is also zero, and its first one
occurs at the same time. The shield makes subsequent bits irrelevant for
this particular eventual-failure predicate. It does not claim they are
irrelevant for the entire state trajectory.

For rotated coordinates, with defect at index L=F-1, x_s(k) is the physical
cell at index L-s and time s+k. The pure-sea value one update before first
contact is the stated Pascal transform. At s=L, age k is exactly the number
of appended zero seeds beyond the original defect. Thus the acceptance
condition in the finite automaton counts the actual first bad apex, with
no time-index or probability-factor mismatch.

An independent implementation builds the 23-state age-two automaton and a
separate eight-state recognizer for (0000)*11(ε|10*). Every reachable state
of their product has matching acceptance. This proves language equality at
all lengths, not just the enumerated control range. Direct word enumeration
through length16 additionally matches floor((L+2)/4). The count equals
floor((F+1)/4), and every accepted word has conditional failure probability
1/4. These words are disjoint from the F immediate-failure words. Hence
c6(F)>=floor((F+1)/4)/2 and the asymptotic ratio lower bound17/16 are valid.
They decisively refute the inherited logarithmic excess and ratio-to-one
claims. This is a random-sea theorem, not a statement about prime gaps.

A second independent script encodes each age vector as a base-four integer,
constructs all reachable transitions, finds the coaccessible states, and
propagates counts backward from acceptance. For each of ages2,4,6,8,10,12,14,
it agrees with the proposed rational generating function through degree
m+2h, where m is the coaccessible state count. The actual denominator has
degree at most m, numerator at most m-1; after combining with the proposed
denominator of degree2h+1, the difference numerator has degree at most
m+2h. Since the denominator has constant term1, these exact coefficient
checks prove equality identically. This validates the all-length claims.

The independently recomputed slope is23533/131072, giving ratio lower bound
285677/262144. The proof does not justify exchanging the infinite age sum
with a limit, so the complete asymptotic and all-age recurrence correctly
remain open. Likewise the seven odd-age tests are not an all-odd-age proof.

Reproduce from this directory:

    python3 verify_age_two_language.py > age_two_language.json
    python3 -O verify_age_two_language.py > age_two_language_optimized.json
    python3 verify_fixed_age_series.py > fixed_age_series.json
    python3 -O verify_fixed_age_series.py > fixed_age_series_optimized.json

Both normal/optimized pairs match byte-for-byte. These independent controls
use no imports from the producer scripts.
