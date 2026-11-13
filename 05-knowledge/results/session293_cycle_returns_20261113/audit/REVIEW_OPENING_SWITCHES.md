# Independent review: renewable openings, activation, and two-switch carry

**PASS, with test-harness correction applied before final replay.**
No all-target coverage or unrestricted return theorem is inferred.

## Neutral opening and renewal

The classification in opening/PROOF.md §1 is correct under its four-distinct-
old-vertices hypothesis. A neutral birth is at an SCC cut. The inverse opening
support has one cycle on the cut's prefix and one path beginning in its suffix;
neither can traverse the other side. A spanning cycle makes that prefix strong,
so it is precisely the first SCC. This also explains why the arc must actually
lie on a Hamilton cycle; strong connectivity alone does not supply that arc-
specific property.

The all-length conveyor returns the full H and all old-old endpoint counts.
A Hamilton path must traverse the SCCs in order and every singleton in the
middle chain. The old endpoint matrix is unchanged, and middle seam vertices
have zero global starts and ends.

The renewed b-port is a DIFFERENT physical vertex. Its first source suffix
counts h_bt(B), while subsequent suffixes count e_t(B). The condition
`e_t(B)=h_bt(B)` is necessary for the nonzero conditioned sector itself to
remain constant, and is sufficient. A source b is a valid sufficient
physical condition. Fixed names cannot span the preceding inserted singleton
on a second neutral opening. These distinctions are all stated correctly.

Our direct permutation hostile uses B rows `(8,1,3,6)`, b=0,t=1. The first
opening sector is1 and the renewed sector2; the full count stays constant.
This independently validates the need for the endpoint-column guard.

## Feedback activation and target budget

Both proofs of the exact feedback formula are sound. In the direct path
proof, a path using the unique backedge has exactly the two displayed A
segments, two B segments and two increasing middle subsequences. Every middle
vertex independently chooses its side of the backedge. The parts beginning
a and ending b must be nonempty; the other pieces can be empty. This gives
`2^m L_a R_b` without an unnoticed multiplicity.

For m>=1 every old path traverses the nonempty seam, so none uses the removed
cross arc. At m=0 the lost `e_a(A)s_b(B)` term is indispensable. The proof
retains this singleton/empty boundary correctly.

For a fixed baseline h, an activated target must have c>h. The positive
integer response product then forces `m<=v_2(c-h)`. This is sharp across the
allowed class, e.g. singleton A and B have L=R=1. It is not a restriction on
all tournament constructions or on a baseline that is allowed to change.

## Two-switch identity and carry classification

The inherited one-apex odd-cycle formula gives the displayed invoice. Directly
subtracting a single adjacent10→01 swap yields its `-1` internal term and
the two geometric tails. For the two disjoint swaps, the cross-term is
`-2^(m-1)` in H/2, canceling the remaining edge-end contribution exactly.
The full difference is therefore twice the binary integer difference B−C.

The carry argument is complete. A nonzero carry-stop j away from the midpoint
contradicts its reflected digit; at the midpoint for odd m>=5, the forced zero
pair at positions1,m−2 contradicts the earlier anti-pair relation. The only
extra positive words occur at m=3. Bit complement supplies the negative case.
Zero is exactly the anti-palindrome language.

The converse-core explanation for zero repairs is correct: after deleting
X's final sink and Y's initial source, the words are `10w0` and `1w10`.
Reversing and complementing the first gives
`1 complement(reverse(w)) 10`, equal to `1w10` for anti-palindromic w.
An initial reviewer concern about this identity was a mistaken suffix
complement and was withdrawn before this report; no source correction was
needed.

The explicit dyadic values, affine rearming laws, and O(log X) sparse positive
repair ray follow exactly. Swapping the same two pairs again is the inverse,
not another +2 step. The text does not mistake the infinite family for a
renewable arbitrary-target selector.

## Independent computation

`verify_opening_switches_independent.py` imports no primary lane code. Its
subset path DP independently checks actual tournament counts; its complete
word check covers all32,764 words of lengths2–14. It performs:

- 60 word-pair comparisons against direct Hamilton path DP (lengths2–5);
- 2,523 feedback cases with arbitrary outer tournaments of orders1–3 and
  seam lengths0–2;
- full old-old endpoint-matrix checks along a neutral conveyor;
- direct permutation verification of the sector1→2 renewal hostile;
- complete zero and ±1 carry-language checks through length14.

All68,575 explicit checks pass. Normal and optimized output are byte-identical.
The primary opening checker initially used bare assertions, which disappear
under python -O; its author replaced them with explicit raises before the
final replay. This was a reproducibility issue, not a mathematical failure.

## Arbitrary-core slot compiler and exact repair tariff

The all-middle compiler in opening/MULTI_FEEDBACK.md passes review after two
scope clarifications: A,B are explicitly nonempty, and the k-specific
logarithmic estimate is asserted only for m>=k. The both-empty exterior needs
its own slot convention and is no longer included accidentally.

For a fixed exterior order, deleting middle runs leaves precisely the
stated possible gaps. Every missing A→B adjacency must be repaired; no
other missing adjacency can be repaired by a uniformly oriented middle.
Choosing occupied legal slots and an ordered path cover is a bijection.
The constant coefficient must be retained separately for the empty middle.
The upper slot bound and nonnegative coefficients follow immediately from
the two-side word. Rédei positivity on each block of each ordered partition
gives the claimed transitive-middle lower envelope for arbitrary middle
structure, not just transitive middle lengths.

The single-backedge formula is correct with nonempty M: the avoiding sector
has h_A h_B H(M), and the using sector has L_a R_b times the two-possibly-empty
path partition count `2H(M)+C_2(M)`. Both coefficients are positive, so equality
in its transitive lower envelope indeed forces H(M)=1, hence transitivity.

For a general feedback set E, the proposed two-slot skeleton contains exactly
one B→A step and exactly two legal slots, including boundary slots when an
outer prefix/suffix is empty. It uniquely recovers both exterior paths and
the selected edge. Thus `d_2>=|E|h_A h_B` and `d_1>=h_A h_B` are valid.
The stronger OCF lower bound also has no double counting: its one cross-module
cycle contains exactly one A vertex and one B vertex, hence recovers its
selected reversed edge and ordered odd middle path. All remaining cycles
are confined to the three residual interiors. The resulting sharp universal
bound `H>=h_A h_B+2^n W_E>=1+|E|2^n` is correct. The exponent n counts middle
vertices; this is not a bound for unconstrained tournament operations.

The two seven-vertex exterior polynomials in switches/PROOF.md §5 have been
independently reconstructed by explicit exterior permutations and inserted
run symbols, checking actual adjacencies rather than using the primary
binomial formula. They are exactly `[15,31,15]` and `[17,33,15]`.
For nonempty C, coefficient k multiplies the ORDERED k-path cover count,
so the identities `H_X=h+15Z2`, `H_Y=3h+15Z2` follow for arbitrary C.
The repair difference is exactly2h, not2; transitivity is exactly the
condition retaining unit repair. Multiplicativity of Z2 under ordinal sum
follows from the induced-subset product formula, so singleton insertion
fixes h and doubles Z2 as claimed.

`verify_slots_independent.py` checks 21 middle cores by direct full Hamilton
path DP, 296 arbitrary-exterior/middle cases (including empty seams), all
slot/OCF lower bounds in that universe, and the exact exterior polynomials.
All1,228 explicit checks pass in identical normal/optimized runs. These
finite checks support the all-order bijections; they are not the proof of
arbitrary-middle scope by extrapolation.
