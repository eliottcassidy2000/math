# Independent review of moving mark repairs and balanced children

**PASS.** The sibling counting agent checked the proof independently and
constructed the small graphs directly from named arcs. No producer
mathematical code is imported by the independent checker. The conclusions
remain scoped source constructions and obstructions; prime/7p coverage is
OPEN.

The cone C2 follows from a genuine injection. Each path of E-y has a last
in-neighbor of y; inserting y after it preserves directedness and creates a
path not starting at y. Restriction recovers the old path, so h-s>=q. With
the inherited nonnegative even defect, this gives A>=1, A+B>=2, A>=B and odd B.
The assumption that y has an in-neighbor is indispensable and is stated.

For the numerically admissible pair(A,B)=(4,-1), the equations force
q=1,delta=2,h-s=2. The transitive deletion then makes s either0 or1, and
oddness forces h=3,s=1. Positive defect supplies two distinct in-neighbors.
Together with the first deletion-chain vertex and y, they give a directed
four-cycle and hence an induced strong four-tournament with five paths.
Insertion monotonicity contradicts h=3. This closes the exceptional M8 case
in C4; the proof does not mistake the numerical cone for sufficient data.

The moving-mark construction keeps every old vertex and arc. Its new apex
pair(c_m,x) is an actual directed edge, and its remaining vertices form the
claimed transitive chain. The appended32 letters give exactly the stated
new arcs. The first birth is a sink and the completed graph is strong.
Removing the other clones gives the actual five-vertex exterior232, so the
inherited clone identity proves C5 for all m>=2. The m1 exception is real:
232 has count9;220 has count11 and supplies the separately required base.

The C6 polynomial is uniform in the clone size. A two-bin cover has at most
two starts inside the clone block and one later entry through y, hence at
most three clone runs. For a fixed exterior skeleton, assigning the k
ordered clones onto r nonempty runs has exactly r!S(k,r) possibilities;
the internal order of each run is forced. This gives coefficients independent
of k. The three exact values21,45,93 therefore determine21,12,0, and the
formula holds for every k>=1. No extrapolation from three samples is used
without this bounded-run theorem.

The inherited two-birth reset at the old physical x yields the claimed
balanced profile. The arithmetic is exact:

    5t-b=26,  t odd,
    13|t iff 2^m=7 mod13 iff m=11 mod12.

The order of2 modulo13 is12, checked directly. For m=3j the primitive guard
is automatic;7 divides b and the quotient obeys q_(j+1)=8q_j+1. Modulo49,
15*8^j-1=(1+14)(1+7)^j-1=7(j+2), proving the stated j=5 mod7 guard.
Neither this congruence nor the three initial prime quotients licenses an
infinitude-of-primes assertion or arbitrary target selection.

The finite all-suffix plus-four claim is correctly bounded source by source.
Trailing3 letters are sinks; after removing them, a non3 final letter reaches
an old apex while the original strong source reaches every added chain
vertex. Thus the completion is strong, and the cited strong-order floor
bounds its entire suffix. The bound is not silently extended to all m.

## Independent finite controls

[independent_audit.py](independent_audit.py) uses literal vertex permutations
for Hamilton counts, root starts, deletion counts and two-bin covers cut
immediately after the named root. It verifies:

* Exterior232 at y has profile(9,3,5,11).
* B_m for m2,3,4 is strong and has the stated H and E_x(2) values.
* The first two reset children are strong and have actual profiles
  (59,17,17,59) and(119,29,29,119), on8 and9 vertices respectively.
* The separate m1 construction has H11, whereas generic232 has H9.
* The gcd, divisibility-by7 and divisibility-by49 guards hold for every
  integer m2..97, supplementing the all-order arithmetic proof.

The normal and optimized receipts agree byte for byte. The producer's
complete bounded suffix census is not independently repeated here; its
logical strongness and floor reduction are reviewed above. The additional
permutation controls target precisely the finite exterior and first
reusable balanced children required by the all-order formulas.
