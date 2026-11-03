# Exact coverage as an explicit axiom: what becomes constructive

**CONDITIONAL lens, with PROVED reductions and inherited decidability.**
Assume only for this section the external standard-finite statement

    E: for every positive odd h outside {7,21},
       a finite tournament T exists with H(T)=h.

E remains OPEN in the repository. It is not made a proved dependency by
calling it foundational. Rédei parity and the two exclusions remain in
force: literal raw H cannot cover 0, any even integer, 7, or 21.
The relevant canon is [THM-4094, two-prime-lane reduction](../../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md),
not the historical reflections that asserted coverage. See the October27
correction “two forbidden values do not prove cofinite coverage” in
`01-canon/MISTAKES.md`.

## 1. The finite decision problem is already unconditional

Delete singleton SCCs without changing H. The inherited exact budget is

    K(h)=max_{b in {0,1,2}, 5^b<=h}
             [4b+3 floor(log_3(h/5^b))].              (2)

For h>1, h is attained iff it is attained by a core of order at most K(h).
This is [Session133, bounded-H kernel, Theorem A and section5](../../session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md),
sharpening [rotation-owners spectrum sections3-4](../../rotation_owners_spectrum_20261020e/spectrum.md).
The classical input is **CITED** Busch's exact strong minimum,
[Theorem1/Corollary1, primary PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf),
reopened in this session. Formula (2) is the repository's elementary SCC
budget consequence, not a new theorem here. Compute its floors by integer
powers; its h=1 value is zero, with a singleton used for a nonempty witness.

Thus enumeration of the finitely many labelled tournaments up to K(h),
with exact Hamilton DP, already decides membership and returns a witness
on acceptance. The work is 2^{O((log h)^2)}, including exact integer
arithmetic. E asserts that this terminating decision procedure accepts
every allowed input; it is not needed for termination. This makes the
remaining question a universal positivity theorem for a finite predicate.
No larger census is required to state it, and no finite prefix proves it.

## 2. E supplies an effective witness section, not local source closure

Fix an ordering of the finite tournament codes. Under E select the first
accepted witness. This is a total computable section h -> T_h of the map
T -> H(T) on the allowed integers. No nonconstructive choice principle is
needed: the inherited decision procedure gives the section explicitly.
It can be evaluated by any faithful general machine implementation.

The section may choose unrelated graphs at successive values. It does not
supply a legal birth word from T_h to T_k, a primitive two-port interface,
an H-neutral native route, or a reusable marked child. The actual hostile
[47 -> 73 -> 117 but 47 not ->117](../../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md)
already proves that even existential induced count arrows do not compose.
This remains true under E: the obstruction's target values already exist.
At H31 the [source-selection branch packet](../../source_selection_closure_20261103/BRANCH.md)
gives a finer obstruction in the retained-port language.

## 3. A conditional strong-lane constructor using native edits

There is nevertheless a precise constructive consequence combining E with
the inherited [stabilized native compiler, METAGRAPH section5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md):

**Conditional theorem.** For every target h=p with odd prime p!=7, or
h=7p with odd prime p!=3, E gives an algorithm from a named C3 to an actual
strong count-h tournament, preserving all three seed names. If the chosen
core has order n, the algorithm uses at most n(n-1) native moves, at most
binomial(n,2) temporary births and deletions, n-3 permanent sink births,
and peak order at most n+1, where n<=K(h).

Proof: select a core witness by section1. It is strong. Indeed its SCC
counts are integers >=3; for h prime there is only one factor, and for
h=7p more than one factor would require a factor 7, which is impossible.
This also covers p=7, where h=49. Hence n>=3. Append n-3 sinks to the
named seed and assign the other target names explicitly. For each pair
whose arc differs from the target, apply the inherited single-arc compiler.
Its two legal native moves reverse that pair and preserve every other old
arc and name; deleting its temporary vertex restores order n. At most
binomial(n,2) pairs need correction, proving the stated clock bounds.
The final graph is exactly the chosen strong witness. QED under E.

Every temporary cut and path can be found by finite exact search and
recorded in the program. Native moves are the actual path-complement
operations, not arbitrary relabellings. The permanent sink births and
temporary at-most-two-port exception cuts enlarge the original single-root
adaptive language. Intermediate sources need not be strong, H-neutral,
primitive, or represented by the old seventeen coordinates. The two chosen
seed ports persist as physical names, but the full labelled adjacency and
the transient exceptional names remain necessary sidecars.

This is a conditional fallback compiler with an explicit native clock;
it does not discharge E. It separates two tasks cleanly: reachability from
an actual seed is available in the enlarged language, while proving that
the finite target search always accepts is the remaining arithmetic issue.
Starting from an arbitrarily larger old source would require target padding
or old-vertex deletion; target padding can lose strongness, so the theorem
is deliberately stated from C3.

## 4. Every integer as an address, with the two holes retained

The increasing bijection from N={0,1,...} onto the allowed odd values is

    sigma(n)=2n+1+2*[n>=3]+2*[n>=9].                   (3)

For allowed h its inverse is

    rho(h)=(h-1)/2-[h>7]-[h>21].

Under E, T_sigma(n) is an exact decoded encoding of every nonnegative n.
Read 2rho(H) or 2rho(H)+1 for every nonnegative even or positive odd
integer. The output H is always sigma(n), not n itself. Equivalently,
(H-1)/2 has exactly the missing addresses 3 and10 under E.

The unconditional [strong counter carrier](COUNTER_CARRIER.md) already
encodes every n through H=3+12*2^n, and executes two-register programs.
Its entire count image misses the prime and allowed 7p lanes. This is a
concrete reason why machine expressiveness, integer coding, and exact H
surjectivity must remain separate predicates.

## 5. What a proof may legitimately import from this lens

Use E to choose target witnesses and inspect their possible inverse
programs, but mark every claim depending on that choice CONDITIONAL.
To remove E, prove acceptance of all required prime and 7p instances, or
replace the global search by a well-founded reduction with actual source
entry and reusable terminal ports. A residue controller, a singular
predecessor solution, or a cyclic cover vector is sufficient only after its
actual-source and next-step guards are supplied. The counter collision
shows that even a complete forward observer can forget a required test.

The distinction is operational: forward evaluation and finite decision
already terminate; target acceptance and a compatible local induction are
the open obligations. No cellular-automaton or Turing-machine simulation
changes the truth of that acceptance predicate by itself.
