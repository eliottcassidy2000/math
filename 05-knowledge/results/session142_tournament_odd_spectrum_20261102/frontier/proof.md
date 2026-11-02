# Session 142: tournament spectrum frontier and a 272-case exclusion core

Status: the all-order omissions 7 and 21 are proved; completeness is OPEN.
This packet gives a smaller independently checked finite core for the known
omissions, explicit carry witnesses, and an arithmetic minimal-obstruction
description. It does not claim a new infinite realization theorem or priority
over the cited results.

## 1. Object and literature status

A tournament is an orientation of every edge of a finite complete graph.
H(T) counts vertex words (v1,...,vn) following directed arcs at each step.
Paths are not quotiented by automorphisms. Relabelling T leaves H unchanged.
Let S be the set of such counts over all finite orders; Redei's theorem
places S in the positive odd integers. The target is

    A = {positive odd integers} minus {7,21}.

[Chan--Pak, Computational complexity of counting coincidences,
Theoretical Computer Science 1015 (2024), 114776](https://doi.org/10.1016/j.tcs.2024.114776),
author [PDF](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
dated February 9, 2024, section 6.5, Conjecture 6.2, states S=A as a conjecture
and attributes it to bof and Gordon Royle in 2016. The original
[Royle computation report](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament)
reports allowed odd values **strictly below 80557**, on at most 11 vertices.
The safe inclusive endpoint is 80555, not 80557. Chan--Pak's informal endpoint
wording does not license the latter. [Grinberg's lecture12](https://darijgrinberg.gitlab.io/t/22s/lec12.pdf),
version April 27, 2026, Remark 1.1.11, explicitly retains the open status and
the 80555 endpoint. No later proof of completeness was located in this
targeted November 2, 2026 search; that search is not an exhaustive literature
certificate.

The repo's selected-order 11 theorem THM-4104 certifies allowed ordinary
counts through 80405 and a strong interval 429:2:80265. It is not the full
order 11 census and is distinct from Royle's stronger reported prefix.
These large prefixes are sourced, not freshly replayed by this packet.

## 2. Two published structural inputs, correctly scoped

[Busch, A Note on the Number of Hamiltonian Paths in Strong Tournaments,
EJC 13(1) (2006), N3](https://doi.org/10.37236/1141),
[PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf):

* Lemma 1, pages 2--3: a strong tournament of order at least 5 is either
  nearly transitive, or its vertex set is A union B, where both induced
  pieces are strong, each has at least 3 vertices, and |A intersection B|=1.
* Corollary 1, page 4: the exact strong minimum at order n>=3 is

      f(n)=3*5^((n-3)/3)       if n=0 mod3,
           5^((n-1)/3)        if n=1 mod3,
           9*5^((n-5)/3)      if n=2 mod3.

The minimum theorem does not imply spectrum completeness. Its first
values are 3,5,9,15,25,45; f(n+3)=5f(n), and f is increasing.

For a nearly transitive tournament of order n>=3, orient i->j for i<j, except n->1.
Exactly one Hamiltonian path avoids n->1. Every path using n->1 is
uniquely determined by assigning each intermediate vertex to its increasing
segment before n or after 1. Thus H=1+2^(n-2).

## 3. An exact structural cover much smaller than a census

For order 4, enumerate its 64 labelled orientations. Exactly 24 are strong;
all 24 are relabellings of one canonical tournament (upper-pair code 4),
and each has H=5. The checker verifies orbit equality, not just the two
sets' cardinalities. Strong order 3 is the directed triangle with H=3.

At order 5, Busch's overlap alternative necessarily has |A|=|B|=3.
Label their common vertex 0 and orient the two triangles

    0->1->2->0,     0->3->4->0.

There are exactly 4 still unspecified cross edges, hence 16 orientations.
Every one is strong because both pieces are strong and share 0. Together
with the nearly transitive tournament, this is an isomorphism-cover of
all strong order 5 tournaments, not an enumeration of all labelled ones.

At order 6, the overlap sizes are3 and4. Use0->1->2->0, choose any of the
four vertices of the canonical strong4 tournament as its common root,
relabel the other three vertices 3,4,5, and orient all 6 cross edges between
{1,2} and{3,4,5}. This gives 4*64=256 cases. Root positions are not silently
identified: the four labelled choices guarantee every rooted strong4
type is represented. Add the nearly transitive tournament.

For every one of these 272 overlap cases, the checker verifies strongness,
and compares Held--Karp subset DP with literal permutation enumeration.
It obtains exactly:

| Order | Strong H values |
|---:|---|
|3|3|
|4|5|
|5|9,11,13,15|
|6|15,17,19,23,25,27,29,31,33,37,41,43,45|

Consequently, no strong tournament has H=7: orders 3,4 miss7, and f(5)=9
excludes every larger order. No strong tournament has H=21: the displayed
small spectra miss21, and f(7)=25 excludes every larger order. Singletons
give1, and no order 2 tournament is strong.

The SCC condensation of any tournament is a total dominance order. A
Hamiltonian path cannot return to an earlier component, so it consists
of one Hamiltonian path in each component concatenated in that order.
Therefore H(T) is the product of the SCC counts. A product of positive
integers equal to 7 needs a factor 7. A product equal to 21 needs a factor 21
or a factor 7. These strong factors were excluded above. Hence 7,21 never
occur, at any order.

## 4. Explicit finite carry certificates

Rows encode outgoing neighbours as bitsets; vertices are 0,...,n-1.
The checker independently verifies validity, strongness, and H by both
DP and literal permutations.

| H | n | Outgoing bitset rows |
|---:|---:|---|
|49|7|(110,44,32,116,7,16,54)|
|63|8|(162,252,249,17,1,152,185,24)|
|343|8|(236,161,74,130,143,220,26,68)|

The 63 witness is inherited from
`04-computation/h63_counterexample_audit_s8.py` at 9e1ae5a.
The 49 and 343 representatives were found with deterministic seed 142
and then frozen; their validity no longer depends on the search.

## 5. Exact minimal obstruction and target-monoid factor length

This section is arithmetic, not an assumption that A is already the
tournament spectrum. A is closed under multiplication. Its multiplicative
atoms, meaning elements that cannot split into two nonunits of A, are

    p              for odd prime p !=7;
    7p             for odd prime p !=3 (including49);
    63 and343.

Proof: the listed elements have no allowed proper factor pair. Conversely,
write h=7^a*m with gcd(m,7)=1. If a=0, factor m into ordinary primes. If
a is positive even, extract 49 factors. If a>=3 is odd, extract 343 and 49
factors. If a=1 and m has a prime factor p!=3, extract 7p. The only remaining
allowed case is h=7*3^b with b>=2, which factors as 63*3^(b-2). This expresses
every allowed h as a product of listed atoms and proves exhaustiveness.

The source's two-prime-lane reduction now follows without unverified carry
hypotheses: S=A iff every p and7p in the two displayed prime lanes is
attainable. Necessity is immediate; sufficiency uses the atom factorization
and the three verified carries. Moreover every witness of an arithmetic
atom has exactly one non-singleton SCC: otherwise SCC multiplication
would split that atom inside A. In particular each attainable value in these two lanes
admits a strong witness after its singleton SCCs are removed.

If any allowed count is missing, its numerically least missing value is
an atom: two smaller allowed factors would both be attained and their
order join would attain the product. Since 49,63,343 are attained, the
least missing value must be p or7p in the two lanes. Using Royle's reported
prefix it exceeds 80555; using only the cited local THM-4104 prefix it
exceeds 80405. These are reductions of the unresolved existence problem,
not existence proofs for its remaining prime targets.

A useful closed arithmetic bound is also available. Let Omega(h) count
prime factors with multiplicity and let a=v7(h). The maximum number L_A(h)
of nonunit factors from A in a factorization of h is

    L_A(h) = Omega(h)-ceil(a/2),

except when h=7*3^b, b>=2, in which case L_A(h)=b-1=Omega(h)-2.
Set L_A(1)=0. Indeed, any allowed factor with k copies of 7 costs at least
ceil(k/2) mergers of prime factors. Summing gives the upper bound. The
factorizations using 49,343, or 7p above attain it. The exceptional h has
only one 7 and needs at least two 3 factors joined to it; 63 attains this
extra cost. The checker compares this formula and the atom classification
with an independent divisor dynamic program on all 4999 allowed values
through 10001 (including 1). For any actual tournament with H=h, the number
of nontrivial SCCs is at most L_A(h), unconditionally. Attaining this bound
for every h would require the missing prime-lane realizations; no such
claim is made here.

## 6. What not to infer from the inherited constructions

* THM-1370's tournament statement proves exclusions and a finite prefix,
  not exact completeness. Always use its full slug: a distinct Jacobian
  theorem uses the same numerical identifier.
* THM-4749's Fermat-family construction is valid independently of
  completeness: H(NT_(m+2))=2^m+1 for m>=1, so products of these values are attained.
  But 2 has multiplicative order 3 modulo 7, so no 2^m+1 is divisible by7.
  Even the infinite monoid generated by this family misses every multiple
  of 7, including the already attained63. It cannot close the 7p lane.
* SCC multiplication and padding never create a new prime outside their
  input bank. A finite verified bank plus those operations cannot settle
  the prime lane.
* THM-4118 correctly distinguishes the additive response lattice from the
  feasible image: gcd 2 is necessary, not sufficient for a solid odd interval.
  A selected insertion matching proves monotonicity only; unmatched/orphan
  paths must remain in any exact increment calculation (THM-4094).
* Ordinary completeness should not be silently replaced by the stronger
  assertion that every allowed composite has a strong witness. The proved
  reduction needs strong witnesses only for the arithmetic atoms. Any
  proposed equivalence with full strong completeness needs its own proof.

## 7. Reproduction and scope

    PACKET_DIR=/absolute/path/to/frontier
    python -B "$PACKET_DIR/verify_frontier.py" > "$PACKET_DIR/verify.normal.json"
    python -B -O "$PACKET_DIR/verify_frontier.py" > "$PACKET_DIR/verify.optimized.json"
    cmp "$PACKET_DIR/verify.normal.json" "$PACKET_DIR/verify.optimized.json"
    python -B "$PACKET_DIR/verify_sources.py"
    python -B -O "$PACKET_DIR/verify_sources.py"
    (cd "$PACKET_DIR" && sha256sum -c SHA256SUMS)

All gates use explicit exceptions. Injecting 21 into the expected order 5
spectrum fails in both modes, recorded in `negative_controls.json`.
`repository_sources.json` pins exact full paths, commits, Git blob IDs,
SHA256 hashes and lengths for the seven inherited theorem/frontier/carry sources
at both 9e1ae5a and 91cc1210. Exact original bytes are preserved in `sources/`;
each manifest row supplies its portable `local_snapshot` selector.
`verify_sources.py` checks all 14 source pins offline, including Git blob
identities, in both normal and optimized modes. No giant mistakes ledger was read, no old cache
was deleted, no tracked repository file was changed, and no publication was attempted
by this subtask. The environment restart cleared the old `/tmp` externally;
the current proof has no dependence on that lost cache.
