# Fixed-port arithmetic traps invisible to symmetry, and a root-change sidecar

Status: PROVED elementary consequences of the inherited fixed-port transfer;
explicit finite witnesses checked independently below. This is an operation-level
obstruction and information-loss theorem, not Chan–Pak completeness.

## Inheritance and concept board

Anchor: odd Hamilton-path spectrum, still OPEN beyond its finite certificates.
Niche: rooted stabilizers versus arithmetic content. Wildcard: whether a
sufficient statistic for every fixed-port continuation also predicts the benefit
of changing port. Live board: subgroup orbit gaps; fixed-port O/Z transfer;
root deletion; scalar interface mimics; root change; converse.

Recovered mechanisms, with honest provenance:

* THM-4888, `01-canon/theorems/THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement.md`,
  and THM-5085, `01-canon/theorems/THM-5085-fixed-port-tails-in-every-base-and-the-ternary-bank-to-ten-million.md`:
  the old Sessions 195–206 fixed-port transfer is KNOWN at base 2. It is not new here.
* THM-4987, `01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`:
  every actual automorphism subgroup has Hamilton-path orbit count outside 7,21.
* THM-4990, `01-canon/theorems/THM-4990-stabilizer-corrected-hamilton-prefix-orbit-transfer.md`:
  marked automorphisms divide ordered cover coordinates. Rigidity alone does
  not imply primitive marked content. This warning is inherited, not new.
* THM-4879, `01-canon/theorems/THM-4879-two-port-transitive-seam-gadgets-closed-transfer-matrices-with-a-three-power-entry-and-the-seam-parity-r-equals-s-mod-2.md`:
  rooted seam parity and integer transfer matrices; converse requires the
  relevant endpoint sidecar, rather than an informal symmetry claim.
* Hostile from THM-4988: endpoint orbit counts can equal 7 or 21, so the global
  exclusion must not be applied to every coordinate. We retain deletion counts
  separately.

The new synthesis here is the exact four-test content certificate, explicit
rigid arithmetic trapping, and a pair of physically realizable equal states
with opposite answers to primitive re-rootability. No external novelty claim.

Equal-interface collisions themselves are already known: Session227
[multiroot section5](../../session227_chan_pak_constructions_20261108/multiroot/PROOF.md)
constructs binom(2m,m) nonisomorphic strong tournaments at the same order7m+5
with identical four-coordinate rooted states and all common fixed-root graft
responses. Its convention (D,E,U,d) translates here to (h,R,e,d)=(D,U,E,d).
That theorem expressly excludes arbitrary rerooting from its equivalence.
The advance here is a witnessed difference in primitive reroot escape and
converse states, not discovery that the four coordinates have collisions.
The larger Session227 family is a useful future root-deck stress test; no
reroot difference within that family is established in this session.

## Definitions

For a rooted nonempty tournament (C,u), let h=H(C), d=H(C-u), e the number of
Hamilton paths ending at u, and

    R = sum over nonempty directed paths P starting at u of H(C-V(P)).

Use H(empty)=1. Thus the rooted state is z=(h,R,e,d). A letter S adds a sink.
A letter B adds w beaten by all existing vertices except u, with w→u. Keep u
fixed. In chronological word order, the inherited formulas are

    S(h,R,e,d) = (h,2R,0,d),
    B(h,R,e,d) = (h+R-e,2R-d,d,d).

Let c(z)=gcd(h,R,e,d)>0, called fixed-port content here. It is odd because h
is odd. This is arithmetic content, NOT an automorphism order.

## Theorem 1: an exact finite certificate for an infinite family

For every physical rooted state z,

    gcd{H(C_w): w in {S,B}*}
      = gcd(H(C),H(C_B),H(C_SB),H(C_BB))
      = c(z).

Proof. Every transfer is an integer homogeneous linear map, so c(z) divides
every descendant coordinate. The four displayed Hamilton counts are

    h, h+R-e, h+2R, h+3R-e-2d.

If g divides these four, it is odd since it divides h. Their differences
show g divides 2R, hence R, then e, then2d, hence d. This proves both
inclusions. There are no bounded-word or numerical assumptions. QED.

Consequently, content can never decrease under a fixed-port continuation:
c(z) divides c(wz). More generally this divisibility holds for any integer
homogeneous linear rooted-interface map, including any seam transfer for which
these exact coordinates close. Do not extend this to arbitrary root changes.

The prime-by-prime onset of a new trap is completely explicit. For an odd prime
p not dividing c(z):

* p divides c(Sz) iff p divides h,R,d and does not divide e.
* p divides c(Bz) iff p divides R,d,h-e and does not divide h.

For S this follows immediately from the formula. For B, divisibility of its
last coordinate and second coordinate first forces p|d,R, and then p|(h-e).
The hypothesis p∤c(z) forces p∤h and p∤e. These criteria are only for a newly
acquired prime, not a complete statement about increases in its valuation.

For a finite bank of roots with contents c_i>1, every prime Hamilton value
represented by this bank must itself equal one of the c_i that is prime.
Indeed c_i divides the prime output. Therefore only content-one seeds can
supply any primes above max_i c_i. This removes all other seeds from the
large-prime part of a coverage proof, but does not prove sufficiency for the
remaining seeds, nor a positive-density miss theorem.

## Theorem 2: rooted symmetry is a strictly weaker obstruction

Let A=Aut(C)_u, of order a. Then a divides all four coordinates and hence c(z).
For h and e, ordered spanning paths are free A-sets. For d, an automorphism
fixing u and an ordered path on C-u fixes every vertex. For R, use pairs
(P,Q), where P starts at u, Q is a spanning Hamilton path on the complement,
and Q is allowed to be empty. Fixing an ordered pair fixes every vertex.
This proves divisibility directly, without assuming anything about newly
created automorphisms of the larger tournament.

Moreover h/a and d/a both avoid 7,21. The h assertion is THM-4987. The action
of A restricts faithfully to C-u; if C-u is nonempty, apply THM-4987 there.
If it is empty, then a=d=1. In contrast e/a need not avoid 7,21, by THM-4988.
The distinction between actual subgroup quotients and formal division by
c(z) is essential: this proof gives no gap theorem for h/c(z) or d/c(z).

Strictness holds even for a strongly connected fully asymmetric tournament.
Use vertices 0,…,4 with out-neighbour bit rows

    T = [24,1,19,22,2].

At root 1 its state is (9,9,3,3), with content 3. Its outdegree classes are {0}, {1,4}, {2,3}. Every automorphism preserves
these classes and cannot swap the two endpoints of an arc. Hence its full
automorphism group is trivial, not just its rooted stabilizer. The spanning
directed cycle 0→3→2→4→1→0 proves strong connectivity. All fixed-port tails have counts
divisible 3. Re-root at 0 instead: the state becomes (9,11,1,1), with content 1,
and one B produces H=9+11-1=19. Thus port motion can release a genuine
arithmetic trap even when there is no symmetry left to break.

The converse tournament at the SAME root 1 has state (9,21,1,3), so content 1;
a B there gives 29. Thus even converse can release the arithmetic trap while
preserving h and d. This is not an automorphism-breaking claim.

This trap is created by one B on the four-vertex root

    C=[8,1,3,6], u=1, z=(5,6,2,3), c(z)=1.

It satisfies the prime-onset criterion at p=3. It also lies on the endpoint
equality h=2e+1 (5=2·2+1); the bridge leaves that equality and enters the trap.
This connects the endpoint-minimality lane with the arithmetic-content lane
without claiming the equality forces trapping in general. The child is exactly T.
A sink-triggered counterpart is C=[20,1,2,7,14], u=4,
z=(9,21,1,3), whose S-child has state (9,42,0,3). Both roots are rigid.

Minimality here refers only to the tested labelled universe: no S-triggered
content increase occurs below order 5, no B-triggered increase below order 4,
and no rooted-rigid state with content>1 occurs below order 5. Do not read
this finite result as a classification of all rigid arithmetic traps.

## Theorem 3: the fixed-port interface loses root-change escape information

Let U be the ordinal sum C3⇒C3, and root it at any vertex of the second C3.
Its state is again (9,9,3,3): every rooted path/cover computation factors by
H(C3)=3, and a rooted C3 has state (3,3,1,1). Consequently T and U have
identical Hamilton counts after EVERY common fixed-port S/B word, by induction
on the transfer formulas.

Nevertheless T has a content-one re-root, while U has none. At every vertex
of U, the rotation of the other C3 is a rooted automorphism of order 3, so
Theorem 2 gives content divisible 3 for every re-root. This proves that no
function of the four fixed-port coordinates can determine whether a
primitive re-root exists. This is an actual pair of physical tournaments,
not a non-realizable integer-state counterexample.

There is also NO universal converse-state map on these four coordinates.
The same equal-state pair has converse-root states (9,21,1,3) for T and
(9,36,0,3) for U. In the latter, after reversal the root is in the first
triangle: e=0 and R=36 (checked by literal path enumeration). Therefore
converse releases T's trap but not U's. This does not contradict known seam
adjoint formulas, which retain their own two-endpoint and transfer-matrix data.

The information loss persists even if the order is supplied. Prepend a single
source to T, obtaining X with rows[62,48,2,38,44,4], and use root 2 (the old
root 1). Prepending a source leaves every old-root coordinate h,R,e,d unchanged:
root-starting paths cannot visit the source, and every complement count retains
its ordinal source factor1. Thus X and U both have order6 and state(9,9,3,3),
but only X admits a content-one re-root (vertex1). X is asymmetric because its
unique source is fixed and T is asymmetric. Their converse states are
(9,42,0,3) and(9,36,0,3), respectively. The former follows by appending a sink
to T's converse. So even (order,h,R,e,d) cannot decide primitive re-rootability
or determine the converse state. The equal-order strengthening is a padding
argument, not a separate search claim.

Connection ledger: the map is (C,u)→(h,R,e,d); its preserved predicate is the
entire fixed-port tail count function; it destroys the root-resolved deletion,
endpoint and response data. An adequate sidecar is the set of rooted states
{z(C,v):v∈V(C)}, retaining vertex incidence if subsequent operations need it.
The cheapest decisive test is the pair T and U above. A complete bank proof
that permits changing ports cannot quotient only by the four coordinates
without retaining at least information about this sidecar.

## Scope and verification

Theorem 1 is an elementary consequence of known transfer matrices, not a new
transfer law. Theorems 2–3 connect this arithmetic consequence with inherited
symmetry restrictions and exhibit the precise loss under state compression.
Nothing here resolves baseline placement, positive density of finite-bank
misses, or the all-odd-except7,21 conjecture.

Reproduce the finite checks from repository root:

    python3 04-computation/experiments/session237_symmetry_grafts/probe.py
    python3 04-computation/experiments/session237_symmetry_grafts/witness_audit.py

The first uses independent subset path DP to compute physical states and
checks short literal tail constructions. The second uses literal permutation
and path enumeration for the witnesses, strong connectivity, full groups,
root spectra, equal-state mimic, and the 19 escape. The complete state universe is all 202,013 labelled rooted tournaments
of orders 1 through 6 (no inherited filters). Physical short-tail replay covers
all roots of all labelled tournaments through order 4 and all words through
length 3, totaling 4,275 comparisons. The separate literal audit checks seven
word pairs through length 2 on T and U, plus their full automorphism groups
and all re-roots. The first appearances claimed above therefore have exhaustive
lower-order exclusions. Exact counts and witnesses are in the output files; no sampled computation is called exhaustive.
