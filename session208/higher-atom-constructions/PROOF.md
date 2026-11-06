# Beyond the order-eight atom dictionary

Status: PROVED constructions and elementary global lower bounds; VERIFIED explicit witnesses. Searches are not exhaustive. This is a continuation of the audited rational-atom-mixtures packet, commit 58264985c7f8a1590d1b5a934473068d4869888d at `session208/rational-atom-mixtures/`.

## 1. The old dual cannot extend

The preceding packet proved optimality of rate40/11 only within the complete dictionary of strong tournaments through order8. It explicitly did not claim global optimality. That distinction is decisive, not merely caution: an order9 strong tournament has H=1875=3*5^4, with adjacency code14359238086. The proposed prices (16/11,24/11) value this atom at112/11>9. Thus they fail the global atom inequality. This is not a retraction of the finite-dictionary theorem or of its global upper bound.

Bits of all codes below correspond to lexicographically ordered pairs i<j, beginning with bit0. A one means i->j. Codes wider than64 bits are stored as decimal strings.

## 2. Explicit improved blocks

Each listed atom is strongly connected and its Hamiltonian-path count is checked by an independent reverse-endpoint subset recurrence. Exact counts multiply under dominance joins, while orders add. The resulting constructions are:

- Order9 atoms1875=3*5^4 and2025=3^4*5^2 give1875^2*2025^3=15^14 on45 vertices, rate45/14.
- Order10 atoms10125=3^4*5^3 and9375=3*5^5 give10125^4*9375=15^17 on50 vertices, rate50/17.
- Order11 atoms54675=3^7*5^2 and46875=3*5^6 give54675*46875=15^8 on22 vertices, rate11/4.
- Order12 atoms151875=3^5*5^4 and253125=3^4*5^5 give151875*253125=15^9 on24 vertices, rate8/3.

The final atoms have respective codes21193152190910216198 and151165614966068813. `WITNESSES.json` contains every count, exponent pair, code, search iteration and strong-connectivity check. All code conventions are fully specified above.

For every k>=1, writing k=9q+s with0<=s<9 and appending s copies of an order5 H=15 atom proves

 m(15^k) <= 24q+5s.

Hence the existing subadditivity limit satisfies

 lim_k m(15^k)/k <= 8/3.

The displayed rates strictly improve40/11. They are upper bounds only. No completeness of the order9–12 dictionaries, exact minimum for any large count, or global optimality of8/3 is asserted. Finding still larger efficient atoms could improve the bound again. The relevant paper mechanism remains rational mixtures and exact exponent balance, independently proved in the preceding packet; its logarithmic tree-embedding theorem is not used.

## 3. A universal lower bound and its limitation

For every n-vertex tournament,

 H(T) <= n! / 2^floor(n/2).

Proof. Partition all vertex words into orbits under independently swapping the entries in the disjoint position pairs (1,2),(3,4),... . Every orbit has2^floor(n/2) words, since all labels are distinct. In a directed Hamiltonian path, the order inside each adjacent pair is forced by its unique tournament arc. Thus at most one word per orbit can be a path. The quotient count proves the bound. This is an elementary general bound; no novelty claim is made.

Consequently an exact computable global lower bound is

 m(15^k) >= min{n>=1: n!/2^floor(n/2) >= 15^k}.

The simpler asymptotic consequence is

 m(15^k) >= (log15-o(1)) k/log k.

Indeed, the order5 H15 atom gives n=m(15^k)<=5k; the factorial bound, or just H<=n!, gives k log15<=n log n<=n log(5k). Thus the explicit lower bound n>=ceil(k log15/log(5k)) holds. Its asymptotic leading constant is log15. The exact factorial bound can be computed using integers and is stronger than this simplified expression.

This is sublinear in k. It does not prove that the linear efficiency limit is positive, much less that it equals8/3. The interval between a sublinear lower bound and a linear construction remains substantial. The nine explicit atom witnesses prove the old global price extension false, not an impossibility theorem for all future price certificates.

## 4. Search versus proof

Discovery used deterministic pseudorandom tournaments with mt19937_64 seed208, at most100000 candidates at orders9,11,12 and200000 at order10. Counts were computed by forward subset DP. Discovery completeness is neither claimed nor needed: the proof uses only the exact listed witnesses. Verification reconstructs each tournament from its stored code, independently computes reverse-endpoint counts and reachability, and checks every exponent/product identity. The additional order9 H2187 witness is a retained control, not needed for the final8/3 construction.

This packet includes standalone search sources and raw outputs so the discovery budget and provenance are explicit. Search sources guard their supported orders and use integer arithmetic throughout. The64-bit code sources support their stated orders only; the order12 source uses unsigned128-bit codes. Verification uses Python arbitrary-size integers. No paper theorem or external numeric optimizer is a computational dependency.
