# An initial complete symmetric 7p spectrum and a mod18 port rule

**PROVED from inherited all-order reductions, Moon's marked-cycle theorem,
and FINITE-EXACT complete quotient tables.** This classifies a symmetry
branch; it does not extend the inherited unrestricted positive-count prefix
or prove completeness of the unrestricted spectrum.

## 1. The exact conclusion

For an odd prime p<4375, there exists a tournament T with nontrivial
Aut(T) and raw H(T)=7p **if and only if**

    p in {593,827,2113,2297,2531,3433,3517,3533,
          3767,3851,3967,4051}.                         (1)

Every such automorphism group has order7. The following orders are minimum
**among symmetric witnesses**, not asserted minimum among all witnesses:

| p | H=7p | least symmetric order |
|---:|---:|---:|
|593|4151|10|
|827|5789|10|
|2113|14791|12|
|2297|16079|12|
|2531|17717|12|
|3433|24031|12|
|3517|24619|12|
|3533|24731|13|
|3767|26369|12|
|3851|26957|12|
|3967|27769|13|
|4051|28357|12|

The boundary593 and its order10 witness are **inherited** from Session142
`minimal_gap/PROOF.md`, newly received at public commit0d2b365955. The present
extension is the complete initial list and its marked-module mechanism.

## 2. Exact homogeneous-module response

Let C_j(B) count **ordered**, nonempty spanning j-path covers of B. Thus
C_1=H and C_2=E-2H. For a tournament Q with marked vertex v, let b_j(Q,v)
count directed words using each exterior vertex once and j indistinguishable,
nonadjacent copies of v, with adjacency checked in Q. Then

    H(Q[v<-B]) = sum_(j>=1) b_j(Q,v) C_j(B).             (2)

This maximal-run kernel is inherited from THM-3121 (the complete slug is
in SOURCE_PINS.json), Session136 gluing and the incoming Session142
`path_cover_amplification/PROOF.md`. Its use here is the regular-seven
response calculation and the complete symmetric-prime branch.

A full path decomposes uniquely into its maximal visits to B; these ordered
runs give C_j, and contracting them gives the quotient word b_j. Conversely
expand each copy by its prescribed run. Homogeneous adjacency makes every
joining arc valid. This proves a bijection; there is no additional j!.
In particular b_1(Q,v)=H(Q), independently of v.

The source is the actual marked quotient and the inner adjacency. The map
forgets labels within a run only after recording its ordered path. The
formula does **not** apply to arbitrary nonhomogeneous attachments. It is the
same boundary-state discipline as the repaired square-sum cyclic sewing:
multiple runs must not be replaced by one locally Hamiltonian cycle.

## 3. The regular seven-vertex response line

The three inherited regular seven-vertex classes have raw counts171,175,189.
Complete permutation/cut enumeration gives, respectively,

    C(B)=(H,1638,6678,15120,20160,15120,5040).            (3)

Thus all responses of two or more runs agree. Circular R7 with steps{1,2,3}
has H175; Paley P7 with steps{1,2,4} has H189. Equation(2) gives the all-context
identity

    H(Q[v<-P7]) - H(Q[v<-R7]) = 14 H(Q).                 (4)

Both objects still have seven labelled module vertices and identical exterior
ports. The quotient in the right-hand coefficient is a counting sidecar,
not a physical deletion allowed to change the source clock.

There is a more informative reason for the agreement. The three classes have
(C3,C5,M33,C7), counting literal directed cycles up to rotation and unordered
pairs of disjoint triangles,

    (14,28,14,17), (14,36,10,15), (14,42,7,24).

They all satisfy C5+2 M33=56. REGULAR_BALANCE.md proves its all-order
regular-tournament formula and gives the two-coordinate order9 repair.
The inherited odd-cycle-collection formula,
applied separately to each colour class of an r-colouring, gives

    E_r(B)=sum_collections 2^(#cycles) r^(7-sum(|cycle|-1))
          =r^7+28r^5+112r^3+2 C7 r.                    (5)

Here E_r counts ordered r-slot path covers, allowing empty slots; every
selected cycle must be monochromatic, and each unused vertex has a free
colour. Also E_r=sum_j binom(r,j) C_j. Five-cycles and two disjoint triangles
have the same deficit4 and hence the same r^3 exponent. The collapse erases
their separate counts. It is an exact instance of a lost structural coordinate,
not an identification of those two cycle configurations.

For any two of these three regular modules, replacing one by the other
changes the contextual H by (H_new-H_old)H(Q). Several independent such
replacements yield a multilinear expansion: a term for a chosen set of
replacements is the product of their H differences times the count with
those modules contracted to singletons and all remaining modules unchanged.
This follows by iterating(4). The forward graphs keep all original labels;
contractions belong only to the coefficient calculation.

**Failure boundaries.** At order9, cyclic regular tournaments already have
different C_2 values42606,43074,43290. The seven-vertex response line is not
an all-regular-order law. Nor does(4) survive a nonhomogeneous exterior
attachment: if one new mark beats precisely vertices0 and1, direct counts
are561 for R7 and579 for P7, a difference18 rather than14. The cut must remain
part of the interface.

## 4. The mod18 connection is an exact map

From(3), gcd(C_2,...,C_7)=126. For the circular module,

    H(Q[v<-R7])/7 = 25 H(Q) +234 b_2+954 b_3
                     +2160 b_4+2880 b_5+2160 b_6+720 b_7.

Consequently, at every quotient size and every homogeneous marked port,

    H(Q[v<-R7])/7 = 7 H(Q) (mod18).                     (6)

On the units modulo18, multiplication by7 has the two cycles

    1 -> 7 -> 13 -> 1,       5 -> 17 -> 11 -> 5.

This recovers the earlier {1,7,13} motif through a specified source, map and
preserved residue. If the left side is prime greater than3, H(Q) must be
coprime to3 and must lie in the residue13p modulo18. These are necessary
source filters, not a converse: the b_j must come from the **same actual
marked quotient**, not arbitrary nonnegative solutions.

For Paley7, all C_j are divisible by21, so every homogeneous substitution
has H divisible by21. Its module has an explicit21-element affine symmetry
group. It cannot supply H=7p with p>3. The near equality of the two regular
cover profiles therefore hides a decisive difference in the prime branch.

## 5. Why the finite table gives an all-order classification

We inherit from Session142's fully stated proof:

1. For H=7p with prime p>=11, the tournament has one nontrivial strong core,
   and its automorphism group has order1 or7. Singleton SCCs are fixed.
2. A nontrivial orbit of an order7 automorphism is a circulant on7 vertices,
   with H175 or189. Two disjoint orbits force H>=175^2=30625, by the injective
   order-preserving merge of induced Hamiltonian paths.
3. The small prime cases p=3,5,7 give no symmetric witness.

For H<30625 there is therefore exactly one nontrivial orbit. All exterior
vertices are fixed, so it is a homogeneous module. Contract it to the marked
vertex of a strong Q. Paley7 is excluded by its21-divisibility, leaving R7.

Moon's *On Subtournaments of a Tournament* (1966), Theorem1, states that every
vertex of a strong tournament lies on a directed cycle of each length3..n.
Thus if |Q|>=8, it contains a marked strong eight-vertex subquotient. Lift its
marked vertex to the **entire original module**. This is an induced
subtournament of the original source.

The complete exact marked-quotient minima are:

| quotient order | strong classes | minimum with R7 | minimum with P7 |
|---:|---:|---:|---:|
|3|1|2163|2205|
|4|1|4151|4221|
|5|6|6489|6615|
|6|35|10815|11025|
|7|353|19467|19845|
|8|6008|32445|33075|

In particular |Q|>=8 would imply H>=32445>30625. Hence Q has at most7
vertices. The complete marked table through7 yields exactly(1), with actual
adjacency witnesses. Q of size1 gives quotients25 or27, not primes; size2
cannot be strong. This proves both directions. Neutral singleton padding
cannot lower order, and two nontrivial7-orbits already need14 vertices, so
the first occurrences in the table give the stated symmetric minimum orders.

The cutoff is tied to the two-orbit capacity30625. Beyond it, two or more
moving orbits may interact nonhomogeneously, so the single-module argument
is not a complete carrier. It would be incorrect to extrapolate the12-value
list to the entire symmetric7p spectrum.

## 6. Exact universe and audit

`research.py` uses all nauty gentourng classes at orders3..8, all marked
vertices and both seven-cycle profiles. This is50782 marked strong inputs.
Its endpoint-oriented skeleton DP is compared with `module_audit.cpp`'s
remaining-vertex recursion. The latter computes inner C_j independently by
all vertex permutations and all legal break positions. Both produce the same
minima and12 prime targets. Full-adjacency DP separately checks each displayed
witness and small quotient controls. Both routes depend on the same trusted
complete nauty generator; they are not independent class-generation proofs.
The input banks, commands and hashes are retained.

Primary sources reopened for this session:
- J. W. Moon, Theorem1, DOI10.4153/CMB-1966-038-7,
  https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on_subtournaments_of_a_tournament.pdf
- A. H. Busch, Corollary1 (inherited strong floor),
  https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/
- The odd-cycle-collection identity is inherited through
  `01-canon/theorems/THM-002-ocf.md`; its canonical literature attribution is
  Grinberg--Stanley, not a new identity of this packet.
