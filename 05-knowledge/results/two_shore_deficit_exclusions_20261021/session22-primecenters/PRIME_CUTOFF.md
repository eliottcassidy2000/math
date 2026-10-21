# Prime-circulant connectivity and the final degree45 cutoff

Status: PROVED; independent analytic and finite-control review passed. This is a separate final refinement; earlier degree47, q2, and24-coefficient support packets remain unchanged as historical checkpoints.

## 1. Exact connectivity of a prime-order circulant tournament

Let T be the circulant tournament on the additive group Z/pZ, p an odd prime. Write D for its connection set, so x->y iff y-x belongs to D, with D disjoint from -D and D union(-D) the nonzero residues. Set d=|D|=(p-1)/2. Then

    kappa(T)=d.

For the lower bound, suppose deleting a set S leaves a nonstrong tournament. Its surviving SCC order admits a nonempty split X->Y. No arc can go from Y to X. Hence

    Y+(D union{0}) is contained in Y union S.

The **CITED Cauchy–Davenport theorem**, for nonempty subsets of Z/pZ, gives

    |Y+(D union{0})| >= min(p, |Y|+d).

The left side is at most |Y|+|S|=p-|X|<p. Therefore the minimum cannot be p, and |Y|+d<=|Y|+|S|, giving |S|>=d.

For the upper bound, delete all d outneighbors of any vertex v. The remainder has v as a sink and has at least one other vertex, so it is nonstrong. Thus kappa(T)<=d, proving equality.

The prime hypothesis is essential. On Z/9Z take D={1,3,4,7}. Its regular circulant tournament has a separator S={1,4,7} of size3<4=|D|. The remaining shores X={2,5,8},Y={0,3,6} satisfy X->Y. The additive step fails precisely: Y+(D union{0}) has size6, smaller than min(9,3+5-1)=7. This is a hostile against applying the prime theorem at composite center orders such as39 or45.

Cauchy–Davenport source: Terence Tao, *An uncertainty principle for cyclic groups of prime order*, arXiv:math/0308286, printed PDF page5, displayed inequality under “Cauchy-Davenport inequality”: https://arxiv.org/pdf/math/0308286 . The primary source was checked during this session. The displayed theorem uses MIN; a later proof-ending line contains a MAX typo and is not the statement being invoked. No literature priority is claimed for the additive theorem or its elementary separator application.

## 2. Exclusion of the remaining large prime orders

Let T be a strong self-converse deficit-six center of prime order p>8. The inherited cycle-cut theorem gives q<=8 and p=qm with m odd. Primality forces q=1. Therefore T is a circulant tournament on Z/pZ, and section1 gives kappa(T)=(p-1)/2.

If p>=37, then kappa(T)>=18>16. Strong pancyclicity supplies a directed16-cycle on a vertex set C, and deleting C leaves a strong tournament. Also p>=32, so the initial16 positions and their reversals are disjoint. Any external vertex has at least8 of its16 arcs to C in one direction. The already proved half-supported prefix lemma, using the entry dual if needed, gives

    Delta(T)>=2(8-1)=14,

contradicting Delta=6. Thus every prime-order center with p>8 has p<=31.

The only additional prime orders within the previously proved cutoff47 are37,41,43,47, all now excluded. No assertion about composite orders follows from section1.

## 3. Updated polynomial support and recurrence

Combining this lemma with the degree47 result and its support reduction (which already excludes order46) removes the top possible degree47 and leaves

    deg R6 <=45.

The previous24-coefficient support reduction, with its explicitly inherited low-order data through9, becomes

    R6(z)=z^5+z^6+z^7+sum_{n in V} r_n z^n,

where the twenty unknown nonnegative integer coefficients have indices

    V={11,13,15,17,19,21,23,25,27,29,31,33,35,39,45}
      union {12,18,20,24,28}.

This is a necessary support reduction, not existence or determination of any r_n. The small coefficients and order9 zero remain inherited, not freshly censused by this lemma.

For the inherited degree-three generating function, the common numerator degree is at most deg(R6)+12<=57. Therefore the already proved order18 recurrence is homogeneous for every n>=58. Its exact minimal eventual order18 remains inherited algebra; n=58 is only a safe onset, not a claimed earliest onset.

## 4. Finite controls and boundaries

The tiny checker exhausts the Cauchy–Davenport inequality for all nonempty subset pairs at primes3,5,7, checks minimum vertex cuts for every cyclic orientation at those primes, verifies the explicit composite9 hostile, and checks the final support arithmetic. These checks support the proof and guard its prime scope; they do not replace the cited additive theorem or enumerate center tournaments at37..47. Literal connection sets and the composite witness are retained in the receipt, with run timestamps and a source hash.
