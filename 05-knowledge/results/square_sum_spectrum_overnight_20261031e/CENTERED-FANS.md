# Midpoint fans, fixed deficits, and lossless palindromic energy

**Status: PROVED scoped all-size theorems, with exact controls. Universal raw-H completeness OPEN.**

The inherited input is THM-4748 / `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`, especially its known one-apex formula and Mersenne fan. We do not claim the Mersenne family or the general one-apex formula as new.

## 1. Exact translated-cluster formula

Take a finite nonempty set S of integers with minimum0 and maximum d. Define

    A(S)=sum_(s in S) 2^s,
    R(S)=sum_(s in S) 2^(d-s),
    W(g)=2^max(g-1,1),
    D(S)=sum_(s<t in S) W(t-s).

For an integer a>=2 set K=2a+d. On vertices0,1,...,K let1,...,K be a transitive chain, and let the apex0 lose precisely to

    (a+S) union {K}.

The tournament is strong because0,1,...,K,0 is a directed cycle. Its raw Hamiltonian-path count is

    H = 2^(K-1)+1 + 2^(a-1)(A(S)-R(S)) - D(S).       (1)

**Proof.** The one-apex formula gives a singleton term2^(K-1), terms2^(a+s-1) for the interior losing vertices, negative pair terms D(S) within that cluster, and negative terms2^(K-a-s-1) for its pairs with K. Since a>=2, none of the last gaps hits the adjacent-edge correction. Summing gives(1). The convention W(1)=2, not1, remains essential for pairs inside S.

For a fixed S, the difference H-2^(K-1) is constant as a grows **if and only if** A(S)=R(S). Binary expansion is unique, so this is equivalent to

    S = d-S.                                           (2)

Thus midpoint symmetry is exactly the cancellation condition in this translated-cluster grammar. It is not merely a visual resemblance. For a non-symmetric S the missing term grows like2^a. The controls S={0,1,3}, {0,2,3}, {0,1,4} exhibit both signs of this term.

## 2. A uniform family below powers of two

If S=d-S and |S|>=2, put c=D(S)-1, a positive odd integer. For every K>=d+4 with K congruent d modulo2, the above explicit strong tournament has

    H=2^(K-1)-c.                                       (3)

It has exactly |S|+1 back arcs in the displayed order, hence feedback arc number at most that many; no minimality is claimed. The actual source labels and K determine the tournament.

Examples:

| Normalized S | d | c | Exponent parity k=K-1 |
|---|---:|---:|---|
|{0,1}|1|1|even|
|{0,2}|2|1|odd|
|{0,3}|3|3|even|
|{0,1,2}|2|5|odd|
|{0,4}|4|7|odd|
|{0,2,4}|4|11|odd|
|{0,1,2,3}|3|13|even|
|{0,1,3,4}|4|21|odd|

The first two recover the earlier Mersenne fan, now as the smallest symmetric clusters. A deficit equal to7 or21 is not a forbidden raw count: it is the amount subtracted from a power of two.

## 3. An exact outer-pair recurrence

Given symmetric S in[0,d] and h>=1, form

    S'={0,d+2h} union (h+S).

Its complete controller update is

    d'=d+2h,
    A'=1+2^h A+2^(d+2h),
    D'=D+2^h A+2^(d+2h-1)+2*[h=1].                  (4)

The old internal pair energy is D. Symmetry makes the cross terms from the two new endpoints equal. Their sum is2^h A, with an extra2 when the two nearest gaps are1. The new endpoint pair contributes2^(d+2h-1). This proves(4), including the exceptional adjacency correction. The singleton seed S={0} also satisfies the formula.

Equation(4) is an actual family recurrence: it constructs new labelled tournaments through their losing sets. It retains the binary mask A and diameter d; it is not a formal substitution in H alone.

## 4. At fixed diameter, the energy itself encodes the whole symmetric pattern

**Theorem.** For symmetric subsets S of[0,d] containing0,d, D(S) is strictly increasing with their binary masks A(S). In particular the pair(d,D) determines S exactly. It can be decoded by binary search over the independent half of the palindrome.

**Proof by an adjacent binary carry.** Order the independent reflection-orbit bits from the outside toward the center. Between consecutive masks, an absent orbit at positions a,d-a is turned on; all inner orbits, previously on, turn off. The common outside set E remains. Put w=d-2a.

If w=0, the operation simply adds the center; the energy increases by the sum of its positive interactions with E. If w=1, it adds the two adjacent central points; the increase is2 plus their positive interactions with E.

For w>=2, the old inner set is the full interval a+1,...,d-a-1, whereas the new inner set is{a,d-a}. The energy of r consecutive points, r>=1, is2^r-2: sum(r-g)W(g) over gaps g. Therefore the internal energy increases by exactly2. For each outside point x, let b be its distance to the nearer new endpoint. Its cross-energy increase is

    W(b)+W(b+w)-sum_(j=1)^(w-1) W(b+j)
      =3*2^(b-1)+[b=1] >0.                            (5)

The same formula applies on either side. All increments are positive. This proves strict increase for every carry and every diameter.

The theorem gives an unexpectedly strong scoped encoding: this weighted distance statistic, which usually loses an entire configuration, retains the complete palindromic bit pattern when its diameter is marked. Its inverse need not enumerate every pattern: an energy evaluation takes O(d^2) arithmetic operations, and O(d) binary-search comparisons suffice.

**The diameter is essential.** S={0,1} and S={0,2} both have D=2. Adding one new outer pair at gap h=1 gives D'=14 and22 respectively. The erased clock is testable by a later legal operation. Conversely, once d is fixed, there is no same-energy ambiguity: the carry theorem restores the mask.

## 5. Why this still cannot prove the complete raw spectrum

These fixed deficits form a sparse set. If c=D(S)-1<=X and the diameter d>=2, the endpoint pair alone gives

    2^(d-1)<=X+1,

so d<=floor(log_2(X+1))+1. At diameter d>=1 there are exactly2^floor(d/2) symmetric subsets with both endpoints. Summing over all allowed d gives only O(sqrt(X)) possible deficits. More explicitly, up to diameter2r there are3*2^r-3 patterns, and up to2r+1 there are4*2^r-3. Duplicate deficits across diameters only reduce this number.

In particular deficit9 is impossible in this grammar: D would be10, forcing d<=4; the complete energies through d=4 are{2,4,6,8,12,14,22,30}. Thus a fixed translated midpoint cluster cannot generate a family2^k-9, even though particular values on that line can be realized by other tournament words. No global exclusion for2^k-9 is inferred.

`centered_fans.json` completely enumerates this grammar's deficits through1,000,000:2399 distinct deficits, not all positive odds. It independently checks the outer recurrence, every adjacent palindrome carry through diameter20, and actual adjacency Hamilton counts for183 small fan instances, with subset DP where the size is small. The all-size proofs above use the formulas, not this finite count.

The next useful expansion must change more than the translation parameter: allow growing cluster width, more than one moving cluster, a second feedback vertex, or an authenticated re-marking with a larger boundary state. The seven-parent bank is still a live candidate. Its100-million finite coverage is stronger evidence of numerical breadth than this sparse fan lane, while the fan lane supplies an exact explanation of a class of cancellations.

## 6. An exact information bridge, with its limits

The branch code for any one square-sum family in README has a fixed length b. Put d=2b+1, include endpoints0,d, and use its b bits to choose the b interior reflection pairs. This gives a symmetric S. Set K=d+4 and construct its fan. Its raw count is2^(K-1)+1-D(S). The strict carry theorem recovers S, then all b bits, from that integer **together with d**. With the original square-sum size and schedule retained, those bits reconstruct the original cycle.

This is a reversible map of the two authenticated code families, not a graph isomorphism or a preservation of their edge relations. Source: a cycle in the designated square-sum grammar. Target: a strong one-apex tournament and its raw count. Preserved predicate: its exact branch history. Lost information without sidecars: source size/schedule and the fan diameter. `decode_fan_energy.py` tests this using actual adjacency counts and also proves the inverse implementation against every pattern through diameter20,10,010 complete small positive/negative membership queries, and exact256-bit addresses. The inverse itself uses binary search, not enumeration.

The arithmetic image is sparse, as Section5 proves. Efficient lossless encoding and universal numerical coverage are different obligations.
