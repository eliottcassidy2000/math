# Independent audit of neighboring, disjoint two-one swaps

20 October 2026. Status: **PROVED for the stated content theorem; FINITE-EXACT for the computational census.** No all-length integer-cycle exclusion or general disjoint-edit theorem is claimed. This audit did not modify the math repository or either inherited frozen packet.

## 1. Verdict, inheritance, and scope

The proof in `../NEIGHBOR_DISJOINT_CONTENT.md` is sound. For every primitive mechanical source

    p>=2, 1<=r<p, gcd(p,r)=1, q=p+r, n=p-r,

moving two consecutive source ones through disjoint adjacent edges in one common direction is parametrized exactly by

    plus:  {u,u+n}, 0<=u<r-n;
    minus: {u,u+n}, n<=u<r.

Every legal parameter has two commuting ordered paths. With the ordinary carry B and D=2^(p+r)-3^p, the actual content satisfies

    gcd(|D|,B) <= 121^(4 p^(2/3)) in the plus family,
    gcd(|D|,B) <=  65^(4 p^(2/3)) in the minus family.

This is an upper bound on the entire integer gcd, not only its radical. Both clock signs and original-cut crossings are included. The family is empty when r<=p/2. There is no p=2 instance; the first is p=3,r=2,q=5.

Closest inherited mechanism: the moving-anchor identity in `../../collatz_two_swap_norm_barrier_20261020/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md` and the prime-5 lattice proof in `../../collatz_two_swap_overlap_lattice_content_20261020/OVERLAPPING_TWO_ONE_CONTENT.md`. These inherited relative paths intentionally use the eventual publication layout; in the staging workspace, the corresponding source folders are `two-edit` and `two-edit-overlap-lattice`. The new valuation prime is 7. The canonical hostile is a sparse resultant that contains clock factors absent from actual content. The repaired near miss is treating resultant divisibility as an iff gate. The least-used sidecar is the whole-modulus common root, retaining prime powers and cut transport.

The concept board was physical disjointness, source-neighbor residues, the whole clock root, three-dimensional exponent relations, valuation noncancellation, and exact content. The bridge is

    physical carry -> common modular root -> two integral relation polynomials
                   -> a nonzero integer resultant divisible by actual content.

It preserves divisibility at every prime-power depth. It loses the original clock-root identity after taking the final resultant, so the resultant cannot replace the iff root gate.

## 2. Exact legality and boundary audit

The source labels have lifted positions and residues

    s_j=j+ceil(jr/p)-1, a_j=-jr mod p.

A gap-two consecutive pair has residue u on the earlier source one and u+n on the next source one, where 0<=u<r. This is equivalent to the standard recurrence a_(j+1)=a_j-r mod p: in the gap-two case the residue wraps by adding p-r=n.

A right move is initially legal precisely at a residue below r. Therefore both are legal iff u+n<r, which is 0<=u<r-n. A left move is initially legal precisely at a residue at least n. Therefore both are legal iff n<=u<r.

For physical positions x and x+2, the two plus edges are {x,x+1} and {x+2,x+3}; the two minus edges are {x-1,x} and {x+1,x+2}. The nonempty legal parameter intervals force p>=3 and q>=5. Hence those four physical bit positions are distinct modulo q. A swap on one edge does not alter either bit on the other edge, proving commutation and legality in both orders. This argument works across the original cut. The boundary q=3 has p=2,r=1 and both intervals are empty. There is no primitive q=4 source in this universe.

A source gap of one cannot contribute to this disjoint family. One member of a same-direction adjacent pair is blocked initially; after the leading member moves, the newly unlocked move uses an overlapping edge. That is the separate prime-5 family, not a missing ordering here.

## 3. Whole-gcd equivalence and excluded primes

For every h>1 dividing |D|, 2 and 3 are units. Choose v,t with pv+nt=1 and set

    z=2^v(4/3)^t mod h.

The clock gives (4/3)^p=2^n. Consequently z^p=2 and z^n=4/3. Any other simultaneous root equals this one by the same Bezout relation. Also

    (z-1)(1+z+...+z^(p-1))=1,

so z-1 is a unit over the entire possibly nonsquarefree modulus.

The moving-anchor formula says h divides the actual carry iff the following sparse quantity vanishes:

    plus:  1+(z-1)(z^u+z^(u+n));
    minus: 1-(z-1)(z^u+z^(u+n))/2.

Substituting z^n=4/3 gives the exact cleared gates

    plus:  3+7z^u(z-1)=0;
    minus: 6-7z^u(z-1)=0.

No nonunit was canceled. In particular, computing z modulo the full |D| gives the exact identities

    g_plus =gcd(|D|,3+7z^u(z-1)),
    g_minus=gcd(|D|,6-7z^u(z-1)),

using any integer representatives. Reduction modulo 7 excludes 7 from g. Since D is coprime to 6, g is coprime to 42. The equations therefore become

    z^u(z-1)=epsilon*2^delta*3/7,

with (epsilon,delta)=(-1,0) for plus and (+1,1) for minus. The prime-7 valuation used below does not assert that 7 divides content; it proves characteristic-zero noncancellation precisely after 7 has been excluded from content.

## 4. Lattice and polynomial audit

For L={(b,c,m):m-bn-cu is divisible by p}, the coefficient 1 on m makes the congruence map onto Z/pZ. Thus det L=p. Define

    a=(m-bn-cu)/p, H=max(|b|,|c|,|m|).

Every legal u satisfies n+u<=p-1. Hence |a|<=(1+n+u)H/p<=H. Negating a relation is harmless and permits c>=0. Put

    e=a+2b+delta*c, f=c-b,
    U=2^max(-e,0)3^max(-f,0)7^c,
    V=epsilon^c 2^max(e,0)3^max(f,0),
    F=U X^max(m,0)(X-1)^c - V X^max(-m,0).

All denominator clearing uses units modulo g, so F vanishes at the actual z modulo the full g.

The following degeneracies were checked separately.

- F is never identically zero for a nonzero relation.
- Neither 0 nor 1 is a root. With c>0, a cancellation requiring equality of the two rational sides contradicts their different 7-adic valuations. With c=0, equality would imply 2^(a+2b)3^(-b)=1 and hence a=b=m=0.
- F is constant exactly when m=c=0. Two independent relation vectors cannot both lie on this one-dimensional axis.
- If m<0 and -m=c, cancellation of the leading terms would require U=V; the 7-adic valuation again excludes this.

For two independent relations, a hypothetical common algebraic root x is neither 0 nor 1. A valuation nu above 7 gives

    m_i nu(x)+c_i(nu(x-1)+1)=0.

Independent projected pairs (m_i,c_i) force nu(x)=0 and nu(x-1)=-1, impossible by the ultrametric inequality. Dependent projected pairs admit a nonzero integer combination eliminating m and c. Independence of the full relations ensures that the resulting (a,b) is not zero. Since ap+bn=0 and gcd(p,n)=1, it is (nt,-pt) with t!=0. The corresponding multiplicative relation would be

    1=2^(-qt)3^(pt),

contradicting unique factorization. This establishes coprimality even for nonprimitive relation vectors and dependent projections.

## 5. Sharp box-based coefficient majorants

The degree is at most |m|+c<=2H. The sum of coefficient magnitudes is at most

    T1+T2,
    T1=2^(max(-e,0)+c)3^max(-f,0)7^c,
    T2=2^max(e,0)3^max(f,0).

Each exponent logarithm is a convex piecewise-linear function of (a,b,c), and exponentiation preserves convexity here. Thus T1+T2 is convex on the box |a|,|b|<=H, 0<=c<=H, and it is bounded by its largest vertex value.

In lexicographic vertex order, the plus values are

    8^H+3^H, 112^H+9^H, 3^H+2^H, 14^H+2^H,
    2^H+3^H,  28^H+9^H, 3^H+8^H, 14^H+8^H.

The minus values are

    8^H+3^H, 56^H+9^H, 3^H+2^H, 14^H+4^H,
    2^H+3^H, 14^H+9^H, 3^H+8^H, 14^H+16^H.

For every H>=1 the maxima are 112^H+9^H and 56^H+9^H, respectively. Therefore

    ||F||1<=112^H+9^H<=121^H,
    ||F||1<= 56^H+9^H<= 65^H.

The second constant improves the separate-maxima bound 56+16=72. The constants 121 and 65 are sharp for this unrestricted-box majorant at H=1, not a claim that some actual shortest lattice vector attains them. The relation a=(m-bn-cu)/p can yield further instance-specific improvements.

The inherited Minkowski second-theorem bound for the cube is lambda1*lambda2*lambda3<=p. Since lambda3>=sqrt(lambda1*lambda2), the first two independent vectors have H1*H2<=p^(2/3). This is a product estimate, not two separate O(p^(1/3)) estimates.

The polynomials are coprime, so their resultant is a nonzero integer. The integral Bezout identity, or direct constant-polynomial case, gives g dividing that resultant at full multiplicity. Sylvester-Hadamard then gives

    g<=C^(4H1H2)<=C^(4 p^(2/3)),

which is the audited content theorem.

## 6. Independent exact computation

The independent producer `audit_neighbor_disjoint.py` uses only the Python standard library and imports no repository mathematical implementation. It constructs the bit words directly, discovers actual source gap-two pairs, tests both edge orders physically, compares all resulting ordered parameters against the claimed intervals, and checks the ordinary carry against an independent position sum. It evaluates the common root modulo the full |D| and compares the two integer gcds exactly.

For all coprime shapes 2<=p<=120 and 1<=r<p it verifies:

- 4,385 source shapes and legality comparisons;
- 351,012 ordered actual paths, 175,506 in each direction;
- 340,824 positive-clock and 10,188 negative-clock paths;
- 4,384 ordered original-cut-crossing paths;
- 351,012 full-gcd iff comparisons.

No integer output occurs in this finite universe. This observation is finite evidence only. The largest actual content is 104233 at (p,r,u)=(119,72,20), plus.

For p<=36, exact sup-norm shell enumeration verifies 3,176 distinct lattices and 4,756 resultants. The selected vectors attain the first two successive minima; the product bound is checked without floating point as (H1*H2)^3<=p^2. Polynomials use exact binomial coefficients; resultants use a fraction-free Bareiss determinant. These checks verify nonzero resultants and divisibility by full g.

Additional depth controls retain minus (19,14,6) with g=25, minus (37,29,25) with g=121, minus (39,34,26) with g=125, and the largest-content plus control. Full g divides each corresponding exact resultant.

### A minimal partial-factor and dependent-projection hostile

At p=3,r=2,n=1,u=0, plus, the legal source 01011 becomes 10110 after both swaps, in either order. Its B=29,D=5,g=1. The independent lattice vectors

    (b,c,m)=(0,1,0), (3,1,0)

have equal projected pairs (m,c). Their polynomials are 7X-4 and 63X-31, whose resultant is 35. The constant vector (3,0,0) gives -5, with resultant of magnitude 5 against 7X-4. Thus gcd(D,R)=5>g in both controls. Dropping the dependent-projection case or interpreting the resultant as an iff content formula is invalid.

Reproduce from this audit directory:

    python3 audit_neighbor_disjoint.py > receipt.replay.json
    python3 -O audit_neighbor_disjoint.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

The frozen normal and optimized receipts are byte-identical. `SOURCE_MANIFEST.json` records the owned audit artifact hashes.

## 7. Stopping boundary

The theorem concerns only primitive k=1 sources, neighboring source ones, one common direction, and disjoint physical edges. It makes no assertion for nonneighbors, mixed signs, repeated sources, general mechanical clocks, or arbitrary two edits. The known negative p=2 overlap cycles are absent here because this family's legal intervals are empty, not because negative cycles were discarded by a sign filter.

The parent proof was audited for all stated quantifiers, iff directions, integer prime-power depth, resultant degeneracies, the exact Minkowski constant, and coefficient maxima. No defect was found. A later signed-descent or Diophantine large-weight closure would be a separate theorem requiring its own hypotheses and audit.
