# Actual-content bounds for same-direction disjoint edits of neighboring ones

20 October 2026. Status: **PROVED**, with cited Minkowski input and an independent mathematical and exact-computational audit. This is a specified primitive disjoint-edit subclass, not a claim about arbitrary two edits.

## 1. Exact family and theorem

Let p>=2, 1<=r<p, gcd(p,r)=1, q=p+r, and let c be the zero-intercept lower mechanical word of slope p/q. Label its ones by

    s_j=j+ceil(jr/p)-1 modulo q,   a_j=-jr modulo p.

Put n=p-r. Consider two consecutive source ones, separated by a source gap of two, and move both one step in the same direction using disjoint physical adjacent edges. Both swaps must be actual exchanges of unequal bits. They are already legal in the source, commute, and either order is allowed.

The complete parameter ranges are

    plus:  residue pair {u,u+n},  0<=u<r-n;
    minus: residue pair {u,u+n},  n<=u<r.             (1.1)

An empty interval contributes no path. Thus r>p/2 is necessary, p=2 contributes none, and the first nonempty shape is p=3,r=2,q=5. Every path in (1.1) has two orders; original-cut crossings are included.

For the final word w define

    D=2^q-3^p,
    B(w)=sum_{i:w_i=1}2^i 3^(number of ones after i),
    g=gcd(|D|,B(w)).

**Content theorem.** Every path in this family satisfies

    g_plus <=121^(4 p^(2/3)),
    g_minus<= 65^(4 p^(2/3)).                        (1.2)

The bounds retain every prime-power depth and apply to both clock signs. Their mechanism is the actual common clock root, not a raw sparse norm. They do not cover mixed directions, nonneighboring source ones, repeated sources, or mechanical clocks outside q=p+r with 0<r<p.

## 2. Legality and the whole-modulus root gate

Consecutive source ones have gap two exactly when their residues are u followed by u+n with 0<=u<r. A right swap is source-legal exactly when its residue is below r; a left swap is source-legal exactly when its residue is at least n. Applying these criteria to both residues gives (1.1).

The gap-two source positions are x and x+2. The right-moving affected edges are {x,x+1} and {x+2,x+3}; the left-moving edges are {x-1,x} and {x+1,x+2}. All nonempty legal shapes have q>=5, so these bit sets are disjoint even cyclically. Both orders therefore remain actual and yield the same output. The gap-one same-direction case has overlapping edges and belongs instead to the previously classified overlap family.

For any h>1 dividing |D|, 2 and 3 are units. Since gcd(p,n)=1, choose integers v,t with pv+nt=1 and define

    z=2^v(4/3)^t modulo h.

The clock relation gives the unique simultaneous root

    z^p=2,  z^n=4/3.                                (2.1)

For example, raising the definition to p and using (4/3)^p=2^n proves its first equation; the other follows similarly. Also z-1 is a unit because (z-1)(1+...+z^(p-1))=1.

The [moving-anchor identity](../collatz_two_swap_norm_barrier_20261020/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md) gives carry vanishing iff

    plus:  1+(z-1)(z^u+z^(u+n))=0;
    minus: 1-(z-1)(z^u+z^(u+n))/2=0.

Using (2.1), these are respectively

    3+7 z^u(z-1)=0,  6-7 z^u(z-1)=0.                (2.2)

All transport factors are units over the entire modulus h, not just its prime divisors. Equivalently the full actual g equals the gcd of |D| with the appropriate integer representative in (2.2), evaluated at the canonical root modulo |D|. When |D|=1 this equality is interpreted directly.

Equation (2.2) excludes 7 from g. Thus gcd(g,42)=1 and, modulo g, the gate can be written

    z^u(z-1)=K,
    K=epsilon*2^delta*3/7,
    (epsilon,delta)=(-1,0) plus, (+1,1) minus.        (2.3)

The prime 7 does not divide content; its characteristic-zero valuation will certify that two relation polynomials cannot cancel identically.

## 3. Relation lattice and noncancellation

Use the determinant-p lattice

    L={(b,c,m) in Z^3 : m-bn-cu is divisible by p}.

For a nonzero vector set a=(m-bn-cu)/p and H=max(|b|,|c|,|m|). Both legal ranges imply n+u<=p-1, so |a|<=H. Negate the vector if necessary to arrange c>=0. The common-root equations imply

    z^m(z-1)^c = 2^(a+2b+delta*c) 3^(c-b) 7^(-c) epsilon^c.

Put e=a+2b+delta*c and f=c-b. Clearing only units modulo g gives

    F_v(X)=U X^max(m,0)(X-1)^c - V X^max(-m,0),
    U=2^max(-e,0) 3^max(-f,0) 7^c,
    V=epsilon^c 2^max(e,0) 3^max(f,0).              (3.1)

These polynomials vanish at the actual z modulo g, including all prime powers.

No nonzero vector gives the zero polynomial. Neither 0 nor 1 is a root: if c>0, any apparent cancellation at 0 or between leading terms would force the rational right side to have 7-adic valuation zero, contradicting its valuation -c; evaluation at 1 leaves the nonzero second term. For c=0, cancellation would require 2^(a+2b)3^(-b)=1, forcing a=b=m=0. Exactly m=c=0 gives a nonzero constant. Two independent vectors cannot both be of that form.

**Coprimality lemma.** Any two independent lattice vectors give coprime polynomials over Q.

If one is constant this is immediate. Otherwise a common algebraic root x is neither 0 nor 1. At a valuation nu above 7, normalized by nu(7)=1, their rational equations give

    m_i nu(x)+c_i(nu(x-1)+1)=0,  i=1,2.

If the pairs (m_i,c_i) are independent, then nu(x)=0 and nu(x-1)=-1, contradicting nu(x-1)>=min(nu(x),0)=0. If the pairs are dependent, a nonzero integer combination kills m and c. Independence of the full vectors leaves (a,b)!=(0,0), with ap+bn=0. Write a=nt,b=-pt for nonzero integer t. The same combination of rational equations gives

    1=2^((n-2p)t)3^(pt)=2^(-qt)3^(pt),

impossible by unique factorization. This handles the projected-degeneracy case rather than assuming it away.

## 4. Explicit height and resultant bounds

For either family, deg F_v<=|m|+c<=2H. The coefficient norm satisfies

    ||F_v||_1 <= T_1+T_2,
    T_1=2^(max(-e,0)+c)3^max(-f,0)7^c,
    T_2=2^max(e,0)3^max(f,0).

On the normalized box -1<=a/H,b/H<=1 and 0<=c/H<=1, each T_i is log-convex and their sum is convex. Its maximum is therefore bounded by a maximum over the eight vertices. At H=1 the vertex sums, in lexicographic order of (a,b,c), are

    plus:  11,121,5,16,5,37,11,22;
    minus: 11, 65,5,18,5,23,11,30.

At general H>=1 each vertex has the form A^H+B^H<=(A+B)^H. Therefore

    ||F_v||_1<=121^H plus,  <=65^H minus.             (4.1)

The lattice determinant is p. Minkowski's second theorem for the sup-norm cube gives successive minima lambda_1 lambda_2 lambda_3<=p. Choose two independent vectors of heights H_1=lambda_1,H_2=lambda_2. Since lambda_3>=sqrt(lambda_1 lambda_2),

    H_1 H_2<=p^(2/3).

The exact theorem and constant are retained from the [audited overlap lattice proof](../collatz_two_swap_overlap_lattice_content_20261020/OVERLAPPING_TWO_ONE_CONTENT.md), Section 5, citing Barvinok's lattice notes, Theorem 22.5. This is a product bound; it does not assert two individually bounded O(p^(1/3)) vectors.

Their resultant R is a nonzero integer by the coprimality lemma. The integral resultant Bezout identity evaluated at the common modular root gives g|R. If one polynomial is constant, its resultant with the other is a positive power of that constant, and the same divisibility follows directly. The case of two constants is impossible.

Sylvester-Hadamard now gives

    g<=|R|<=||F_1||_1^deg(F_2) ||F_2||_1^deg(F_1)
           <=C^(4H_1H_2)<=C^(4p^(2/3)),

with C=121 or 65. This proves (1.2). The resultant is only a necessary gate, not an iff formula for content.

## 5. Exact controls and remaining scope

The producer `verify_neighbor_disjoint.py` tests every coprime p<=48, 1<=r<p. It independently enumerates same-direction disjoint physical edges moving neighboring source labels, compares the entire ordered path set with (1.1), checks the two commuting orders, and computes carries in two ways. It reconstructs the unique clock root over |D| and checks the whole-gcd iff in (2.2). Exact shell enumeration finds the first two independent sup-norm vectors, and fraction-free Sylvester determinants check nonzero resultants and full-g divisibility.

The universe contains 711 shapes, 11,440 unordered edit pairs, and 22,880 ordered actual paths. There are 710 ordered original-cut crossings. No integral output occurs in this finite universe; that observation alone is not an all-length theorem.

A nontrivial prime-power control is the minus path p=19,r=14,u=6, where actual g=25. The largest content in this finite universe is 2615 at p=23,r=22,u=17, plus. These illustrate why radical-only testing would not suffice.

At p=5,r=3,u=0, plus, the independent vectors (b,c,m)=(0,1,0) and (5,1,0) have dependent (m,c) projections. Their polynomials are 7X-4 and 567X-311, with resultant 91. The constant vector (5,0,0) instead gives -13. Here D=13 but actual g=1, so gcd(D,R) can exceed g. Both projected-degeneracy and constant-polynomial cases are retained in the test script.

Reproduce from this directory:

    python3 verify_neighbor_disjoint.py > receipt.replay.json
    python3 -O verify_neighbor_disjoint.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

No all-length exclusion, general disjoint classification, or historical novelty is claimed in this content theorem. The separate mixed-direction or nonneighboring patterns are not imported into the prime-7 argument merely because their affected edges are disjoint.
