# Independent audit of the exact O/Z integer index

**PASS, in the stated four-carrier scope.** This is an all-order parameterization and a sharper upper-density bound for four already physical O/Z orbits. It does not prove arbitrary-base occupancy, F(8,137), prime-track coverage, or an eventual spectrum threshold.

The reviewed proof is `reviewed-inputs/PROOF.md`, SHA256 `7e1ebb8fe9229bd4604f915c9ac23c8b1eacca6ffeb977b14957f4cdb9fe555b`. The producer verifier is SHA256 `6c27a517f59e1ad6798d34b2f0d848cf3b9da3a6c93014f8cdfa97dedc657212`. All seven entries of the producer manifest, plus that manifest, are preserved under `reviewed-inputs/`.

## Exact result and conventions

For the row (X,e,V,d), the physical letters are

    O: (X,e,V,d) -> (X,0,2V,d),
    Z: (X,e,V,d) -> (X+V-e,d,2V-d,d).

With 0<=e0<=d, d>=1, and A=V0-d>0, all finite words produce exactly

    X(u)=X0-V0+A*u+d*g(u)+(d-e0)*(1-(u mod 2)), u>=1.

The g convention is g(1)=1, g(2h)=g(h) for h>=1, and g(2h-1)=g(h)+h-(h mod 2) for h>=2. The proof correctly matches this inherited convention.

Chronological letters are the binary digits of k=u-1 from least to most significant, with O=0 and Z=1. The empty word has u=1. Any trailing O letters add leading zero digits and leave X unchanged. This use of u=1 denotes the empty continuation of an existing physical carrier; it does not insert a singleton into the different central source family.

The q recurrence, its relation q(k)=k+1-(k mod 2)-g(k+1), and substitution into the displayed X formula are correct. In particular the parity correction has the stated sign and is present on even u. The four expanded formulas and all four empty-word values agree with their initial rows.

The inverse bound follows from the exact nonnegative remainder

    X(u)-X0-A*(u-1)
      = d*(g(u)-1)+(d-e0)*(1-(u mod 2)) >= 0.

Thus every representation of a fixed target X has u<=1+floor((X-X0)/A), and X<X0 has no representation. This is a complete finite inverse for membership in these orbits, not an existence assertion for arbitrary X. No monotonicity or injectivity of X(u) is needed. Replacing A by V0 is unsafe: a single Z at each of the four seeds already violates that stronger proposed inverse.

## Independent derivation

For a chronological bit word b_0,...,b_(m-1), put k=sum b_i*2^i. Expanding the row dynamics directly gives

    V_i=2^i V0-d*sum_(j<i) b_j*2^(i-1-j),
    e_i=d*b_(i-1) for i>=1,
    X_final=X0+V0*k-e0*b_0-d*Q(k),

where the binary-pair polynomial is

    Q(k)=sum_(0<=i<j<m) b_i*b_j*(2^(j-i-1)+[j=i+1]).

The bracket is 1 precisely for adjacent bits. This provides an independent expansion without importing the producer's q implementation. Splitting off the lowest bit recovers its q recurrence. Combining with the specified g recurrence gives the claimed formula at every length. Leading zero padding contributes no binary pair, which also verifies trailing-O neutrality directly.

The independent symbolic checker expands the complete coefficient vector in (X0,V0,e0,d), so its word tests cover every numerical initial row rather than just the four seeds. Its separate g routine is iterative, and its correction is computed from binary pairs.

## Physical inheritance and density scope

The immutable paired-carrier dependency is the published proof at [commit 9ca780a3fe19d11846e7afa54c4f974f4567344d](https://github.com/eliottcassidy2000/math/blob/9ca780a3fe19d11846e7afa54c4f974f4567344d/session203/staggered-paired-carriers/PROOF.md), with SHA256 `e54792c206d6945a11dda2c4e6ca1e8d20140bbf3aaec91081987753c4a16f38`. These are the same bytes independently audited in the paired-carrier review. It supplies the eight explicit graphs, their roots, fixed starts and j values, all-order physical transfers, and the actual four-core readout. The g convention is frozen from the published central-arithmetic proof at commit `c8411ae1cc3748391d8ad2c4142cc68201b75f84`, SHA256 `49ad1123e80fbf751eac3e72ed0f5ad2ec39e7705f9ee89c77d4826ef5e1066c`. Both immutable URLs are recorded in the input provenance.

Every u therefore gives two actual rooted sources with the same X and the inherited fixed admissible shifts. Their four-core outputs have H=3X+2 and H=3X+4. No freely selected state coordinates are treated as physical.

Counting canonical integer indices, rather than binary-word lengths, gives

    C(N)<=sum_i max(0,1+floor((N-X_i)/(V_i-d_i))).

Distinctness of orbit values is unnecessary because this is an upper bound. All physical bases are odd. Dividing by the number of positive odd integers through N gives upper relative density at most

    2*(1/1012+1/898+1/1582+1/1046)
      = 693609731/93988567442 < 3/400.

This is less than 0.75 percent and is half the previously audited bound for this same four-orbit bank. It is a limsup bound, not a claim of a density limit or a matching lower bound. It excludes an odd tail from this bank alone, not from general tournaments or from other profile constructors.

The ancillary prime-atom discussion was removed before the final freeze. No such reduction or prime existence claim is needed for, or promoted by, this audit.

## Independent controls

Normal and optimized Python receipts agree on:

- 131,071 symbolic words of lengths 0 through 16.
- 65,536 binary-pair/g identities and positivity/parity checks.
- 880 direct long-word controls, up to length 256, including all-zero, all-one, alternating, and deterministic random words, with boundary and interior e0 choices.
- 80,016 target-by-target finite-inverse checks.
- 262,144 odd-base checks for the four formulas.
- Explicit hostile controls for reversed chronology, the unsafe V0 denominator, and dropping e0<=d.

The separate C++ graph checker uses literal directed-path enumeration on all subsets, independently of the producer's arithmetic verifier. Optimized and debug receipts agree on 120 actual source graphs, all common words of length at most 3 from all eight sources, 60 paired-row comparisons, and 24 literal four-core output graphs. Source order reaches 11 and output order reaches 12. It checks X, e, V=R+b, d, fixed starts, and H=3X+2r directly.

The producer verifier was also replayed in normal and optimized Python from its frozen snapshot. Both outputs match the original receipt exactly.

## Reproduction

From this audit directory:

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json
    g++ -O3 -std=c++17 check_physical.cpp -o check_physical
    ./check_physical > PHYSICAL_NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check_physical.cpp -o check_physical_debug
    ./check_physical_debug > PHYSICAL_DEBUG.json
    cmp PHYSICAL_NORMAL.json PHYSICAL_DEBUG.json
    python3 reviewed-inputs/verify.py > PRODUCER_NORMAL.json
    python3 -O reviewed-inputs/verify.py > PRODUCER_OPTIMIZED.json
    cmp PRODUCER_NORMAL.json PRODUCER_OPTIMIZED.json
    cmp PRODUCER_NORMAL.json reviewed-inputs/NORMAL.json

Independent symbolic checker SHA256: `aa014c1d2a8baa002990574c39cc3843250c001cdb1a1c7fae229f9d8b872b03`.

Symbolic receipt SHA256: `15d18edfaedfff6f1fc38e5edbf2082cb4948977a54730f3cc341d233e5cfc4a`.

Physical checker SHA256: `e89ada99762f1a9c31401e23eae18540bf35a16f0e8585bbbcc8c3552c495ca2`.

Physical receipt SHA256: `a41856a78e473b8c82f5269c9a1437b860edbc4f7a0d1ed2ef72e377bcff6dbf`.

The manifest freezes the reviewed snapshots, independent sources and receipts. It excludes local binaries. No producer input or shared checkout was edited.
