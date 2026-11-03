# Digit rules, legal division and the price of preserving a source

**PROVED elementary interface theorems DM1–DM6.** Finite controls are separate
from the proofs. No literature-priority or universal Collatz claim. Domain:
ordinary integer bases B>=2, moduli positive; digit suffixes may have leading
zeros. Sources are positive integers unless a complete residue universe is
explicitly used. An invalid partial operation is an observed failure, not an
arbitrary transition. Quotients below depend on the declared operations.

## DM1. A common engine for positional divisibility rules

Let m=s u, where s contains exactly the prime-power factors of m whose primes
divide B; thus gcd(u,B)=gcd(s,u)=1. Let L be the least nonnegative integer
with s dividing B^L. Then m divides n iff both of these tests accept:

* s divides the value of the last L base-B digits;
* u divides the weighted digit sum sum_i d_i B^i modulo u.

For u>1 the second weights are periodic, with least period ord_u(B).
If B^ell=1 mod u, use sums of ell-digit blocks. If B^ell=-1 mod u, use
alternating sums of those blocks. The latter option need not exist. For a
right-trimmed digit n=B a+d, the equivalent coprime test is

    u | n  iff  u | a+(B inverse mod u)d.

The full remainder update, reading from left to right, is r -> B r+d mod m.
All statements follow from positional expansion, invertibility modulo u,
and the Chinese remainder theorem. In particular no inverse of B modulo s
is used. Decimal examples: modulo 7 trim to a-2d; modulo 13 to a+4d;
modulo 11 alternate digits; modulo 16 the last four digits suffice.

For a precise completeness statement about *positional congruence identities*,
evaluation at B gives

    (Z/mZ)[X]/(X-B) = Z/mZ.

Indeed polynomial division by the monic X-B gives remainder P(B). Thus any
replacement of digit polynomials preserving the residue differs by a multiple
of X-B in this ring. This organizes all such congruence rules, including
block and cyclic-weight versions. It is not a classification of every
conceivable nonlinear divisibility algorithm or its optimal running time.

For positive n, the base-B digital root is

    dr_B(n)=1+((n-1) mod (B-1)).

Zero has its separate usual convention. A finite bank of digital roots in
bases B_i is equivalent, on positive integers, to n modulo lcm_i(B_i-1).
This follows by CRT compatibility, including noncoprime moduli. It does not
retain the unbounded quotient above that modulus. With all bases the bank
does distinguish integers: choose B>max(n,n')+1.

### The multiplication cycles behind the digital-root diagram

For multiplier a>0 on Z/m, use CRT prime-power components p^e. If p divides a,
the capped valuation of a^j x is min(e,v_p(x)+j v_p(a)); this component reaches
zero. If p does not divide a and v=v_p(x)<e, its exact period is
ord_(p^(e-v))(a). A zero component has period one. The global transient is
the maximum time to zero in the nonunit components; its eventual period is
the least common multiple of the unit-component periods. This is proved
by cancelling p^v and applying CRT.

For decimal doubling, this gives modulo 9 the cycles

    (1,2,4,8,7,5), (3,6), (0).

The displayed digital root 9 stands for residue zero. In base 8 the modulus
is 7 and doubling instead gives (1,2,4),(3,6,5),(0); in base 12 its nonzero
modulo-11 cycle has length ten. The six-cycle is therefore a property of
the chosen multiplier and modulus, not a base-independent number partition.
The repo's four-vertex tournament “vortices” are a different definition.

## DM2. The exact smallest divisibility observer

Consider prefixes n, followed by arbitrary base-B suffixes, with only the
Boolean response m|result. All future responses depend on r=n mod m, but
the m residues need not be distinguishable.

Use s,u,L from DM1, and choose K>=L with B^K>=m. For each 0<=k<K put

    a_k(r)=(-B^k r mod m),
    c_k(r)=a_k(r) if a_k(r)<B^k, and bottom otherwise.

Then the following signature is exactly the future-response equivalence:

    Sigma(r)=(r mod u, c_0(r),...,c_(K-1)(r)).

**Proof.** A k-digit suffix has value t in [0,B^k). It accepts precisely
when t=a_k(r) mod m. This set is empty if a_k(r)>=B^k, otherwise is the
truncated progression starting at a_k(r). Thus c_k records exactly its
acceptance set. For k>=K, B^k kills s and is invertible modulo u, so r mod u
determines B^k r mod m and every suffix response. Conversely, if r mod u
differs, at k=K the two least admissible suffixes differ and both are <B^K;
one separates the sources. If a short coordinate differs, its two acceptance
sets differ and a short suffix separates them. QED.

Every residue has a digit representation, even allowing positive prefixes
only. Thus distinct signatures are all reachable, and their separating
suffixes prove that no deterministic exact divisibility machine has fewer
states. Appending a digit maps equivalent prefixes to equivalent prefixes,
so the signature classes themselves attain the bound. This proves minimality
without assuming all m residues are necessary. For gcd(B,m)=1 there are m
classes; in binary, divisibility by 4 needs three classes {0},{2},{1,3}.
The exact finite checker also finds nine for decimal divisibility by 16,
and four for base-12 divisibility by 18. These latter numbers are finite
instances, not a general state-count formula.

## DM3. Appending digits and legal prime division: a sharp dichotomy

Fix a prime p. Allowed actions are appending any base-B digit and replacing
n by n/p only when p divides n. Observe whether a proposed division is legal.

**Theorem.** On positive integers:

1. If B=p^a for some a>=1, two sources have identical responses to every
   context iff their p-adic valuations agree.
2. If B is not a power of p, distinct sources are always separated by a
   finite context. The minimal faithful quotient is the full source set.

**Proof of 1.** Division decrements a positive valuation and is illegal at
valuation zero. Appending digit zero increases valuation by a; appending
0<d<p^a changes it to v_p(d), independently of the old source, since
v_p(d)<a<=v_p(p^a n). These rules close on valuation alone. Repeated legal
division separates different valuations, proving both directions.

**Proof of 2, with an explicit separator.** Write B=p^a C with gcd(C,p)=1
and C>1. For n!=n' set r=v_p(n'-n). Choose k>=1 with C^k>=p^(r+1), and set

    t=(-C^k n mod p^(r+1)),   d=p^(ak)t.

Then 0<=d<B^k, so d is a genuine k-digit suffix, leading zeros permitted.
After appending that same suffix to both sources,

    v_p(B^k n+d)>=ak+r+1,
    v_p(B^k n'+d)=ak+r.

The first is construction; their difference has exactly valuation ak+r,
giving the second. The (ak+r+1)st proposed division separates them. QED.

Example: 1 and 3 cannot be separated by binary appends and legal halving,
but decimal suffix 6 sends them to 16 and 36, with halving clocks 4 and 2.
The older decimal-halving theorem is the case (B,p)=(10,2); its proof did
not imply the pure-prime-base case, which has a genuinely smaller quotient.
Allowing legal division by two distinct primes makes the quotient injective
in every base, because a base cannot be a power of both primes.

These are operation-relative statements. Adding an ordinary numerical-value
observer, subtraction, or a zero test is a different interface. Counting the
unbounded valuation register as “one state” would be a resource error.

## DM4. Exact modulus required by a guarded affine quotient

For integers a,b and d,M>0 define

    F(n)=(a n+b)/d mod M if d divides a n+b, and bottom otherwise.

If gcd(a,d) does not divide b, F is constant bottom, of period one.
Otherwise its least positive translation period is

    P=dM/gcd(a,dM).

**Proof.** There is a legal source n. A translation P preserves legality
and response there iff dM divides aP. This condition also suffices at every
legal source. It preserves legality of every source because d divides aP;
therefore it preserves the illegal set too. The least positive solution is
the displayed formula. Positive sources alone give the same necessity by
choosing a positive legal representative. This also covers a=0. QED.

For a half-step even branch, the required period is 2M. For the odd Collatz
branch (3n+1)/2, whose integrality guard is precisely oddness, it is
2M/gcd(3,M). Division spends a binary digit; multiplication can conceal
ternary source phase. A translation period is not automatically the number
of contextual states (DM2 is a counterexample to that identification).

## DM5. Exact joint parity-history and endpoint-phase size

Use T(n)=n/2 on even n and (3n+1)/2 on odd n. This is the shortcut map.
It continues on the 1,2 cycle; no strict first-ROOT convention is implicit.
The inherited parity-cylinder theorem says that each word w of length k
with o odd symbols has exactly one source class r_w mod 2^k and

    T^k(n)=(3^o n+b_w)/2^k.

The affine coefficients are recovered from (A,C,D)=(1,0,1) by
(A,C,D)->((1+2e)A,(1+2e)C+eD,2D). The residue is -b_w (3^o)^(-1) mod 2^k.
The two lifts of a length-j source class have endpoints differing by the
odd integer 3^(number of odd steps), hence give opposite next parities;
this proves the inherited uniqueness and existence directly.

Observe the pair (whole word w, T^k(n) mod M). Its exact number of responses
over all positive integers is

    N_k(M)=sum_(o=0)^k binom(k,o) M/gcd(M,3^o).              (1)

Each response in a word of odd weight o has exactly gcd(M,3^o) source
residues in the universal period 2^k M. To prove this, write

    n=r_w+2^k t,  T^k(n)=d_w+3^o t,
    d_w=(3^o r_w+b_w)/2^k.

As t varies modulo M, this affine phase has M/gcd(M,3^o) images, each
with the stated multiplicity. Every source residue has a positive lift.
Summing over words proves (1). The ordinary Collatz clock is exactly k+o:
each odd shortcut combines one 3n+1 step and one halving. The clock and
full word are retained, not replaced by the endpoint alone.

Writing M=3^a u, gcd(u,3)=1, gives the closed form

    N_k(M)=u [2^k + sum_(o=0)^(a-1) binom(k,o)(3^(a-o)-1)]. (2)

Out-of-range binomials are zero. In particular

    N_k(3)=2^k+2,
    N_k(9)=2^k+2k+8,
    N_k(27)=2^k+2 binom(k,2)+8k+26.

If k<=a then N_k(M)=M(4/3)^k, an integer. If 3 does not divide M,
N_k(M)=2^k M. These are cardinalities of the declared observation, not
counts of convergent sources and not independent bits supplied for free.
At k=6,M=9, the 576 source residues give exactly 84 responses; the
nonuniform fibre sizes are 1,3,9 according to the word's odd weight.

There is an exact weighted statement if, and only if, a probability model
is specified. Give the 2^k M source residues equal weight. A response in
word w then has probability gcd(M,3^o)/(2^k M), not 1/N_k(M). For
M=3^a u the conditional entropy of the original residue given the response
is exactly

    log_2(3) * 2^(-k) sum_o binom(k,o) min(a,o).

This follows because each conditional source fibre is uniform of size
3^min(a,o), while words themselves have probability 2^(-k). At k=6,M=9
the lost entropy is (15/8)log_2(3). This is finite uniform-residue sampling,
not a probability distribution on all positive integers or an independence
claim along a marked Collatz orbit.

The period 2^k M is least for this combined observation: any period must
preserve the word, hence is a multiple of 2^k, and on the all-even word it
must change the endpoint by a multiple of M. This period can be larger
than N_k(M). The incoming positive-realization theorem implies that an
exact nonnegative hidden-state representation of every pure-source response
needs at least N_k(M) states; a signed-linear rank has a different meaning.

## DM6. A changing terminal modulus repairs nonnested refinement

Increasing k at a fixed endpoint modulus M need not refine the old source
partition. For M=9, sources 3 and 9 have the same one-step response
(odd,5 mod9), despite distinct original phases 3 and 0. More elapsed time
can discard source information. Cardinality growth is not partition refinement.

Fix a source-tail modulus M0. On a word w of length k and odd weight o,
we want to retain n modulo 2^k M0, equivalently t modulo M0 in the DM5
source parameterization. The smallest positive endpoint modulus accomplishing
this, given w, is

    M_w = M0                  if 3 does not divide M0,
    M_w = 3^o M0              if 3 divides M0.             (3)

**Proof.** Endpoint equality modulo M identifies parameters whose difference
is a multiple of M/gcd(M,3^o). It retains t mod M0 iff

    M0 divides M/gcd(M,3^o).

Compare valuations prime by prime. Away from 3, M must contain M0's full
prime powers. If v3(M0)=a>0, it additionally needs v3(M)>=a+o; when a=0
no factor 3 is necessary. This proves the divisibility-minimal and numerical-
minimal solution (3). There are exactly M0 attainable terminal phases per
word, even when the ambient terminal modulus is much larger. QED.

Consequently the adaptive response (w,T^k(n) mod M_w) is *exactly equivalent*
to n mod 2^k M0. It has 2^k M0 classes, and as k grows these form a genuinely
nested source partition. Every new level has exactly two source lifts per
class. This gives a precise changing-modulus rule, preserving original
source residues and integer clocks. For a decimal digital-root bank M0=9,
the needed terminal modulus is 3^(o+2); for M0=7 it remains 7. The rule
preserves finite source information, not the whole unbounded source integer.
Only M0 terminal phases are attainable in each word even when M_w grows:
increasing the ambient modulus does not supply that many independent states.

An associated exact precision ledger is useful. On a length-k word with o
odd steps, a source difference 2^A 3^B t, A>=k, becomes

    2^(A-k) 3^(B+o) t.

Thus the 2-adic digit loss and 3-adic gain are simultaneous. With the usual
prime-adic norms, the slope has magnitudes 3^o/2^k, 2^k, 3^(-o) at the
real, 2-adic and 3-adic places, whose product is one. This is exact rational
arithmetic bookkeeping; it is neither an independence assumption nor a
negative-real-drift theorem.

### The surviving Collatz obligation

The exact source descent test remains

    (2^k-3^o)n > b_w.

The changing-modulus controller does not imply this inequality. The inherited
hostile family n=q2^k-1 has all-odd prefix and endpoint q3^k-1, and every
intermediate step grows for n>1. Every finite depth therefore has infinitely
many positive hostile sources. The new result repairs source-phase refinement;
it does not settle termination, uniform ROOT debt, or the tournament {7,21}
spectrum. A useful proof must additionally pay the actual source's height
comparison or provide another well-founded certificate.
