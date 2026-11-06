# Session 207: exact inverse sizes and the strong-atom reduction

Status: elementary PROVED consequences of SCC multiplication, with CITED Busch bound and FINITE-EXACT controls through order five. This does not solve global range completeness or give a simple closed form for every fixed-order spectrum.

Use nonempty finite tournaments; H counts directed Hamiltonian vertex words, not automorphism orbits. Put S_n={H(T): |T|=n}, and use S_0={1} only as an algebraic empty-product convention. Put A_r={H(C): C strongly connected, |C|=r}; A_1={1}, A_2=empty.

## 1. Every inverse size set is an upper ray or empty

Define m(h)=min{n>=1:h in S_n}, with infinity if no witness. Then, for every integer h,

    {n>=1:h in S_n} = {n>=1:n>=m(h)}.

Proof: append a new universal source. It must be first in any Hamiltonian path, so deleting it is a bijection with paths of T. Thus S_n is contained in S_(n+1), and every finite first occurrence pads to every later size. Conversely, no witness precedes the minimum by definition. In particular m(1)=1; nonpositive/even h and h=7,21 have m(h)=infinity. No other all-order omissions are claimed.

The inverse question therefore has no sporadic missing orders after first occurrence. The remaining problem is determining the single threshold m(h).

## 2. Exact fixed-order recurrence

    S_n = union_(1<=r<=n) A_r S_(n-r),
    AB = {ab:a in A,b in B}.

Proof: every tournament's SCCs form a linear dominance order. Any Hamiltonian path traverses the components in that order, independently choosing one internal path in each. Remove the first component to obtain the right-hand side. Conversely, take the dominance join of any witnesses from the two factors. This is an equality, not a heuristic interpolation from minimum and maximum counts.

## 3. Exact divisor/min-plus inversion

For d>1 let a(d) be the smallest order of a nontrivial strong tournament with H=d, or infinity. Define b(1)=0, and for h>1 define

    b(h) = min_{d|h, d>1} (a(d)+b(h/d)).

Then m(h)=b(h) for h>1, while m(1)=1. All infinities obey the usual extended arithmetic. The recursion decreases h because d>1. Singleton SCCs are never present in a smallest witness for h>1: deleting them preserves H and leaves a nonempty tournament. The first nontrivial SCC gives the lower inequality, and the dominance join of minimizing witnesses gives the upper inequality. Equivalently m(h) is the minimum sum a(d_i) over multiplicative factorizations h=product d_i with d_i>1.

One must allow composite d_i. Prime factorization alone is insufficient: m(9)=5 from a strong five-vertex tournament, while two triangle components cost six vertices. Prime h must have a single nontrivial strong component, so m(p)=a(p). More generally this holds for any h with no proper factorization into globally attainable values, including every attainable 7p for prime p!=3 (proper factors require the forbidden value 7).

## 4. Effective finite computation

First handle the boundary inputs: return m(h)=infinity for h<=0 or even h, and return m(1)=1. The remaining algorithm applies only to odd integers h>=3, for which both cutoff sets below are nonempty.

The inherited current spectrum note proves the reduced-witness bound B(h)=max{r:3^r<=h^3}. Therefore enumerate all tournaments of orders 1 through B(h), count H exactly, and take the first matching order; if absent, certify m(h)=infinity. Use integer powers, not floating-point logs. This is a terminating exact algorithm for every h, although not practical for large h.

The strong-atom route needs only tournaments through R(h)=max{r>=3:f(r)<=h}, where Busch's exact strong minimum is

    f(3+3k)=3*5^k; f(4+3k)=5*5^k; f(5+3k)=9*5^k, k>=0.

Enumerate strong tournaments through R(h), retain their H values dividing h, record a(d), then apply the divisor recursion. All required atoms divide h and have H<=h, so no larger strong tournament can matter. This can be smaller than enumerating complete reduced witnesses through B(h). It is still an exact finite algorithm, not a new infinite existence theorem.

## 5. Controls and provenance

Run `python3 check.py` from this directory. `check.py` exhaustively enumerates all labelled tournaments through n=5, counts paths by direct permutation enumeration, determines strong connectivity by transitive closure, and independently verifies the fixed-order SCC recurrence. It gives:

- S_1=S_2={1}.
- S_3={1,3}.
- S_4={1,3,5}.
- S_5={1,3,5,9,11,13,15}.
- First occurrence: 1 at n=1; 3 at n=3; 5 at n=4; 9,11,13,15 at n=5.

The finite cutoff is inherited, not claimed new: `05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md` section 4, fetched from GitHub main on this session, blob 56a5b0b8b0a3a264cef4e95eb947fbfb62dbac4f. That note also preserves the open all-odd-except-7-and-21 completeness status and its citation to Busch, *Electronic Journal of Combinatorics* 13 (2006), N3. The present note's contributions are the explicit inverse upper-ray formulation and its strong-atom min-plus recurrence, without a novelty claim.

Additional inherited provenance: `05-knowledge/results/spectrum_basins_20261020e/spectrum.md` sections 3–4 contains the same logarithmic bound and a positive rational series whose tropicalization yields the min-plus formulation above. Thus this session packages an explicit inverse-size characterization and algorithm from established mechanisms rather than claiming a new bound or independent discovery of that algebra.
