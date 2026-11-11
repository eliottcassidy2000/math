# Fast exact queries, a million-target atlas, and the state between them

**Session270, 11 November 2026. Research complete. Chan–Pak remains OPEN.**

The session delivers two practical query tools and new exact generating-function
bridges. It does not replace A000568 by a finite scalar recurrence or turn finite
Hamiltonian coverage into an all-target theorem. Earlier mechanisms were
inherited and extended, rather than subjected to another audit.

## 1. What “instantaneous A000568” can now mean

A000568(n) counts ordinary tournament isomorphism classes on n vertices. It is
neither a count of Hamiltonian-path values nor a census after identifying
converse. The exact Burnside input sums over odd permutation cycle types.

The new portable query tool has two regimes:

- **Saved exact queries, 0<=n<=250.** All 251 coefficients are included and were
  independently compared using the degree-local odd-factorial and factorial
  response recurrences. A fresh Python process returns a(250), including its
  full **8,878 decimal digits**, in a recorded median **0.0616 seconds**.
  Loading, hash-checking and parsing the whole table takes about 0.017–0.023s.
- **Fresh recurrence.** The C++/GMP engine computes a(120) in **0.0306s** and
  a(200) in **1.096s**, interleaved three-run medians including process launch
  and output capture. A single additional run computes a(250) in **8.713s**.

The table required 17.38s to build plus 45.87s for independent full-prefix
comparison. Filesystem caches were not flushed. Warm indexing excludes decimal
formatting and transmission; its nanosecond loop figure is not end-to-end
latency. Outside the saved range, the Python query falls back to fresh exact
recurrence and may be slow. Exact arbitrary-n output has Theta(n^2) bits, so
constant-time printing is impossible.

The implementation combines Session98's **already proved** fixed-mass
projective response with degree-local odd-factorial scaling, exact signed
power-of-two shifts, packed half-defects, late gauge switching and closed
terminal evaluation. The Python engine additionally contracts the remaining
{1,3} alphabet by a Laurent polynomial. The integrated invariant and signed
shift proof are in [enumeration/PROOF.md](enumeration/PROOF.md).

The comparisons are deliberately separated:

- At120 the new compiled engine is about25 times faster than the recovered
  historical C/GMP partition code (0.765s). That combines inherited response
  compression with implementation; it is not a new asymptotic theorem.
- The strongest recovered Python prefix engine takes2.502s at 200, but returns
  **all** coefficients rather than one. This is a practical query comparison,
  not identical output work.
- The same-language C++ raw/projective comparison at 200 is1.107s versus1.096s.
  State compression alone gives no robust wall-time win and loses slightly
  at 160. The isolated Python {1,3} terminal saves about2.5% at 200.

This is an honest practical advance, not a claim of polynomial state growth.
All full counts use exact integer, rational or GMP arithmetic. No small
relative tail is rounded to an exact integer; the inherited floating/absolute
error examples at 14 and 15 are explicitly covered by exact controls.

[Query tool, complete benchmarks and cache](enumeration/README.md).

## 2. Every eligible odd target through a million, with literal witnesses

The other query tool returns an actual tournament with a requested number H
of directed Hamiltonian paths. It uses the inherited two-core/transitive-spine
family G(2,0), not a newly invented family. A four-letter word chooses the
orientations from each ordered spine vertex to two fixed core vertices.

**FINITE-EXACT:** every positive odd H<=1,000,000 except7 and 21 has a saved
witness on at most16 vertices. There are **499,998 witnesses**.

All shorter words were exhausted first, so witness order is minimal within
this fixed marked language. It is not globally minimal; H=1 returns the
language's two-vertex core, although a singleton suffices. Strongness is not
required. Distinct words need not represent distinct unmarked classes.

The exhaustive run has477,218,583 prefix visits, including iterative-deepening
repeats. A deterministic compressed2MB packed atlas retains one word for each
odd target. After one decompression, direct-address lookup selects the word;
the compiler emits literal adjacency rows and recounts H with the transfer.
All499,998 saved words replay correctly. Additional independent subset-DP
controls cover626 atlas witnesses and 2,205 exhaustive small-core examples.

This converts finite existence into a usable compiler. It does not prove the
infinite Chan–Pak conjecture, global completeness of G(2,0), or an asymptotic
logarithmic witness-size bound.

## 3. The returning state is now explicit and exactly minimal over Q

Scalar H cannot safely be iterated through arbitrary contexts. The compiler
instead retains endpoint statuses of core-separated spine slots. This gives
a transparent **36-state nonnegative integer representation** for all words.
The new work is the explicit smaller realization and reusable compiler;
THM-4878 already supplied the family, the finite-state principle, a 200-state
realization and empirical rank17.

Exact reachable and observable invariant closures have dimensions28 and 20.
Their pairing has rank17, with a stored17-by-17 minor of determinant 150994944.
This certifies a **minimal 17-dimensional rational representation for the
entire infinite word series**, and the four matrices are exported. It upgrades
a sampled rank observation to an all-word proof. It does not prove that a
17-state nonnegative integer automaton exists.

[Transfer, rank certificate, million-target atlas and compiler](chanpak/README.md).

## 4. A new averaged generating function exposes the path-cover vector

For any fixed k-vertex core C, let q_r(C) count ordered covers by r nonempty
directed paths, and let A_m(C) sum H over all 2^(km) orientations between C and
an m-vertex transitive spine. With y=exp(2^k z),

    sum_m A_m(C) z^m/m! = (y+1)^2/4 * sum_(r=1)^k q_r(C)((y-1)/4)^(r-1).

Deleting the spine from a path leaves core runs. Required internal spine slots
and optional exterior slots give a Stirling-surjection transform; that is the
mechanism behind the formula. It converts the full cover polynomial into a
finite exponential sum, hence a rational ordinary GF and a constant-coefficient
recurrence in spine length. The zero-frequency term only contributes at m=0,
so the recurrence's starting index matters.

For the two-core family,

    A_0=1,
    A_m=(12^m+3*8^m+3*4^m)/8  for m>=1,
    A_m=24A_(m-1)-176A_(m-2)+384A_(m-3)  for m>=4.

The transform is triangular: A_0,...,A_(k-1) recover every q_r. A first
q_r difference is visible with r-1 extra spine vertices per graph, rather
than Session268's2r-2 auxiliary vertices in one homogeneous context. The
price is summing2^(k(r-1)) orientations. This is a precise vertex/query tradeoff,
not a violation of the earlier homogeneous-module lower bound.

## 5. A genuine junction with enumeration, with the markings retained

The same arbitrary-core/ordered-spine carrier has a class-counting observer.
Let B(k,m) count classes where isomorphisms preserve the core and fix every
ordered spine position. For odd cycle partitions lambda of k,

    B(k,m)=sum_lambda 2^e(lambda) 2^(m ell(lambda))/z(lambda).

Thus B(k,0)=A000568(k). For fixed k its spine-length OGF is rational, with
simple poles at 2^(-t), t=k,k-2,...,1 or2, and exact minimal recurrence order
ceil(k/2). For example B(3,m)=(2*2^m+4*8^m)/3.

This is a marked census, not the global unmarked census at order k+m. The
hard k-dependent inventory has not disappeared. It explains why the
construction language can support simple extension recurrences without
contradicting the known non-P-recursiveness of unrestricted A000568.

Averaging the path formula over all labelled cores gives another exact
simplification. In the uniform labelled ensemble with the spine fixed,

    E H = k! 2^(1-k) E[X^m],  X~Binomial(k+1,1/2).

The class-count recurrence and path-mass recurrence are different observers
of the same physical graphs. Neither determines H-value support. Uniform
ordinary-class averaging would require the missing automorphism weights.

[Common carrier, proofs and direct small-universe controls](carrier/PROOF.md).

## 6. Two further exact acceleration regimes

**Arithmetic low bits.** For even n>=4 set U_n=a(n)/2^(n/2),
c=phi(n)/2, A=sum_(1<=a<n/2,gcd(a,n)=1) a, and
D=[3|n,n>=12]phi(n/3)/2. Then

    U_n = nA-c+4D+6[n=4]+2[n=6] (mod8).

A Möbius divisor sum computes A using only squarefree divisors. Given the
factorization, the query takes O(2^omega(n)) arithmetic operations, excluding
factoring. It answers huge indices such as2^1000 and 2*3^1000 without forming
the full count. Repeated-cycle exceptions at 4 and 6 are essential. This extends
the inherited odd-order precision engine; the O(sqrt(b)) cycle cutoff and
odd-order congruences were already THM-5098.

**Constant-gcd sectors.** If feasible distinct cycle lengths share gcd g,
strip the common energy g*binom(K,2), convolve small bivariate kernels in mass
and cycle number K, then restore that energy. The method has polynomially
many states in this restricted sector. More generally the residual gcd graph
factors over components; edges whose endpoint masses exceed the degree budget
can be deleted formally. The Gaussian common mode becomes an exact deterministic
convolution, not a merely suggestive analogy.

For the composite support3*{1,3,5,7,11,13,17,19,23,29,31,37,41,43}, the prefix
through degree900 took15.40s versus34.49s in an identical-state conventional
response control: **2.24 times faster** with the same1,694,066 transfers.
This is a weighted Burnside sector, not the complete A000568 sequence or the
count of a selected tournament subclass. The universal version gave only
modest, inconsistent speed gains. Scope is the result, not an afterthought.

[Huge-index residue proof](bridges/PROOF.md) · [Sector convolution](terminal/PROOF.md).

## 7. What remains worth doing

The practical next step for full counting is dynamic detection of cheap residual
gcd components and choosing terminal contraction versus transfer by an actual
cost model. A new asymptotic state-count bound is still missing.

The constructive next step is to invert the17-state semigroup or find a
returning positive cone with a decreasing digit algorithm. The million-target
atlas is a source of possible structure, not evidence sufficient for an
all-target claim. The averaged rational GF controls mass, not range.

[Open questions](OPEN_QUESTIONS.md) · [Reproduction and validation scope](REPRODUCE.md).
