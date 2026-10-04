# Golden digits, arithmetic memory, and two additional return families

Session begun 2026-10-03, completed 2026-10-04 (America/Denver).

**PROVED:** the elementary reader, carry-memory, Fourier-transfer and guarded
first-descent statements below, with proofs in this note or its linked
companions. **FINITE-EXACT:** the declared digit, ring, head-ranking and full
route audits. **CITED:** classical golden numeration and the historical
cyclotomic fact. **OPEN:** primality classifications beyond the stated ring
criteria, universal successful certification, and Collatz convergence. No
novelty claim is made for the classical ingredients.

## Inheritance, scope and live board

The user's anchor is a structured family retaining its route home; the niche
is exact base-phi arithmetic and prime quotients; the wildcard is Fourier
transport of digit states. The live board is: literal digits, marked
Zeckendorf expansions, reversible carries, residue clocks, timed routes, and
new coverage. The inherited LRC(14) anchor remains OPEN; these results supply
no new LRC bound.

The closest proved mechanism is source-preserving affine composition in
[the recursive entry grammar](entry_20260927_recursive.md), together with
[THM-4528, golden parity and holonomy](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md)
and its [audited proof](collatz_golden_holonomy_20261001.md). The recent
[modular/Fourier/Zeckendorf synthesis](modular_fourier_zeckendorf_routes_20261003.md)
supplies the exact guards and the corrected ordered extra-unit grammar.
The canonical hostiles are forgetting the original descent threshold and
forgetting a digit's radix or elapsed time. The corrected near miss is the
signed norm in the [September 30 circulant note](collatz_circulant_20260930_circulants_lucas_cubic_monotile.md):
N(phi^4-1)=-5, although its quotient has size 5. That error is repaired in
place and recorded in MISTAKES. The least-used sidecar was the integer
polynomial discarded by carrying.

Companion results and reproducible implementations:

- [Literal digits and reversible carries](golden_digit_carry_20261003.md).
- [Prime quotients and exact residue clocks](golden_prime_clocks_20261003.md).
- [A sufficient head/cycle compiler and added coverage](collatz_mixed_return_compiler_20261003.md).
- [Cross-interface experiment](../../04-computation/experiments/golden_route_interface_20261003.py)
  and its [saved output](golden_route_interface_20261003.out).

## 1. The source identities fit a general law, with a primality boundary

Put O=Z[phi], phi^2=phi+1, and L_k=phi^k+(-1)^k phi^-k. Then

    phi^(2k)+(-1)^k=L_k phi^k.

The user's identities are exactly k=4,5, with L_4=7 and L_5=11. The next
instance has L_6=18: sparse base-phi formulas also occur for composites.
The useful common reader represents A+Bphi by the pair (A,B):

    (A,B) --append d--> (B+d,A+B),
    M=[[0,1],[1,1]],  det M=-1.

The carry 11_phi=100_phi is the kernel relation f(t)=t^2-t-1. It preserves
value but forgets which expression supplied it. Bergman's primary
[1957 paper](https://math.berkeley.edu/~gbergman/papers/base_tau.pdf) supplies
the classical numeration setting; [Shallit's automata paper](https://arxiv.org/abs/2305.02672)
is a relevant existing formalization direction. Our claims concern the
explicit exact readers and certificates here, without asserting a new
global normalization theorem.

Three readings of binary digits must retain their types:

| Input | Exact object read | Essential extra coordinate |
| --- | --- | --- |
| Positional base-phi digits of n | A+Bphi=n, hence (n,0) | Radix and carry/factor provenance |
| Ordinary Zeckendorf digits of n, read as phi powers | n-2b(n)+b(n)phi, b(n)=floor((n+1)/phi^2) | The second register and the extra-unit marker |
| Ordinary Collatz parity itinerary of n | Theta(n)=sum b_j phi^(-j-1) | Time, legal edges, endpoint and certified suffix |

For the second row, (X,Y)=(A+2B,2A+3B)=M^3(A,B) is an invertible integral
change of coordinates, exactly recovering the earlier Fibonacci reader.
Its projection to X is not multiplicative: the lift of 2 is phi, whose
square is the lift of 3. Literal positional digits preserve ordinary
multiplication; the Fibonacci lift preserves a different structure.

## 2. What 105 retains after multiplication and normalization

Exact arithmetic gives

    105 = 1001010101.0101001001_phi,

with occupied exponents 9,6,4,2,0,-2,-4,-7,-10. Let sigma(n) be its selected
normal Laurent polynomial, C=sigma(3)sigma(5)sigma(7), and B=sigma(105).
Then

    t^10(C-B)=f(t)Q(t),
    Q=-t^4 Phi_3(t) Phi_6(t)^2 Phi_12(t).

This is a checked factorization of the carry record. The tuple (B,Q,10)
recovers the raw convolution C exactly, including repeated coefficients;
it does not recover a chosen factor tree. Indeed Q(phi)=-16phi^8 has norm
256. Its specialization is not a recovery of the factors 3,5,7.

**PROVED arithmetic memory law.** For n,a,b>=0, write a decorated state
(n,q), with q an integer Laurent polynomial, for the raw expression
sigma(n)+f q. Define

    A(a,b)=(sigma(a)+sigma(b)-sigma(a+b))/f,
    C(a,b)=(sigma(a)sigma(b)-sigma(ab))/f.

The divisions are exact in Z[t,t^-1], by the minimal polynomial of phi.
Arithmetic on decorated states is

    (a,q)+(b,r) = (a+b, A(a,b)+q+r),
    (a,q)*(b,r) = (ab, C(a,b)+q sigma(b)+sigma(a)r+fqr).

Expanding sigma(a)+fq and sigma(b)+fr proves both identities and shows that
they preserve the complete raw sum or product. This provides composable
aggregate memory of normalization. If chronology or factorization matters,
retain its operation DAG as well; q alone intentionally merges those histories.
The companion proves the additive and multiplicative associativity laws.

There are two other precisely typed roles for 105. Modulo 105, the shift M
has order lcm(8,20,16)=80, from its orders modulo 3,5,7. Separately,
Phi_105 is the first cyclotomic polynomial with a coefficient outside
{-1,0,1}, as explained in [Garrett's primary note](https://www-users.cse.umn.edu/~garrett/m/algebra/notes_2023-24/105th_cyclotomic_poly.pdf)
and independently checked through index 105. These are a modulus and a
polynomial index; neither is identified with the preceding carry polynomial.

The clock note also gives a concrete obstruction to prime-index reasoning:

    Phi_23(phi)=phi^12(11+2phi)(21+phi),
    N(11+2phi)=139,  N(21+phi)=461.

A binary repunit with prime length can therefore factor in O. Conversely,
Phi_3(phi)=2phi^2 has composite norm 4 but prime ideal (2), with quotient F_4.
Integer primes, prime indices, and prime elements are separate predicates.

**Concurrent exact connection.** The incoming [signed golden-carrier note,
section 4](collatz_golden_carriers_20261004.md) identifies the denominator
ideal of Theta(-5)=(-1+7phi)/11. It is precisely

    (Phi_5(phi))=(phi^5-1)=(phi-4)=(11,phi-4).

Indeed Phi_5(phi)=phi(phi^5-1), phi^5-1=-phi^3(phi-4), and both displayed
multipliers are units. The quotient is F_11 with phi mapped to 4, of order
5. This connects the binary repunit, monodromy quotient and golden
denominator module by an actual ideal equality. The rational 1/13 cycle
has coordinate (1+4phi)/11 and the same ideal: multiplying these two
coordinates by phi-4 gives respectively 1-2phi and -phi. Neither coordinate
is integral, since its coefficient pair is not divisible by 11. Its
denominator ideal is therefore proper and contains the maximal index-11
ideal (phi-4), forcing equality. The numerator pairs themselves are not
units modulo 11: their norms are -55 and -11. An ideal and its clock
therefore do not certify integer realization
of the ordered word. Only these elementary identities are independently
imported here, not the incoming larger denominator-gate census.

## 3. The three-colour connection has a faithful register and a boundary

The user's ordered red/black/blue construction keeps its marked extra copies
of 1, as specified in the [guarded Zeckendorf note](zeckendorf_guard_automaton_20261003.md).
The M^3 coordinate change proves its arithmetic-register connection to O.
It does not identify that ordered grammar with every three-state colouring.

Modulo 2, O/(2)=F_4, phi^3=1, and the three nonzero states cycle under M.
Zero remains a fourth state. An additional useful necessary condition follows
for literal integer words. Let c_j be the parity of the number of occupied
exponents congruent to j modulo 3, including negative exponents. Since
phi^2=phi+1 in characteristic 2, their represented value is

    (c_0+c_2)+(c_1+c_2)phi mod2.

Thus an exact integer n necessarily has

    c_1=c_2,  c_0+c_2=n mod2.

Odd integers have phase signatures (1,0,0) or (0,1,1). These are only
necessary congruences: phi^3=1+2phi passes the same test without being an
integer. They also do not detect primality: 3 and 9 have signature (0,1,1),
while 5 and 105 have (1,0,0). The exact second coordinate, marker, and digit
positions each repair a different information loss.

## 4. Fourier modes now transport through the digit reader exactly

For v=(A,B) modulo m define

    chi_(r,s)(v)=exp(2*pi*i*(rA+sB)/m),
    D_d(v)=Mv+d(1,0).

**PROVED:** direct substitution gives

    chi_(r,s)(D_d v)=exp(2*pi*i*r*d/m) chi_(s,r+s)(v).

The digit reader permutes modes and adds a known phase. In real coordinates,
that phase rotates cosine and sine into one another. This is the concrete
connection to the earlier sine/cosine discussion: each output character is
a known phase times its specified input character. Retain the full mode
array or the required frequency orbit, rather than one fixed-mode scalar.

**Hostile to cosine-only storage.** Modulo 3, states (1,0) and (-1,0) have
the same cosine at every mode. After appending digit 1 they become (1,1)
and (1,2); mode (1,1) sees phases 2 and 0, so its cosine differs. A
deterministic digit update cannot recover this information from cosines
alone. Full phase, or equivalently the signed sine component, resolves it.

For exact binary/ternary guards, the companion proves

    ord_(2^a 3^p)(phi)=2^max(a-1,3) 3^max(1,p-1),  a,p>=1.

This lets the modular reader store the radix correction as a finite phase.
The exact source still retains its actual radix and magnitude. A full clock
shift preserves residues but can destroy literal rational integrality.

## 5. The creative connection that added proved coverage

The [compiler theorem](collatz_mixed_return_compiler_20261003.md) joins an
expanding legal head to repetitions of an actual negative cycle, then pays
for descent below the original positive source. Its guards, positivity
domain and final extra division are explicit. This implements a restricted
but useful operation on certificates instead of relying on visual similarity.

For the -17 valuation block w=(1,1,1,2,1,1,4), put P=2187,Q=2048.
The bounded search inspected all 13 expanding heads of length 1..4 with
letters 1..4, at repetition counts 1..3. Head (1) gave the largest added
density in that declared universe. Its all-height theorem is:

For m>=1 choose the least t>=1 such that

    2^t(2Q^m-35)>3P^m-51.

For every positive odd u satisfying

    uQ^m=1 mod3,  uP^m=17 mod2^t,

the source n=(2uQ^m-35)/3 has exact first odd descent after 7m+1 steps,
to oddpart(uP^m-17). Each m gives a full dyadic cylinder; distinct m are
disjoint, and v2(3n+35)=11m+1 recovers the repetition count from n.
The first cylinder is n=6815 mod8192. The second retained head (1,2)
begins at n=27291 mod32768 with exact first odd descent at step 9.

Both infinite families miss the named old 171-row/65-cylinder bank and
its -5 and -17 extensions, and they miss one another. Their additional
densities among all positive integers are respectively approximately
0.00012212994626282017 and 0.00003053248656570591. The proof retains exact
partial sums, tail bounds and shrinking parent cylinders; these are source
densities, not orbit visit frequencies. This establishes new coverage
relative to those named certificates, with no literature novelty claim.

A deliberately weakened division budget produces explicit failures, despite
preserving the proposed legal word. This is why the source threshold belongs
inside the stored certificate.

## 6. Complete routes, tournament carriers, and a concrete storage test

For a legal finite ordinary Collatz prefix with bits b_0,...,b_(L-1), set
G=sum b_j phi^(-j-1). At its actual endpoint r,

    Theta(n)=G+phi^-L Theta(r).

The pair (G,L) composes by (G,L)*(H,K)=(G+phi^-L H,L+K), provided endpoints
match. A certified suffix to 1 gives a complete home route. Since
Theta(1)=phi/2, twice the whole certificate has integer coordinates in O.
The pair (G,L) verifies a projection of a route; retain the actual legal
word/guard and endpoint because evaluated G alone does not reconstruct it.
For example, prefixes [1] and [1,0] from source 3 both give G=phi^-1 but
end at 10 and 5 after different lengths.
This loss statement concerns the finite prefix G with its length discarded.
It does not assert that the full infinite Theta, equipped with its boundary
convention, loses the admissible itinerary; the incoming carrier note makes
that distinction explicit.

The cross-interface script independently tests both retained heads, m=1..6,
and coefficient lifts 0,1,2: 36 sources. It reads their literal digits,
verifies the dyadic guard and radix phase, rejects neighbouring sources,
checks the numerator carry, and splices an independently iterated suffix.
Every one reaches 1 within the declared 10,000-step suffix bound. This is
a FINITE-EXACT full-route result, distinct from the all-height first-descent
theorem.

| Source | Certified endpoint | Ordinary prefix steps | Suffix steps to 1 | Total |
| ---: | ---: | ---: | ---: | ---: |
| 6815 | 5459 | 21 | 160 | 181 |
| 27291 | 8197 | 25 | 158 | 183 |

The [inherited route tournament](collatz_atom_route_memory_20261003.md)
uses each compressed inverse index j as a cyclic strong component of order
2j+3, followed by the terminal singleton. The script recovers all 36
first-hit routes from their ordered component sizes, excluding root padding.
The largest carrier has 942 vertices. This uses an existing intrinsic
tournament relation; symmetric product-table entries are not reoriented.

Sharing exact suffixes reduces 4,607 odd-edge occurrences in those 36 routes
to 3,315 distinct stored odd edges. This measures one concrete storage gain,
without claiming a complexity bound or identifying equal-valued colour
markers. A small certificate object can retain a head/cycle/count, actual
terminal division, source guard/budget, and a pointer to the verified suffix.
Arithmetic carry records and marked source expansions are additional views
when the provenance of the number itself matters.

## 7. Next focused steps and decisive tests

1. Extend the compiler by certificate composition while retaining the original
   source threshold. Rank candidates by exact coverage gain and description
   size, after subtracting the named existing union. Reject legal-but-unpaid
   candidates using the stored weak-budget hostiles. The smaller winning head
   already demonstrates why this ranking is more informative than long traces.
2. Build the suffix library as a directed acyclic proof graph rooted at 1.
   Add an entire family only with a proved constructor or a verified suffix
   schema. An unmatched endpoint stays a proof obligation. Measure coverage
   and shared-edge savings separately; neither density nor compression proves
   complete coverage.
3. Minimize the synchronized digit/guard automaton using full phases and exact
   acceptance predicates. Quotient states only after checking that future
   digit continuations preserve acceptance. The mod-3 cosine collision and
   the mod-2 integrality hostile are inexpensive rejection controls. Keep the
   marked extra-unit grammar as a separate language unless a lossless map is
   proved.

These steps emerge from the successful operations above. Searching for a
prime-only visual pattern or a route rank derived solely from an untimed
finite parity sum has met explicit obstructions and is not the next priority.

## Reproduction and independent review

Run the four scripts with the common stem date 20261003:

    python -X utf8 -B 04-computation/experiments/golden_digit_carry_20261003.py
    python -X utf8 -B 04-computation/experiments/golden_prime_clocks_20261003.py
    python -X utf8 -B 04-computation/experiments/collatz_mixed_return_compiler_20261003.py
    python -X utf8 -B 04-computation/experiments/golden_route_interface_20261003.py

The carry script uses assertions and must be run without -O. The other three
also pass with -O and retain explicit checks. Their companion notes declare
each finite universe and hostile controls. The interface adds 79,946 exact
Fourier checks, 1,001 exponent-signature checks, and 128 decorated arithmetic
checks. All full routes are independently replayed by the ordinary map;
no floating point is used for acceptance or proof conditions.

Independent readers audited the clock, carry, compiler and interface claims.
The all-height conclusions rely on the written proofs, not on extrapolation
from the tested instances. The exact outputs are retained beside the notes.
