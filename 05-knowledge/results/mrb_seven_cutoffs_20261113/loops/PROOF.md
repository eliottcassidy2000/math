# Exact finite tails and loops of transitive tournament context profiles

**PROVED IN NOTE; independently FINITE-EXACT checked on 270 moduli.**
The result classifies one infinite physical tournament family under *every
uniform exterior test modulo a fixed integer*. It does not classify the
whole tournament metagraph, preserve integer Hamiltonian counts across a
loop, or prove Chan–Pak. No literature-priority claim is made.

## Inheritance and the connection contract

The closest proved mechanism is the [Session295 exact contextual quotient](../../session295_construction_atlas_20261113/bridges.md):
write `F_k(C)` for the number of ordered spanning covers of a nonempty
tournament C by k nonempty directed paths, and

    D_k(C)=(k-1)! F_k(C).

For every uniform exterior E (each exterior vertex sees the entire module
C in one direction),

    H(E[C]) = sum_k b_(E,k) D_k(C),    b_(E,k) nonnegative integers.

Alternating-separator probes give an integer unitriangular recovery of all
D_k. Thus modulo q, equality in *every* such exterior is equivalent to
agreement of D_k for k<K_q, where

    K_q = min{k>=1 : q divides (k-1)! k!}.

The quotient is a congruence for uniform substitutions. Selective attachment
is outside its scope: Session295 gives two full-profile-equivalent cores
whose same marked selective addition has H=9 versus H=5. This is the
canonical hostile and the missing port coordinate.

The older loop mechanisms are [THM-3359, modular C-finite supports](../../../../01-canon/theorems/THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar.md)
and [THM-3364, cyclotomic Boolean clocks](../../../../01-canon/theorems/THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase.md).
They preserve eventual periodic index supports, with Fourier phase retained;
they do not classify arbitrary value sets or reconstruct physical sources.
The present result gives an exact tail length, exact loop length and exact
state count for a concrete all-order tournament family inside the newer
contextual quotient.

| Field | Exact contract |
|---|---|
| Source | Actual transitive tournaments TT_n, n>=1, under adjoining a terminal vertex |
| Target | Their complete uniform-context response modulo q |
| Map | TT_n -> (D_1,...,D_(K_q-1)) modulo q |
| Preserved | H modulo q in every uniform exterior, including nested uniform substitutions |
| Lost | Integer counts, induced-subset ports and selective incidences |
| Sidecar | q and absolute n; selective operations need additional marked endpoint data |
| Decisive hostile | TT_1 and TT_7 agree in every uniform context modulo14, but gating gives integer H=3 and129 |

The working board is: full context observer; factorial channel cutoff;
prime-power nonunit decay; surviving unit clocks; CRT composition; exact
count lift. The first five close below; the sixth remains a separate proof
obligation.

## L1. The complete state of TT_n

Each nonempty subset of a transitive tournament has exactly one directed
spanning path. Hence its ordered k-path covers are the ordered set
partitions, so with S(n,k) the Stirling number of the second kind,

    D_k(n)=(k-1)! k! S(n,k)
          =(k-1)! sum_(j=0)^k (-1)^(k-j) binom(k,j) j^n.       (1)

Here n>=1 and 0^n=0. In particular D_1=1,

    D_2=2^n-2,
    D_3=2*3^n-6*2^n+6.                                     (2)

Stirling recursion gives the closed physical transition

    D_k(n+1)=k D_k(n)+k(k-1)D_(k-1)(n),                     (3)

with the second term zero at k=1. Initialize at n=1 by (1,0,...,0).
The truncation at K_q is legitimate for every n, since (k-1)!k!
divides D_k(n). This is a finite deterministic state map, but its exact
orbit is considerably smaller than its formal state space.

## L2. Complete orbit classification for every q>=2

Factor q=product_p p^a_p. Define

    tau(q)=max_p a_p,
    ell(p^a)=p^(a-1)(p-1)                  for odd p,
    ell(2^a)=1                            for a=1,2,
    ell(2^a)=2^max(1,a-3)                 for a>=3,
    L(q)=lcm_p ell(p^a_p).                                  (4)

Then the first periodic state occurs at **n=tau(q)** and the least
period is **L(q)**. There are exactly

    tau(q)-1+L(q)                                           (5)

different states. For arbitrary m,n>=1, all uniform exterior counts of
TT_m and TT_n agree modulo q if and only if

    m=n, OR [m,n>=tau(q) AND m=n modulo L(q)].               (6)

Thus this metagraph orbit is a path of tau(q)-1 states feeding a cycle
of L(q) states. Both the transient and the cycle are exact, not upper
bounds from finite-state pigeonhole arguments.

### Odd prime-power proof

For q=p^a with p odd, every nonunit base j in (1) has j^n=0 modulo p^a
once n>=a. Every unit base has period dividing phi(p^a)=p^(a-1)(p-1).
This gives the asserted upper bounds. When a=1, all surviving k satisfy
k<p, so every nonzero j in (1) is already a unit, and the orbit is periodic
from n=1.

To show sharpness, the triangular finite-difference transform (1) recovers
j^n modulo p^a for j<=p-1, since (k-1)! is a unit for k<=p-1. If a>=2,
K_(p^a)>p, so it also recovers p^n using k=p, where (p-1)! is still a unit.
Consequently a periodic state at n forces

    p^n(p^L-1)=0 modulo p^a,

and therefore n>=a. This proves the exact transient for a>=2.

The recovered unit clocks with j=1,...,p-1 force the loop length to be
multiple of p-1: their reductions modulo p run through the nonzero field,
whose exponent is p-1. For completeness, if their exponent were e<p-1,
all p-1 nonzero elements would be roots of X^e-1, impossible by the field
root bound. The exponent divides p-1 by the elementary finite-group law.
Furthermore j=p-1 has order 2p^(a-1) modulo p^a. Indeed
v_p((p-1)^2-1)=1, and binomial induction gives
v_p(((p-1)^2)^(p^r)-1)=r+1. Combining this p-power factor with p-1 gives
the full phi(p^a). No assumption about a particular small primitive-root
lift is used.

### Two-power proof, including small boundaries

For a=1, K_2=2 and only D_1=1 survives: period1 from n=1.
For a=2, K_4=3 and D_2=2^n-2: period1 from n=2.
For every a>=2, D_2 forces the first periodic state to have n>=a,
because 2^n(2^L-1)=0 modulo2^a is equivalent to n>=a.

Suppose now a>=3. Then k=3 survives and, once n>=a, equation (2) reads

    D_3=2*3^n+6 modulo2^a.

Thus a period must be a multiple of ord_(2^(a-1))(3), namely
2^max(1,a-3). The boundary values are ord_4(3)=ord_8(3)=2;
for b>=3, binomial induction gives v_2(3^(2^r)-1)=r+2 for r>=1,
and hence ord_(2^b)(3)=2^(b-2).

For the upper bound, even j vanish once n>=a. For every k>=3 the prefactor
(k-1)! in (1) is even, so the odd-base powers need only be periodic modulo
2^(a-1). The exponent of its unit group is exactly the same
2^max(1,a-3). One can check this upper bound directly by repeated squaring:
every odd j has j^2=1 modulo8, and squaring raises the guaranteed 2-adic
valuation by at least one. The remaining k=1,2 coordinates are constant.
The single factor2 in the observed odd-base channels explains why using
phi(2^a), or even the full Carmichael exponent at 2^a, overestimates the
true observer period.

### CRT and the exact graph

Reducing the q-profile to p^a gives precisely its prime-power profile;
coordinates k>=K_(p^a) are zero there. Thus the first simultaneous periodic
index is the maximum of the prime-power indices and the least simultaneous
period is their lcm. If any state repeated before tau(q), deterministic
transition (3) would make the sequence periodic from that earlier index,
contradicting sharpness. Inside the cycle, any earlier repetition would
contradict the least period. This proves (5) and (6).

## L3. Fourteen, seventy-seven and the first25 primes

| q | Surviving run coordinates | First periodic n | Least period | Different states |
|---:|---:|---:|---:|---:|
|2|1|1|1|1|
|4|2|2|1|2|
|8|3|3|2|4|
|14|6|1|6|6|
|16|3|4|2|5|
|49|7|2|42|43|
|77|10|1|30|30|
|121|11|2|110|111|
|196=14^2|7|2|42|43|
|5929=77^2|11|2|2310|2311|

For every a>=1 the two requested bases belong to explicit infinite types:

    tau(14^a)=tau(77^a)=a,
    L(14^a)=3*2^max(1,a-3)*7^(a-1),
    L(77^a)=30*77^(a-1).                                  (7)

Thus increasing prime-power precision exposes longer transients and longer
loops. It does not reveal a single bounded loop underlying all exact counts.

Let P be the product of the first25 primes (through97). This squarefree
observer has tau(P)=1, K_P=97, and

    L(P)=lcm_(p<=97)(p-1)=39,419,059,680.                  (8)

The [99*2^k+23 first25-prime sieve](../../triples_ports_sieve_20261112/arithmetic/WHEEL.md)
has the much smaller exact Boolean period61,209,720, with

    39,419,059,680 = 644 * 61,209,720.

The two periods observe different data. The sieve retains only whether
one affine power-of-two value avoids those25 primes, with redundant
exclusions removed. The full context observer retains every surviving
power-base response from (1), which is sufficient for all uniform
exteriors. The 644 factor is an exact comparison, not an identification of
the underlying states, and sieve survival still does not imply primality.

## L4. Exact returns versus modular returns

Replacing TT_n by TT_(n+L(q)) for n>=tau(q) is an actual uniform module
substitution which preserves *all exterior counts modulo q*. Any nested
uniform construction transports the equivalence. This gives a repeatable
modular return with a complete state certificate, rather than a scalar
congruence tested on one exterior.

It does not give an integer target-preserving return. Already D_2 is
strictly increasing, and

    D_2(n+L)-D_2(n)=2^n(2^L-1)>0                         (9)

for every L>0. Gating TT_n gives the strong tournament Gamma(TT_n) with
H=1+2^n, so the same witness is physically visible:

    q14: n1 and n7 -> H3 and129;
    q77: n1 and n31 -> H3 and2,147,483,649.

The respective differences are divisible by14 and77, but are nonzero.
The physical full count is the lost magnitude coordinate. A bounded
observer becoming periodic cannot by itself control a growing defect or
prove a source-return construction at a fixed exact target.

At bounded order one can restore exactness by retaining a sufficient size
bound and using a modulus larger than the possible count difference. But
an all-order construction needs that bound or a separate positive identity;
letting the modulus grow without controlling construction size does not
supply it. The new theorem gives a safe rejection and caching layer for
uniform module searches, not a Chan–Pak coverage theorem.

## L5. An exact finite search box for infinite transitive inflations

Fix any quotient tournament R with r>=1 vertices and replace its vertices
by nonempty transitive modules TT_(n_1),...,TT_(n_r), with arbitrary positive
sizes. Put B(q)=tau(q)-1+L(q) and

    red_q(n)=n                              if n<tau(q),
    red_q(n)=tau(q)+(n-tau(q) mod L(q))       otherwise.

Then **every uniform-context response modulo q**, hence the immediate full
Hamiltonian count modulo q, of R[TT_(n_1),...,TT_(n_r)] is unchanged when
all n_i are replaced by red_q(n_i). Indeed replace the modules one at a
time; each complement is a uniform exterior, and the contextual equivalence
is itself preserved under any later uniform exterior. Therefore the entire
unbounded-size residue spectrum of this inflation family is exactly its
spectrum on the finite box

    1<=n_i<=B(q) for i=1,...,r.                            (10)

This gives an all-order, source-realized residue search with at most B(q)^r
parameter tuples, without assuming every formal state is realizable.
For q14 the box has side6; for q77 side30; for77^2 side2311. The ambient
quotient R and its orientations are fixed. This is a finite residue-spectrum
reduction for an infinite family, not a finite enumeration of ordinary
isomorphism types or exact integer H-values. It can safely reject an
impossible target residue before a source-return search, while positive
residue membership still requires an exact-target lift.

The companion verifies this consequence independently on every C3 inflation
with module sizes1..4, using actual subset-DP Hamiltonian counts, modulo3
and4. Its reduced representatives lie in the respective two-element boxes.

### A concrete all-order three-core obstruction

Take R=C3, and write e for the number of module sizes greater than1.
The eight representatives with sizes in{1,2} have H=3,5,11,45 when the
number of twos is respectively0,1,2,3. Cyclic rotation covers the three
placements for each middle case. Direct permutation enumeration verifies
all eight tiny counts independently of the profile recurrence.

Modulo3, reduction maps odd sizes to1 and even sizes to2. Therefore, for
all positive a,b,c,

    H(C3[TT_a,TT_b,TT_c]) = 0 mod3  if a,b,c have the same parity,
                            2 mod3  otherwise.             (11)

Modulo4, reduction keeps1 and maps every larger size to2, giving

    H(C3[TT_a,TT_b,TT_c]) = 3-2e mod4.                     (12)

These are all-order consequences of the finite-box theorem, with exact
small representatives as a certificate. They do not come from observing
several large sizes. In particular, this entire triangle-inflation family
misses every exact target H=1 mod3, including7 but also the realizable
count13. It cannot by itself prove the Chan–Pak coverage statement.
Changing the quotient core or allowing nontransitive/selectively attached
modules is a necessary way out of this particular source family. No claim
of literature novelty is made for the congruences (11)-(12).

## L6. Recovering the exact-target size sidecar

The finite logarithmic target box is **inherited**, not new here:
[Session175 clone obstruction](../../session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md)
and [Session181 root-module multiports, §§3-4](../../session181_multiport_defect_20261104/older/PROOF.md)
prove it, with Session181 extending the size bound to arbitrary module
interiors. For a fixed strong quotient R and target H=N,

    |B_i| <= floor(log_2(N-1)) for every module B_i.         (13)

The mechanism is the directed triangle through every quotient vertex:
the corresponding induced C3[B_i,1,1] has H(B_i)+P_(B_i)(2)>=1+2^|B_i|,
and Hamilton-path counts increase under induced extension by vertex
insertion. For transitive modules this yields only O((log X)^r) possible
counts up to X for a fixed r-vertex core. Session181 already applies it to
multiple independent transitive-module clocks and proves that no fixed
finite bank can close both branches H->3H+2 and H->3H+4 at all depths.

The new use of L5 is to combine this **exact-size** box with **complete
modular context** classes. Enumerate the allowed module sizes within(13),
reduce each modulo its exact tail-plus-loop class for cheap context tests,
and evaluate the surviving physical candidates integrally. This is a
terminating exact-target procedure for a fixed core with certified modular
pruning. It supplies no escape from the inherited fixed-bank sparsity
barrier. One must increase core structure, allow richer internal modules,
or use selective incidences with their full sidecars.

The independent-base automaticity criterion in incoming
[Session297 CUTOFF, §3](../../session297_parity_cutoff_20261113/CUTOFF.md)
is a different two-clock claim: it needs recognizability of the *same
physical target spectrum* in two multiplicatively independent numeral
bases. Our size-index loops give neither such recognizability theorem.
Conflating module-size periodicity with a digit language for exact H would
remove precisely the needed physical lift.

## Verification and scope

Run `python -B verify.py` and `python -O -B verify.py` in this directory.
The programs use explicit runtime checks (no assert dependence), compare
recurrence orbits against independent finite-difference evaluations, and
recover the first repetition rather than injecting the proposed cutoff.
The exact finite universe is q=2..256, 2^a for a<=14, 3^a for a<=7,
5^a for a<=5, 7^a for a<=4, 11^a for a<=3, and14^2,14^3,77^2:
270 distinct moduli. Direct permutation counts verify gated TT_n for n1..6.
The first25 modulus itself is not enumerated; its period follows from the
proved prime-power/CRT formula and exact lcm arithmetic.

The proof classifies this modular context orbit for all q and all n.
The receipt verifies a stated finite universe. Neither statement asserts
that every formal context state is realizable, every tournament lies on
this orbit, or an unmarked quotient supports selective port moves.
