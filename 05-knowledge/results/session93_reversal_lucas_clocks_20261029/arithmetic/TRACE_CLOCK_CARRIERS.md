# A trace is a quotient; two clocks can repair its lost orientation

**Session93, 2026-10-29. Status: elementary proofs and finite-exact tests.**
No claim of external novelty for trace recurrences, Smith normal forms,
finite-field norm groups, or Fibonacci strong divisibility. The assembled
clock-menu loss theorem and hostile controls are the contribution of this note.
No prime infinitude, density, primitive-root conjecture, or tournament
classification follows from the recurrence alone.

## Inheritance and the decisive map

Source pin: `59cb22d0e93b6c364c4ec5d2954d9fe14ed468d8`.
Relevant source paths:

- `05-knowledge/results/golden_prime_clocks_20261003.md`: two-coordinate
  golden ring, norm/unit guards, split/inert/ramified clocks; explicitly
  distinguishes the ideal `(phi^5-1)` from `(11)`.
- `05-knowledge/results/session91_packing_golden_carriers_20261029/zeckendorf/NORM_CONTROLLED_OBSERVERS.md`:
  exact observer-kernel and Smith-form mechanism. Our trace calculation is
  a clock-dependent specialization, not a replacement theorem.
- `05-knowledge/results/lucas_mirror_gap_spaces_20261028/LUCAS_MIRRORS.md`:
  named Lucas words, reflection, and the eleven-point arithmetic reader.
- `01-canon/theorems/THM-4763-the-cat-map-on-the-logarithmic-side-lucas-traces-and-spectral-radius-phi-squared.md`:
  classical FGST decomposition, Lucas traces of the norm-one cat-map unit.
- `04-computation/fib_lucas_primes_s90bo.py`: historical prime-sieve exploration
  with a later warning that its probabilistic/geometric interpretations are
  heuristics. Its title is not a primality theorem.

Anchor: retained conjugate choice. Niche: decimation changes the observation
kernel. Wildcard: two individually lossy coprime clocks jointly recover the
state. The objects are ring states and named clock operations, not a cosmetic
orientation of a graph.

Let `O_m=(Z/m)[phi]/(phi²-phi-1)`, `a=phi²=1+phi`, and `x=u+v phi`.
Conjugation is `bar(x)=(u+v)-v phi`; trace is `2u+v`, norm `u²+uv-v²`.
The unit a has norm one. Multiplication and conjugation are retained exactly.

## 1. Exact full-history carrier for a sampled clock

For stride `d>=1`, put `b=a^d=r+s phi`, so `s=F_(2d)` and
`Tr(b)=L_(2d)`. The observation history is

`y_k(x)=Tr(b^k x) mod m`, for `k>=0`.

**Theorem.** Equality of the entire infinite histories is equivalent to
agreement of the first two observations. Every equivalence class has size
`g=gcd(m,5F_(2d))`, and the minimal number of exact history states is `m²/g`.
This statement concerns all `m²` ring states, not only the norm-one subset.

**Proof.** Since `N(b)=1`, `b²-Tr(b)b+1=0`, whence
`y_(k+2)=L_(2d)y_(k+1)-y_k`. The first two observations are the linear map

```
[ y_0 ]   [ 2       1       ] [u]
[ y_1 ] = [ 2r+s    r+3s    ] [v]  mod m.
```

Its determinant is `5s`; the entry 1 makes its integer Smith form
`diag(1,5s)`. More directly, a kernel point has `v=-2u`, and the second
condition is `5s*u=0 mod m`. There are exactly g such u. This proves the
kernel and image counts. The recurrence proves sufficiency for every future
observation, and reading times 0 and 1 proves necessity. The quotient is
closed under advancing the sampled clock:
`(y_0,y_1) -> (y_1,L_(2d)y_1-y_0)`.

At stride one the matrix is `[[2,1],[3,4]]`. For `gcd(m,5)=1`, reconstruction
is explicit:

`u=(4y_0-y_1)/5`, `v=(-3y_0+2y_1)/5 mod m`.

The denominator is guarded by an actual unit test. This is a pair of ordered
observations: an unordered pair is not the same carrier.

### The recovered state lies on a conic, not just a recurrence

The same two-coordinate map retains the norm through the exact identity

`y_0²+y_1²-L_(2d)y_0*y_1 = -5F_(2d)² N(x)`.

This follows by substituting the displayed linear forms and using
`L_(2d)²-5F_(2d)²=4`. At stride one, on norm-one states it is
`y_0²+y_1²-3y_0*y_1=-5`. When 5 is invertible modulo m, the explicit inverse
shows a bijection between norm-one ring states and points of this conic.
Clock advance `(t,z)->(z,3z-t)` preserves the form exactly. Thus the forward
recurrence, the chosen state and a quadratic conservation law are one exact
model, not merely analogous integer sequences. At m=5 the conic degenerates
to `(t+z)²=0`; this is precisely where the source recovery fails. It does not
mean that a real Pell orbit or every integral point has been enumerated.

## 2. Changing the clock can hide or restore a state

For finitely many positive strides `d_1,...,d_t`, include the common
observation `Tr x` and the observation `Tr(a^(d_i)x)` for every i.

**Theorem.** The joint kernel has cardinality

`gcd(m,5F_(2d_1),...,5F_(2d_t))`
`=gcd(m,5F_(2 gcd(d_1,...,d_t)))`.

Proof: substitute `v=-2u` as above. The simultaneous congruences have exactly
the stated gcd number of solutions. The final identity is Fibonacci strong
divisibility, obtained by the Euclidean algorithm from the addition formula
`F_(r+s)=F_(r-1)F_s+F_r F_(s+1)` and `gcd(F_r,F_(r-1))=1`.
Each selected clock's entire future adds nothing beyond its first pair.

Thus coprime strides recover everything recoverable by stride one. In
particular modulo 6, stride 2 loses a factor 3 and stride 3 loses a factor 2,
but together they recover all 36 ring states. This is genuinely free *joint*
information from correlated observations of the same source. It is not a
license to combine separately sampled states without retaining their identity.

**Hostile decimation.** Modulo 3, take `x=a=(1,1)` and
`x'=a^-1=(2,2)`. Both traces are zero. Since `a²=-1`, every stride-two
observation is zero for both. But the next stride-one traces are 1 and 2.
The underlying source and recurrence are unchanged; only the clock menu
changed. A closed quotient for one operation need not be closed for another.

## 3. The missing conjugate is not universally one bit

For an odd prime `p!=5`, the norm-one group U has order
`p-(5/p)`. In the inert case it is the norm-one subgroup of `F_(p²)^*`;
in the split case it is `{(t,t^-1):t in F_p^*}`.
On U, the fibre of a **single** trace is exactly `{x,x^-1}`. In the inert
case, elements with trace t and norm one are roots of `X²-tX+1`, a polynomial
of degree two over a field. In the split description it is the same
quadratic in t's coordinate. The branch collapses at x=x^-1, namely ±1.
This is a reflection quotient of a cyclic group; it is not yet a tournament
anti-automorphism or the Lucas-word reflection.

There are two independent caveats.

1. **CRT permits independent local signs.** For squarefree odd m coprime to 5,
   the single-trace fibre on norm-one states is the product of the local
   inversion fibres. It can have `2^omega(m)` points, not just two. Modulo 33,
   trace 3 has four norm-one preimages:
   `(1,1),(2,32),(13,10),(23,23)`.
   Global conjugation only pairs these; it does not connect the full fibre.
   Stride-one adjacent traces still recover all four states because 5 is a
   unit modulo 33.
2. **Ramification has a nilpotent phase, not just inversion.** Modulo 5 write
   `eta=phi-3`, so `eta²=0`, `bar(eta)=-eta`, and `a=-1+eta`.
   The ten norm-one states are `epsilon(1+c eta)`, epsilon=±1 and c in F5.
   For every k,
   `Tr(a^k epsilon(1+c eta))=2epsilon(-1)^k`.
   Each whole-history class therefore has five states. No number of later
   trace observations, and no alternative stride, can find c. A coefficient
   or equivalent source tag must be retained. This is a sharp five-valued
   hostile to an unconditional one-bit sidecar claim.

The full-history theorem holds also at m divisible by 2 and at prime powers;
the simple single-trace inversion-fibre assertion above was deliberately
restricted to odd unramified primes. No unsupported extension across
zero-divisors is being made.

## 4. Exact order, Frobenius, and the prime-distribution boundary

For odd p!=5:

- Split: `phi^(p-1)=1`, and `ord(a)` divides p-1.
- Inert: `phi^p=1-phi=-phi^-1`; hence `phi^(p+1)=-1`, while
  `a^(p+1)=1`. The norm-minus-one clock phi and norm-one clock a are distinct.
- At p=5, `ord(phi)=20` and `ord(a)=10`; reduction to the scalar residue
  discards the nilpotent factor 5.
- At p=2, O_2=F4 and both phi and a have order 3.

These are order **bounds**, not generator assertions. For example, at p=47
(inert), U has 48 elements but a has order 16. At p=29 (split), U has order
28 but a has order 7. At p=11, a has order 5; phi in O/(11) has order 10,
whereas the selected root phi=4 in O/(phi^5-1) has order 5. A root choice and
an ideal choice cannot be reconstructed merely from a sequence name.

Algorithm: given a proved multiple B of the order, factor B; for each prime
q dividing B, replace B by B/q as long as the actual ring power `a^(B/q)=1`.
The returned D is certified by `a^D=1` and `a^(D/q)!=1` for each prime q|D.
This uses O(sum valuation(B)) repeated-squaring tests after factorization;
no enumeration of the entire field is needed. It does not avoid factoring B.
The checker compares these certified orders with all structural bounds.
Nothing here proves that maximal orders occur infinitely often or with any
particular density. Likewise a Lucas or surd trace recurrence does not prove
infinitely many prime trace values.

## 5. What transfers to reversal carriers, and what does not

Loss ledger:

| Source | Quotient | Preserved | Lost | Minimal repair stated here |
|---|---|---|---|---|
| norm-one state over odd unramified prime | one trace | quadratic polynomial, inversion orbit | conjugate choice | at most one branch bit, or next named trace |
| ring state O_m | sampled trace history | all selected-clock scalar predictions | kernel of explicit matrix | coset tag of size gcd(m,5F_(2d)) |
| same state with several clocks | joint traces | correlations sharing source x | intersection kernel | gcd formula, often no tag |
| O_5 norm-one state | all trace clocks | sign | nilpotent phase | five-valued phase tag |

The mathematical bridge is an operation-relative observable algebra, already
present in the repo's carriers: the quotient must remember precisely what
allowed future operations distinguish. Here its kernel is computable by
Smith form, and coprime clocks force complementary hidden directions into
view. In tournament round trips, the retained centralizer and triangle errors
are noncommutative source data; a trace here is a commutative scalar quotient.
We do **not** identify the two kernels or infer centralizer generation from
this calculation. Their meaningful common template is: retain the marking
until a demonstrated joint observation separates its fibre.

A useful new-definition proposal is a **clock-resolving menu**: named
operations whose joint first-response map has trivial kernel on the intended
source class. This definition includes the explicit coprime-stride theorem
above and suggests checking determinant ideals for arithmetic carriers,
versus stabilizer intersections for tournament port carriers. It is a program,
not a theorem equating arithmetic and tournaments.

## Evidence and sources

Run `python3 check_trace_clock.py` and `python3 -O check_trace_clock.py`.
Both return **2,806,465 checks** with identical output. Universe: all ring
states for m=1..70 and strides 1..12; joint stride menus (2,3),(2,4),(3,6),
(4,6),(5,7) for m=1..50; complete norm-one groups for odd unramified primes
through 100; explicit ramification, decimation and CRT hostile controls.
The arbitrary-m theorems are proved above, not extrapolated from this census.

Primary author/institutional expositions checked 2026-10-29:
- Keith Conrad, [Finite Fields](https://kconrad.math.uconn.edu/blurbs/galoistheory/finitefields.pdf),
  cyclic multiplicative groups, Frobenius conjugacy, trace and norm. These
  are classical inputs; the concrete golden equations above are derived.
- NIST DLMF [18.9, recurrence relations](https://dlmf.nist.gov/18.9),
  classical three-term recurrence context. The trace recurrence here is
  independently proved by the quadratic identity, so no Chebyshev theorem
  beyond that identity is required.

The exact recent surd family referred to conversationally was not identifiable
by filename in this bounded source pass. No asserted surd-family instance
has therefore been silently substituted: phi² and the source-pinned cat-map
Lucas traces are the explicit concrete map used here.
