# Independent audit: finite lattice intervals and coupled prime flow

**Verdict: PASS after a minor termination-schedule clarification.** The
finite-lattice interval theorem, both primitive-prime selections, explicit
N5 witness, coupled contradiction, and parity obstruction are correctly
scoped. They do not assert an odd-perfect-number support or a global bound.

## 1. Universal interval realization

The labels `d_x=product_(x not<=t) ell_t` form an injective join embedding
into divisibility. If x not<=y, the label prime ell_y distinguishes d_x
from d_y; y cannot be the omitted top in this case. Also d_0=1 and every
nonzero d_x contains ell_0, so every primitive-divisor invocation has an
odd exponent greater than one.

I checked the classical source directly: [Gyory and Smyth, Theorem 13,
page 325](https://math.colgate.edu/~integers/k27/k27.pdf). Its positive-base
exceptions concern exponent two with a special sum, or exponent six for
2^n-1. The additional negative-b exception in that statement cannot occur
here because both uses have b=1. Odd exponents greater than one exclude
all exceptional cases in both applications.

The first selected primitive prime p has ord_p(2)=M, hence p=1 mod M.
Each target r_x has ord_(r_x)(p)=d_x>1. Different x give different orders
and therefore different target primes. None of the targets is p, is even,
or equals a label prime ell_t: every ell_t divides p-1 and would have
order one. Consequently no target divides any d_y. This is precisely the
condition preventing unintended LTE valuation increments in the image.

Thus `Phi_p(d_y)_x=c_x*1[x<=y]` is exact. The ambient coordinatewise meet
is the lattice meet. Closure of the coordinatewise maximum is its join,
because its demand rank is the LCM of the active d_x. More importantly,
every realizable profile between zero and c has each coordinate either
zero or c_x; applying closure to its active labels makes it exactly the
profile of their join. This proves that the image is the whole ambient
interval, not merely an embedded sublattice.

The zero profile uses length one only as a compiler representative. A fresh
odd prime length greater than every label avoids every nontrivial d_x and
has positive even exponent. The one-element lattice with no targets is
also valid. No arithmetic size bound or global perfect-number realization
is inferred.

## 2. Explicit arithmetic and independent method

The independent checker imports no producer code. For N5 it verifies all
five primes, proves each claimed order using its prime-factor tests, and
checks that the order does not lift to the target square. It then enumerates
all sixteen possible capped profiles and evaluates the actual geometric
sum at each demand LCM. Exactly the five advertised N5 profiles are fixed.
Length three independently realizes the zero corner.

For the coupled and parity examples, the independent path does not use the
producer's LTE/rank compiler. It computes complete modular periods at the
valuation upper-bound guard prime powers. The denominators p-1 are units
at all these targets, so the geometric sum modulo each guard is recovered
directly from p^n. Testing a complete period together with ordinary/Euler
length parity gives the entire capped row domain, not a sample of lengths.

The coupled exact domains are

    source 3: {(0,0),(2,0),(0,1)}, period 2530;
    source 5: {(0,0),(1,0),(0,1)}, period 5060.

Each period has 1,265 role-allowed classes. Each column total two is separately
attainable, but no pair of rows gives totals (2,2). The forced source-3 demand
(1,1) has least length 55 and actual profile (3,1), exceeding its first cap.
Both physical exponent parities in the displayed separate witnesses are
correct: ordinary n is odd and Euler n is two modulo four.

All four parity sources have complete capped domain {0,2}; their guard
period is 78, with 39 odd length classes. Their possible sumset is exactly
{0,2,4,6,8}. The requested total five is excluded although every interval
endpoint is attainable and the residual column inequalities leave [0,2]
unchanged. This is a genuine modular obstruction after interval propagation,
not merely an untested failure of the local closure routine.

## 3. Propagation and the corrected scheduling statement

Each lower closure is the least feasible row completion, and the column
residual updates follow from the exact column equality. In Step 5, a threshold
k is feasible exactly when its least row completion remains inside the box.
The feasible thresholds are downward closed. Their greatest member is the
actual greatest attainable coordinate: if its completion had a larger
coordinate, that larger threshold would also be feasible. Thus tightening
to it loses no allocation.

Each strict update decreases the integer sum of box widths. The original
phrase "repeatedly apply any of" did not literally rule out an infinite
schedule of no-op updates. The producer clarified the theorem to repeated
full sweeps stopping at an unchanged sweep, and separately noted that an
arbitrary schedule has only finitely many strict updates. The supplied code
already implements full sweeps, so no numerical receipt changed.

The checker deliberately excludes a source from its own target tuple.
The theorem's diagonal source=target entries are permanently zero and can
be omitted before using that routine; the tested examples have no such
diagonals. No unsupported diagonal arithmetic was used in a certificate.

## 4. Replay and scope

I read the full proof and checker, and freshly ran the producer normally
and under Python optimization. Both outputs equal the supplied normal.json
byte-for-byte. The independent modular-period checker also passes identically
in both modes. Source and proof hashes are in its receipt.

The finite lattice census tests symbolic encodings of naturally ordered
finite lattice relations, not actual enormous arithmetic witnesses for
each relation. The general proof supplies their existence. The random
propagation tests corroborate soundness but are not an exhaustive claim
about all branches. Feasible partial allocations do not certify global
exponent compatibility or exclusion of outside factors.

No further mathematical correction or blocker was found.
