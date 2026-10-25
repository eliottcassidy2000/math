# Independent review of the integral adaptive lattice

2026-10-25. **PASS / PROOF AUDIT.** Reviewed the complete local
[README.md](README.md) and [check.py](check.py), together with the inherited
[adaptive moment proof](../../unprotected_ports_cancellation_20261025/counting/README.md)
and the exact clone-cut interpretation. No mathematical gap was found in
the claimed integer coordinate change. The all-order proof and the finite
computational universe are distinguished below.

## 1. Inheritance and review scope

The closest mechanism is the finite-horizon moment quotient; the hostile
is its genuine past-5-cycle prefix with a negative power of two in the
old dense transfer. The corrected near miss is to clear only odd
factorials: that preserves actual integral answers but does not make the
transition integral on arbitrary integer coordinate vectors. The
least-used sidecar is gcd(c,d)>=1, which bounds the cycle count using
the state already retained.

The review board is cycle count / retained defect mask / odd factorial /
full product lattice / exact clone inflation. The audit checks the
coordinate operation, not a claim of positivity, a smaller state graph,
or a complete description of realizable moment vectors.

## 2. Independent algebraic derivation

Write O(u)=u!/2^v2(u!), with u including excess. Every actual prefix with
m cycles satisfies

    0 <= D(d) <= (d-1)m,

so b=max ceil(D(d)/(d-1)) over retained d>=3 is an integer lower bound
for m. An empty maximum is zero. Consequently

    Y_r = O(u) 2^(-rb) W_r
        = sum_m x_m 2^(r(m-b)),     x_m=O(u)w_m,

is integral whenever the sparse normalized weights are integral.
This embeds every reachable bank into the stated full product lattice;
it does not assert that every lattice vector is reachable.

The unnormalized moment transition is independently

    W'_r = beta/(ell^k k!) * 2^(A-kD(ell)+rk) * W_(ell k+r).

Substituting W_a=2^(ab)Y_a/O(u), multiplying by O(u')2^(-rb'), and
writing k!=2^v O(k) gives exactly

    Y'_r = T * 2^J * Y_(ell k+r),
    T = beta O(u')/(O(u) ell^k O(k)),
    J = A-v+k(ell b-D(ell))+r(k+b-b').

The odd quotient is integral before multiplying by beta. Indeed
u'!/(u! ell^k k!) is an integer: choose the ell*k cycle positions,
partition them into k unordered directed ell-cycles, and multiply by
(u'-ell*k)!/u!, an integer because excess is nonnegative. Taking odd
parts proves the asserted exact division.

All three contributions to J are nonnegative. First
v2(k!)<=binom(k,2)<=A. Next, for an affordable ell>=3,
D(ell)<=(ell-1)b; for ell=1 both b and D(1) vanish. Finally the new
mask is a subset of the old mask and

    D'(d) <= (d-1)(b+k),

so b'<=b+k, even when coordinates disappear. Thus each displayed
transition is integral on the entire source product lattice. A terminal
channel has r=0; it is covered by the same calculation.

The source channel is valid: ell*k+r<=S-ell*e<=S. When S>0 it cannot
be zero in an allowed completion. The k=0 convention forbids e>0;
therefore no unrepresented zero moment is read. At S=0 the separate
Y_0 coordinate is correctly retained.

## 3. Conversion, merging, and exact scope

Sparse-to-dense conversion uses only nonnegative shifts because m>=b.
All histories sharing a retained state use the same b and O(u), so
conversion commutes with both transfer and merging. The coordinate
change is invertible over Q, and therefore preserves the inherited
decorated continuation readouts, including signed clone excess.

The adaptive rank min(S,h) is the **nonterminal** rank; the terminal
readout is one coordinate. The support-monotonicity argument applies to
actual reachable banks: at fixed u,m every contribution has the same
sign (-1)^(u-m), hence cannot cancel. It need not apply to arbitrary
signed lattice vectors, and the proof does not require it there.

The exact clone interpretation remains the inherited one: weighted
quotients are automorphism orbits of cuts inside the canonical maximal
transitive modules of the inflated tournament. This normalization
preserves that signed inventory. It does not recover literal ports,
individual core identities, or automorphism actions from a moment bank.

## 4. Reproduced finite controls

Both commands passed, using explicit exceptions that survive -O:

    python 05-knowledge/results/port_residual_integral_20261025/counting/check.py
    python -O 05-knowledge/results/port_residual_integral_20261025/counting/check.py

The receipts cover ordinary and core targets 0..26, the independent
unmerged core oracle through 22, and 2,642 intermediate state audits.
The full order-44 core control has 719 state keys, 1,090 scalar cells,
and two genuine negative ungauged dense-channel exponents. There are
374 seeded relaxed-lattice transfer squares, containing 430 negative
ungauged channel exponents. Both the single-5-cycle hostile and the
genuine dense u=41 hostile pass. These are finite controls of the
all-order proof, not evidence of a new runtime bound.

## 5. Independent 105-cycle hostile and changing-mask control

A past 105-cycle, at stage ell=9 with N=114, has u=105,m=1,S=9 and

    (D(1),D(3),D(5),D(7),D(9)) = (0,0,0,0,6).

The correct ceiling bound is b=1. Replacing ceilings by floors would
give b=0. Adding one 9-cycle and terminating gives A=4, D(9)=6.
The floor variant has shift -2 and sends an arbitrary unit coordinate
to T/4 with T odd, so fails the full-lattice claim. The actual bound
gives shift 7 and the integral answer 128*T. This identifies the first
failed implication in the tempting weaker normalization.

There is also a coordinate-deletion control. Taking k=e=0 at ell=9
deletes D(9), dropping b from 1 to 0. The raw moments are unchanged,
but Y'_r=2^r Y_r. Hence a zero batch is not generally the identity in
normalized coordinates. The producer correctly claims identity only
when ell>S, where the retained affordable mask does not change.

The following standalone control was run successfully from the
repository root; its gates are explicit and survive optimized Python:

```python
from fractions import Fraction
from math import factorial, gcd
import importlib.util

p = '05-knowledge/results/port_residual_integral_20261025/counting/check.py'
spec = importlib.util.spec_from_file_location('candidate', p)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x

D = tuple(d-gcd(105, d) for d in range(1, 10, 2))
b = max((D[i]+2*i-1)//(2*i) for i in range(1, len(D)))
floor_bound = max(D[i]//(2*i) for i in range(1, len(D)))
need(D == (0, 0, 0, 0, 6) and b == 1 and floor_bound == 0,
     'ceiling versus floor hostile')
scale, remainder = divmod(oddpart(factorial(114)),
                         oddpart(factorial(105))*9)
need(remainder == 0 and scale % 2 == 1, 'odd quotient')
need(Fraction(scale, 4).denominator == 4, 'floor gauge must fail')
key, correct = c.dense_contribution(114, 9, 9, D, [0]+[1]*9, 1, 0)
need(key == (0, ()) and correct == [scale*128], 'terminal shift')
key, skip = c.dense_contribution(114, 9, 9, D, [0]+[1]*9, 0, 0)
need(key == (9, (0, 0, 0, 0)), 'deleted defect coordinate')
need(skip == [0]+[2**r for r in range(1, 10)], 'gauge-changing skip')
x = Fraction(oddpart(factorial(105))*2**52, 105)
need(x.denominator == 1, 'actual prefix sparse weight')
_, actual = c.dense_contribution(114, 9, 9, D, [0]+[int(x)]*9, 1, 0)
need(actual[0] == scale*128*x, 'actual prefix completion')
print('PASS: independent 105-cycle ceiling and mask-deletion controls')
```

**Disposition:** the producer's integer-transfer and product-lattice
claims pass independent review. No producer correction is required.
The ceiling hostile and changing-mask behavior are added review
evidence; they do not enlarge the theorem's claimed scope.
