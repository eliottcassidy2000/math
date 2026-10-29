# Pruning partial cycle tuples with a future multiplicity bound

**PROVED elementary extension / FINITE-EXACT benchmarks.** This extends the
incoming `f8846701cd` precision theorem, not the full-count formula. The
baseline is pinned as `baseline_precision.py` in `source_pins.json`.

## Incoming mechanism and the extra coordinate

The concurrent [precision packet](../metagraph_scheduled_recursion_20261029/PRECISION.md)
proves that only O(sqrt(b)+1) odd permutation cycles can affect b normalized
low bits of F_n(r,s), where r>=1 tournament layers and s>=0 binary attributes
share the vertex relabeling. It also gives a nonnegative gcd-collision
budget, but its implementation checks completed tuples. We retain its
normalization, odd/even boundaries, allowed cycle counts and odd denominators.

For a nondecreasing odd partition lambda of length k define

    Delta(lambda)=sum_(i<j) (gcd(lambda_i,lambda_j)−1),
    V(lambda)=sum_a v_2(m_a!).

For n>=3 odd the normalized summand valuation is

    T+r Delta−V,  T=r(k−1)^2/2+s(k−1).

For n>=4 even it is the same expression with
T=r k(k−2)/2+s(k−2). Cases n=0,1,2 remain explicit boundary cases.

Consider an unfinished nondecreasing tuple alpha. Its last value occurs
m times, t positions remain, and all earlier value groups are closed.
Let F be the sum of v_2(m_a!) for those closed groups. Define

    U(alpha,t)=F+v_2((m+t)!).                            (1)

At the empty prefix put U=v_2(k!). Every completion satisfies

    V(lambda)<=U(alpha,t).                              (2)

Indeed, if x of the remaining positions repeat the last value and new
groups have sizes b_1,...,b_j with x+sum b_i=t, the multinomial coefficient
(m+t)! / ((m+x)! product b_i!) is an integer. Taking its 2-adic valuation
gives (2). Earlier groups cannot change because the tuple is nondecreasing.
The bound is sharp when all remaining positions may repeat the last value;
fixed-sum restrictions can make it loose, but cannot invalidate it.

Since all further collision terms are nonnegative, every completion has
normalized valuation at least

    T+r Delta(alpha)−U(alpha,t).                         (3)

**Pruning theorem.** If (3) is at least b, discard the entire unfinished
branch when computing the normalized count modulo2^b. Equality is safe.
Each omitted term is separately zero modulo2^b in the odd-denominator
local ring; cancellation is not assumed. This retains all genuinely
contributing summands and preserves the exact requested residue.

The bound never weakens the coarser factorial ceiling v_2(k!), since a
second multinomial comparison gives U<=v_2(k!). It can strengthen it
substantially once a prefix has several closed multiplicity groups.
At a completed tuple it is exact, so every completed tuple emitted by
the implementation genuinely contributes at the requested precision.

The extra information is tiny but necessary: the closed-group valuation
F, the size m of the current repeated group, and t remaining slots.
Using only current multiplicities as a ceiling would be invalid. Already
(1,1,1) has valuation1 in the ordinary odd normalization: at precision2,
forgetting the future factorial refund would discard its nonzero term.

## A second bound from the remaining sum

Let a be the last chosen odd part and R the remaining sum. Put
E=(R−t a)/2, a nonnegative integer. If E=0 every remaining part must equal
a, so (1) is exact. If E>0, let v be the number of new parts strictly larger
than a. Then 1<=v<=min(t,E): every such odd part consumes at least two
units of excess. Their multiplicity groups together have factorial
valuation at most v_2(v!). A sharper safe ceiling is therefore

    U_sum = F + max_(1<=v<=min(t,E))
                    [v_2((m+t−v)!)+v_2(v!)].             (4)

This is at most U. Apply (3) with U_sum in place of U. Prefix maxima
precomputed for each m+t make (4) a constant-time lookup per state.
It remains an upper bound: the maximizing choice may require all new
parts to be equal when the residual sum does not allow that.

`check_refund.py` independently enumerates every odd partition through30
and every prefix, using literal factorials to check
V(actual completion)<=U_sum<=U<=v_2(k!). This tests the bounds separately
from modular outputs, where accidentally omitted terms could cancel.

## Implementation and matched tests

`pruned_precision.py` retains the original denominator and compares its
completed valuation with an independent incremental calculation. Four
modes separate the inherited completed-tuple evaluator, the global factorial
ceiling, the completion ceiling (1), and the remaining-sum refinement (4).
All actual odd cycle lengths and multiplicities remain allowed.

The checker compares four channels (r,s)=(1,0),(1,1),(2,1),(3,0), every
n=0..30 and precisions1,2,4,8,12, against full Fraction sums over all odd
partitions. No precision truncation is used in that independent reference.
An additional repeated-unit hostile tests the refund boundary.

Three-run alternating-order benchmark medians on this machine:

| n, bits | inherited | global ceiling | completion ceiling | remaining sum |
|---|---:|---:|---:|---:|
| 201,8 | 0.15057 s | 0.08077 s | 0.03501 s | 0.02985 s |
| 401,8 | 2.16933 s | 1.26321 s | 0.44525 s | 0.34475 s |
| 101,16 | 0.07175 s | 0.01798 s | 0.01261 s | 0.01325 s |

At n=401,b=8 the inherited evaluator inspects622235 completed tuples;
the new evaluator emits4983, exactly the number of contributing tuples.
The completion ceiling visits565149 prefix nodes and rejects336925 prefix
branches; the sum refinement lowers those to363825 and334876. Completed
tuples and prefix nodes are different measures and are not equated.
The measured runtime improvement is about6.3-fold for this request, not
a uniform or full-count speedup. All raw timings and exact statistics
are saved in `precision_checks.json`. At101,16 the extra sum lookup saves
nodes but is slightly slower than the simpler completion ceiling; no
uniform timing dominance is claimed.

The matched exact outputs are

    a_201 / 2^100 = 67    mod 256,
    a_401 / 2^200 = 103   mod 256,
    a_101 / 2^50  = 53763 mod 65536.

The101 and201 residues also agree with complete counts from the independent
full-prefix engine (`precision_full_controls.json`). These are requested
**partial outputs**. In particular this is not a
complete A000568 count at401 or a 401-vertex class enumeration. The
incoming cycle-cap asymptotic bound remains unchanged; the present result
adds a proved early-exclusion rule and its implementation.

## Connection back to reset debt and metagraph ports

The successful move is to bound every possible future negative adjustment
before discarding a prefix. Here future gcd energy is nonnegative, and
the largest future factorial refund has the explicit ceiling (1).
For Collatz reset debt, the analogous missing ingredient is a bound on
the actual later payment with the original source and clocks retained.
The existence of (1) in the Burnside sum does not supply that arithmetic
bound. Likewise these global sum valuations do not bound entries of an
individual native port-response operator.

The next useful extension is an exact feasible-completion refund: each
pattern of new multiplicities imposes a numerical-semigroup constraint on
the remaining excess, with generators given by suffix multiplicity sums.
This is not implemented or benchmarked here. Every refinement must still
range over all legal completions; substituting one cheap example would
repeat the source-versus-average error.
