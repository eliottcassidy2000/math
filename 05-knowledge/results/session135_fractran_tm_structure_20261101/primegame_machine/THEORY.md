# A total basin decider for the fixed, nonhalting PRIMEGAME

**PROVED elementary, using the pinned proved S126 divisor-loop and S129
complete-cycle theorems; FINITE-EXACT independent controls.** This concerns
Conway's exact ordered fourteen-fraction program with `15/2,1/7,55` at the end.
It is not the `15/14` variant, a Collatz result, or a universality result.
No claim of literature priority is made.

## Result

Every positive integer belongs to exactly one of three decidable classes:

1. it eventually enters a least-period-two orbit;
2. it eventually enters a least-period-six orbit;
3. after removing its invariant spectator factor, its trajectory eventually
   **coalesces exactly with the standard orbit from seed 2**.

The third class tends to infinity, but not monotonically. The basin decision
procedure below takes polynomial time in the binary input length. All three
basins have polynomial-time computable natural densities (as reals, in precision
bits); no closed rational values for those
basin densities are asserted. S129's exact rational densities count periodic
points, not these eventual basins.

The original FRACTRAN program **never halts** on any positive input. It is the
new **decision procedure** that always halts with one of these three answers.

## 1. Inheritance and the new coordinate

Write

`N = S*2^a*3^b*5^c*7^d*11^e*13^f*17^g*19^h*23^i*29^j`,

where `S` is coprime to the ten displayed primes and never changes. Full
factorization is unnecessary: repeated division by these ten fixed primes
extracts all relevant coordinates and leaves `S` untouched.

The closest mechanisms are S126's ordered ten-register semantics and
four-counter/seven-control onehot machine, and S129's complete classification
of cycles (only least periods two and six), height grading, and periodic-point
density. The hostile is S129's three-control six-cycle at `4301`; discarding
extra control tokens is unsound. The corrected near miss is treating general
FRACTRAN universality as a property of this fixed rack. S133's semantics make
that language/program distinction explicit. The least-used sidecar here is a
nonnegative budget that decreases during *positive-control phases*.

Use the six-control-token count `k`, the inherited nondecreasing coordinate
`T`, and the new phase budget `Q`:

```
k = e+f+g+h+i+j,
T = a+c+i,
Q = b+d+g+h+i+j.
```

The exact rule increments are:

| Rule | Fraction | delta T | delta Q | delta k |
|---|---|---:|---:|---:|
| A | 17/91 | 0 | 0 | 0 |
| B | 78/85 | 0 | 0 | 0 |
| C | 19/51 | 0 | -1 | 0 |
| D | 23/38 | 0 | 0 | 0 |
| E | 29/33 | 0 | 0 | 0 |
| F | 77/29 | 0 | 0 | 0 |
| G | 95/23 | 0 | 0 | 0 |
| H | 77/19 | 0 | 0 | 0 |
| I | 1/17 | 0 | -1 | -1 |
| J | 11/13 | 0 | 0 | 0 |
| K | 13/11 | 0 | 0 | 0 |
| L | 15/2 | 0 | +1 | 0 |
| M | 1/7 | 0 | -1 | 0 |
| N | 55 | +1 | 0 | +1 |

These are exact linear substitutions, not observations extrapolated from runs.
**Q is not globally decreasing:** the seed `2` fires L and increases Q.

## 2. Finite positive-control phases

If `k>0`, at least one of the six lone-control denominator rules F,G,H,I,J,K
is enabled. Thus priority excludes L,M,N. Until the first arrival at `k=0`,
T is constant, Q is nonincreasing, and k is nonincreasing. In particular

```
a+c <= T_initial,  b+d <= Q_initial,  e+f+g+h+i+j <= k_initial.
```

The number of possible exponent states in this box is at most

```
B(T,Q,k) = (T+1) * binom(Q+2,2) * binom(k+6,6).
```

The last factor counts six-control vectors with bounded sum. Conditional on
such a vector, i is fixed, so the invariant equation a+c=T_initial-i gives at
most T_initial+1 possibilities for (a,c), or none when i>T_initial. The middle
factor counts pairs (b,d) with bounded sum. Overlap between budgets can only
make the true carrier smaller. This is a bound on full exponent states, not
on a lossy observer.

Consequently every positive-control phase either reaches `k=0` or repeats a
state. By S129, a repeat enters a two- or six-cycle. No trajectory can escape
to infinity while retaining positive control forever. The same argument,
stopping earlier at `k<=1`, bounds the multi-control prefix used by our decider.
The pigeonhole bound concerns actual FRACTRAN updates, not only abstract
decider work: a positive-control phase reaches k=0 or a cycle within B steps.

At a positive-control phase's first exit to `k=0`, the last rule is I. Just
before I there is only one 17 token; B and C must be disabled, so `c=b=0`.
The resulting state is exactly `S*2^T*7^d`. If T>0, the inherited explicit word

```
L^T M^d N (EF)^T K
```

carries it to `S*5^(T+1)*7^T*13 = S*C(T+1,T)`. S126 proves the seed-2 orbit
visits C(n,n-1) for **every candidate n>=2**, not only primes. Hence this is
literal coalescence with a spectator-scaled seed-2 trajectory. If T=0, the
same word instead leads to the `55 <-> 65` cycle, scaled by S.
Since d at phase exit is at most Q_initial, the original orbit reaches either
a periodic point or its stated standard-orbit target within
`B+3*T_initial+Q_initial+2` actual FRACTRAN updates. This bound makes no
assumption about prime gaps: the target is a candidate state, not the next
pure-power prime-output event.

**Extra closed domain.** If initially k>0 and T=0, the orbit is in a two-cycle
basin: a positive-control phase cannot diverge; S129's six-cycles have T>=1;
and an exit to k=0 leads to `55 <-> 65`.

## 3. The exact onehot escape bit

When `k<=1`, put `q=0` if no control token is present and otherwise let q be
the unique control prime. This domain is invariant under all fourteen rules.
The following Boolean function is an exact invariant of its dynamics:

| q | Escape bit P(a,b,c,d,q) |
|---|---|
| 0 | a+b > 0 |
| 11 or 13 | a+c > 0 **and** b+d > 0 |
| 17,19,29 | a+c > 0 |
| 23 | true |

**Guarded invariance proof.** In control 0, L replaces `(a,b)` by `(a-1,b+1)`;
M leaves them alone; and N is reached only when a=d=0, making the new
11-control condition `(c+1>0 and b>0)`. Thus all preserve P. E takes control 11
with b>0 to 29, where both predicates reduce to a+c>0; K at b=0 and J at d=0
preserve the displayed conjunction. A at d>0 takes13 to 17, preserving a+c>0.
B at c>0 preserves a+c and supplies b+1>0 in control 13. C requires c=0,b>0,
so its 17 and 19 conditions both reduce to a>0. I requires c=b=0, so its17
and 0 conditions both reduce to a>0. D requires a>0, making both 19 and 23 true.
H requires a=0 and supplies a positive d in control 11; both conditions reduce
to c>0. G supplies c+1 in control 19, making true remain true. F supplies a
positive d in control 11 and preserves a+c. All priority guards are retained.

**P=true implies standard-orbit coalescence.** For k=1, P excludes every
onehot JK cycle, because those cycles have b=d=0. A positive-control phase is
finite unless periodic, and S129 permits no other onehot cycle. It must
therefore reach I and then k=0. Invariance of P gives a=T>0 at that exit;
Section 2 supplies the exact standard-orbit reset. If initially q=0 and P is
true, the finite word

```
L^a M^d N (EF)^(a+b) K
```

gives control 13 with `(a',b',c',d')=(0,0,a+c+1,a+b)`. Its P is true, so the
previous argument applies. This also handles a=c=0,b>0, where the first trial
candidate is 1, not a prime-output event.

**P=false implies eventual period two.** If k=1 and T>0, the table forces
q=11/13 with b=d=0, already a JK cycle. If k=1 and T=0, Section 2's closed-domain
corollary applies. If q=0 and P=false, a=b=0, and `M^d N` directly enters the
cycle `S*5^(c+1)*11 <-> S*5^(c+1)*13`.

Thus P is not merely a predictor seen in a finite sample. It is a proved,
guarded Boolean invariant whose two values are the two possible onehot basins.
It does **not** preserve instruction words, timings, unbounded counters, or
prime-output positions; S129's no-finite-trace-observer theorem remains intact.

## 4. The all-input algorithm and complexity

`primegame_basin.py` supplies the actual ordered fraction list, ten-register
stepper, integer input stripping, exact period test, and the following decider:

1. Extract the ten valuations and invariant spectator S.
2. If k<=1, return `prime` or `cycle2` by P.
3. If F^6(x)=x, return `cycle2` if F^2(x)=x, otherwise `cycle6`.
4. Execute one actual first-enabled step and repeat from 2.

The period test uses S129's complete-cycle theorem, not an arbitrary six-step
cutoff. If neither terminal test has succeeded, no previously encountered state
can recur: a repeat would lie on a cycle, and every cycle point is detected by
F^6(x)=x. Before first entry to k<=1 the full state lies in the box of Section 2.
Therefore there are at most B nonterminal iterations. This proves totality;
the executable's explicit bound check is a consistency gate, not a timeout
interpreted as a mathematical outcome.

For input binary length ell, T,Q,k <= ell, since each is a sum over disjoint
subsets of the input's prime exponents. Thus B=O(ell^9). Each iteration makes
at most seven PRIMEGAME updates on ten O(log ell)-bit exponent registers.
Elementary fixed-prime valuation extraction takes polynomial time as well.
A conservative bit-operation bound is O(ell^2 + ell^9 log(ell+1)); no factoring
oracle, undisclosed search, sequence table, or unbounded run-until-it-looks-done
is used. This is an upper bound, not an optimized complexity claim.

The result is a **halting decision-machine bridge for a fixed nonhalting
FRACTRAN program**. Ordinary implementation of the finite algorithm on a Turing
machine proves decidability; we do not claim to deliver a newly minimized
Turing transition table, a universal simulation by PRIMEGAME, or a general
FRACTRAN halting decider.

For any diverging initial state with k>0, the guarantee in Section 2 is more
specific: it reaches the candidate state C(T_initial+1,T_initial), scaled by S.
The three classes are disjoint because the standard seed-2 orbit is unbounded
(S126's output theorem) and the other two have different least periods.
As a structural corollary, T is bounded along an orbit iff the orbit is
eventually periodic: cycles never use N, whereas the standard orbit executes
N once per successful candidate and T grows without bound. Thus the inherited
monotone coordinate T becomes an exact escape criterion only after the new
basin proof; it is not by itself an immediately decidable threshold test.

## 5. Computable natural densities of eventual basins

For any fixed valuation vector x, the set of positive integers with exactly
those ten valuations has natural density

```
w(x) = product_p (1-1/p)*p^(-x_p).
```

This follows by inclusion-exclusion for exact fixed-prime valuations. Let D be
any of our three basin classes of vectors. For caps h_p, the finite sum

```
L_D(h) = sum_{x in D, 0<=x_p<=h_p} w(x)
```

is exactly computable using the total decider. The omitted inputs have density

```
epsilon(h) = 1 - product_p (1-p^(-(h_p+1)))
           <= sum_p p^(-(h_p+1)).
```

The finite-box lower density and its lower-plus-epsilon upper bound sandwich
both liminf and limsup. Sending every h_p to infinity proves existence of
`delta_D=sum_{x in D}w(x)`. Choosing caps for any requested error gives an
algorithm computing this real to that accuracy. The three densities sum to 1.
All are positive: initial valuation cylinders for 11, 4301, and 2 give positive
measures in the respective classes.

In fact these are **polynomial-time computable reals** in requested precision
bits m. Take the uniform cap H=m+4. Then the error is at most
`10*2^(-m-5) < 2^(-m)`. There are `(H+1)^10` vectors, each with an O(m)-bit
encoded integer, so the preceding polynomial-time decider makes this a
polynomial-size outer enumeration. All weights can share the exact denominator
`D=product_p p^(H+1)`, of O(m) bits, and the numerator of w(x) is exactly
`product_p (p-1)*p^(H-x_p)`. Summation and error computation therefore also have
polynomial bit cost. A very loose resulting bound is O((m+1)^19 log(m+2));
this is a complexity classification, not a practical-speed claim.
`density_approximation.py --bits m` is the actual finite algorithm, while
`--cap 1` is the small normal/optimized smoke test. Even m=0 requests 5^10
states: the exponent 10 is material. Density existence needs only finite-prime
valuation dependence and the tail bound; **computability and this complexity
bound** additionally use the total basin decider.

There are simple exact lower bounds without a large valuation enumeration.
Put `C0=product_(q in {11,13,17,19,23,29})(1-1/q)`. The total onehot density is
`C0*(1+sum_q 1/q)`. Independence of the four data valuations and Section 3 give
its exact prime-basin contribution

```
C0 * [2/3 + (1/11+1/13)*9/35
            + (1/17+1/19+1/29)*3/5 + 1/23]
  = 50136286961664/86285158710179 > 0.58105.
```

Here `Pr(a+c>0)=3/5`, `Pr(b+d>0)=3/7`, and `Pr(a+b>0)=2/3`.
The remaining onehot mass is in two-cycle basins. The disjoint closed domain
`k>=2,T=0` contributes a further

```
(2/5)*[(22/23)-C0*(1+sum_(q!=23)1/q)]
  = 365260463044/18757643197865
```

to those two-cycle basins. Combining gives the rigorous two-cycle-basin lower
bound `753557212/2002524095 > 0.37630`. S129's period-six-point density is also
a six-cycle-basin lower bound. `density_bounds.py` emits these exact rationals
and the resulting complementary upper bounds; its decimals are explicitly for
display only. These are densities of **eventual basins**;
S129's exact periodic-point density remains unchanged. We assert neither that
these new densities are rational nor that their tight numerical evaluation is
cheap. The tiny finite integer census is not used as an asymptotic estimate.

## 6. Connection contract and reproducibility

* **Source:** the exact prioritized arithmetic rack on all positive integers.
* **Target:** a total three-output basin decision procedure.
* **Map:** extract ten valuations, retain spectator S, then use positive-phase
  budgets, periodic-point detection, and the guarded onehot escape bit.
* **Preserved predicate:** eventual least period2, eventual least period6, or
  eventual coalescence with the S-scaled standard prime-generating orbit.
* **Lost:** exact trajectory, first-arrival time, instruction word and outputs.
* **Restoration sidecar:** full exponent vector plus original priority and S.
* **Hostiles:** 4301 requires three controls and enters a genuine six-cycle;
  17 versus23 share zero data counters but have different basins; 2*11 has
  positive T but is trapped because b=d=0; Q increases on 2→15.

Run from this directory:

```sh
python3 -B verify.py > VERIFIED.normal.json
python3 -B -O verify.py > VERIFIED.optimized.json
cmp VERIFIED.normal.json VERIFIED.optimized.json
python3 -B density_bounds.py > DENSITY.normal.json
python3 -B -O density_bounds.py > DENSITY.optimized.json
cmp DENSITY.normal.json DENSITY.optimized.json
python3 -B density_approximation.py --cap 1 > APPROXIMATION.normal.json
python3 -B -O density_approximation.py --cap 1 > APPROXIMATION.optimized.json
cmp APPROXIMATION.normal.json APPROXIMATION.optimized.json
python3 -B verify_entry_bound.py > ENTRY.normal.json
python3 -B -O verify_entry_bound.py > ENTRY.optimized.json
cmp ENTRY.normal.json ENTRY.optimized.json
python3 -B primegame_basin.py 4301
```

The checker has explicit exception-based gates, active under optimization. It
compares the producer against separately printed literal integer fractions,
uses actual candidate states C(n,n-1) as coalescence witnesses, and detects raw
repeats independently of the producer's period-six test. It checks all 59,049
vectors in `{0,1,2}^10`, all 9,072 onehot states with data 0..5, integer inputs
1..30,000, and 300 seeded higher-exponent cases. It also verifies every budget
delta and pinned input hash. The separate raw-entry checker tests 2,760 inputs
(2,744 with positive control), compares each claimed reset letter to literal
first-eligible execution, and checks the exact candidate endpoint and bound. These are finite controls, not the all-state proof.

`SOURCE_PINS.json` identifies the exact input commit, paths, and SHA-256 bytes.
`inputs/` makes the inherited proof dependencies locally inspectable. The
parent's independent audit supplies a separate implementation and proof review.
