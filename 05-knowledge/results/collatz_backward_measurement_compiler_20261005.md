# Working backward from a positive Collatz measurement

2026-10-05. **PROVED:** an explicit conversion from a source-specific odd-step
deadline to a positive weight floor and a signed-measurement margin.
The conversion is **CONDITIONAL** on the deadline being independently proved.
The mixed-family application below discharges that premise by a proper rank.
**FINITE-EXACT:** the declared arithmetic checks and failed candidate deadlines.
**OPEN:** a proved deadline or other positive measurement certificate for
every positive odd source. Unique certificate selection does not discharge
that universal existence obligation.

[Script](../../04-computation/experiments/collatz_backward_measurement_compiler_20261005.py)
and [output](collatz_backward_measurement_compiler_20261005.out).

## 1. The proposed reversal, with its missing premise exposed

The user's proposal is useful when an assumed positive measurement becomes a
finite proof obligation with a separately checkable source of truth. The
distinction is between using an assumption to **design** a certificate and
using it as an unproved premise of the claimed conclusion.

The closest mechanism is the
[source-floor deadline theorem, F1-F3](collatz_floor_transport_deadlines_20261005.md):
an independently proved `W(n)>=epsilon>0` forces a bounded actual ROOT word.
The earlier [localized kernel](collatz_localized_resolvent_floor_20261005.md)
turns independently certified measurements into such a floor. The present
reverse direction supplies an explicit floor from a stopping-time bound,
without needing that source's completed valuation word.

The incoming [refinement-floor note, T2](collatz_refinement_floor_shadows_20261005.md)
identifies the stopping-time obligation, and
[source-refinement SF1/SF3](collatz_source_refinement_floor_20261005.md)
transports existing floors across guarded edges or common futures. None starts
a floor at an arbitrary unknown source. The new step here is a quantitative,
source-only reverse compiler. Its finite tests do not upgrade its premise.

The portfolio is: **anchor**, close an independent measurement obligation;
**niche**, a localizing-matrix certificate; **wildcard**, a global Poisson
subsolution. The board is **deadline / hidden counter budget / signed
measurement / unique selected proof / well-founded dependency / residual**.
The least-used sidecar is the number of even valuations. The hostile is a
unique solution with a zero target atom, or a deterministic dependency cycle
with no justified seed. The corrected near miss is treating a unique proof
format as a proof that a witness exists.

The META-PATTERNS cards used are **"One item left" requires a typed residual**
and **Canonicalize mechanics; retain the predicate**. The residual below is
universal source coverage by independently valid premises, not a finite
missing numerical calculation.

## 2. A deadline bounds the hidden valuation budget

Let n>1 be odd. Write U(n)=oddpart(3n+1), and suppose an independent argument
proves that its first visit to1 occurs after tau<=T odd steps. This is the
premise; the compiler does not infer it from a positive numerical guess.

For the actual strict word a_1,...,a_tau, put

    A=sum_i a_i,       e=#{i:a_i is even},
    L=tau-1,          K=sum_i floor((a_i-1)/2),
    N=L+K.

Every preterminal odd state n_i is at least3. Telescoping the actual relations
`2^(a_i)n_(i+1)=3n_i+1` gives

    2^A = n product_(i=0)^(tau-1) (3+1/n_i)
        <= n(10/3)^tau <= n(10/3)^T.                  (1)

The last valuation is even and at least4, so e>=1. The exact parity identity is

    A=2K+tau+e,       N=(A+tau-e)/2-1.                (2)

All needed quantities can therefore be bounded with integers alone:

    Q=floor(n 10^T/3^T),
    a_max=bitlength(Q)-1,
    C=floor((a_max+T-1)/2)-1.

Under the premise, A<=a_max and N<=C. In particular C>=1. If the formula gives
C<1, the proposed nonroot deadline is impossible already by these necessary
constraints; n=3,T=1 is a small example. ROOT n=1 is handled separately.

This is sharper than bounding every intermediate integer independently and
summing the resulting logarithmic bounds: telescoping removes all internal
states. The preserved data are n and a proved T. The discarded word costs a
larger counter envelope, not a change of logical direction.

## 3. A finite counter envelope gives a positive source floor

For the inherited beta weight,

    W(n)=w(L,K)=2/[(N+2) binom(N+1,K)].

For fixed N>=1 and 1<=K<=N, the largest binomial denominator is the central
one. Also `(N+2)binom(N+1,floor((N+1)/2))` increases with N. Hence N<=C implies

    W(n)>=eta(n,T):=
       2/[(C+2) binom(C+1,floor((C+1)/2))]>0.         (3)

This is the exact minimum over the formal counter simplex N<=C, K>=1.
It need not be attained by an actual integer with the proposed n,T.
The formula is evaluated without reading a ROOT word or a certificate bank.
Its validity still depends on the proved deadline in section2.

For example, the conditional premise tau(27)<=41 yields C=56 and
`eta=1/435975364243345080`. The script uses the known route only as an
independent arithmetic control of this conversion; it does not claim an
independent new proof of that particular deadline.

## 4. Compile the floor into a positive measurement

If n is divisible by3, set z=n and T_z=T. Otherwise use the least
three-divisible predecessor z=rho(n). Its exact valuation table at
`n mod9=(1,2,4,5,7,8)` is `(6,5,4,1,2,3)`, so

    z=(2^a n-1)/3,       U(z)=n,       T_z=T+1.

This is a guarded, specified leaf, not an unspecified positive residue class.
Apply (3) at (z,T_z), obtaining eta>0 and index m=(z-3)/6.

For the source-local kernel h_m(j) of the earlier note, write

    H_(m,d)=sum_j lambda(6j+3) h_m(j)^d,
    A_(m,d)=9H_(m,d+1)-8H_(m,d).

Its unconditional error theorem says

    A_(m,d)>=lambda(z)-8(16/25)^d.

Choose the **least** nonnegative integer d with

    8(16/25)^d <= eta/2.

Then the independent deadline proves, without evaluating H numerically,

    boxed: A_(m,d)>=eta/2>0.                          (4)

The exact-rational loop terminates geometrically. Thus the backward design can
be reversed into a proof:

    independently proved deadline
      -> hidden counter bound
      -> source weight floor
      -> explicit positive signed-measurement certificate.

The canonical choice of rho, the integer budget, and least d make the compiled
plan unique **for the supplied pair (n,T)**. Different valid deadlines or
different proofs can still yield different plans. No uniqueness assertion
supplies the missing deadline.

## 5. Three ways to discharge a premise without a stored target word

### A. A proper induction domain, with an actual infinite positive result

The companion [inductive floor receipts](collatz_inductive_floor_receipts_20261005.md)
constructs a guarded two-step tree rooted at5. Each new unit source n has an
actual first-descent word `(1,a)` leading to a smaller parent y, with
`3<=a<=20`. Its first step rises; its two-step parent map decreases. The
membership decoder and the proof therefore use a well-founded ordinary-size
rank, rather than a stored ROOT entry for n.

The tree has two selected children at every vertex. Each child is at least
3/2 times its parent, so for b=bitlength(n) its depth is at most2b and
its odd ROOT time is at most4b+1. This is a proved family theorem. Equations
(1)-(4) therefore give a positive measurement certificate for the designated
leaf above every member of this infinite family.

The companion also proves a simpler counter-only floor
`lambda(rho(n))>=w(4b+1,18b+3)` and chooses d from it directly. Both routes
avoid using target ROOT words as production inputs. They are new explicit
measurement guarantees on a constructive rooted domain; they do not assert
a new general convergence basin or cover arbitrary integers.

The two constructions can be composed, rather than merely taking their
union. Permit both the source-refinement note's two one-step children and
the companion's two two-step children at **every** parent. The first valuation
selects the rule uniquely: at least2 for the former,1 for the latter.
The resulting four-child tree has a common strictly decreasing recognition
rank and a source-only deadline `tau(n)<=6bitlength(n)+1`.

It strictly enlarges the union of the two pure trees. For example,

    739 --(1,8)--> 13 --(3)--> 5 --(4)--> 1

uses both rules;739 is in neither pure tree. The production API
`hybrid_measurement_plan` discharges its deadline premise through that proved
membership criterion, without receiving a stored ROOT word or moment values.
At739 it uses the family deadline61, maps to the designated leaf15765 with
deadline62, and gives C=90 and

    eta=1/9448295337167360434640071920,
    m=2627,       d=151,
    A_(2627,151)>=1/18896590674334720869280143840>0.

This is an analytic positive measurement guarantee for the actual injection
measure. No H-value was evaluated. The separate rectangular-counter family
bound requires d=371 at this example; telescoping the valuation budget reduces
that sufficient degree to151. These are sufficient degrees, not claimed
optimal ones. Membership is itself a structural convergence certificate;
the construction does not claim to prove positivity without proving some
form of convergence.

The extension problem is precise: supply additional guarded parent rules,
prove a well-founded rank for every dependency, and prove that their domains
cover every source. More than one construction may discharge different
domains. A finite sample of covered sources does not prove that union is full.

### B. A finite measurement certificate found by exact optimization

The [moment-localizer companion](collatz_moment_localizer_feasibility_20261005.md)
tests a true source-sensitive moment packet against the zero-target model.
All off-target kernel values lie in[0,8/9]. Therefore, if the target mass is
zero, the matrix

    C_ij=(8/9)H_(i+j)-H_(i+j+1)

must be positive semidefinite. A rational polynomial P with P(1)=1 and
negative quadratic value gives a signed minorant

    q(h)=9(h-8/9)P(h)^2

and a positive atom floor. Exact interval errors must be included.
This is a finite algebraic way to **find** a stronger measurement witness;
it can outperform the single monomial selector on the same moment orders.
The moment packet must still be independently justified. A unique optimized
solution with value zero is an explicit hostile to inferring positivity
from uniqueness.

The same-order synthetic improvement has coefficient norm28577/81 instead
of17, so its gain in selectivity has an explicit precision cost. The companion
keeps that cost when verifying interval measurements. Its signed matrix has
at most one negative eigenvalue. Moreover the actual Collatz law has
infinitely many already-proved off-target atoms, making the finite normalized
optimization strictly convex and its minimizer unique. The missing step is
the **sign of that unique optimum**, not uniqueness of its definition.

### C. A globally proved Poisson inequality

The [Poisson source-dual companion](collatz_poisson_source_dual_20261005.md)
turns the measurement into a globally checkable functional inequality for the
adjoint of a specified contractive inverse operator. A bounded subsolution
with positive ROOT evaluation gives a source lower bound. This leaves a
concrete analytic or symbolically compressed proof obligation, with its norm
and tail control retained.

Finite locally supported positive duals force a literal connecting path;
they cannot manufacture support on a disconnected component. This scopes
certificate recycling precisely while leaving compressed global inequalities
available as a research target. A disconnected contracted self-loop has a
unique zero solution, so uniqueness of the Green solution is not positivity.

## 6. What assuming a unique certificate does and does not buy

For deterministic U, a rooted source has one strict first-hit valuation word.
One may also choose the least sufficient kernel degree or the first rational
dual in a fixed enumeration. These conventions remove redundant descriptions.
They do not change the logical distinction

    at most one canonical witness
    versus
    at least one valid witness for every source.

Backward reasoning is useful for synthesizing candidate deadlines, parent
rules, and inequalities. Each must be checked by an independent proof rule
before its consequence is promoted from CONDITIONAL to PROVED at a source.
A cycle of assumptions such as "source A is paid if B is paid, and B if A"
has no seed and no decreasing rank. A unique cyclic dependency is still unpaid.

The script checks three tempting universal deadline guesses. In the declared
odd-source universe below1024, `tau(n)<=bitlength(n)` first fails at7;
`tau(n)<=2bitlength(n)` and `tau(n)<=4bitlength(n)+1` first fail at27.
The latter formula is valid on the companion's guarded tree, demonstrating
why the domain guard cannot be dropped from a successful family proof.

This is the current recommended division of labor: generate candidate
well-founded domains, use exact moment optimization to discover useful
inequalities, and attempt global Poisson inequalities on residual domains.
For every route retain its source label, proof premises, finite checker,
and the exact set still uncovered. No unproved premise is erased by choosing
a unique certificate format.

## 7. Reproduction and finite scope

```text
python -B 04-computation/experiments/collatz_backward_measurement_compiler_20261005.py
python -B -O 04-computation/experiments/collatz_backward_measurement_compiler_20261005.py
```

There are29,343 exact checks: all511 nonroot positive odd sources below1024
with independently completed control routes, three deadline margins,
formal counter simplices through C=50, five compiled measurement controls,
three failed unguarded deadline guesses, malformed/inconsistent inputs,
and the structurally grounded mixed-family measurement plan.
Checks remain active under optimization. The universal conditional formulas
are proved above; bounded success is not a proof of all-source coverage.
