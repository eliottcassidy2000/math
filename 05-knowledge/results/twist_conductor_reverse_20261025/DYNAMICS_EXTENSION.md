# The golden tilt reverses the barrier drift, within its measure scope

2026-10-25. **PROVED conditional-law specialization / FINITE-EXACT.**
This extends the incoming Session64 harmonic calculation, not its open
marked-integer conclusion. It is the elementary Doob transform, explicitly
computed here; no novelty claim for conditioning a Markov chain is made.

The current stronger source is
[Session64, sections 1--7](../session64_lca_reverse_valuation_20261026/session64-reverse-valuation/REPORT.md),
with its [independent dynamics audit](../session64_lca_reverse_valuation_20261026/session64-audit-dynamics/REPORT.md).
It derives H(x)=3^(k-1)(3x-1)/2^k from Antihydra's parameter macros,
k=v2(3x-1), and the resource update a'=a+2k-3. It proves the punctured
2-adic coding, its exact LCA isometry, integer inverse uniqueness, and
the Haar ruin probability q^a, q=(sqrt(5)-1)/2. We use those results rather
than re-proving weaker coding statements. Its native-Turing-machine clock
is expressly outside scope; one odd return uses k parameter macros.

## Exact conditional valuation law

Under normalized odd 2-adic Haar measure, valuation letters have independent
probabilities 2^-k. Start a>=1 and stop on first arrival at zero. The source
proves that the probability of eventual ruin is h(a)=q^a. Until stopping,
conditioning on eventual ruin gives the transition probabilities

    P_*(k)=2^-k h(a+2k-3)/h(a)
          =2^-k q^(2k-3)
          =(1-r) r^(k-1),       r=q^2/2.              (D1)

They do not depend on the positive current counter a. The k=1 move can
land at zero when a=1; h(0)=1 includes that boundary. They sum to one
because q(2-q^2)=1. The conditional mean letter and counter drift are

    E_* k=1/(1-r)=2q,
    E_* (2k-3)=4q-3<0.                              (D2)

Thus the original positive mean drift +1 becomes negative under conditioning.
No negative or signed probability was introduced. The change is an exact
positive weighting in Q(sqrt(5)). For a surviving prefix w of m letters and
sum S, the likelihood ratio is

    P_*(w)/P(w)=q^(2S-3m)=h(a_m)/h(a_0).           (D3)

Every complete ruin path has the same endpoint ratio q^(-a_0). Therefore
the conditional probability of ruin by a finite horizon t is exactly
the original finite-horizon probability divided by q^a_0. Taking the
source's proved limit gives eventual ruin probability one under P_*.
This also proves the conditional description directly on cylinder events.

## What the tilt transports, and what it does not

Source: Haar-distributed odd 2-adic initial states with the marked initial
counter. Target: that probability measure conditioned on eventual barrier
hit. Map: retain each exact valuation word and multiply its prefix mass by
(D3). Preserved: word identity, all counter prefixes, odd-return and macro
clocks, and the barrier event. Changed: the probability weights and their
drift. Lost if only the counter is stored: the actual initial integer and
the full arithmetic endpoint. Keep the affine inverse branch as a sidecar.

The conditional law is not Haar on the original source space. It is not a
change of the deterministic orbit of x=27, nor a proof of its survival or
halting. A specific integer still has mass zero before and after conditioning.
Even positive density and a complete full-shift model cannot supply the
missing fixed-source lower bound.

There is a useful structural comparison with reverse tournament potentials.
There the factor 2^[V(w)-V(v)] moves valuation cost between endpoints of a
finite DAG edge while retaining exact weighted paths. Here the factor
h(a')/h(a) moves probability between paths and needs harmonicity to normalize
each outgoing row. In both cases the factor telescopes, but they preserve
different predicates. Neither produces a well-founded rank for an arbitrary
infinite deterministic orbit.

The conductor comparison is equally precise: a one-step compatibility test
must use the correct ambient operators. In I4 the missing operator is the
integral ratio w; here it is the full row-normalization identity, which fails
if a Haar expectation is replaced by one deterministic successor. One cannot
erase the operator/measure and keep only the scalar golden number.

## Hostiles and next question

The words (1,2) and (2,1) have the same final change and likelihood ratio,
but different least initial safe resource, two and one. Keep prefix minima.
Under the 2-adic conjugacy Collatz's integer fixed point 1 goes to 3/5,
not to a native Antihydra state. Keep the integer locus as well as the code.

**PROVED refinement obstruction.** Let C_m be the nested valuation cylinders
of one fixed source, with letter sum S_m. Under Haar measure,
mu(C_m)=2^(-S_m). Every fixed finite measure absolutely continuous with mu
assigns C_m mass tending to zero, even if its density is unbounded: continuity
from above reduces the limit to the singleton, which has mu-measure zero.
The same conclusion holds for a varying uniformly integrable family of
densities f_m. Indeed split f_m at height L; the bounded part contributes
at most L*mu(C_m), and uniform integrability makes the remaining integral
uniformly small. Thus a fixed positive mass on every refinement forces
failure of uniform integrability, not merely a large bounded constant.

Quantitatively, if p>1 and ||f_m||_p<=M, Holder gives

    integral_(C_m) f_m dmu <= M*2^[-S_m(1-1/p)].       (D4)

Conversely a lower bound c>0 requires ||f_m||_p >=
c*2^[S_m(1-1/p)]. The exact conditional density 1_(C_m)/mu(C_m)
attains the exponent, so this concentration cost is sharp. It stores a
finite source cylinder, not a terminal certificate. The conditional ruin
law (D1) is absolutely continuous and cannot avoid the obstruction.

**OPEN helper question:** can this explicit source-concentration cost be
charged against a genuine independently bounded transport quantity, while
the event being certified remains terminal entry? A candidate must retain
the integer source, exact cylinder and native clock. Alternatively use an
atomic source measure from the outset, but then re-prove the transition law:
Haar independence and (D1) cannot simply be inherited. These are concrete
routes around the obstruction, not a no-go for Collatz itself.

The checker computes both finite-horizon ruin recursions independently over
Q and Q(sqrt(5)), retaining the infinite-letter tail as a zero-ruin tail at
the stated finite horizon. No truncated alphabet is renormalized.
