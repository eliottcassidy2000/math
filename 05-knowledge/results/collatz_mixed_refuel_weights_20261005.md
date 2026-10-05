# Mixed refuel weights have finite certificate approximants with uniform tails

2026-10-05. **PROVED:** the beta mixtures, source-aware weight formulas,
incoming-flow identities, polynomial comparisons, finite backward compiler,
and unconditional numerical tail bounds. **FINITE-EXACT:** the stated
implementation controls. **OPEN:** positivity at every positive odd input,
equivalently universal Collatz. A small numerical tail of these weights is
not a bound on the atomic mass of still-unproved sources.

Artifacts: [script](../../04-computation/experiments/collatz_mixed_refuel_weights_20261005.py)
and [output](collatz_mixed_refuel_weights_20261005.out).

## 1. Inherit the rooted flow, not a convergence assumption

Use U(n)=oddpart(3n+1), S(n)=4n+1, and the primitive bases
B={b positive odd:v2(3b+1) in{1,2}}. Every positive odd source uniquely
decomposes as n=S^j(b). A nonroot base has

    U(b)=S^k(G(b)),  G(b) in B, k>=0.

These decoders terminate by arithmetic division; following G indefinitely
is a different operation. The exact definitions and finite ROOT-kernel
verifier are inherited from
[three-bit sibling flow](collatz_three_bit_sibling_flow_20261005.md), §§5–7.
The closest current existence result is incoming commit69f6b3903,
[critical flow](collatz_three_bits_critical_flow_20261005.md), C8: for each
fixed computable r in(0,1), the rooted chain weight is already an
unconditional computable nonnegative summable real-valued function, even
at rho=1. Its support is exactly the rooted component, and it is maximal
under its specified constraints. We do not claim that existence result anew.

Here take rational0<r,rho<1. If b reaches1 under G in T steps, and K is
the sum of the T sibling depths k, put

    f_(r,rho)(S^j b)=rho^T(1-r)^T r^(K+j).             (1)

If its G-chain never reaches1, define the value to be0. ROOT has T=K=j=0
and weight1. The full known ROOT sibling ray has T=K=0, arbitrary j.
Throughout write L=K+j. T is a **base-chain length**, not the number of
actual odd U-steps. This distinction matters even on the ROOT sibling ray.

The inherited killed incoming operator K_U excludes the ROOT target.
Its inequality is K_U f<=rho f. Equality holds at nonroot targets with
nonempty inverse fibre, namely3∤y; at multiples of3 the incoming sum is0.
The target3 is already a hostile to unqualified off-root equality.

Our board is **verified backward addresses / original source / refuel depth /
latent parameters / summable weight tails / support positivity**. The niche
is a mixture over parameter choices; the wildcard is the finite coefficient
alphabet after dividing by the earlier source prior. The corrected near miss
is reading numerical convergence of weights as universal source coverage.

## 2. Integrate the refuel parameter exactly

Give r the probability density2(1-r) on(0,1). On a rooted source with
parameters(T,L), the resulting weight is

    F_rho(n)=2rho^T B(L+1,T+2)
            =2rho^T L!(T+1)!/(L+T+2)!.                (2)

It is0 on all other sources. This follows by integrating
2rho^T r^L(1-r)^(T+1). All displayed values on a supplied finite ROOT
certificate are exact rationals for rational rho.

In particular

    F_rho(S^j(1))=2/((j+1)(j+2)),
    sum_(j>=J) F_rho(S^j(1))=2/(J+1).                 (3)

The entire rooted sibling ray has total weight2. The first source is ROOT,
with the empty strict first-hit word; j>0 has the one-letter word(2j+2).
Using a uniform density in r instead would give1/(j+1), a nonsummable
harmonic root ray. Thus the boundary weight near r=1 is a real obligation.

Mixing preserves K_U F_rho<=rho F_rho by nonnegativity. More explicitly,
for a nonroot3-unit target y, its inverse base b0(y) is unique and every
predecessor is S^j(b0(y)). If that base has parameters(T,K), then

    sum_(j>=0) F_rho(S^j b0)=2rho^T B(K+1,T+1)
                            =rho F_rho(y).            (4)

The last equality uses the actual base edge, including its k. At multiples
of3 there is no inverse base and the sum is0. At ROOT the nonroot incoming
sum is1; imposing the discounted inequality there would be false. The
operator must remain killed at ROOT.

There is a pathwise comparison with every fixed r. Put H=T+L. Since
binom(H,T) r^L(1-r)^T<=1,

    F_rho(n) >= [2(T+1)/((H+1)(H+2))] f_(r,rho)(n).   (5)

The mixture pays only this polynomial comparison factor relative to the
best choice of a fixed refuel parameter on a specified certificate. It does
not prove that every input has such a certificate.

## 3. An optional second mixture removes the fixed discount choice

Also give rho density2(1-rho). The double mixture is

    D(n)=4 L!T!/((T+2)(L+T+2)!)                      (6)

on the rooted component, and0 elsewhere. For each nonroot3-unit target
with base length T=T(y), its incoming ratio is exactly

    K_U D(y)/D(y)=(T+1)/(T+3)<1,                     (7)

when D(y)>0. This is the posterior mean of rho under density proportional
to rho^T(1-rho). Rooted multiples of3 have ratio0. On unrooted targets
both sides vanish, so no positive-ratio claim is made there.

For every fixed r,rho in(0,1), the same binomial comparison gives

    D(n) >= [4/((T+2)(H+1)(H+2))] f_(r,rho)(n).       (8)

Strict subinvariance plus summability would suffice for universal ROOT if
D were positive everywhere: each surviving edge would strictly increase
the weight. A nonroot cycle is then impossible; a distinct infinite orbit
would contain infinitely many weights at least its positive starting
weight, contradicting summability. The unproved part is positivity at
every source, not the strict inequality on the already positive component.

## 4. Structural recursion retains two integer parameters and a certificate

The r posterior for a verified address is Beta(L+1,T+2). Appending one
sibling step has the exact update

    F_rho(T,L+1)/F_rho(T,L)=(L+1)/(L+T+3).             (9)

It approaches1 for large refuel depth, unlike a fixed exponential penalty.
Adding a backwards base edge of sibling depth k updates(T,L) to(T+1,L+k):

    ratio =rho(T+2) (L+1)^(overline k)
                   /(L+T+3)^(overline(k+1)).          (10)

For D, replace rho in(10) by(T+1)/(T+3); the sibling update(9) is unchanged.
These formulas follow by cancelling factorials in(2),(6). They are an
adaptive weight recursion without a search of an unknown forward orbit.

The counters determine the weight, not the source or a legal route.
Reordering k-values or valuation letters cannot be justified by their
commutative product. Retain the actual parent, source, depth and checked
ROOT word. The exported receipt_weight API first replays the supplied
strict ROOT certificate and independently decodes its actual(T,K,j).
An unvalidated two-counter record is not accepted as a source certificate.

## 5. Finite backward boxes have unconditional global error bounds

Start at the only seed1. For a verified parent c and k>=0 form y=S^k(c).
If3|y, there is no inverse base. Otherwise form

    b0(y)=(2y-1)/3 if y=5 mod6,
          (4y-1)/3 if y=1 mod6.

Omit b0=1, which is the root self-loop. This gives all inverse base edges,
each with an actual valuation1 or2. The parent's strict ROOT word transports
to S^k(c) by adding2k to its first valuation; when c=1,k>0, use(2k+2).
Prepending the actual b0 edge yields a ROOT receipt for the new base.
Finally generate its source siblings j the same way.

Thus all sources with T<N and K+j<M are generated by a **finite** backward
construction. No supplied child beyond the literal ROOT seed is assumed.
Unique primitive sibling decomposition and the deterministic G-parent
show that the parameter addresses are unique; a repeated base would force
a cycle in a previously verified root chain. Every rooted input appears
in a sufficiently large box. No assertion about other inputs is required.

At fixed(T,L), there are at most binom(L+T,T) sources: record the T
nonnegative inverse depths and the final j, whose sum is L. Some weak
compositions are forbidden by the residue or root guard, so this is an
upper bound. Multiplying by(2) gives the simple majorant

    total_(T,L) F_rho
      <= 2rho^T(T+1)/((L+T+1)(L+T+2)).                (11)

Summing its telescoping L-tail proves both summability and finite-box error:

    ||F_rho||_1 <=2/(1-rho),
    error(T<N,L<M)
      <=2rho^N/(1-rho)+2/((M+1)(1-rho)^2).            (12)

For D the corresponding majorant is

    total_(T,L) D <=4/((T+2)(L+T+1)(L+T+2)).          (13)

It gives ||D||_1<=4 and the rational bounds, for N>=1,M>=2,

    tail(T>=N)<=4/(N+1),
    tail(L>=M)<=4(H_M-1)/(M-1),
    error(box)<=4/(N+1)+4(H_M-1)/(M-1),               (14)

where H_M=sum_(i=1)^M1/i. The latter sum follows by partial fractions in
sum_(T>=0)4/((T+2)(M+T+1)). Both error bounds tend to0 along explicit
choices of N,M. The finite approximants are rational, supported on finite
sets, and include independently checkable first-hit ROOT words.

These proofs also give computable real values at every supplied source:
evaluate successive finite boxes with their uniform error bounds. A value
converging numerically to0 is not a decision that it is0. A positive value
is eventually exhibited by a literal certificate. The uniform L1 error
controls omitted **weight**, not omitted mu/nu source-prior mass; an
unrooted source, if any exists, has zero weight from the beginning.

A slightly sharper norm bound is available after summability is established.
Let M3=sum_(3|n)F_rho(n). The complete ROOT predecessor ray has total2,
so summing(4) over killed target rows gives

    ||F_rho||_1-2=rho(||F_rho||_1-1-M3).

Hence ||F_rho||_1<=(2-rho)/(1-rho); integrating this bound gives ||D||_1<=3.
This numerical bound is unrelated to the earlier three-bit source entropy.

## 6. The finite-alphabet reduction still needs unbounded address information

For the earlier prior nu(n)=8/(3*4^ell(n)), write g(b)=nu(b)h(b) at
r=1/16,rho=1/2. Its exact sibling scaling cancels the k factor, leaving

    h(b)<=a_b h(G(b)),
    a_b=(15/32)4^(ell(b)-ell(U(b))) in{15/128,15/32,15/8}. (15)

This is a three-element coefficient alphabet, not a three-state dynamical
system. For any fixed modulus m, sources b=16mt-1 satisfy
G(b)=24mt-1>b, k(b)=0, and both lie in the same residue class-1 modulo m.
Here a_b<=15/32<1. Consequently a positive residue-only correction h
cannot satisfy(15). More generally, no correction taking finitely many
positive values can work: arbitrarily long actual G-runs with valuation1
strictly increase the source, retain k=0, and repeat one of those values;
their cumulative a-factor is at most(15/32)^length<1. The run
3^i2^(R+3-i)-1,0<=i<=R, supplies exact witnesses.

This is the precise finite-projection obstruction. It does not rule out
unbounded source-dependent corrections or the certificate-dependent
mixtures above. It also does not license replacing pointwise constraints
by average constraints on a projected graph.

## 7. Exact checks and boundaries

Run normally and optimized:

    python -B 04-computation/experiments/collatz_mixed_refuel_weights_20261005.py
    python -B -O 04-computation/experiments/collatz_mixed_refuel_weights_20261005.py

The code uses no external packages and has no import-time computation.
Controls include exact polynomial integration for T0..8,L0..12; five
fixed-r and three fixed-rho comparisons; all posterior updates at
T0..6,L0..9,k0..6;40 root-ray partial sums and exact tails; complete
backward boxes T<6,L<12 and T<7,L<14; independent actual routes on all2048
odd inputs below4096 to check finite box membership; and exact full-fibre
incoming sums at1332 targets in the smaller box. Its1990 sources include
688 primitive bases, with a strict ROOT word for every source. The larger
box contains8292 sources. These counts are finite, not a global census.

The smaller box has total F_(1/2) weight427967041/196035840 and total D
weight890213479/428828400. Invalid source types, valuations, root padding,
rates and box parameters are explicitly rejected. The independently searched
finite literal controls are used only to audit the backward generator;
the generator itself consumes no forward-orbit oracle or unproved seed.

Preserved: exact source, actual ROOT receipt, sibling/base counters, mixture
weight and a global numerical remainder. Lost by storing just(T,L): address,
word order, actual guard and source identity. Still missing for universal
coverage: a positive lower bound or another proof for every supplied source.
