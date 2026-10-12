# Terminal congruence completion gives independent, refinement-stable ROOT floors

2026-10-11. **PROVED, scoped and self-audited:** an explicit two-parameter
family of strict five-step ROOT words; its exact quarter-price and beta
scores; a positive marked measurement margin at every stated refinement.
**FINITE-EXACT:** literal controls on36 family members. **OPEN:** positive
independent readouts for arbitrary prescribed sources and universal Collatz.
This is a new application of inherited word arithmetic, not a claim that
bounded-length inverse families are new in the literature.

## Inheritance, hostile, and the reroute

The closest mechanism is the authenticated native bridge in
[atomic run bridges](collatz_atomic_run_bridge_20261011b.md), with the
refinement-uniform source-weight transport in
[bridge weight transport](collatz_bridge_weight_transport_20261011.md).
The measurement and its sharp residual are inherited from
[quarter-price ROOT readouts](collatz_quarter_price_source_20261011.md).

The canonical hostile is an ungrounded child: a paid edge to a smaller
integer has not yet proved that the child reaches ROOT. The corrected near
miss is replacing the selected source by a rooted representative in each
finer residue class. In fact
[THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md)
already proves density of the ROOT basin in every binary prefix. Proving
that density again would not improve the marked-source problem.

The least-used sidecar here is the **terminal exponent congruence**, which
can ground a proper infinite subfamily without an orbit census or a stored
certificate for any target member. The board is **actual source / common
future / terminal phase / exact score / alias bill / persistent floor**.
The anchor is this ROOT target; the orthogonal field-context and clock-divisor
lanes are linked in the [geometry packet](cyclotomic_unit_component_20261011.md).

Let U(n)=oddpart(3n+1). For a>=18 divisible by18 and t>=0, the inherited
bridge is

    N=(2264*2^a-23)/27+96*2^a*t
       --(1,2,a,1)--> 425+486t
       <--(1,1,1,1,3)-- 223+256t.                       (1)

Instead of repeatedly asking whether the child belongs to a known controller,
choose t so the common endpoint has a single authenticated ROOT step.

## 1. Solve the terminal equation exactly

Require

    425+486t=(4^k-1)/3,
    t=(4^k-1276)/1458.                                  (2)

For k>=1, divisibility by2 is automatic. Modulo729, multiplying by4 gives

    4^k=1276 mod729 iff 4^(k+1)=1 mod729,

because4*1276=7*729+1. The order of4 modulo3^6 is3^5=243. One elementary
proof is v3(4^b-1)=1+v3(b), obtained by the standard cubic lifting induction.
Thus (2) is integral exactly when k=-1 mod243. Its least positive k is242,
which also makes t positive.

For integers j,r>=1 put

    a=18j,      k=243r-1,
    N_(j,r)=[8*2^a*(2*4^k-5)-207]/243.                 (3)

Then N_(j,r) is a positive odd integer with strict first-hit word

    boxed: (1,2,18j,1,2k).                             (4)

Indeed (3) is exactly (1) after substituting (2), and the endpoint in (2)
is positive odd with3z+1=4^k. Every preterminal state exceeds1. All valuation
guards in (1) remain inherited, rather than just verifying an affine formula.
The child itself has the strict six-step word (1,1,1,1,3,2k).

This grounds an infinite family in both parameters. It does not cover every
t in (1), and it does not alter a user's prescribed source to fit (3).
The first member has a499-bit source while requiring just five odd
steps. The [exact output](source_context_localization_20261011.json) records
its integer and bit length. Large numerical size is not a complexity claim:
this family was designed to have a short certificate.

### Decode the parameters from the actual source

Membership is decidable without following a Collatz orbit. Given positive
odd N>1, form D=243N+207, v=v2(D), a=v-3, and

    z=(D/2^v+5)/2.

Accept precisely when a is a positive multiple of18, z=4^k is a power of4,
and k>=242 with k=-1 mod243. Then j=a/18 and r=(k+1)/243 are the unique
parameters in (3). Rearranging the defining equation proves both directions.
The accompanying `terminal_parameters` function implements this decoder.
Rejecting membership says nothing about whether N reaches ROOT by another
route. A receipt includes N and these decoded parameters; neighbouring
integers cannot silently reuse it.

## 2. Exact weights and a source-height floor

Use the inherited counters K=sum floor((a_i-1)/2), L=tau-1 and scores

    f(n)=4^(-K)(3/4)^L,
    W(n)=2 K! tau!/(K+tau+1)!.

For (4), tau=5 and K=9j+k-2, so

    f(N)=81/4^(9j+k+2),
    W(N)=240/[(K+1)(K+2)(K+3)(K+4)(K+5)(K+6)].        (5)

The source-normalized quarter score has the exact expression

    Z(N)=(3N+1)f(N)/4
        =1/4 - 5/(8*4^k) - 63/(32*2^a*4^k) > 1/5.    (6)

The last inequality follows already from a>=18 and k>=242; no approximation
to a logarithm or a simulated trajectory is used. Consequently

    boxed: f(N)>4/[5(3N+1)].                            (7)

This source-height bound is compatible with the earlier Mersenne hostile:
that hostile rules out such a floor uniformly over all positive sources.
Here the exact guard (3) is essential. Neither the score nor the guard is
inferred from the approximate size of N.

## 3. The independent lower readout passes the alias bill

Fix one actual source N from (3); it is unchanged throughout this section.
Set m=(N-1)/2 and

    nu(m)=(m+1)^(-6)-(m+2)^(-6),
    q0=3m+2=(3N+1)/2,
    q_d=q0(4^d-1)/3,               d>=1.

The first alias N+2q_d is the actual d-th sibling S^d(N), S(n)=4n+1.
For the true quarter-score readout define

    M_d=sum_{b>=0} nu(m+b*q_d) f(N+2b*q_d),
    A_d=nu(m)+4^(-d)nu(m+q_d),
    E_d=9/(448*q_d^8).

The inherited theorem gives

    A_d f(N) <= M_d <= A_d f(N)+E_d.                    (8)

Equations (3)--(5), proved symbolically, provide the independently justified
lower measurement L_d=A_d f(N). This does not estimate M_d from a finite
ROOT bank and does not count ROOT positivity as an unproved input. The
convergence information enters through the explicit family proof (4).

The elementary integral inequality

    nu(m)=integral_(m+1)^(m+2) 6x^(-7) dx
          >6/(m+2)^7 >= 6/q0^7

and (7) imply

    nu(m) f(N) > 12/(5q0^8) > 9/(448q0^8).             (9)

Thus the readout clears the residual already at d=1. More strongly, q_d>=q0
and A_d>=nu(m), so the exact signed extraction obeys

    f(N) >= (L_d-E_d)/A_d
           =f(N)-E_d/A_d
           > (1777/1792) f(N) >0                       (10)

for every d>=1. In other words the certified floor loses less than15/1792
of the exact family score, uniformly over these source-preserving readings.
For nested congruence refinement one may take d=2^h, because q_d divides
q_(2d). The larger all-d statement does not falsely declare every consecutive
q_d pair nested. Retaining the maximum of valid prior floors is also lawful.

The rightmost quantity in (10) is an explicit rational formula (5), not an
unknown value queried from a score oracle. Applying the inherited inequality
W >= (9/16)f^2 gives a positive beta-mixture floor as well; (5) gives the
stronger exact beta score on this family.

## 4. What this changes, and what remains unproved

The completion move replaces the terminal obligation, not the marked source
during measurement. Its general reusable form is:

1. Retain an actual guarded word whose endpoint is c+d*t.
2. Solve 3(c+d*t)+1=4^k, including ordinary positivity and every modulus.
3. Substitute that solution into the source parametrization, not into the
   definition of an already fixed target.
4. Compile the resulting word into an explicit score and pay the source's
   own residual bill.

A residue equation can have no solution, or it can have a solution only on
a sparse exponential parameter set. Here the large phase243 is the actual
cost. It would be invalid to use this result to certify t=0 through the new
five-step word; t=0 has its separate inherited223 certificate. Wrong terminal
phases are rejected by the arithmetic before any score is issued.

The productive remaining question is whether **changing connector words**
can solve such an endpoint completion while leaving an arbitrary specified
source fixed, with a proved finite selection rule. The fixed connector (1)
does not do that. Universal positivity remains open, and the older
[Poisson source dual](collatz_poisson_source_dual_20261005.md) and
[localized resolvent](collatz_localized_resolvent_floor_20261005.md) remain
stronger general interfaces than this proper-family application.

The graph analogue is similarly precise: generic Cartesian refinement keeps
all old boundary colourings, while an intentional authenticated seam can
change the relation. In both cases the new constraint must be supplied;
repackaging or increasing the representation's size does not supply it.

## Reproduction

```text
python3 -B 04-computation/experiments/source_context_localization_20261011.py
python3 -O -B 04-computation/experiments/source_context_localization_20261011.py
```

The Collatz universe is j,r=1..6 with four measurement depths per source.
The verifier independently executes the specified words and checks their
literal valuations, strict first hit, child/common endpoint, exact counters,
weights and margins. It does not discover a word by running an unbounded
orbit. Adjacent invalid phases and the second valid phase485 are retained
as hostile and positive controls. The universal family proof is (1)--(10).

The aggregate package has **2,473** always-active checks. Normal and optimized
runs agree byte for byte; the source hash was frozen before both runs and
checked again before saving their output. Script SHA256:
`9dc3590775b92c27463b1ee910df88925787457086686e4de63bcc628b5d8f19`.
Output SHA256:
`879073a6c804e3a22e6424eb8ff8571e944903dfa0151990e6c309748544167e`.
