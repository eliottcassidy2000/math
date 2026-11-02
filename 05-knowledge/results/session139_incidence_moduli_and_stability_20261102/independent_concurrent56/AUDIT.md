# Independent audit of the concurrent n >= 56 frontier

**Conclusion: PASS.** The pinned concurrent proof establishes exact weighted
stability `K_n=kappa_n` for every integer n>=56, with precisely the stated
positive-variance ordinal fair-pair equality family and transitive equality
at zero variance. It leaves orders 10 through 55 unresolved. No mathematical
correction is required.

This independent review supplements the source's historical same-agent
review statement without modifying its scientific text. The source's
attribution to distinct concurrent and Session139 contributions is preserved;
no sole-priority or external novelty claim is made.

## Immutable reviewed sources

Both files were read with targeted `git show` from commit
`c1e79aee0e226c2c03409b96b17a92ba88bfcca1`:

- `05-knowledge/results/stability_priority_spectrum_20261102/STABILITY_56.md`,
  SHA-256 `c07cdd8fbb97d1d5848c32768c8decb23dad66df470de62751f33742c4610e0a`.
- `05-knowledge/results/stability_priority_spectrum_20261102/verify_stability_56.py`,
  SHA-256 `3c9d29138c26eba4fc73a83ed4a46226128ae582b8fc6056bbcbc558b782c6e0`.

The factor-four graph bound was independently reviewed in the preceding
`independent_concurrent68/AUDIT.md`, including its exact discarded-degree
identity and disjoint selected/discarded payment. It is credited to the
concurrent `stability_priority_moduli_20261102/STABILITY_N500.md` source at
`c16226249f1f27a335cf7c84afdeaa90879b9580`. The retained negative-BR mechanism
and flexible local dual are credited as in the reviewed spectrum packet.
The tensor reserve and Session139 bootstrap are accepted inputs, as are the
Session138 binary, retained-matching, and equality theorems.

The present audit does not depend on the n>=68 theorem being the strongest
result. That earlier proof remains a valid alternate proof on its own range.
All upstream scientific files and all earlier audited packets are unchanged.

## 1. Nonlinear input and stronger square completion

For the full positive-weight greedy matching, round its complement first.
The nonnegative omitted corrections and the accepted matching-only tensor
reserve imply

    Delta>=B+1-product_i(1-q_i)-1_(n=2r)product_i q_i
         >=B+1-exp(-M)-4^(-R).

The product correction is nonzero only for a perfect chosen matching;
then r=R and its product is at most 4^(-R). For an imperfect matching
subtracting 4^(-R) is merely a weakening. The matching weights do not
change on complement-rounding branches. Their terminal expectations average
to the original independent binary-rounding expectation B exactly once.
This is not a false conditional single-edge factor inequality.

Assume V>0 and Delta<=kappa V. Since kappa<=4/R<=1/7 and g(V)/V decreases,
V<=1 is excluded by `g(1)>49/200>1/7`. Thus V>1 and
`kappa V>B+49/200`. The partial Taylor sum `16319/315>50` justifies the
initial strict bound on g(1).

For any established positive reserve c, the accepted stronger graph bound
`B>=W>=2D^2/n` gives

    M=V-D > c/kappa+2D^2/(n kappa)-D
          =c/kappa-n kappa/8
             +[2/(n kappa)](D-n kappa/4)^2
          >=cR/4-1-1/(2R)=m_c(R).

All denominators are positive. The last step uses kappa<=4/R and
n<=2R+1, not an erroneous parity replacement. The derivative of m_c is
`c/4+1/(2R^2)>0`.

## 2. Independently certified bootstrap

At R=28 the three mass bounds are exactly 122/175, 139/56, and 7493/1400.
The upstream script's degree-12 rational checks pass. Independently, this
audit proves the same advances using lower Taylor degrees 4,5,5 and the
strictly weaker penalty bound `4^(-R)<1/10000`:

| Advance | Taylor degree | `(1-c'-1/10000)T(m)-1` |
|---|---:|---|
| 49/200 to 1/2 | 4 | 28404006581/9378906250000 |
| 1/2 to 91/100 | 5 | 21073575206681/660878131200000 |
| 91/100 to 99/100 | 5 | 340875870324523013869/2151296000000000000000 |

All three rational slacks are positive. The masses and remaining
probability factors are positive, and exp(m) strictly exceeds each finite
positive Taylor sum. Increasing m_c(R), together with decreasing perfect
penalty, extends each advance to every R>=28. Hence all candidate equality
or counterexample configurations satisfy `kappa V>B+99/100`.

No final reserve is assumed before it is established. The same B is retained
through the bootstrap and the later matching argument.

## 3. Heavy matching, exact removed mass, and tail

Take `c=99/100`, `t=1/7`, `d=381/500`, and `lambda=t^2/d=500/18669`.
If L>=1 selected weights are at most t and their total is U, then U<=tR.
The retained heavy matching satisfies

    h>=4M_h>(c-4t)R-4-2/R.

At R=28 the right side is 2677/350, exceeding 7 by 227/350. Its derivative
is `c-4t+2/R^2>0`, so integer h>=8 uniformly. Thus the retained matching
is nonempty before its coefficient is used. With m=R-L>=h>=8, L>=1 also
ensures n>2m, allowing the accepted nonperfect coefficient a_m.

Retaining B and U, the matching inequality expands to

    L(kappa V-B)-4(1-epsilon)U
      <=(kappa R-4)V+4epsilon V+4(1-epsilon)D-BR
      <=4epsilon V+4D-BR,

where epsilon=(3/4)^(R-L). Since kappa V-B>c and U>=0, this implies
`cL-4U<4epsilon V+4D-BR`. No negative-BR term is discarded, and U is not
replaced by the weaker estimate tL at this stage.

The variance bound gives `4V<=Q_R=R+33/2+8/R`. One can independently obtain
this by `4V<=n/2+2n kappa`: even n gives R+16 and odd n gives Q_R. It is
valid at n=56 without importing any large-n bound such as V<=n/7.

For fixed R, the second derivative of
`log((3/4)^(R-L)/L)` is 1/L^2>0. Therefore its maximum on `[1,R-8]` is at
an endpoint. Both endpoint bounds decrease for R>=28:

- `Q_R(3/4)^(R-1)` has successive ratio below `(3/4)(1+1/R)<1`.
- `Q_R/(R-8)=1+(49/2+8/R)/(R-8)` has positive decreasing excess.

At R=28 both are strictly below c-d=57/250. The tighter endpoint has the
claimed exact slack `1751781/458752000>0`. Thus

    dL-4U<4D-BR.

Every light summand satisfies `d-4w>=d-4t=667/3500>0`, so D>0 before the
strict contradiction is used.

## 4. Flexible dual and endpoint strictness

The exact coefficient identity is

    wS/2+S^2/4-tS+lambda(d-4w)
      =(S/2-t+w/2)^2+w(t-4lambda-w/4).

Its constant terms agree because lambda*d=t^2. On the entire light domain
`0<=w<=t`, the second factor is at least
`3t/4-4lambda=1/74676>0`. Thus the local bound holds for every S>=0,
including S=0 and w=0.

Greedy groups count disjoint wedge pairs, so the light local inequalities
and the heavy inequalities W_i>=tS_i can be summed simultaneously. Using
the strict retained deficit budget gives

    B>=W>(t-4lambda)D+lambda BR,
    t-4lambda=667/18669>0.

If `1-lambda R<=0`, then B>=0 makes this impossible because D>0. Otherwise
multiplication of `B<=kappa D<=4D/R` by that positive coefficient is valid:

    (t-4lambda)D < B(1-lambda R)
                  <=(4/R-4lambda)D
                  <=(t-4lambda)D.

The final comparison is `Rt>=4`, with equality exactly at R=28. The earlier
strict inequality closes this equality endpoint; no strict Rt>4 assumption
is hidden. At the endpoint `1-lambda R=667/2667>0`, but the proof correctly
handles the opposite coefficient sign at larger orders as well.

## 5. Boundaries, scope, and equality

- V=0 means every edge is binary. The accepted binary theorem supplies
  nonnegative gap and transitive equality for n>=56.
- Positive V ensures a nonempty greedy matching. Under L>=1, the heavy
  matching is proved nonempty before using its coefficient.
- The positive deficit forces D>0. Neither B nor D is divided by anywhere;
  B=0 cannot escape the strict positive-left-side contradiction.
- Individual light discarded masses, or their entire total, may be zero.
  The pointwise dual and disjoint summation still apply.
- Once L=0, every selected weight exceeds t>=4/R>=kappa. Any D>0 gives
  `W>=sum w_i S_i>tD>=kappa D>=B>=W`; hence D=0 and B=0.
- The remaining uncertainty is a matching. The accepted retained-matching
  bound gives Delta>=kappa V. For equality, every positive-probability
  binary completion has zero nonnegative gap and is transitive; its
  Hamilton-path expectation is H=1. The accepted scalar classification
  gives exactly maximal ordinal fair pairs and one singleton in odd order.
  Those configurations attain the displayed coefficient.

The newly covered interval within 10..1499 is 56..1499, containing 1,444
orders. The remaining 10..55 comprises 46 orders. This is not a claim that
55 is an obstruction or that the threshold is optimal.

## Exact executable evidence

`audit_n56.py` is independently written and emits 40 active checks using
explicit exceptions. Its polynomial identity is checked coefficient-wise,
and its bootstrap uses the distinct degrees 4,5,5 and coarse penalty above.

- `normal.json` and `optimized.json`: PASS and byte-identical.
- `negative_threshold55{,_optimized}.json`: threshold 55 fails the first
  reserve advance, exit status 1 in both modes. This rejects this
  certificate there, not the theorem itself.
- `negative_dual{,_optimized}.json`: replacing d by 3/4 while retaining
  lambda=t^2/d preserves the square identity but makes the local remainder
  `3t/4-4lambda=-1/588`. It is rejected in both modes, exit status 1.
  This is a genuine local counterexample: w=S=t makes the square zero
  and the proposed lower bound's difference equal -1/4116.
- `producer_normal.json` and `producer_optimized.json`: the unmodified
  upstream script, materialized into an isolated audit-only copy, passes
  4,351 gates in both modes; reports are byte-identical. Its 56..1499 replay
  is finite checking evidence, not the justification for all larger orders.

The producer's local output was generated only beside the isolated audit
copy; no upstream scientific file or previously frozen result was modified.
No Git or remote write was performed. No additional threshold optimization
was undertaken.
