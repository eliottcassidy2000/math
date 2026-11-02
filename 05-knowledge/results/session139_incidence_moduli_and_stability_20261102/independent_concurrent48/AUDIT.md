# Independent audit of the concurrent n >= 48 theorem

**Conclusion: PASS.** The pinned proof establishes exact weighted stability
`K_n=kappa_n` for every integer n>=48, with the stated ordinal fair-pair
positive-variance equality family and transitive equality at V=0. No
mathematical correction is required for this theorem. Orders 10 through
47 remain unresolved.

The optional general scalar "if and only if" statement is valid in the
intended positive-parameter domain **d>0**, as assumed by its beta>0 proof.
This domain is made explicit here; the theorem's concrete d=5/6 is positive.
There is no effect on the n>=48 conclusion.

## Immutable sources and attribution

The following files were read by targeted `git show` at commit
`dde756b2c2d4506156d637a536910fc80a35450d`:

- `05-knowledge/results/stability_748_guarded_encodings_20261102/STABILITY_48.md`,
  SHA-256 `711863194c6440d153315f4d8b5153f02ca57750f590922e0be316ebb72a493b`.
- `05-knowledge/results/stability_748_guarded_encodings_20261102/verify_48.py`,
  SHA-256 `40f9acda78dbba05639e7b735309c93f0bd3c9dcdf0aafe23c62904478ffb7a8`.

The reviewed source credits the concurrent n56 flexible-dual framing, the
concurrent stronger discarded-degree estimate, Session139's accepted tensor
reserve, and the previously accepted matching and equality inputs. Those
credits are preserved. This audit does not claim discovery or sole priority.
The original source's same-agent-audit statement is left unchanged; this
separate artifact supplies an independent review of the pinned version.

Earlier n56, n68, n140, and other scientific proofs/audits remain unchanged.
No source/canonical file, Git state, or remote state was modified.

## 1. The maximum is a valid single-group bound

For a greedy selected edge w and its assigned discarded mass S, the two
accepted inequalities apply to the same counted group energy W_i:

    W_i>=wS,       W_i>=wS/2+S^2/4.

Therefore W_i is at least their maximum G, not their sum. The stronger
quadratic bound follows from the selected/discarded pairs and the two arm
energies, using every assigned edge weight <=w. Extra pairs only increase
the actual energy. Greedy groups partition the edges, so the within-group
unordered wedge pairs are disjoint even if groups share vertices.
Consequently summing the group maxima is legitimate and never spends the
same wedge twice.

## 2. Full-continuum two-regime minimization

For t>0, d=5t, lambda=t/5, put u=w/t in [0,1] and v=S/t>=0. The desired
normalized difference is

    F(u,v)=max(uv,uv/2+v^2/4)-v+1-(4/5)u.

The maximum switches at v=2u. More generally write beta=4t/d and define

    P1(u)=u(1-beta-u/4),
    P2(u)=2u^2-(2+beta)u+1.

The independent checker verifies these polynomial identities coefficient
by coefficient, not by sample interpolation:

    F_linear-P2 = (1-u)(2u-v),
    F_quadratic-P1 = (v-2+u)^2/4,
    P2-P1 = (3u-2)^2/4,
    F_quadratic-P2 = (v-2u)(v/4+u-1).

For 0<=u<=2/3, the quadratic vertex v=2-u is feasible (v>=2u).
Both branches are at least P1, and that vertex attains P1. For
2/3<=u<=1, the linear branch is minimized at v=2u, while the final
identity shows the quadratic branch cannot fall below that boundary
value. Thus the claimed minima are exactly P1 and P2, respectively.

At beta=4/5 these are

    P1=u(1/5-u/4)>=u/30>=0,
    P2=2(u-7/10)^2+1/50>0.

This proves the local dual for every real 0<=w<=t, S>=0, including w=0
and zero discarded masses. It is not inferred from the producer's grid.

For the optional generalization with d>0, completing the second square
shows `P2=2(u-(2+beta)/4)^2+1-(2+beta)^2/8`. Nonnegativity therefore
holds for `0<beta<=2sqrt(2)-2`; that range also lies below 5/6, which
makes P1 nonnegative on its interval. Above the boundary and beta<=2,
the center lies in [2/3,1] and gives a negative value. For beta>2,
u=1 gives a negative value. Hence the scalar dual's exact boundary is
`d/t>=2+2sqrt(2)` in the positive-d domain.

The hostile d/t=24/5 has beta=5/6. At u=17/24, v=17/12 its normalized
deficit is exactly -1/288. This disproves that **scalar relaxation** for
that ratio, not the weighted tournament inequality. No actual tournament
counterexample is asserted.

## 3. Nonlinear bootstrap with the same binary expectation

Under V>0 and Delta<=kappa V, kappa<=4/R<=1/6 and decreasing g(V)/V
exclude V<=1 because `g(1)>49/200>1/6`. The accepted partial exponential
sum justifies the strict first reserve.

The stronger discarded-degree input is `B>=W>=2D^2/n`. Completing its
square gives, for every already established c>0,

    M>c/kappa-n kappa/8
      >=cR/4-1-1/(2R)=m_c(R).

The first expression follows by retaining B, not bounding M and D
independently. The second uses n<=2R+1 and kappa<=4/R. Its derivative in
R is `c/4+1/(2R^2)>0`.

The accepted tensor reserve is

    Delta>=B+1-exp(-M)-4^(-R).

The product subtraction occurs only when the chosen matching is perfect,
when its size is R and its variance product is at most 4^(-R). On every
other complement-rounding branch the subtraction is zero. The original
binary expectation B is retained once; matching weights are unchanged.
No false conditional single-edge multiplicative factor is used.

The upstream degree-12 exact bootstrap passes. Independently, this audit
uses the weaker penalty `4^(-R)<1/10000` and Taylor degrees **2,3,6,11**.
At R=24 the four masses are exactly

    539/1200, 259/240, 691/240, 5543/1200.

For the advances `49/200 -> 7/20 -> 13/20 -> 47/50 -> 99/100`, the exact
slacks `(1-c_next-1/10000)T(m)-1` are, respectively,

    212322379/28800000000,
    3765319801/829440000000,
    50314354522091675519/1375941427200000000000,
    4211957060596623741919072456319700832001/
      4279728215162880000000000000000000000000000.

They are all strictly positive. Every mass and remaining probability factor
is positive. Finite positive Taylor sums lie strictly below exp(m), and
increasing mass bounds/decreasing perfect penalties extend every advance
to every R>=24. Thus all candidates satisfy `kappa V>B+99/100`.

## 4. Exact matching slots force twelve retained edges

Let L>=1 selected weights be light (w<=t), U their total mass, and r_h the
number of heavy selected edges. There are at most R selected edges in
all, so the exact slot constraint is

    U<=tL<=t(R-r_h).

Together with r_h>=4M_h and M_h=M-U, this gives

    r_h > (c-4t)R-4-2/R+4t r_h,
    r_h(1-4t)>g(R)=(c-4t)R-4-2/R.

This implication does not assume a perfect greedy matching or assume
r_h>0 beforehand. With t=1/6 the coefficient 1-4t=1/3 is positive.
At R=24, `g(24)-11/3=1/100>0`, and g increases. Hence integer r_h>=12.
In particular the retained matching is proved nonempty before using its
coefficient. With m=R-L one has 12<=r_h<=m<=R-1 and n>2m.

## 5. Retained light mass, tail payment, and strict contradiction

The decreasing nonperfect matching coefficient a_r, retaining B and U,
gives the exact expanded upper bound

    L(kappa V-B)-4(1-epsilon)U
      <=(kappa R-4)V+4epsilon V+4(1-epsilon)D-RB
      <=4epsilon V+4D-RB.

Since kappa V-B>c, replacing the left side by cL-4U is a strict weakening.
Here epsilon=(3/4)^(R-L). The negative-RB term and actual removed mass U
are both retained with their correct signs.

The variance-wedge bound gives `4V<=Q_R=R+33/2+8/R` for both parities:
the direct even envelope is R+16 and the odd envelope is Q_R. Convexity
in L of `log((3/4)^(R-L)/L)` places the maximum at L=1 or L=R-12.
The first endpoint decreases because its successive ratio is less than
`(3/4)(1+1/R)<1`. The second decreases because

    Q_R/(R-12)=1+(57/2+8/R)/(R-12)

has positive decreasing excess. At R=24 the exact endpoints are

    7688359604205/140737488355328,
    14467005/134217728,

both below c-d=47/300. Their exact positive slacks are

    1077038517859729/10555311626649600,
    492032929/10066329600.

Therefore `dL-4U<4D-RB`. Every light summand is at least d-4t=1/6>0,
so D>0. Applying the local maximum dual to light groups and W_i>=tS_i
to heavy groups, with disjoint summation, yields

    (1-lambda R)B>hD,
    lambda=1/30,    h=t-4lambda=1/30>0.

If 1-lambda R<=0, B>=0 contradicts D>0. Otherwise the positive coefficient
allows `B<=kappa D<=4D/R`, giving

    h<(1-lambda R)kappa<=4/R-4lambda<=t-4lambda=h.

This is contradictory even at R=24, where Rt=4 exactly. The strict
inequality was established earlier; strict Rt>4 is not required.

## 6. Boundary cases, equality, and scope

- V=0 is binary and is covered by the accepted binary nonnegative gap and
  transitive equality classification.
- Positive V ensures an initially nonempty positive-weight greedy matching.
- Under L>=1 the slot argument proves r_h>=12, L>0 permits per-light-edge
  tail division, and the positive deficit proves D>0 before division by D.
- B=0, individual zero discarded masses, and zero total discarded light
  mass cannot escape the proof. The local continuum inequality covers them.
- Once L=0, selected weights exceed t>=4/R>=kappa. Any D>0 gives
  `W>=sum w_i S_i>kappa D>=B>=W`; hence D=B=0.
- The remaining uncertainty is a matching. The retained-matching theorem
  gives Delta>=kappa V. Equality makes every supported binary completion
  zero-gap and transitive; its Hamilton-path expectation is H=1. The
  accepted scalar equality theorem gives exactly maximal ordinal fair
  pairs, plus a singleton in odd order. These constructions attain kappa.

The new interval beyond the incoming n56 result is precisely 48..55.
The remaining unresolved interval is 10..47. At R=23 the present slot
bound no longer forces twelve heavy edges and Rt<4; its failure there is
not evidence against the stability theorem itself.

## Independent and upstream executable evidence

`audit_n48.py` is independently written, with 47 explicit active checks.
It verifies the continuum polynomial identities coefficient-wise and uses
the distinct bootstrap degrees and penalty described above. In addition,
it exhausts 760 weighted graphs on 1..4 vertices with edge weights in
`{0,1/6,1/4}`, a different grid from the producer. It checks each group's
actual energy against the maximum invoice, disjoint summation, and the
whole light/heavy dual. These finite controls supplement the uniform proof.

- `normal.json` and `optimized.json`: PASS and byte-identical.
- `negative_threshold47{,_optimized}.json`: threshold 47 fails the first
  reserve advance, exit status 1 in both modes.
- `negative_dual{,_optimized}.json`: ratio 24/5 makes the second continuum
  minimum -1/288; rejected in both modes, exit status 1. This is the
  formal scalar hostile, not a tournament counterexample.
- `producer_normal.json` and `producer_optimized.json`: the unmodified
  upstream script, run in an isolated audit-only copy, passes in both
  modes and produces identical reports: 25,351 local grid checks, 729
  complete four-vertex graph controls, and 277 finite-R checks.

The upstream verifier uses only the Python standard library and has no
sibling-import dependency. Its finite controls are not used to infer the
continuum local lemma or the all-order conclusion. No new threshold search
or optimization was undertaken during this audit.
