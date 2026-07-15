---
id: THM-806
title: Scale-one Hamming-three collar closure
status: PROVED (oriented collar-handoff reduction, using settled lower-runner LRC) + FINITE-EXACT (5,713,539 component-containment rows)
source: codex-2026-07-15-S10 (continuation of THM-804)
depends_on:
  - LRCUpTo13  # only the 10- and 11-speed bounds
  - THM-804    # only for the arbitrary-AP-scale corollary
related: [THM-770, THM-795, THM-800, HYP-6775, HYP-6800, HYP-6820]
verification:
  - 04-computation/lrc13_scale_one_hamming_three_collar_closure_codex_S10.py
  - 05-knowledge/results/lrc13_scale_one_hamming_three_collar_closure_codex_S10.out
---

# THM-806 — Scale-one Hamming-three collar closure

Put

```text
delta=1/13,                         [12]={1,...,12},
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

## Theorem

### A. Every proper scale-one triple lift is loose

Let `r,s,t` be distinct members of `[12]` and let `i,j,k>=1`.  Then

```text
B=([12]\{r,s,t}) union {r+13i,s+13j,t+13k}              (1)
```

satisfies

```text
M(B)>1/13.                                               (2)
```

Thus the scale-one chart left open by THM-804 is empty.

### B. Uniform shallow rigidity through Hamming radius three

Let `13` not divide `c`.  Replace `m in {1,2,3}` distinct members `cr_a` of
`c[12]` by positive integers `w_a` satisfying

```text
w_a=cr_a (mod 13),                  w_a!=cr_a.            (3)
```

Every resulting packet is loose at `1/13`.  The cases `m=1,2` are THM-795
and THM-800.  For `m=3`, THM-804 forces `c|w_a` under a hypothetical equality;
division by `c` gives (1), contradicting Part A.

If a three-labelled presentation allows one or two zero-height labels
`w_a=cr_a`, it is simply a Hamming-one or Hamming-two packet and is covered by
the corresponding earlier theorem.  The all-zero row is the original tight
AP and is, of course, excluded by the word "proper."

## 1. A universal owner collar

Write

```text
R={r,s,t},                         P=[12]\R.
```

For `r in R`, choose `a_r in [12]` with `a_r r=1 (mod 13)` and put

```text
tau_r=a_r/13.                                            (4)
```

At `tau_r`, a retained core speed `p in P` has clearance at least `2/13`,
except possibly `p=13-r`, whose signed phase is `-1/13`.  For

```text
0<h<1/156,                                               (5)
```

every noncomplementary speed remains strict-safe because

```text
||p(tau_r-h)|| >= 2/13-p h > 2/13-12/156=1/13.          (6)
```

The complementary phase moves from `-1/13` farther into the negative safe
side; throughout (5) its norm lies strictly between `1/13` and `2/13`.
If the complement is also missing, all retained speeds fall under (6).
Consequently

```text
(tau_r-1/156,tau_r)
  subset {x:min_(p in P)||px||>1/13}.                    (7)
```

This `1/156` collar is uniform in the missing residue triple.

## 2. The half-open cross-handoff relation

Let `u_r=r+13i` be the replacement belonging to `r`.  At `tau_r` its signed
phase is `+1/13`, so its closed danger tooth immediately to the left is

```text
[tau_r-2/(13u_r),tau_r].                                 (8)
```

If `u_r>24`, then `2/(13u_r)<1/156`.  The left endpoint

```text
x_r=tau_r-2/(13u_r)                                     (9)
```

lies strictly inside the core-safe collar (7).  At (9), the own replacement
has phase `-1/13` and becomes safe immediately to the left.  If the full
packet were tight, at least one other replacement `u_q` would therefore have
to cover the left germ at `x_r`.

For a positive speed, a closed danger tooth covers immediately to the left of
an endpoint phase exactly on the half-open arc

```text
(-1/13,1/13].                                            (10)
```

The `+` endpoint is included because moving left enters its tooth; the `-`
endpoint is excluded because moving left exits.  Put

```text
z=q r^(-1) in F_13^*,             lambda=u_q/u_r.        (11)
```

Multiplying the phase at (9) by `13`, the exact handoff condition is

```text
-1 < z-2lambda-13m <= 1            for some m in Z.      (12)
```

Draw an arrow `q -> r` when (12) holds.  If all three replacements exceed
`24`, tightness and (7)--(10) give positive indegree at every one of the three
owner vertices.  Such a digraph contains a directed 2- or 3-cycle.

## 3. Handoff cycles are impossible

First observe that if `0<lambda<1`, then (12) forces

```text
z=2,                         1/2<=lambda<1.              (13)
```

Indeed `0<2lambda<2`, distinct labels exclude `z=1`, and every translate
`m!=0` lies outside `(-1,1]`.  The remaining inequality has the unique
solution `z=2` and the stated ratio band.

### No directed 2-cycle

Let `u_q<u_r` and suppose both arrows exist.  The arrow `q->r` gives `z=2`
and `u_r/u_q<=2`.  The reverse residue is `2^(-1)=7`, but a handoff with
residue `7` requires provider/owner ratio in

```text
[3,4)                                                      (14)
```

in its first positive band.  It cannot occur at a ratio at most two.

### No directed 3-cycle

Name the three speeds `A>B>C`.  There are two cyclic orders.

1. If `B->A`, `C->B`, and `A->C`, the two descending arrows give residue
   ratios `2,2` and speed ratios at least `1/2`.  Hence the residue of `A/C`
   is `4^(-1)=10`, while `A/C<=4`.  But residue `10` first becomes eligible
   only for ratio in `[9/2,11/2)`, a contradiction.

2. If `C->A`, `A->B`, and `B->C`, then `C/A>=1/2`, so `A/C<=2`.  Both ratios
   `A/B` and `B/C` lie strictly between one and two.  In that metric band,
   (12) permits only residues

   ```text
   {2,3,4}.                                               (15)
   ```

   The residue product around the cycle is one.  After the descending
   residue `2`, the other two would therefore have to multiply to
   `2^(-1)=7`.  Yet

   ```text
   {ab mod 13:a,b in {2,3,4}}={3,4,6,8,9,12},            (16)
   ```

   which omits seven.

Both possible cycle lengths are excluded.  We have proved the first decisive
reduction:

```text
hypothetical tightness => min(u_r,u_s,u_t)<=24.          (17)
```

Because every lift is proper, the small replacement in (17) lies in
`{14,...,24}` and has lift height one.

## 4. Lower-runner bounds make the residual finite

Choose `x` to be the least replacement satisfying `x<=24`, and order the
other two replacements as

```text
x<v<w.                                                    (18)
```

The ordering is automatic: a replacement smaller than `x` would itself be at
most `24`, contrary to the choice of `x`.

We use twice the elementary periodic-danger estimate.  For

```text
D_u(gamma)={y:||uy||<=gamma}
```

and an interval `I` of length `L`, partition the scaled interval into
unit-length pieces plus one residual interval.  Periodicity gives

```text
meas(I intersect D_u(gamma)) <= 2gamma L+2gamma/u.       (19)
```

Adjoin `x` to the nine-speed core `P`.  The settled ten-speed LRC bound gives
a point of clearance at least `1/11`; the maximum core speed is at most `24`.
Its `delta`-safe Lipschitz interval has radius

```text
rho_10=(1/11-1/13)/24=1/1716                            (20)
```

and length `L=1/858`.  If the `v`- and `w`-combs covered this interval, (19)
would imply

```text
L <= 4delta L+2delta(1/v+1/w),
1/v+1/w >= 3/572.                                       (21)
```

Since the left side is at most `2/v`,

```text
v <= 1144/3,                 hence v<=381.              (22)
```

Now adjoin `v`.  The resulting eleven-speed core has maximum speed `v` and,
by settled eleven-speed LRC, clearance at least `1/12`.  Its `delta`-safe
interval has length

```text
2(1/12-1/13)/v=1/(78v).                                 (23)
```

The last danger comb must cover this connected interval.  It must therefore
fit inside one tooth of length `2/(13w)`, which forces

```text
1/(78v)<=2/(13w),                 w<=12v.                (24)
```

Equations (17), (22), and (24) give a practical exact box.

## 5. The exact component-containment predicate

For each row in the box, put

```text
Q=P union {x,v},
E_Q={y in R/Z:min_(q in Q)||qy||>1/13}.                 (25)
```

The strict-safe bands of one speed `q` are

```text
((13n+1)/(13q),(13(n+1)-1)/(13q)),       0<=n<q.         (26)
```

Intersecting the eleven finite unions (26) gives every component of `E_Q`
with exact rational endpoints.  If `(l,h)` is one component, let

```text
c=(l+h)/2,                       eta=(h-l)/2.             (27)
```

Because `D_w(delta)` is closed,

```text
(l,h) subset D_w(delta)
  iff [l,h] subset D_w(delta)
  iff ||wc||+w eta <= delta.                             (28)
```

Thus the twelve-speed packet is tight exactly when (28) holds for every
component.  A failure of even one inequality leaves a nonempty open strict
witness interval.  Indeed the packet is a complete nonzero residue transversal
modulo `13`, so every nonzero thirteenth already gives clearance exactly
`1/13`; absence of a strict witness is therefore equivalent to equality.
The verifier clears denominators in (28); its hot loop is integer-only.

## 6. Exhaustion of the finite box

The enumeration is canonical.

1. Choose the missing label triple `R` (`220` choices).
2. Choose its unique least replacement `x<=24`; necessarily
   `x=a+13` with `a in R\{12}`.
3. Assign the other two residue labels to the numerical order `v<w`.
4. Range over every proper `v<=381` and every proper `v<w<=12v`.

The guard `x<v` removes exactly the duplicate presentations in which another
height-one lift should have been the anchor.  An independent arithmetic count
and the literal loop both give

```text
5,713,539 rows.                                           (29)
```

The exact result is

```text
distinct eleven-speed cores        33,730
rows                            5,713,539
rows satisfying every (28)                0.             (30)
```

The SHA-256 digest of every row, its first failed component, and the
denominator-cleared positive surplus is

```text
d85978aceb2f460036b3d2b1919d78381948e458c44fd7e43d23db4a1829479a.
```

As an independent low-height audit, a full pair-crossing/self-cusp maximin
scan of all `220` height-one rows finds minimum

```text
2/21,
```

attained for missing labels `(2,5,6)` and `(4,5,6)`.  There are no finite-box
tight rows, proving Part A.

## 7. Descent back to arbitrary AP scale

Suppose a proper three-coordinate packet around `c[12]` were tight.  THM-804
gives `c|w_r,w_s,w_t`.  Write `w_a=cu_a`.  Since `c` is a unit modulo `13`,

```text
u_a=r_a (mod 13).                                       (31)
```

Positivity and properness leave no negative lift: necessarily
`u_a=r_a+13h_a` with `h_a>=1`.  The packet is `cB` for a set (1).  Circle
multiplication by `c` is onto, so `M(cB)=M(B)>1/13`, contradicting tightness.
Together with THM-795 and THM-800, this proves Part B. ∎

## 8. Tournament Analysis and assumption challenge

There are two exact carriers, used at different recursive depths.

- In the unbounded reduction, vertices are the three **owner-collar exit
  obligations** and the three replacement colours.  The pair observable is
  (12), with speed ratio and the half-open left-germ flag retained.
- In the finite closure, vertices are strict-safe **components** and last-speed
  danger teeth.  The exact incidence is (28), which retains component width.

A runner tournament loses simultaneous owner obligations.  A residue graph
retains `z` but loses `lambda`.  Tooth vertices without owner labels lose the
aligned exit (9).  Fixed circle sections lose the moving tooth scale.  Fourier
modes retain average density but not the connected interval used in (24).
Matroid circuits and bare proof-obligation vertices lose the metric width in
(28).  The proof therefore changes carrier rather than forcing one quotient
to encode both stages.

For telemetry, orient a residue pair provider-to-owner when the subunit
handoff residue is `2`; use increasing residue on silent pairs.  Reverse the
twelve live arrows as the switch.  Both gauges have

```text
live edges                         12
score histogram                    {1:2,3:2,5:2,6:2,8:2,10:2}
directed triangles                 18
SCC sizes                          [12]
Hamiltonian paths                  4,029
edge flips under the switch        12.
```

The strongly connected tournament shadow advertises cyclic residue pressure,
but the theorem is precisely that no cycle survives after the ratio bands and
oriented boundary bit are restored.

THM-806 closes the first unramified three-colour base exposed by THM-804.  It
does not blur the next frontier into a generic "radius at least four": the
first live shallow ramification is the four-colour order-three quartic-coset
interface, which descends to an `s=3` deep packet rather than to another
scale-one triple chart.

## Reproduction

```bash
python3 04-computation/lrc13_scale_one_hamming_three_collar_closure_codex_S10.py
```

The replay uses only Python integers and `fractions.Fraction`.  It checks the
half-open handoff predicate through lift height twelve, the symbolic cycle
alphabets, the analytic constants, all rows in (29), an independent height-one
piecewise-linear maximin census, and the tournament fingerprints.
