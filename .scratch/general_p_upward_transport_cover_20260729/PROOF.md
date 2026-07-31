# Upward one-marginal transport is exactly a fractional cover

Status: **PROVED** (finite-dimensional real LP, arbitrary `p`).  The exact
small-`p` program in this directory is a hostile control, not a dependency of
the proof.

## Statement

Let `[p]={1,...,p}` and let `A` be an upward-closed family of subsets of
`[p]`.  Prescribe one-marginals `r=(r_1,...,r_p) in [0,1]^p`, and write

```text
F_A(r) = max_mu mu(A),
```

where the maximum is over all probability laws `mu` on `2^[p]` satisfying
`mu{i in E}=r_i`.  Let `H=min(A)` be the clutter of inclusion-minimal true
sets.  Its fractional-cover polyhedron and weighted cover value are

```text
C(H) = {w in R_{\ge 0}^p : sum_{i in H} w_i >= 1 for every H in H},
tau_H(r) = min_{w in C(H)} sum_i r_i w_i.
```

Then

```text
                 F_A(r) = min(1,tau_H(r)).                 (*)
```

The conventions handle the two constant events: if `A` is empty, `H` is empty
and `tau=0`; if `A=2^[p]`, then `H={empty set}`, `C(H)` is infeasible,
`tau=+infinity`, and the right side is one.

## Upper bound

For every fractional cover `w` and every state `E`,

```text
1_A(E) <= sum_{i in E} w_i.
```

Indeed this is trivial off `A`, while a true `E` contains some minimal true
set `H` and hence has weight at least one.  Taking expectations gives
`mu(A)<=r.w`; also `mu(A)<=1`.  Thus `F_A(r)<=min(1,tau_H(r))`.

## Attainment when `q=tau_H(r)<1`

The dual of the cover LP is the fractional packing LP

```text
max  sum_H lambda_H
s.t. sum_{H containing i} lambda_H <= r_i       (all i),
     lambda_H >= 0.
```

Choose optimal primal and dual solutions `w,lambda`.  Put mass `lambda_H` on
each minimal true state `H`.  This gives an `A`-supported subprobability
measure `nu` of total mass

```text
q=sum_H lambda_H
```

and coordinate loads

```text
a_i=sum_{H containing i} lambda_H <= r_i.
```

Let `s_i=r_i-a_i`, `b=1-q`, and `Z={i:w_i=0}`.  Complementary slackness says

```text
w_i s_i=0,
```

so all residual marginal demand lies on `Z`.  Moreover no true set is
contained in `Z`: otherwise it would contain a minimal true `H subseteq Z`,
contradicting `w(H)>=1`.  Consequently every subset of `Z` lies in `A^c`.

Some `s_i` may exceed the residual mass `b`.  For each such coordinate, set

```text
theta_i=s_i-b.
```

There is enough `nu`-mass currently lacking `i` to add `i` to mass `theta_i`,
because that available mass is `q-a_i` and

```text
theta_i = r_i-a_i-1+q <= q-a_i
```

is exactly `r_i<=1`.  Split atoms if necessary and push the selected mass from
`E` to `E union {i}`.  This does not change total mass, does not affect any
other marginal, and preserves membership in `A` by upward closure.  Processing
the coordinates successively therefore produces an `A`-supported measure
`nu'` of mass `q` whose remaining marginal demands are

```text
s'_i=min(s_i,b).
```

They are still supported on `Z`.  Realize them by a measure `rho` of mass `b`
on `2^Z`: under the normalized probability `rho/b`, include each `i in Z`
independently with probability `s'_i/b` and never include coordinates outside
`Z`.  This is legal because `0<=s'_i<=b`, and every state in its support lies
in `A^c`.  Now `mu=nu'+rho` has the prescribed marginals and
`mu(A)=q=tau_H(r)`.

This argument also covers `q=0`.

## Attainment when `tau_H(r)>=1`

Define the upward blocking polyhedron

```text
Q = conv{1_H : H in H} + R_{\ge 0}^p.
```

We claim `r in Q`.  Otherwise strong separation of `r` from the closed
polyhedron `Q` gives a vector `u` with

```text
u.r < inf_{x in Q} u.x = min_{H in H} u(H).
```

The recession cone of `Q` is `R_{\ge 0}^p`, so any separating vector with a
negative coordinate would make the displayed infimum `-infinity`; hence
`u>=0`.  Since `r>=0`, the strict inequality forces
`m=min_H u(H)>0`.  Then `w=u/m` is a fractional cover with `w.r<1`,
contradicting `tau_H(r)>=1`.  Thus

```text
r >= y=sum_H alpha_H 1_H
```

coordinatewise for some `alpha_H>=0` summing to one.  Put probability
`alpha_H` on the true state `H`.  Its mean is `y`.  For each coordinate `i`,
add `i` to mass `r_i-y_i` of the states currently lacking `i`.  The available
mass is `1-y_i`, and `r_i-y_i<=1-y_i` because `r_i<=1`.  Successive additions
preserve `A` and do not alter other coordinates.  The resulting law is
supported entirely on `A` and has mean `r`, so `F_A(r)=1`.

## Unnormalized fibre form

Suppose a fibre has `t` outer positions and the `p` needles have exact
activity counts `R_i`.  Relax the status-cell counts to reals:

```text
n_E >= 0,
sum_E n_E=t,
sum_{E containing i} n_E=R_i.
```

For an upward heavy-status event `A`, the maximum **real-relaxed** number of
heavy positions is

```text
t F_A(R/t)
  = min(t, min_{w>=0, w(H)>=1} sum_i R_i w_i).             (**)
```

Thus the `2^p`-cell threshold-transport LP can be replaced exactly by a
fractional vertex-cover LP with `p` variables and one constraint for each
minimal heavy status.  For integer status tables, (**) remains a valid upper
bound (and the integer maximum is at most its floor); equality need not be
asserted without a separate integrality argument.

In the LRC application the heavy event is upward whenever the status capacity
`c_E` is monotone under adding active needles.  Its minimal true statuses are
the minimal `E` with `c_E` at least the target threshold.  This monotonicity is
the exact hypothesis that makes the cover formula sharp.

