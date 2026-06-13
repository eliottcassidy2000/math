# Faulhaber Anchors, Square Pyramids, And Bernoulli Addresses

codex-2026-06-13

The cleanest thing from this session is that the midpoint balance is not a
generic high-degree polynomial mess.  Around the midpoint `c=a+n`, every even
Faulhaber moment cancels:

```text
D_p(c,n)=c^p - 2 * sum_{r odd} binom(p,r)c^(p-r)S_r(n).
```

So the power tower is an odd-moment object.  The visible row center is the
triangular carrier `u=n(n+1)`, but the proof-bearing address is the sequence of
odd Faulhaber moments `S_1,S_3,S_5,...`.

That makes the exactness of `p=1,2` feel much less mysterious.  The real root
has

```text
c_p(n)=p*u + alpha_p + beta_p/u + gamma_p/u^2 + ...
```

and the first three corrections all carry `(p-1)(p-2)`.  Linear and square
towers are not merely the first two cases; they are the two cases where the
Bernoulli address is invisible.

The square-pyramid picture is the geometric face of the same fact:

```text
1^2+...+n^2 = n(n+1)(2n+1)/6,
6P_2(n)=n(n+1)(2n+1)=2*S_1(n).
```

So the p=2 tower is simultaneously:

```text
an exact consecutive-square balance,
a square-pyramidal cuboid packing,
and the last integer anchor before Bernoulli denominators appear.
```

The transfer lesson for LRC14 is not "use Faulhaber sums" literally.  It is:
when a scalar balance almost works but only two low levels are exact, look for
the hidden address whose first invisible terms vanish in those levels.  For
the triangular towers that address is odd moments.  For LRC14 it is probably a
ledger of odd walls, owner support, shell-27 class, divisor fiber, carry
residue, and endpoint atom.

The practical target is HYP-2454's bracket:

```text
D_p(p*n(n+1),n)<0<D_p(p*n(n+1)+1,n),  p>=3.
```

The asymptotic says the root sits at

```text
p*n(n+1) + (p-1)(p-2)/(12p) + lower corrections,
```

so the bracket should be true if the lower terms can be bounded before they
push the root across either endpoint.  This is now a concrete inequality task,
not just a pattern.
