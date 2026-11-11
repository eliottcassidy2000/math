# The unbounded prefix-rewrite family is still sparse

**PROVED, scoped to the family below.** The new neutral rewrite permits
unbounded prefix length, repair size and root response. Nevertheless its
entire scalar image, allowing every eligible suffix, has only O(sqrt X)
values up to X. This does not obstruct compositions with additional graph
operations, other word languages, or Chan–Pak in general.

## The actual family being counted

The [Berggren return note](../berggren/PROOF.md) constructs the following
generalization of its three-arc repair. For j>=3, t>=1, let z be any
strong one-apex cut word satisfying

```
H(z)=h_j=3*2^(j-2)-1.
```

The source word and its neutral replacement after the usual two births are

```
w  = O^(j+1) I^t z,
w' = O I O^(j-1) I O I^t z.
```

The replacement reverses exactly j arcs incident to the newest vertex,
retains all names, and has H(w')=3H(w)+2. In particular the class is
nonempty at every j: `z=OO I^(j-2)` has the required h_j. The theorem
below counts ALL eligible z, not only that example or the two original
H5 suffixes.

## Two elementary bounds retain the right scale

For a strong cut word z define
`Q(z)=sum_(i:z_i=I)2^i`, with positions starting at0. The exact one-apex
count law gives

```
H(z)-1 <= Q(z) <= 2H(z)-4.                              (1)
```

For the lower bound on Q, each appended I contributes P+e<=2^i, while
appended O contributes nothing to H. Summing from initial H1 gives
H<=1+Q. For the upper bound, the first O paired with every subsequent I
contributes exactly Q/2 to the positive cut formula, and the two endpoint
paths contribute2. No cancellation is used.

The exact source count is

```
H(w)=(2^(j+1)-1)(Q(z)+1)2^t+1-5*2^(j-2).                (2)
```

Combining Q+1>=h_j with j>=3,t>=1 gives

```
H(w) >= 2^(2j+t).                                      (3)
```

For completeness, after subtracting the right side the lower expression
is `2^j[(2^(j-1)-11/4)2^t-5/4]+2^t+1`, strictly positive on this domain.

A strong word of length l has H>=2+2^(l-2), so H(z)=h_j forces l<=j+1.
There are at most `sum_(l=2)^(j+1)2^(l-2)=2^j-1` strong words at these
lengths. This counts pointed words; it is an upper bound even if several
words define the same ordinary graph or H value.

## Summing the complete family

Put L=floor(log2 X). Every source of count at most X must satisfy
2j+t<=L. Therefore its number of presentations is at most

```
sum_(j=3)^floor((L-1)/2) (2^j-1)(L-2j).                 (4)
```

If m=floor((L-1)/2), write L-2m=epsilon in {1,2} and j=m-r. The upper
bound obtained by dropping -1 and extending the sum is

```
2^m sum_(r>=0) (epsilon+2r)/2^r
  =2^m(2epsilon+4) <=8*2^m = O(sqrt X).
```

The target counts are just3H+2, so their counting function satisfies the
same bound. Adding the mirrored one-step continuation only introduces
another fixed affine image, also O(sqrt X).

This is a quantitative boundary on a successful unbounded construction:
escaping a fixed-defect theorem is necessary for some intended programs,
but it does not by itself supply positive density or exact target coverage.
The preserved variable is numerical H, not operation depth or number of
unrestricted binary words. The next useful change must add another source
mechanism or a proved compositional closure, rather than merely letting
j and t grow in this same repair family.

The proof is elementary and uses no fitted exponent. Independent literal
wide-family controls are recorded with the Berggren construction. The
inequalities (1) can also be checked directly on all511 cut words in that
lane's existing positive/hostile test universe, restricting to strong words.
