# Large-divisor first-flag Newton block

**Status: VERIFIED, not used by THM-3161.**

Exact probe:
`04-computation/factorial_large_divisor_first_flag_block_probe_thm3161.py`
(`d3226bbc352731900ed8fdcd4af0b41abeca4f5671bf7b127cb027d5b1ef8d03`),
with byte-identical normal/optimized transcript
`05-knowledge/results/factorial_large_divisor_first_flag_block_probe_thm3161.out`
(`77d45df581d4966248e52b3c53d1fcba86cffb9be325d8407151991d015e5bcd`).

Let `d=mp`, where `p` is prime and `p>2m`, and let

```text
P=A_(d-2)^(d),  Q=A_(d-1)^(d),  R_1=the THM-3152 first full remainder.
```

The exact tests suggest the complete `p`-adic lower polygons

```text
P:   slope 0       capacity 1,
     slope 1/p     capacity (m-1)p,
     slope 1/(p-3) capacity p-3;

Q:   slope 1/p     capacity (m-1)p,
     slope 1/(p-1) capacity p-1;

R_1: slope 0       capacity 2,
     slope 1/p     capacity (m-1)p,
     slope 1/(p-5) capacity p-5,
```

with the zero-capacity last row edge omitted for `p=5`.  Consequently the
sole common finite block is

```text
(slope,capacity,denominator)=(1/p,(m-1)p,p).
```

There were zero failures in all `68` pairs with prime `p<=43`, `2<=m<=8`,
and `p>2m`, and in the high controls `(d,p,m)=(1002,167,6)` and
`(1384,173,8)`.  The hostile point `(999,37,27)` lies outside `p>2m` and
has common capacity `10p`, not `26p`.

For `P` and `Q`, the displayed polygons are accessible directly.  In the
explicit coefficient sum for `A_n^(d)`, the `k=0` term gives the lower bound

```text
v_p([v^j]A_n^(d))
 >=v_p(binomial(n,j))+floor((n+j)/p),
```

because `p^2>2d`; the base-`p` digits of `n=mp-2` and `mp-1` give exactly
the displayed step hulls, and the decisive vertices have a unique minimal
term.  The analogous step hull for `A_(d-3)^(d)` matches the displayed
polygon of `R_1` in every test.  A proof that the Euclidean combination
`R_1` preserves those decisive valuations is the remaining audit step for
promoting the full flag formula.

There is a useful exact reduction for that remaining step.  With `n=d-2`,
the moment recurrence simplifies the full remainder to

```text
R_1=3(2n+1)d^n
   +n(n+1)(2n-1) A_(n-1)^(d)
   +2dn(n-1)(n+1)(1-4dv) A_(n-2)^(d).
```

For `p>5`, the middle scalar is a `p`-adic unit, while the last row has an
extra factor `p`.  This reduces the promotion audit to comparing the explicit
step hull for `A_(mp-4)^(mp)` with the four decisive vertices of the
`A_(mp-3)^(mp)` hull.  The small admissible cases `(p,m)=(5,2),(7,2),(7,3)`
are already exact finite controls.  This proof route has not yet been audited
as a promoted theorem, so the status of this note remains `VERIFIED`.

Reproduce with

```text
python3 04-computation/factorial_large_divisor_first_flag_block_probe_thm3161.py
python3 -O 04-computation/factorial_large_divisor_first_flag_block_probe_thm3161.py
```
