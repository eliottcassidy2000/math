# Independent quantifier audit of THM-4453

**Verdict: VERIFIED, with the scope stated in THM-4453.** This audit derived
the composition directly from the maintained THM-3818, THM-4153, and THM-4450
statements and replayed both THM-4450 exact implementations. It did not use the
producer's enumeration to establish the five factorizations.

For coprime odd 3-units `p<q`, strict failure with `mu(G_H)>=4/91` forces the
primitive ratio into

```text
(1,11), (1,23), (5,11), (1,37), (1,25).
```

Their sums are respectively

```text
12=2^2*3, 24=2^3*3, 16=2^4, 38=2*19, 26=2*13.
```

Each violates THM-3818's requirement that every prime dividing `p+q` be
`2 mod 3` with exponent at most two. Thus the five-ray and 5,855-ray atlases
are disjoint. Compact containment in a proper open cross-comb makes the
`4/91` boundary inclusive.

The six-coprime hypothesis is necessary: in the larger all-odd lane `(1,9)`
has sum `10=2*5`, is inert-admissible, and survives the `4/91` filter.

The original-body transfers were checked independently:

```text
q=4: mu(G_(2C union {r})) >= mu(G_C)/2,
     mu(G_C)>=8/91 gives 4/91;

q=2 one-even: mu(G_(C union {r})) >= mu(G_C)-8/63,
     4/91+8/63=20/117.
```

Both are genuine sufficient gates only inside the inert-pair class. For an
actual THM-3818 `11+2` consequence, `{a,b}` must be the actual two-vertex
decoder component; the theorem does not close every abstract decoder row.

The incoming unit-core result remains complementary below these mass gates.
For `H=2C union {r}` with odd `r`, its literal normalized unit condition is
equivalent to `r` dividing every member of `C`; its component estimate remains
height-dependent and supplies no universal body-width floor.

LRC(14) remains **OPEN**.
