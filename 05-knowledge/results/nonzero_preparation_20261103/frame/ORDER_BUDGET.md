# Permanent names impose an asymptotic birth budget

**PROVED from the cited strong-order bound — [independent audit](ORDER_AUDIT.md) passed.** This applies when
every completed output is strong and every previous physical name survives.
It charges no permanent vertices for temporary helpers that are deleted or for old-arc changes.

The **CITED** bound of [Busch, Theorem 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
gives `H(T)>=5^((n-1)/3)` for every strong tournament of order n>=3.
The sharper congruence-class values in Corollary1 are retained below.
The primary theorem and corollary were inspected during this session.

Suppose a physical chain starts at `(n0,H0)` and satisfies

```
H_(i+1)=3H_i+d_i,       d_i in {2,4}.
```

Then `H_i<=3^i(H0+2)-2`. If each step spends two permanent C2 vertices,
and e_i additional permanent preparation vertices have been spent in total,
name retention gives `n_i>=n0+2i+e_i`. Combining these inequalities yields

```
e_i < 1-n0+3 log_5(H0+2) + (3 log_5(3)-2)i.           (1)
```

Consequently any unbounded all-strong rearming scheme of this form must
satisfy

```
limsup e_i/i <= 3 log_5(3)-2 = 0.04781858345... .     (2)
```

The decimal is explanatory; (1) and (2) are exact logarithmic inequalities.
In particular spending one additional permanent preparation vertex at every
step cannot close an unbounded chain. A finite prefix that does so is not
itself a contradiction: the constant budget in (1) matters.

For the actual eleven-vertex H1429 source, four steps that each add three
permanent names would have order at least23 and count at most
`3^4(1429+2)-2=115909`. But the exact minimum at strong order23 is
`9*5^6=140625`. Thus no such four-step all-strong chain exists. This rules
out repeating that birth expenditure, not the six-flip rearming itself:
the new +4 preparation has no permanent birth, and subsequent preparations
may have different costs.

The coordinate exposed here is retained order, alongside H and the changing
root defect. A returning class must control all three. Source padding is a
separate weak-tournament mechanism and does not inherit this strong bound.
