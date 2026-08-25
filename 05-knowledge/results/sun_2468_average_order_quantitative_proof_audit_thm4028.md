# THM-4028 quantitative lattice-count proof audit

**Status (2026-08-24): PROVED + INDEPENDENTLY PROOF-AUDITED.** This report
checks the error term and fixed-residue extension in
`THM-4028-sun-two-four-six-eight-average-order-criticality.md`. It is a proof
audit, not a numerical extrapolation.

## General lattice error

For fixed positive degrees `d_i` and fixed lower bounds, put

```text
sigma=sum_i 1/d_i,       d_max=max_i d_i.
```

The pure-power body `R_X={u_i>=0: sum u_i^(d_i)/d_i!<=X}` has volume
`V_d X^sigma`. Use disjoint half-open cubes `p+[0,1)^s` based at its integer
points. Monotonicity gives the lower containment `R_X` inside their union.
Across one cube, coordinate `i` changes the defining function by
`O(X^(1-1/d_i))`, so the union lies inside

```text
R_(X+O(X^(1-1/d_max))).
```

Taking the volume difference proves

```text
#(R_X cap Z^s)=V_d X^sigma+O(X^(sigma-1/d_max)).
```

The exact inequalities

```text
(t-d+1)_+^d/d! <= C(t,d) <= t^d/d!
```

squeeze the binomial sublevel set against this count. The translated map is
bijective off coordinate faces; duplicate zero fibres, fixed lower bounds,
and all face strips contribute
`O(sum_i X^(sigma-1/d_i))`, within the same error. The Dirichlet substitution
gives

```text
V_d=product_i((d_i!)^(1/d_i) Gamma(1+1/d_i))/Gamma(1+sigma).
```

For degrees `2,4,6,8`, `sigma=25/24`, `d_max=8`, hence the error exponent is
`25/24-1/8=11/12`.

## Fixed residue classes

Fix `q`. Safe coordinate periods are `P_i=q*d_i!`; the exact periods from
THM-4027 may be used instead. Split the input lattice into the finitely many
product residue classes modulo the `P_i`. The same fixed-fundamental-box
argument tiles a translated orthant and divides the box-union volume by
`product_i P_i`. Only fixed-width coordinate slabs are missed. Away from
zero-coordinate faces, `t_i -> t_i-d_i+1` sends an input coset to a translated
coset with the same leading density. Thus every class has leading term
`V_d X^sigma/product_i P_i` and error
`O_q(X^(sigma-1/d_max))`.

If `C_q(r)` product classes give target residue `r`, then
`sigma_q(r)=q*C_q(r)/product_i P_i`. Summing the accepted cosets proves

```text
sum_(n<=X, n=r mod q) a(n)
  =(sigma_q(r)/q) V X^(25/24)+O_q(X^(11/12)).
```

The implied constant is not uniform for growing `q`.

## Mesoscopic boundary and non-consequences

Writing `J=(25/24)V`, subtraction at `X+H` and `X` is valid when
`H=o(X)` and `H/X^(7/8)->infinity`: the main increment is
`J H X^(1/24)(1+o(1))`, while the inherited error is `O(X^(11/12))`.
The error/main ratio is `O(X^(7/8)/H)`. The fixed-residue statement follows
the same way.

At `H=1` the error is much larger than the formal shell main term. The audit
therefore supports no pointwise asymptotic, Poisson law, density-one claim,
or eventual positivity. THM-4026 remains the exact hostile. **PASS.**
