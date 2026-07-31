# Exact Faber-to-Stieltjes denominator bridge

**Status:** proved algebraic addendum for audit; not canon.

Use the nonsplit response notation

```text
U^2=V,                 R_Q=U G,
q=A/U,                 q^2=T=A^2/V,
K=R_Q/q,               F_resp=R_Q^2=T K^2.
```

Then the Faber quotient already contains the balanced Stieltjes denominator:

```text
A K
 =A(R_Q/q)
 =U R_Q
 =V G.                                                     (1)
```

For the balanced factor data of the Stieltjes--Padé theorem,

```text
V=v S D T_pole^2,      G=g E/(D T_pole),
lambda=vg,             M=S E T_pole,
```

so `(1)` becomes

```text
A K=lambda M.                                             (2)
```

In particular, `A K` must be a nonzero scalar multiple of a squarefree
polynomial of degree `r+1`; its roots are exactly the simple-zero,
double-zero, and pole locations of `F_resp`.  The residues of
`d log(F_resp)` at those roots are respectively `+1,+2,-p_j`.

THM-2245's Keller one-form

```text
A(2T K'+K T')=2 kappa T
```

therefore gives, without introducing an auxiliary Padé denominator,

```text
(F_resp)'/F_resp
 =(T'K+2TK')/(TK)
 =2 kappa/(A K)
 =C/M,

C=2 kappa/lambda.                                        (3)
```

The square-contact identity also becomes

```text
(A K)^2=V F_resp=lambda^2 M^2,

V-vM^2
 =V(1-F_resp/mu),
mu=vg^2,                                                (4)
```

and hence the balanced theorem forces

```text
deg(V-vM^2)=r+2,
LC(V-vM^2)=2 kappa v/(r lambda) !=0.                    (5)
```

Thus the cheapest degree-fourteen continuation test is:

1. form the explicit THM-2245 quotient `K=R_Q/q`;
2. require `A K` to be polynomial, squarefree, and of the required degree;
3. read the residue alphabet from `d log(TK^2)`;
4. impose the exact nonzero first square defect `(5)`;
5. only then solve the remaining source-polynomial and exact-prefix sidecars.

This is an identity internal to the inherited quartic chart.  It does not
derive that chart from an arbitrary Keller pair, and `T_pole` above is the
monic polynomial of finite pole locations, not THM-2245's coefficient
`T=q^2`.
