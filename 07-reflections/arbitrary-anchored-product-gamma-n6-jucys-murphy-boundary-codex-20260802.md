# The first new central layer needs three content statistics

**Exact reflection, 2026-08-02.** This note records a positive structural
signal and three sharp stopping results for the arbitrary anchored
product-Gamma width-three problem after THM-3110. It is reproducible exact
evidence, not a theorem and not an all-degree proof.

## 1. Normalization and the central element

For either cleared THM-3110 response bank, delete its exact common root
alphabet and write

```text
Phi(f)=sum_R c_R f(R),
base_1=a^4 b^2(b-a)^3,        base_2=a^3 b^2(b-a)^4.
```

For a cycle type `rho |- N`, put

```text
q_rho=Phi(p_rho)/base_j,
Z_(j,N)=(1/N!) sum_(sigma in S_N) q_type(sigma) sigma.       (1)
```

If `mu |- N`, the scalar by which `(1)` acts on the irreducible `S_N`
module `V_mu` is

```text
z_mu=(1/f^mu) sum_(rho|-N) q_rho chi^mu(rho)/z_rho.          (2)
```

Dual Cauchy identifies the THM-3110 Schur coefficient indexed by `lambda`
with

```text
A_(j,lambda)=base_j f^lambda z_(lambda').                    (3)
```

Thus Schur positivity is exactly positivity of a finite sequence of central
spectra. The conjugation in `(3)` is load-bearing.

## 2. Degree five is the one-generator wall

At `N=5` the only nonzero class quotients are

```text
q_(2,1,1,1)=30,
q_(1^5)=60 L_j,
L_1=3a+5b,                  L_2=3a+4b.                (4)
```

By the Jucys--Murphy content theorem, `(2)` is

```text
z_mu=L_j/2+p_1(C(mu))/4.                               (5)
```

After `mu=lambda'`, `(5)` is exactly THM-3110's
`f^lambda(L_j-kappa_lambda/4)/2` formula. This explains why the first
surviving degree sees only the transposition sum.

## 3. The exact degree-six Farahat--Higman packet

At `N=6`, `q_rho` vanishes for

```text
(6), (5,1), (4,2), (3,3).                              (6)
```

The universal nontrivial class values are

```text
q_(4,1,1)=36,       q_(3,2,1)=-21,       q_(2,2,2)=48,
q_(3,1,1,1)=9(18a+18b-5).                              (7)
```

The identity, transposition, and double-transposition coefficients carry
the chamber and invariant dependence. Hence `Z_(j,6)` lies in the
Farahat--Higman filtration of defect at most three.

Let `C(mu)` be the multiset of contents `column-row` and put
`p_k=p_k(C(mu))`. Exact character inversion gives the unique filtered
content polynomial

```text
z_mu = c_0+c_1 p_1+c_2 p_2+c_11 p_1^2
       +(7/18)p_3-(31/240)p_1p_2+(1/90)p_1^3,          (8)
```

where the first four coefficients depend on the invariant and on whether
`a<b<2a` or `b>=2a`. The homogeneous cubic part is universal:

```text
H_3=(8p_1^3-93p_1p_2+280p_3)/720.                     (9)
```

The exact THM-3110 Newton certificate proves the full `(8)` is positive on
every relevant degree-six representation in both chambers. Formula `(9)`
identifies what that positivity must control; it does not prove it by
itself.

## 4. Three sharp no-go results

### The transposition sum is no longer enough

The partitions `(4,1,1)` and `(3,3)` have the same first content sum
`p_1=3`, but `(8)` separates them by

```text
z_(4,1,1)-z_(3,3)
 = (4a+18b+31)/10       for I1,
 = (4a+11b+31)/10       for I2.                       (10)
```

Therefore `Z_(j,6)` is not a polynomial of the transposition class sum
alone. Any induction that retains only `kappa_lambda=2p_1` has already
forgotten necessary information.

### There is no scalar multiple of the degree-five factor

The degree-five central factor is `L_j/2+p_1/4`. If it divided the universal
degree-six polynomial in the Farahat--Higman content ring, substitution
`p_1=-2L_j` would kill the result. Instead the coefficient of `p_3` remains
exactly `7/18`. Thus the tempting scalar Jucys raising recurrence fails at
the first next degree.

### The universal cubic is not positive

On the three conjugation test shapes,

```text
H_3(6)=295/16,
H_3(3,2,1)=0,
H_3(1^6)=-295/16.                                    (11)
```

So `(9)` is conjugation-odd and is not a positive central operator or a sum
of Hermitian squares. The chamber-dependent lower filtration terms in `(8)`
are load-bearing.

## 5. Holotopy interpretation and the missing sidecar

THM-3110's labelled zeta current lives on rank-four flats of a coloured
partition lattice. Two natural attempts to regard it as a bare cycle fail:

- the unsigned Hasse coboundary `delta W(kappa)=sum_(lambda covers kappa)
  W(lambda)` is nonzero on `500/1050` rank-three flats for `I1` and
  `900/2646` for `I2`;
- after lifting `W(partition(F))` uniformly to lexicographically oriented
  four-edge forests of `K_8` or `K_9`, the ordinary edge-deletion boundary is
  nonzero on `1754/3220` and `3502/7056` three-edge forests.

For the common three-edge forest

```text
G={(0,1),(0,2),(3,4)},                                  (12)
```

the two oriented boundary coefficients are `-1/10` and `+1/10`. In the
first bank the three nonzero extensions are the three edges among the three
`b`-vertices, each contributing `-1/30`; in the second bank the six such
extensions each contribute `+1/60`.

The first failed implication is therefore precise: the component partition
does not remember the rooted tree, broken-circuit presentation, or edge
order that carries an Orlik--Solomon boundary sign. The strongest live
replacement is a rooted NBC/flag lift whose central shadow retains at least
`p_2` and `p_3`. A useful next test is to construct that lift at rank four
and ask whether insertion of one redundant edge becomes an effective
circuit divisor. Repeating a scalar majorization or `p_1`-only induction
cannot succeed.

## 6. Reproduction and scope

Run

```text
python 04-computation/gmc_product_gamma_n6_jucys_murphy_boundary.py
python -O 04-computation/gmc_product_gamma_n6_jucys_murphy_boundary.py
```

The companion derives the class quotients from the actual THM-3110 banks,
reconstructs the `S_6` characters by Murnaghan--Nakayama, verifies the
seven-dimensional content interpolation, and checks `(4)--(11)` exactly.
The forest-boundary counts above were independently enumerated from the
rank-four zeta table; they diagnose the missing lift and are not used as a
proved dependency.

This reflection does not prove all-degree Schur positivity, arbitrary
anchored product-Gamma goodness, GMC(2), NC2, LRC(14), JC(2), or DC(2).
