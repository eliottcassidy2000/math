---
id: THM-3618
title: "Compiler one-graph-observable fibre separator and no embedding"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
  AUDIT.  The polynomial observables separating every closed fibre of the
  THM-3561 compiler are exactly r+lambda e^(m-1)xq with r in C[b,c,e],
  lambda nonzero, and m>=1.  None generates C[x,q] over C[b,c,e], while the
  two observables xq and x do.  This result is not yet proved canon.
source: root / audit_thm3611 off-diagonal-conic wildcard, 2026-08-21
audit: PENDING -- provisional package frozen for independent hostile audit
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3614-russell-cylinder-collision-free-full-linear-projection-rigidity
script: 04-computation/jc2_compiler_one_graph_observable_separator_no_embedding_thm3618.py
output: 05-knowledge/results/jc2_compiler_one_graph_observable_separator_no_embedding_thm3618.out
script_sha256: 2dfb92206451b68033e0cc1404e4320709ed866334d78f3450eb9138db69385f
output_sha256: 69ea7d635abdce6708e331aabe542e3b3dc602502792b400f1eaec497017741a
hash_basis: raw LF bytes
---

# THM-3618 -- compiler one-graph-observable fibre separator and no embedding

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  This file is not proved canon.  No proved theorem may depend on it
until the pending audit promotes its status.

All rings are over `C`.  Put

```text
D=1+x^2q,
b=(D-1)(D+2)^2,
c=xD(D+2),
e=q(D+3),
z=xq,
A=C[x,q],                    R=C[b,c,e].               (1)
```

Say that `h in A` is a **fibre separator** if

```text
(b,c,e,h)(x_1,q_1)=(b,c,e,h)(x_2,q_2)
                 implies (x_1,q_1)=(x_2,q_2)          (2)
```

for all closed points of `A2`.

## 0. Statement

The fibre separators are exactly

```text
h=r+lambda e^(m-1)z,
r in R,              lambda in C*,        m>=1.        (3)
```

For every separator in `(3)`,

```text
R[h] is strictly contained in A.                             (4)
```

Consequently no single polynomial graph observable makes
`(b,c,e,h):A2 -> A4` a closed embedding, even though the observables in `(3)`
make it injective on closed points.  Two sidecars are sufficient and give the
sharp positive boundary

```text
e=4q+z^2,             q=(e-z^2)/4,
R[z,x]=C[x,q].                                           (5)
```

The classification is literal: `m` is an integer at least one, `lambda` is
nonzero, and the only freedom suppressed by the fibres is addition of an
element of `R`.

## 1. The normal compiler image and its puncture

Write

```text
f(T)=(T-1)(T+2)^2=T^3+3T^2-4,
g(T)=T(T+2).                                             (6)
```

The compiler satisfies

```text
b=f(D),       c=xg(D),       b+4=D^2(D+3),
c^2e=b(b+4).                                             (7)
```

Thus

```text
R = C[B,C,E]/(C^2E-B(B+4)).                              (8)
```

The displayed hypersurface `Y` is smooth.  Indeed, for
`F=C^2E-B(B+4)`, a singular point would have `C=0` and `B=-2`, whereas
`F(-2,0,E)=4`.  In particular `R` is normal.

The image of `Phi=(b,c,e)` is exactly

```text
Y minus {p_infinity},              p_infinity=(-4,0,0). (9)
```

Here is the full fibre inventory.  If `c!=0`, choose a root `t` of
`f(t)=b` with `g(t)!=0`; then

```text
x=c/g(t),
q=(t-1)g(t)^2/c^2.                                     (10)
```

For `b=0` only `t=1` is usable, and for `b=-4` only `t=-3` is usable.
For every other `b`, all three roots are usable, counted without multiplicity
on the generic locus.  When `c=0`, the fibres are

```text
(b,e)=(0,e):
  x=0, q=e/4, and, if e!=0, D=-2, q=e, x^2=-3/e;

(b,e)=(-4,e), e!=0:
  D=0, q=e/3, x^2=-3/e;

(b,e)=(-4,0): no point.                                (11)
```

This proves `(9)` and also shows that `Phi` is quasi-finite.  It gives the
intersection fact used below:

```text
A intersect Frac(R)=R.                                  (12)
```

For completeness, if an element of `Frac(R)` has a pole along a prime divisor
of the normal surface `Y`, that divisor meets the image because the complement
in `(9)` has codimension two.  Quasi-finiteness pulls the pole back to a pole
on `A2`, so the element cannot lie in `A`.  Normality, as the intersection of
the height-one valuation rings, proves `(12)`.

## 2. The off-diagonal conic and the one-weight gate

For two distinct roots `t,s` of the same equation `f(T)=b`, cancellation of
`t-s` gives the irreducible off-diagonal conic

```text
t^2+ts+s^2+3(t+s)=0.                                   (13)
```

Put `p=t+s`.  On this conic,

```text
ts=p^2+3p,
(t-s)^2=-3p(p+4),
g(t)g(s)=p(p+1)(p+3)(p+4).                             (14)
```

The generic usable-pair locus therefore removes precisely
`p=0,-1,-3,-4`.

Give `A` the grading

```text
wt(x)=1,                 wt(q)=-2.                     (15)
```

Then `b` has weight zero, `c` weight one, and `e` weight minus two.
Write `h=sum_w h_w`.  On `c!=0`, every component has the form

```text
h_w=c^w R_w(t).                                         (16)
```

For a fixed generic pair `(t,s)`, its difference is a Laurent polynomial in
the free coordinate `c`.  A Laurent polynomial over `C` with at least two
nonzero terms has a root in `C*`.  Since the conic `(13)` is irreducible, if
two weights are active there is a generic pair on which both corresponding
differences are nonzero; choosing that root of `c` contradicts separation.
Hence exactly one weight can be active.  Every inactive component is constant
on the generic fibres, belongs to `Frac(R)`, and hence belongs to `R` by
`(12)`.

An even active weight cannot distinguish the sign pairs in `(11)`, because
its `x`-parity is even.  The unique active weight must therefore be odd.

## 3. Positive odd weights do not occur

Let the active weight be `w=2a+1>0`.  There is a polynomial `P` such that

```text
h_w=x^wP(D),              R_w(t)=P(t)/g(t)^w.           (17)
```

On `(13)`, define the symmetric polynomial `Q(p)` by

```text
P(t)g(s)^w-P(s)g(t)^w=(t-s)Q(p).                       (18)
```

Generic separation says that every root of `Q` lies among the four deleted
values in `(14)`.  Separation of the `D=0` and `D=-2` sign pairs gives
`P(0)P(-2)!=0`.  Expanding `(18)` at the two diagonal boundary points gives

```text
ord_(p=0) Q=ord_(p=-4) Q=a.                            (19)
```

At the two mixed points, up to the harmless orientation convention,

```text
Q(-3)=3^(w-1)P(0),
Q(-1)=-3^(w-1)P(-2).                                   (20)
```

Thus

```text
Q(p)=kappa [p(p+4)]^a,               kappa!=0.         (21)
```

This forced divided difference is impossible.  Let `t_1,t_2,t_3` be the
generic roots of `f(T)=b`, let `g_i=g(t_i)`, and orient the three pairs
cyclically.  Since `f'(t_i)=3g_i`,

```text
(t_1-t_2)/(g_1g_2)
=(t_2-t_3)/(g_2g_3)
=(t_3-t_1)/(g_3g_1) !=0.                               (22)
```

Equations `(14),(18),(21)` say that the three cyclic differences of
`P(t_i)/g_i^w` are the same nonzero scalar multiple of the `w`-th power of
the common value in `(22)`.  Their telescoping sum is both zero and three
times a nonzero number, a contradiction.  No positive odd active weight is
possible.

## 4. Negative odd weights give exactly one family

Write the active weight as

```text
w=1-2m,                         m>=1.                  (23)
```

Then

```text
h_w=xq^mP(D).                                           (24)
```

On the dense locus `e!=0`, division by the common fibre value `e^m` changes
no collision and gives

```text
h_w/e^m=c P(t)/ell_m(t),
ell_m(t)=g(t)(t+3)^m.                                  (25)
```

Define `Q(p)` on `(13)` by

```text
P(t)ell_m(s)-P(s)ell_m(t)=(t-s)Q(p).                   (26)
```

Again all roots of `Q` lie among the four deleted values.  The two sign
fibres and the actual mixed `b=0` fibre force `Q` to be nonzero at
`p=0,-4,-1`.  Consequently

```text
Q(p)=kappa(p+3)^d                                      (27)
```

for some `d>=0` and `kappa!=0`.

For the three generic roots, if `{i,j,k}={1,2,3}`, then

```text
(t_i+3)(t_j+3)=t_k^2,
t_i+t_j+3=-t_k.                                        (28)
```

The telescoping sum of the three cyclic differences in `(26)`, after the
common nonzero factor from `(22)` is removed, is

```text
sum_(k=1)^3 t_k^(d-2m)=0.                              (29)
```

Among integer powers, this identity holds identically only for exponent
`-1`.  Indeed, `sum 1/t_k=0` is the vanishing second elementary symmetric
function of the roots.  Exponent zero gives `3`; positive exponents fail at
`b=-4`, whose roots are `0,0,-3`; and exponent `-n<=-2` fails at `b=0`, whose
roots are `1,-2,-2`, because

```text
1+2(-2)^(-n) !=0.                                      (30)
```

Thus `d=2m-1`.  There is an exact representative:

```text
P_0(t)=(t+3)^(m-1),
Q_0(p)=-2(p+3)^(2m-1).                                 (31)
```

After scaling `(31)` to match `(27)`, the remainder in `(26)` has zero
divided difference.  It descends to `Frac(R)` and, being polynomial in
`x,q`, lies in `R` by `(12)`.  Finally

```text
xq^m(D+3)^(m-1)=z e^(m-1).                             (32)
```

Together with the one-weight gate, this proves the necessity of `(3)`.

## 5. Converse and fibre-by-fibre sharpness

It remains to check that every observable in `(3)` separates.  Addition of
`r in R` and multiplication by `lambda!=0` do not affect a fibre.  On a
generic fibre with `e!=0`, `(26),(31)` show that `z` distinguishes every
usable pair, so `e^(m-1)z` does as well.

The boundary fibres in `(11)` are equally explicit:

```text
(b,c,e)=(0,0,e), e!=0:
  central point z=0;
  D=-2 pair z^2=-3e, with opposite nonzero signs;

(b,c,e)=(-4,0,e), e!=0:
  D=0 pair z^2=-e/3, with opposite nonzero signs.       (33)
```

Every image fibre with `e=0` is a singleton.  This proves the converse and
the full equivalence `(3)`.

Notice that `z` itself, the case `m=1`, is already a global fibre separator.
The extra powers of `e` are invisible to this property because all
non-singleton fibres have `e!=0`.

## 6. Why one observable never closes the graph

Restrict to the punctured `D=0` curve using `u in C*`:

```text
x=u^(-1),              q=-u^2,
(b,c,e,z)=(-4,0,-3u^2,-u).                             (34)
```

For every separator `(3)`, all generators of `R[h]` restrict to `C[u]`,
whereas `x` restricts to `u^(-1)`.  Therefore `x` is not in `R[h]`, proving
`(4)`.  Geometrically, the graph has a missing limit above
`p_infinity=(-4,0,0)` at `u=0`; pointwise injectivity has not paid the
properness/closure debt.

Two sidecars do pay it.  Since

```text
e=q(D+3)=q(4+x^2q)=4q+z^2,                             (35)
```

`q=(e-z^2)/4` lies in `R[z]`, and adjoining `x` gives `(5)`.

## 7. Connection to THM-3614 and scope ledger

The sole degree-cancellation graph in THM-3614 is

```text
h=-xq+n=-z+n.                                          (36)
```

It is exactly the separator class `(3)` with `m=1`, `lambda=-1`, and
`r=n`.  Thus its stable graph is collision-free but nonclosed.  The residual
coefficient `[q]Jac=8` found there records nonproper-closure debt, not a
surviving fibre collision.

What is preserved:

* the exact closed-fibre equivalence relation of `Phi`;
* the complete one-observable separator classification;
* the codimension-two missing point and its `D=0` valuation witness;
* the sharp two-observable recovery identity.

What is not claimed:

* that an injective polynomial map is automatically closed or proper;
* that every nonseparator has only a finite or prescribed collision set;
* a classification of arbitrary multi-observable generating sets;
* any Jacobian-conjecture counterexample or reduction.

The exact companion checks the compiler and surface identities, complete
special-fibre controls, off-diagonal-conic algebra, weight normal forms,
positive cyclic obstruction, the negative survivor family, integer-power
hostiles, the puncture valuation, and the two-sidecar recovery.  The
all-degree descent and Laurent-unit steps above are proof-driven.
