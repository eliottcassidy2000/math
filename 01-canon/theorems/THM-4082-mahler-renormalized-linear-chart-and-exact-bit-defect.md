---
id: THM-4082
title: "Mahler renormalized linear chart and exact bit defect"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. After the odd-unit
  rescaling x=3^s t, the THM-4077 denominator-19 isometries converge uniformly
  to one linear odd-unit isometry. Their error has the sharp valuation
  s+2+2v_2(x), so the nonlinear and linear carry words have exactly that many
  common initial bits. For every pair of scales s<t, their exact pairwise
  defect is governed by the shallower scale s, so each nonzero point converges
  without ever stabilizing. The induced near-identity isometry transports the
  strict-safe and output-termination fibres exactly, while the separately
  scaled parameter-termination loci form a decreasing dense countable
  filtration with empty intersection. Neither open strict-safe termination
  intersection is decided, and no Mahler Z-number is produced or excluded.
source: codex-frontier-synthesis-creative-20260825e / Mahler renormalization niche
audit: >
  PASS. The primary affine/geometric replay uses 48-bit exact arithmetic for
  0<=s<=12 and 1<=x<=4095: 53,235 sharp-valuation gates, 53,235 exact
  carry-first-divergence gates, 16,382 safe-prefix transport gates, 143
  inverse-fibre gates, 13 zero fibres, 13 cross-scale secants, and 50 exact
  vertical divisibility controls. A genuinely independent 72-bit
  logarithm/exponential replay constructs both series term by term, checks
  1,768 analytic valuations, 1,768 original-map comparisons, 1,768 exact
  exponential remainders, and independently lifts the same 143 fibres at 48
  bits. The audit suite includes s=0, s=1, x=0, even-x, safe/nonoutput, and
  output/unsafe hostiles. Normal and optimized outputs byte-match; both sources have zero
  assert nodes and zero floating literals.
depends_on:
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
  - THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction
  - THM-4077-mahler-denominator19-2adic-tangent-full-shift-isometry
related:
  - THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation
  - THM-4074-mahler-denominator19-postterminal-arbitrary-delay
script: 04-computation/mahler_renormalized_linear_chart_thm4082.py
output: 05-knowledge/results/mahler_renormalized_linear_chart_thm4082.out
script_sha256: 68dbe901a10ca6673cbe3ab9ce0c6582fa6164d7a56aa6e923e72eadedda451e
output_sha256: d7634b499ebabc9b87ff536641494e1ff0efd3456cea4e04c271f696ae4fc8ec
independent_audit_script: 04-computation/mahler_renormalized_linear_chart_thm4082_independent_audit.py
independent_audit_output: 05-knowledge/results/mahler_renormalized_linear_chart_thm4082_independent_audit.out
independent_audit_script_sha256: 03e696bd762f17a972102d593ba915568cc902d228cd9fe55f0a9abafcd512df
independent_audit_output_sha256: a51b2f50388ec6fb7144c4e5cb9837c740e10de55f2c45f48588abbc8126a2ff
hash_basis: raw working-tree bytes (LF)
---

# THM-4082 -- the denominator-19 charts have an exact linear blow-up

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** THM-4077
constructs one onto 2-adic isometry at every reset depth. The present theorem
compares those different scales. One odd-unit change of parameter exposes a
common linear limit, with an exact rather than asymptotic first-error bit. The
same comparison also keeps parameter termination and output-state termination
as different typed loci; it does not close either strict-safe intersection.

## 1. The renormalized analytic normal form

Retain exactly the notation of THM-4077. For `s>=0`, put

```text
L_s=18*2^s,                 k_s=s+3,
g_s=3^(L_s),
U_s=9*3^(k_s)(g_s-1)/(19*2^(k_s)),
Fhat_s(t)=U_s*(g_s^t-1)/(g_s-1).                       (1)
```

Here `Fhat_s:Z_2 -> Z_2` is the onto isometry proved in THM-4077. Let `log`
denote the power-series logarithm on `1+4Z_2`, and set

```text
ell=log(3^18) in 8Z_2,
Lambda=(243/152)ell in Z_2^x.                           (2)
```

The normalization at the apparently exceptional unit `3` is

```text
log_2(3):=log(-3),              ell=18log_2(3),
Lambda=(2187/76)log_2(3).                              (3)
```

Thus no use is made of the power series at `3` itself. Since
`v_2(3^18-1)=3`, the standard logarithm valuation on `1+4Z_2` gives
`v_2(ell)=3`; consequently `Lambda` is an odd unit.

Rescale the parameter by the odd unit `3^s`:

```text
H_s(x)=Fhat_s(3^(-s)x),
z_s(x)=(2/3)^s ell*x,
C_s=9*3^(s+3)/(19*2^(s+3)).                            (4)
```

Because `log(g_s)=2^s ell`, the convergent 2-adic exponential gives the exact
normal form

```text
boxed: H_s(x)=C_s(exp(z_s(x))-1),
       C_s z_s(x)=Lambda*x.                            (5)
```

Both sides of (5) are maps `Z_2 -> Z_2`; the displayed coefficient `C_s`
alone need not be integral. The divisibility in the exponential difference is
load-bearing.

## 2. Sharp nonlinear defect

For every nonzero `x in Z_2`, write `r=v_2(x)`. Then

```text
boxed: v_2(H_s(x)-Lambda*x)=s+2+2r.                    (6)
```

At `x=0`, both maps are zero, so (6) is replaced by exact equality rather
than by a finite valuation.

To prove (6), first note from (2)--(4) that

```text
a:=v_2(z_s(x))=s+3+r>=3.                               (7)
```

For any nonzero `z` with `v_2(z)=a>=3`, the quadratic term in the exponential
remainder is uniquely minimal:

```text
v_2(z^2/2)=2a-1,
v_2(z^n/n!)=na-v_2(n!)
            >=na-(n-1)
            =(2a-1)+(n-2)(a-1)>2a-1       (n>=3).     (8)
```

Hence the infinite tail cannot cancel the quadratic term and

```text
v_2(exp(z)-1-z)=2v_2(z)-1.                             (9)
```

Finally `v_2(C_s)=-(s+3)`, so (5), (7), and (9) give (6). In particular,
`H_s` converges uniformly on `Z_2` to the linear isometry

```text
L_infinity(x)=Lambda*x.                                (10)
```

The equality and even-input boundaries are sharp. The smallest hostile
controls are

```text
                 x=1   x=2
s=0 defect bits    2     4
s=1 defect bits    3     5.                            (11)
```

Thus the odd-unit law `s+2` must not be silently extended unchanged to even
`x`; every extra factor of two in `x` contributes two additional common bits.

### Pairwise cross-scale separation

The common limit makes every pair of finite scales exactly comparable. For
every `x!=0` in `Z_2` and integers `0<=s<t`,

```text
boxed: v_2(H_t(x)-H_s(x))=s+2+2v_2(x).                (11a)
```

Indeed, put `E_j(x)=H_j(x)-Lambda*x`. Equation `(6)` gives

```text
v_2(E_s(x))=s+2+2v_2(x)<t+2+2v_2(x)=v_2(E_t(x)).
```

The two valuations are unequal, so the strong triangle equality makes the
shallower error survive in `E_t-E_s=H_t-H_s`, proving `(11a)`. At `x=0`, all
maps vanish and the finite valuation formula is replaced by exact equality.

Consequently, for fixed nonzero `x`, the sequence `(H_s(x))` is Cauchy but
never eventually constant. Its pairwise ultrametric distance is determined
exactly by the shallower scale. If

```text
R_(s,t)=H_s^(-1) o H_t,                               (11b)
```

then the isometry of `H_s` also gives

```text
v_2(R_(s,t)(x)-x)=s+2+2v_2(x).
```

This is pointwise chart separation. It says nothing about preservation of
ordinary integers, positivity, Mahler time evolution, or either termination
fibre.

## 3. Exact carry agreement and the near-identity transport

Let `Phi:{0,1}^N -> Z_2` be the THM-2228 carry homeomorphism. THM-4077 recalls
that `Phi` is an isometry: two states differ by exact valuation `d` if and
only if their carry words agree at positions `0,...,d-1` and differ at
position `d`. Therefore, for every nonzero `x`,

```text
Phi^(-1)(H_s(x)) and Phi^(-1)(L_infinity(x))
agree for exactly s+2+2v_2(x) initial digits,           (12)
```

and their next digits are opposite. At `x=0` the two entire words coincide.
This is an exact first-divergence statement, not merely a congruence lower
bound.

Likewise, `(11a)` says that for every `0<=s<t`, the carry words
`Phi^(-1)(H_s(x))` and `Phi^(-1)(H_t(x))` agree for exactly
`s+2+2v_2(x)` digits and differ at the next. Thus passing to a deeper chart
cannot repair the first finite-scale defect; it only postpones its own new
defect beyond that position.

Both `H_s` and `L_infinity` are onto isometries: the first is the composition
of THM-4077's `Fhat_s` with multiplication by `3^(-s)`, and the second is
multiplication by the odd unit `Lambda`. Consequently

```text
R_s=H_s^(-1) o L_infinity : Z_2 -> Z_2                 (13)
```

is an onto isometry fixing zero. Applying the isometry `H_s` to a displacement
in (13) and using (6) gives the second sharp law

```text
boxed: v_2(R_s(x)-x)=s+2+2v_2(x)       (x!=0).         (14)
```

Thus `R_s -> id` uniformly. Nothing here says that `R_s` preserves ordinary
integers, positivity, addition, the binary shift, or Mahler time evolution.

## 4. The two termination fibres remain different under transport

Let `S` be the THM-4072 strict safe carry set and `[1]` the words beginning
in one. Define the odd state locus

```text
Y_safe=Phi(S intersect [1]) subset Z_2^odd.             (15)
```

In the original THM-4077 parameter coordinate, write

```text
Acal_s=Fhat_s^(-1)(Y_safe),
G_s=Fhat_s^(-1)(N_odd),       N_odd={1,3,5,...}.        (16)
```

Thus `Acal_s` is THM-4077's strict-safe parameter locus, whereas `G_s` is its
**output-state termination** locus. In the limiting linear chart define

```text
Acal_infinity=L_infinity^(-1)(Y_safe),
G_infinity=L_infinity^(-1)(N_odd).                     (17)
```

The type-correct scaling identity is

```text
H_s(3^s t)=Fhat_s(t).                                  (18)
```

Combining (13), (17), and (18) gives the exact pullbacks

```text
boxed: 3^s Acal_s=R_s(Acal_infinity),
       3^s G_s=R_s(G_infinity),
       3^s(Acal_s intersect G_s)
          =R_s(Acal_infinity intersect G_infinity).    (19)
```

The last line is precisely the output-state/Z-number gate of THM-4077 after
rescaling. It remains **OPEN**.

Parameter termination is a different intersection. Put

```text
P_s=3^s N_odd.                                         (20)
```

Then

```text
boxed: 3^s(Acal_s intersect N_odd)
       =R_s(Acal_infinity) intersect P_s.              (21)
```

This is the fixed denominator-19 family gate, also **OPEN**. Replacing `P_s`
in (21) by `R_s(G_infinity)` would change parameter termination into output
termination and would be a type error.

## 5. Arithmetic, topology, and Archimedean sparsity

The parameter fibres in (20) form a strict filtration

```text
P_(s+1) proper subset P_s,
intersection_(s>=0) P_s=empty.                         (22)
```

Indeed, membership in every `P_s` would make one positive integer divisible
by every power of three. Each `P_s`, however, is countable and dense in the
odd 2-adics, hence Haar-null and meagre there. Thus (22) is not a compactness
argument: none of the `P_s` is closed.

At a fixed scale define the rescaled output fibre

```text
Q_s=3^s G_s=H_s^(-1)(N_odd)=R_s(G_infinity).           (23)
```

THM-4077's inclusion `N_odd subset G_s` gives

```text
P_s subset Q_s.                                        (24)
```

The inclusion is proper, and the three possible classes

```text
P_s,                    parameter and output terminating;
Q_s setminus P_s,       output terminating only;
Z_2^odd setminus Q_s,   neither terminating            (25)
```

are all dense in `Z_2^odd`. The first two are countable, Haar-null, and
meagre; the third is comeagre and has full Haar measure. There is no
parameter-only class because of (24).

Here is the only nonformal density point. Let

```text
E_s=Fhat_s(N_odd)=F_s(N_odd) subset N_odd,              (26)
```

where `F_s` denotes (1) on ordinary nonnegative integers. Since `F_s` is
strictly increasing, for every real `X>=0` put

```text
B_s(X)=floor(log_(g_s)(1+(g_s-1)X/U_s)).               (27)
```

Then the exact height count, followed by its `X -> infinity` asymptotic, is

```text
#(E_s intersect [1,X])=floor((B_s(X)+1)/2)
                      =log X/(2log g_s)+O_s(1).        (28)
```

On the other hand, positive odds in each fixed odd residue class modulo
`2^h` have linear growth. Hence every such class contains a positive odd
outside `E_s`. Since `H_s` maps odd cylinders isometrically to odd cylinders
and

```text
H_s(P_s)=E_s,                                          (29)
```

this proves density of `Q_s setminus P_s`. Density of the complement of
`Q_s` follows because `Q_s` is countable and every 2-adic cylinder is
uncountable. The apparent tension in (28) is exact: modulo `2^h`, the finite
permutation inherited from THM-4077 sends successive odd parameter classes
bijectively to the odd output classes, so their counts up to a parameter
cutoff differ by at most one. `E_s` is simultaneously Archimedean-sparse and
2-adically dense.

## 6. Crossed hostile fibres survive the limit

The two canonical THM-4077 hostile states are

```text
y_*=-9/19,          Phi^(-1)(y_*)=(100)^infinity in S,
y_1=1,              Phi^(-1)(y_1) begins 1011 and is unsafe.   (30)
```

The first is not in `N_odd`; the second is. For nonzero `y in Z_2`, let

```text
x_(s,y)=H_s^(-1)(y),        x_(infinity,y)=Lambda^(-1)y.       (31)
```

Equation (14), applied to `x_(infinity,y)`, yields the exact inverse-fibre
law

```text
v_2(x_(s,y)-x_(infinity,y))=s+2+2v_2(y).               (32)
```

In particular both odd hostile fibres converge with exact valuation `s+2`,
yet at every scale

```text
x_(s,y_*) in 3^s Acal_s setminus 3^s G_s,
x_(s,y_1) in 3^s G_s setminus 3^s Acal_s,              (33)
v_2(x_(s,y_1)-x_(s,y_*))=v_2(1+9/19)=v_2(28/19)=2.    (34)
```

Thus uniform convergence to the linear chart neither merges the two gates
nor supplies a termination theorem. Similarly, agreement with the greedy
equality word through any finite height is only membership in the finite
safe-prefix language; THM-4072 shows that the infinite equality ray is in
the closed boundary `K setminus S`, not in the strict safe set.

## 7. An orthogonal vertical divisibility sidecar

For every positive odd ordinary parameter `t`, direct use of (1) and
`g_(s+1)=g_s^2` gives

```text
boxed: F_(s+1)(t)=3(g_s^t+1)/2 * F_s(t).               (35)
```

Because `g_s=1 mod 4`, the multiplier in (35) is an odd positive integer.
Thus the ordinary outputs at fixed `t` form a divisibility chain across
scales. This arithmetic arrow is not a carry conjugacy, does not preserve the
strict-safe predicate, and gives no inclusion between the open intersections
in (19) and (21).

## 8. Exact replays and scope

The primary replay avoids logarithm and exponential series. It builds the
maps from affine/geometric recurrences, recovers `Lambda` from stabilized
cross-scale secants, checks (6), (12), finite safe-prefix transport, (32), and
(35), and includes zero, even, and crossed-fibre controls.

The independent replay constructs `log(3^18)`, `log(-3)`, and `exp(z)` term
by term with dyadic denominator clearing. It checks the normalization
`log(3^18)=18log(-3)`, the unique quadratic valuation in (8)--(9), equality
with the original geometric map, and independent Hensel lifts of the inverse
fibres.

Replay from the repository root; every command is silent on success:

```bash
python3 -B 04-computation/mahler_renormalized_linear_chart_thm4082.py \
  | diff - 05-knowledge/results/mahler_renormalized_linear_chart_thm4082.out
python3 -B -O 04-computation/mahler_renormalized_linear_chart_thm4082.py \
  | diff - 05-knowledge/results/mahler_renormalized_linear_chart_thm4082.out
python3 -B 04-computation/mahler_renormalized_linear_chart_thm4082_independent_audit.py \
  | diff - 05-knowledge/results/mahler_renormalized_linear_chart_thm4082_independent_audit.out
python3 -B -O 04-computation/mahler_renormalized_linear_chart_thm4082_independent_audit.py \
  | diff - 05-knowledge/results/mahler_renormalized_linear_chart_thm4082_independent_audit.out
```

Both sources contain zero Python `assert` nodes and zero floating literals.
The theorem proves an exact inter-scale normal form, defect valuation, and
typed fibre transport. It proves neither `Acal_s intersect G_s` nor
`Acal_s intersect N_odd` nonempty, does not turn `R_s` into a dynamical or
arithmetic conjugacy, and produces or rules out no Mahler Z-number.
