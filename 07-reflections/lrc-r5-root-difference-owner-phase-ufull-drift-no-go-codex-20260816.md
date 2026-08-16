# Root difference is the right drift coordinate, but the endpoint connection is still missing

**Status: FINITE-EXACT TYPED REPRESENTATION THEOREM CANDIDATE; independent
audit pending.**  The exact companion is
`04-computation/lrc_r5_root_difference_ufull_drift_alignment_probe_20260816.py`.
It starts from THM-2594's common-base joint table and THM-3518's independently
audited endpoint phase convention.  It proves no U_full ancestry relation,
physical current, absolute `H^1`, grouped coefficient, row exclusion, or
LRC(14).

## 1. The corrected coordinate map

The first folded transporter treated the deep coordinate

```text
theta=t-2u                                               (1)
```

as though it were the U_full drift.  THM-3514 identifies the honest
equivariant comparison.  Guard sheets transform by `a->a+tau`, while
THM-2471 collision roots transform by `u->u+tau`.  A common torsor gauge is

```text
u=a+c.                                                  (2)
```

Applied to two legs with the same `c`, equation (2) gives

```text
u_R-u_L=a_R-a_L=d.                                     (3)
```

Therefore U_full endpoint drift is typed like a *difference of collision
roots*.  In THM-2594's current/source convention, put

```text
s=u-q.                                                  (4)
```

Changing the leg orientation sends `s` to `-s`, one of the lawful
`F_13^*` torsor dilations.  The corrected common-base source table is

```text
B(ell,s)=sum_(u,theta)N(u,u-s,ell,theta).              (5)
```

Unlike a post hoc endpoint pairing, (5) is a marginal of the one Boolean
integral already proved in THM-2594.

## 2. The root-difference bank is much richer than the deep marginal

The exact table (5) has

```text
72/91 nonzero physical cells,
B(ell,0)=0 for every ell.                               (6)
```

The zero column is structural: the same-root current/source slice vanishes
in THM-2471.  Every other missing cell is retained in the exact digest rather
than silently filled.

Fourier transform in `ell` and `s` gives

```text
91/91 nonzero C7 x C13 coordinates,
rank of the seven septimal drift rows=6.                (7)
```

Thus the typed drift bank is not the earlier three-window rank-three object.
If the three deep windows are retained separately before their sum, the
twenty-one rows span rank eight:

```text
B_theta(ell,s)=sum_u N(u,u-s,ell,theta),
rank{Fourier_ell B_theta:theta=0,1,2}=8.               (8)
```

The corrected coordinate therefore removes the original rank obstruction
decisively.  It does not by itself create the endpoint connection.

## 3. No common projective convolution exists

For source Fourier vectors `X_k` and four U_full Walsh vectors `Y_k`, one
fixed channel map followed by one common drift convolution requires

```text
MX_k=mu_kY_k                 for every k in F_13.       (9)
```

As in the fixed-root probe, eliminate `mu_k` by the six pairwise wedge
equations and subtract the dimension of maps annihilating the source span.
The exact results are

```text
source bank               rank   equation rank   nullity   excess
-----------------------------------------------------------------
root difference, 7 rows     6         24             4        0
theta resolved, 21 rows     8         32            52        0. (10)
```

For the first row, the four-dimensional nullspace is exactly the space of
`4 x 7` maps annihilating the rank-six source.  For the second, the
fifty-two-dimensional nullspace is exactly

```text
4(21-8)=52.                                             (11)
```

Every formal solution kills the entire source curve.  All twelve nonzero
dilations of `F_13`, including the orientation reversal in (4), reproduce
the same zero-excess verdict for both banks.

The folded `C_7` census is even sharper.  For each of the eight choices of
one member from `{1,6},{2,5},{3,4}`, the four selected rows

```text
b=0 plus three marked nontrivial rows                    (12)
```

have rank four, and the `78 x 16` projective system has rank sixteen.  Its
nullspace is literally zero.  Exhausting all `7^4=2401` allocations finds no
exact or amplitude/shift-gauge common kernel; every allocation needs four
inequivalent channel kernels.

## 4. THM-3518 identifies the next sidecar—and its linear shadow still fails

THM-3518's independent route proves the translated guard identity

```text
sum_tau g_C(a+tau)g_D(b+tau)zeta^(-k tau)
 =zeta^(ka)K_(C,D,b-a)(k).                             (13)
```

The sign in the left-sheet phase `zeta^(ka)` is experimentally hostile: the
opposite sign gives a different endpoint tensor.  The bridge subtracts a
role class with `q_t=0` from one with `q_t=1`.  This suggests retaining the
left collision root `u`, not only `s`.

Accordingly the probe forms the two source owner slices

```text
B_k(ell,s)=sum_(u,theta)zeta^(-ku)N(u,u-s,ell,theta),
k=0,1.                                                  (14)
```

Their exact rank ledger is

```text
rank B_0=6,
rank B_1=6,
rank(B_1-B_0)=6,
rank(B_0 union B_1)=9.                                  (15)
```

The projective systems remain annihilator-only:

```text
bank                  rank   equation rank   nullity   excess
-------------------------------------------------------------
B_0                     6          24            4        0
B_1                     6          24            4        0
B_1-B_0                 6          24            4        0
B_0 union B_1           9          36           20        0. (16)
```

All twelve torsor dilations preserve the last verdict.  Hence merely adding
the correct owner Fourier coordinate to an already marginalized source bank
does not reproduce the U_full endpoint connection.

## 5. What has actually been learned

The sequence of hostiles separates three debts that had been blended:

```text
wrong coordinate debt:
  theta has rank 3;
  root difference has rank 6;                         paid by (3)--(7)

missing basepoint-phase debt:
  root difference forgets u;
  k=0,1 owner slices raise the union rank to 9;        paid by (13)--(15)

missing support/connection debt:
  even the correctly typed rank-nine bank has no
  common projective map to the four endpoint rows.     remains open.       (17)
```

The endpoint response is formed from chamber-dependent atom supports and a
translated guard kernel before its role-character contraction.  A fixed
linear map between already aggregated character rows has forgotten exactly
that support incidence.  The obstruction in (16) says the next object cannot
be another marginal or another frequency-independent channel matrix.

## 6. The next lawful branch transplant

The minimal source object now has coordinates

```text
(u,s,ell,theta),       q=u-s,                          (18)
```

while the endpoint object has

```text
(a,d,C,D),             b=a+d,                          (19)
```

with common gauge `u=a+c`, drift `s=+/-d`, and the owner phase (13).  What is
still missing is a Boolean support relation sending the source word/cell/deep
factors to the endpoint chamber predicates `C,D` *before* either side sums
over `u` or `a`.

The next decisive experiment should therefore test an address relation on

```text
(u,s,ell,theta; a,d,C,D)                               (20)
```

under the affine constraints above.  It should first recover the frozen
`H-q5` bridge bucket and its phase-normalization hostile, then ask about the
two `K4` tree factors.  A nonconvolutional or nonlinear Boolean relation is
now the expected shape; another endpoint spectral fit is not.

## 7. Tournament, tree, and cohomology boundary

The four Walsh channels still form the character basis of the chamber
`K4`, and the six pairwise comparisons are its tetrahedral edges.  Equations
(7), (8), and (15) show that having six, eight, or nine independent source
directions does not determine a `K4` connection.  This is the same
state-versus-connection distinction found in the Fibonacci/Berggren atlas.

THM-3518 also proves that all endpoint edge responses are gradients and all
56,592 cycle pairings vanish.  Nothing in the present source comparison
changes that: nonzero tree determinants remain `B^1` data, not a nonzero
absolute graph `H^1` class.

## 8. Connection ledger

| field | exact content |
|---|---|
| source | THM-2594 common-base `N(u,q,ell,theta)` |
| corrected drift | `s=u-q`, equivariantly `+/-d` under one common root/sheet torsor gauge |
| retained refinements | three theta windows and owner frequencies `k=0,1` |
| target | THM-3518 phase-normalized four-row U_full Walsh bank |
| positive result | full `91/91` spectrum; ranks `6`, `8`, and owner-union `9` |
| failed map | any fixed channel map plus one common drift convolution |
| exact obstruction | every projective nullspace is source-annihilator-only |
| missing sidecar | chamber/guard support relation before owner/root marginalization |
| next test | recover the frozen bridge bucket on a common `(u,s;a,d)` base before tree factors |

## 9. Reproduction

Run

```text
python -B 04-computation/lrc_r5_root_difference_ufull_drift_alignment_probe_20260816.py
python -B -O 04-computation/lrc_r5_root_difference_ufull_drift_alignment_probe_20260816.py
```

The pinned semantic digest is

```text
a954aa6ac96fa0a4e7b77e84473971c8444393409d61ca251b2c31de1b27779f.
```

The result is a typed finite-exact stopping theorem for the shared linear
convolution ansatz, not a proof that no broader branch transplant exists.
