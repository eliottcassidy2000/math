# Independent hostile audit: Sun 2-4-6-8 energy/support theorem

Status: **PASS, with two proof-presentation clarifications; no canonical edit**
(2026-08-24).

Audited source: `.scratch/sun_minimal_classifier/REPORT.md`, Sections 2--4.

## Verdict

The divisor-difference lemma, energy peeling inequality, exponent `13/12`,
global support consequence, and fixed-arithmetic-progression support
consequence are valid under THM-4028's canonical domains

```text
w>=2, x>=3, y>=5, z>=7.
```

I found no missing factor, zero-fibre exception, or direction reversal. The
theorem is ready for a normal independent referee, subject to explicitly
retaining the canonical one-zero-representative convention.

## 1. Strictness at the zero boundary

For degree `d`, the higher-atom domain begins at `L_d=d-1`, where
`f_d(L_d)=0`. For `q>=1`,

```text
g_q(v+1)-g_q(v)=C(v+q,d-1)-C(v,d-1).
```

This remains strictly positive at the endpoint `v=d-1`: the second term is
`C(d-1,d-1)=1`, while the first is strictly larger than one. There is no flat
step caused by the zero atom. For `d=2`, the canonical domain begins at `2`,
so the triangular atom is already positive and the same argument applies.

Exact first boundary values `g_q(L_d)`, `q=1,...,5`, are

```text
d=2: 2,5,9,14,20
d=4: 1,5,15,35,70
d=6: 1,7,28,84,210
d=8: 1,9,45,165,495.
```

Thus `f_d` and every `g_q` are strictly increasing on the declared domains.

For `h>0`, `u>v`, `q=u-v` divides `d!h`, and each positive divisor `q`
supports at most one `v`. For `h<0`, swapping `(u,v)` is a bijection to the
positive case. Therefore the oriented count really satisfies

```text
R_d(h)<=tau(d!|h|)
```

with no factor two. Exhaustive exact checks for every `1<=|h|<=5000` and all
four degrees pass.

## 2. The `h=0` diagonal

When `G(s)=G(t)=m`, strictness forces `u=v`, but `s` and `t` need not be the
same tuple. The exact diagonal contribution is

```text
sum_{0<=m<=X} b_D'(m)^2 N_d(X-m),                      (1)
N_d(Y)=#{u>=L_d:f_d(u)<=Y}.
```

Since `N_d(Y)<=N_d(X)<<X^(1/d)`, (1) is at most

```text
X^(1/d) E_D'(X).
```

This is exactly what the scratch proof's phrase “over each equal pair counted
by `E_D'(X)`” means. It does not collapse equal-sum tail tuples to identical
tuples, and it does not treat `h=0` using the divisor bound.

## 3. Off-diagonal tuple multiplicity

For every **ordered**, role-labelled tail pair `(s,t)`, the difference
`h=G(t)-G(s)` is fixed. The number of oriented outer pairs `(u,v)` is
`R_d(h)`. Hence

```text
offdiag <= max_{0<|h|<=X}R_d(h)
           * #{ordered tail pairs with unequal sums}
        <= O(X^eta) B_D'(X)^2.                         (2)
```

There is no additional target variable: the target sum is determined by the
chosen outer and tail tuples and the condition `<=X` only deletes choices.
There is also no missing factor two. Reversing a tail pair produces the
separate ordered pair and changes `h` to `-h`; both are already in `B^2`, and
`R_d(-h)=R_d(h)`.

Finite exact decompositions confirm the multiplicities. For `X=50`:

```text
E_68=16       = diagonal 14  + off-diagonal 2,
E_468=108     = diagonal 70  + off-diagonal 38,
E_2468=3049   = diagonal 677 + off-diagonal 2372.
```

For the final peel, `B_468=60`, `E_468=108`, so the number of ordered
unequal-sum tail pairs is exactly

```text
B_468^2-E_468=3492,
```

which the independent decomposition reproduces.

## 4. Exponent and epsilon bookkeeping

The exponent ledger is correct:

```text
E_8       << X^(1/8),
E_68      << X^(7/24+epsilon),
E_468     << X^(7/12+epsilon),
E_2468    << X^(13/12+epsilon).
```

Presentation clarification: the `eta` in each divisor bound should not be
read as one error exponent blindly added at every row. Given a requested
final `epsilon`, choose the divisor-bound parameters sufficiently small at
each finite peel (for example `epsilon/4`) and absorb the first strict exponent
gaps. Equivalently, prove the displayed estimates inductively with a fresh
fraction of `epsilon`. This yields every requested final `epsilon>0`.

## 5. Global and fixed-AP support

For a fixed class `r mod q`, define

```text
S_qr(X)=sum_{n<=X,n=r mod q} a(n),
E_qr(X)=sum_{n<=X,n=r mod q} a(n)^2.
```

THM-4028 and THM-4027 give

```text
S_qr(X)=(sigma_q(r)/q)V X^(25/24)+O_q(X^(11/12)),
sigma_q(r)>0.
```

Thus `S_qr(X)>=c_qr X^(25/24)` eventually. Exact Cauchy--Schwarz gives

```text
S_qr(X)^2 <= R_qr(X) E_qr(X)
            <= R_qr(X) E_2468(X).
```

The energy theorem therefore yields, for every fixed `q,r` and every
`epsilon>0`,

```text
R_qr(X)>=c_(q,r,epsilon) X^(1-epsilon).
```

This consequence is valid. The constant and threshold depend on the fixed
class, and nothing is uniform for growing `q`; the scratch report already
states the fixed-modulus scope. Exact finite Cauchy--Schwarz checks pass for
every class modulo `1<=q<=8` at `X=5000`.

The same calculation globally gives `R(X)>=c_epsilon X^(1-epsilon)` and,
using `R(X)<=X`,

```text
lim log R(X)/log X=1.
```

It gives logarithmic support exponent only, not positive density, density one,
or a bound `o(X)` for the holes.

## Reproduction

```text
python -B .scratch/sun_energy_hostile_audit/audit.py
```

The audit checks endpoint strictness, the divisor bound through `|h|=5000`,
exact diagonal/off-diagonal decompositions at `X=50,200,1000`, ordered tail
multiplicity, and finite classwise Cauchy--Schwarz.
