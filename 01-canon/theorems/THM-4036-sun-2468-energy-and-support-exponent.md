---
id: THM-4036
title: "Sun 2-4-6-8 energy and support exponent"
status: >
  PROVED + VERIFIED-EXACT FINITE CONTROLS + INDEPENDENTLY PROOF-AUDITED.
  The canonical representation function has second moment
  O_epsilon(X^(13/12+epsilon)). Its represented support, globally and in
  every fixed arithmetic progression, is at least X^(1-epsilon), and the
  same lower bound holds for targets carrying any fixed positive fraction
  of the relevant mean multiplicity. Thus support and average-scale-rich
  support have logarithmic exponent one. The higher-role shadow
  C(x,4)+C(y,6)+C(z,8) alone represents at least X^(1/2-epsilon) targets;
  the two-role C(y,6)+C(z,8) support has exact logarithmic exponent 7/24.
  This proves neither positive density nor eventual coverage; leastness and
  classification of all holes remain open.
source: codex sun-minimal-sturmian session + independent hostile audit, 2026-08-24
depends_on:
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
  - THM-4028-sun-two-four-six-eight-average-order-criticality
related:
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
script: 04-computation/sun_2468_energy_support_thm4036.py
output: 05-knowledge/results/sun_2468_energy_support_thm4036.out
independent_audit_script: 04-computation/sun_2468_energy_support_independent_audit_thm4036.py
independent_audit_output: 05-knowledge/results/sun_2468_energy_support_independent_audit_thm4036.out
independent_audit_report: 05-knowledge/results/sun_2468_energy_support_independent_audit_thm4036.md
structured_hostile_script: 04-computation/sun_2468_structured_difference_multiplicity_thm4036.py
structured_hostile_output: 05-knowledge/results/sun_2468_structured_difference_multiplicity_thm4036.out
---

# THM-4036 -- Sun 2-4-6-8 energy and support exponent

**PROVED + VERIFIED-EXACT FINITE CONTROLS + INDEPENDENTLY PROOF-AUDITED.**
Let `a(n)` be the canonical role-labelled representation count from
THM-4028:

\[
a(n)=\#\left\{(w,x,y,z):
 \binom w2+\binom x4+\binom y6+\binom z8=n,\quad
 w\ge2, x\ge3, y\ge5, z\ge7\right\}.              \tag{1}
\]

Then, for every `epsilon>0`,

\[
\sum_{n\le X}a(n)^2\ll_\epsilon X^{13/12+\epsilon}.  \tag{2}
\]

If

```text
R(X)=#{1<=n<=X:a(n)>0},
```

then, for every `epsilon>0`, there are constants `c_epsilon>0` and
`X_epsilon` such that

\[
R(X)\ge c_\epsilon X^{1-\epsilon}
\quad\hbox{for every }X\ge X_\epsilon.                \tag{3}
\]

Moreover, for every fixed `q>=1` and residue `r mod q`, if

```text
R_(q,r)(X)=#{1<=n<=X:n=r mod q and a(n)>0},
```

then for every `epsilon>0` there are constants
`c_(q,r,epsilon)>0` and `X_(q,r,epsilon)` such that

\[
R_{q,r}(X)\ge c_{q,r,\epsilon}X^{1-\epsilon}
\quad\hbox{for every }X\ge X_{q,r,\epsilon}.          \tag{4}
\]

Thus the represented set has logarithmic support exponent one globally and
inside every fixed arithmetic progression. Constants in `(4)` are not
uniform when `q` grows.

More sharply, fix `0<theta<1` and put

```text
H_theta(X)=#{1<=n<=X:a(n)>=theta A(X)/X}.              (4a)
```

Then, for every `epsilon>0`,

\[
 H_\theta(X)\gg_{\theta,\epsilon}X^{1-\epsilon}.       \tag{4b}
\]

The analogous statement holds in every fixed class `r mod q`, with the
threshold equal to `theta` times that class's own mean multiplicity. Thus
logarithmically full support persists even after requiring multiplicity of
order `X^(1/24)`, rather than merely positivity.

For the discarded higher-role packets, write

```text
R_D(X)=#{0<=n<=X:b_D(n)>0}.
```

Then

\[
 R_{(6,8)}(X)\asymp X^{7/24},\qquad
 R_{(4,6,8)}(X)\gg_\epsilon X^{1/2-\epsilon}.         \tag{4c}
\]

The upper bound `R_(4,6,8)(X)<<X^(13/24)` is immediate from tuple count, so
its presently certified logarithmic support exponent lies in `[1/2,13/24]`.

## 1. Inheritance, domains, and the divisor-difference lemma

For `d in {2,4,6,8}`, put

```text
L_2=2, L_4=3, L_6=5, L_8=7,
f_d(t)=C(t,d),  t>=L_d.                                (5)
```

These are THM-4028's minimal domains: each higher role has exactly one zero
representative. Pascal's identity gives

\[
f_d(t+1)-f_d(t)=\binom t{d-1}>0                       \tag{6}
\]

throughout its declared domain, so every `f_d` is strictly increasing.

For `h!=0`, define the oriented difference count

\[
R_d(h)=\#\{(u,v):u,v\ge L_d,\ f_d(u)-f_d(v)=h\}.      \tag{7}
\]

Then

\[
R_d(h)\le \tau(d!\,|h|).                              \tag{8}
\]

Here is the exact variable-difference factorization behind `(8)`. Write

\[
S_j(u,v)=\sum_{i=0}^j u^{j-i}v^i.
\]

Expansion of the falling factorials gives

```text
2! (f_2(u)-f_2(v)) = (u-v) Q_2(u,v),
Q_2 = S_1-1;

4! (f_4(u)-f_4(v)) = (u-v) Q_4(u,v),
Q_4 = S_3-6S_2+11S_1-6;

6! (f_6(u)-f_6(v)) = (u-v) Q_6(u,v),
Q_6 = S_5-15S_4+85S_3-225S_2+274S_1-120;

8! (f_8(u)-f_8(v)) = (u-v) Q_8(u,v),
Q_8 = S_7-28S_6+322S_5-1960S_4+6769S_3
      -13132S_2+13068S_1-5040.                        (9)
```

It remains to control how many pairs a divisor can support. By swapping the
two variables, take `h>0`; strictness gives `u>v`. Put `q=u-v>=1`. Equation
`(9)` implies `q | d!h`. For a fixed `q`, define

\[
g_q(v)=f_d(v+q)-f_d(v).
\]

Another use of Pascal gives

\[
g_q(v+1)-g_q(v)
=\binom{v+q}{d-1}-\binom v{d-1}>0.                   \tag{10}
\]

Thus each positive divisor `q` of `d!h` supports at most one `v`, proving
`(8)`. Swapping `(u,v)` is a bijection between the `h` and `-h` counts, so
there is no extra factor two. In particular, for every `eta>0`,

\[
\max_{0<|h|\le X}R_d(h)\ll_{d,\eta}X^\eta,            \tag{11}
\]

by the standard elementary divisor bound. The case `h=0` is deliberately
excluded from `(8)` and is treated exactly below.

## 2. Exact diagonal and the energy-peeling inequality

For an ordered packet `D` of distinct roles from `{2,4,6,8}`, let

```text
b_D(n) = number of canonical D-tuples with atom sum n,
B_D(X) = sum_(0<=n<=X) b_D(n),
E_D(X) = sum_(0<=n<=X) b_D(n)^2.                       (12)
```

Since all atoms are nonnegative, a tuple counted by `B_D(X)` has every atom
at most `X`. The product box therefore gives

\[
B_D(X)\ll_D X^{\sigma_D},\qquad
\sigma_D=\sum_{d\in D}{1\over d}.                    \tag{13}
\]

Suppose `D=(d,D')`. Write the tail sum of a tuple `s` as `G(s)`. An ordered
pair counted by `E_(d,D')(X)` satisfies

\[
f_d(u)+G(s)=f_d(v)+G(t)\le X,
\qquad h=G(t)-G(s)=f_d(u)-f_d(v).                     \tag{14}
\]

For `Y>=0`, put

\[
N_d(Y)=\#\{u\ge L_d:f_d(u)\le Y\}.
\]

On the diagonal `h=0`, strictness forces `u=v`, but it does **not** force
the two role-labelled tail tuples `s,t` to coincide. Grouping them by their
common tail sum gives the exact diagonal formula

\[
\operatorname{Diag}_{d,D'}(X)
=\sum_{0\le m\le X} b_{D'}(m)^2 N_d(X-m).             \tag{15}
\]

Since `N_d(Y)\ll_d 1+Y^(1/d)`, this is

\[
\operatorname{Diag}_{d,D'}(X)
\ll_d X^{1/d}E_{D'}(X)                                \tag{16}
\]

for `X>=1`.

Off the diagonal, an ordered, role-labelled tail pair `(s,t)` fixes `h`.
Both tail sums lie in `[0,X]`, so `0<|h|<=X`, and `(11)` bounds the number
of oriented outer pairs `(u,v)`. There are at most `B_D'(X)^2` ordered tail
pairs. The target sum in `(14)` is then already determined; the terminal
condition `<=X` only deletes choices and introduces no additional
`X`-summation. Hence

\[
\operatorname{OffDiag}_{d,D'}(X)
\ll_{d,\eta}X^\eta B_{D'}(X)^2.                       \tag{17}
\]

Combining `(15)--(17)` proves the peeling inequality

\[
E_{(d,D')}(X)
\ll_{d,\eta}X^{1/d}E_{D'}(X)+X^\eta B_{D'}(X)^2.      \tag{18}
\]

This exact diagonal/off-diagonal split is the load-bearing gate: equal tail
sums retain their full multiplicity `b_D'(m)^2`, while unequal tail sums do
not acquire either a factor two or an extra target variable.

## 3. The exponent ledger

For a one-role packet, strictness gives

\[
E_{(d)}(X)=B_{(d)}(X)=N_d(X)\ll_d X^{1/d}.             \tag{19}
\]

Starting with degree eight and applying `(13)` and `(18)` in the order
`2,4,6,8` from the outside gives

| Packet | diagonal exponent | off-diagonal exponent | resulting bound |
|---|---:|---:|---:|
| `(8)` | -- | -- | `E_8 << X^(1/8)` |
| `(6,8)` | `1/6+1/8=7/24` | `2/8=1/4` | `E_68 << X^(7/24+eta)` |
| `(4,6,8)` | `1/4+7/24=13/24` | `2(1/6+1/8)=7/12` | `E_468 << X^(7/12+eta)` |
| `(2,4,6,8)` | `1/2+7/12=13/12` | `2(1/4+1/6+1/8)=13/12` | `E_2468 << X^(13/12+eta)` |

The epsilon bookkeeping is finite and explicit. It is enough to prove the
claim for `0<epsilon<=1/12`, since that implies every weaker bound with a
larger epsilon. Invoke every divisor estimate `(11)` with
`eta=epsilon/4<1/24`. Inductively retain the same loss `eta`: the strict exponent
gaps at the first two peels absorb the smaller branch, and the two branches
at the final peel both have exponent `13/12+eta`. Since `eta<=epsilon`, this
proves `(2)`. The primary finite control evaluates the recurrence induced by
`(18)` over all 24 peeling orders; within this peeling scheme only
`(2,4,6,8)` and `(2,4,8,6)` attain the minimum `13/12`.

The independent audit checks `(8)` exhaustively for `1<=|h|<=5000`, including
the zero-atom endpoints, and directly decomposes the finite energies. At
`X=50` it obtains

```text
E_68   = 16   = 14 diagonal  +     2 off-diagonal,
E_468  = 108  = 70 diagonal  +    38 off-diagonal,
E_2468 = 3049 = 677 diagonal + 2372 off-diagonal.      (20)
```

These are finite exact controls on the proof's multiplicities, not
extrapolations used to establish `(2)`.

The first intermediate loss can be removed completely. At the `(6,8)` peel,
`(18)` gives

\[
 E_{(6,8)}(X)\ll X^{7/24}+X^{1/4+\eta}.
\]

Choose once and for all `eta=1/48`; then the second exponent is `13/48`,
strictly below `7/24=14/48`. Hence

\[
 E_{(6,8)}(X)\ll X^{7/24}.                            \tag{20b}
\]

### Support ladder for the discarded packets

The product-box estimate `(13)` also has the matching lower bound

\[
 B_D(X)\gg_D X^{\sigma_D}.                            \tag{20a}
\]

Indeed, restrict each coordinate in a nonempty packet `D` to a fixed positive
fraction of its `X^(1/d)` range so that every atom is at most `X/|D|`.
All tuples in that box have total at most `X`, and the box contains a constant
multiple of `X^(sigma_D)` tuples.

Cauchy--Schwarz gives `B_D(X)^2<=R_D(X)E_D(X)`. For `D=(6,8)`, `(20a)` and
`(20b)` yield

\[
 R_{(6,8)}(X)\gg X^{7/24}.
\]

The reverse bound `R_(6,8)(X)<=B_(6,8)(X)<<X^(7/24)` proves the first claim in
`(4c)`. For `D=(4,6,8)`, the same argument
uses `2sigma_D=13/12` and `E_D(X)<<_epsilon X^(7/12+epsilon)`, giving the
second claim in `(4c)`. Thus the exact peeling proof also supplies an
unconditional large-support theorem for the three-higher-binomial sumset;
it supplies no short-interval or maximum-gap control.

### Sharpness and the polylogarithmic-loss frontier

The main power `13/12` is best possible. There are at most `floor(X)` target
fibres, so `(21)` and Cauchy--Schwarz in the opposite direction give

\[
 E_{(2,4,6,8)}(X)\ge {A(X)^2\over\lfloor X\rfloor}
 =V^2X^{13/12}\left(1+O(X^{-1/8})\right).             \tag{20c}
\]

No valid argument can replace the exponent by `13/12-delta`. Removing the
arbitrary `X^epsilon` loss is subtler: the pointwise divisor-multiplicity
input cannot be replaced by any fixed power of `log X`.

Let `P_k` be the product of the first `k` odd primes, put `t=24P_k`, and let
`h=C(t,4)`. For every divisor `a|P_k`, define

```text
b=2h/a,       u=(a+b+1)/2,       v=(b-a+1)/2.         (20d)
```

Then `a` is odd, `b` is even, `u,v>=2`, and

\[
 \binom u2-\binom v2=h.
\]

Distinct divisors give distinct ordered pairs, so `R_2(h)>=2^k`. Moreover

\[
h=\left[\binom t4+\binom56+\binom78\right]
 -\left[\binom34+\binom56+\binom78\right],          \tag{20e}
\]

so this is an actual `(4,6,8)` tail difference, not an artificial input.
Bertrand's postulate gives `log h=O(k^2)`, hence `2^k` eventually exceeds
every fixed power of `log h`. Averaging over tail pairs is therefore
essential.

There is an exact sufficient replacement. For a tail sequence `b(n)` on
`[0,X]`, define

\[
 B_r(q)=\sum_{n\equiv r\ (q)}b(n),\qquad
 K_b(q)=\sum_{r\bmod q}B_r(q)^2.                     \tag{20f}
\]

At an outer degree-`d` peel, summing first over the index gap `j=u-v` gives

\[
 \operatorname{OffDiag}_{d,b}(X)
 \le\sum_{1\le j\le U_d(X)}
 K_b\!\left({j\over\gcd(j,d!)}\right).               \tag{20g}
\]

Indeed `j|d!h` forces `j/gcd(j,d!)|h`; `K_b` counts all such ordered tail
pairs and harmlessly includes the diagonal. Consequently the growing-modulus
gate

\[
 K_{b_D}(q)\ll L_D(X)\left({B_D(X)^2\over q}+E_D(X)\right)             \tag{20h}
\]

for `D=(6,8)`, `q<<X^(1/4)`, and for `D=(4,6,8)`, `q<<X^(1/2)`, would imply

\[
 E_{(2,4,6,8)}(X)
 \ll X^{13/12}(\log X)^{2A+1}                        \tag{20i}
\]

whenever both `L_D(X)<< (log X)^A`. This follows from `(20g)`, the harmonic
sum, `(20b)`, and the two successive peels. The gate `(20h)` is **OPEN**:
THM-4027/4028 control each fixed modulus, not the anisotropic polynomial boxes
uniformly up to these growing ranges. Thus `(20i)` is a proved reduction, not
an unconditional estimate.

## 4. Global and fixed-progression support

THM-4028 proves

\[
A(X):=\sum_{n\le X}a(n)
=V X^{25/24}+O(X^{11/12}),\qquad V>0.                 \tag{21}
\]

Because `11/12<25/24`, eventually
`A(X)>=(V/2)X^(25/24)`. Cauchy--Schwarz on the nonzero fibres gives

\[
A(X)^2\le R(X)\sum_{n\le X}a(n)^2.                   \tag{22}
\]

Apply `(2)` with the same requested `epsilon`. Equations `(21)--(22)` give

\[
R(X)\gg_\epsilon
X^{2(25/24)-(13/12+\epsilon)}=X^{1-\epsilon},        \tag{23}
\]

with the quantified eventual form `(3)`. Since `R(X)<=X`, letting
`epsilon` tend to zero yields

\[
\lim_{X\to\infty}{\log R(X)\over\log X}=1.           \tag{24}
\]

For fixed `q,r`, set

\[
A_{q,r}(X)=\sum_{\substack{n\le X\\n\equiv r\pmod q}}a(n).
\]

THM-4028 and THM-4027 give

\[
A_{q,r}(X)={\sigma_q(r)\over q}V X^{25/24}
+O_q(X^{11/12}),\qquad \sigma_q(r)>0.                 \tag{25}
\]

Thus `A_(q,r)(X)>=c_(q,r)X^(25/24)` eventually. Applying
Cauchy--Schwarz only inside this class and then enlarging its square sum to
the global energy gives

\[
A_{q,r}(X)^2
\le R_{q,r}(X)
   \sum_{\substack{n\le X\\n\equiv r\pmod q}}a(n)^2
\le R_{q,r}(X)\sum_{n\le X}a(n)^2.                   \tag{26}
\]

Equations `(2)`, `(25)`, and `(26)` prove `(4)`. The error in `(25)` is
smaller than its main term by the fixed power `1/8`, so it is absorbed before
Cauchy--Schwarz and causes no exponent loss.

### Average-scale-rich support

Fix `0<theta<1` and abbreviate `T=theta A(X)/X`. The targets outside
`H_theta(X)` contribute strictly less than `XT=theta A(X)`. Hence the targets
inside it carry more than `(1-theta)A(X)` total mass. Cauchy--Schwarz on just
those targets and `(2)` give

\[
 (1-\theta)^2A(X)^2
 \le H_\theta(X)\sum_{n\le X}a(n)^2,
\]

which, together with `(21)`, proves `(4b)`. Moreover

\[
 {A(X)\over X}=V X^{1/24}+O(X^{-1/12}),              \tag{26a}
\]

so every counted target has a fixed positive fraction of the natural
cutoff-scale mean multiplicity.

For a fixed class `r mod q`, let

\[
 M_{q,r}(X)=\#\{1\le n\le X:n\equiv r\pmod q\}
            ={X\over q}+O(1)
\]

and define `H_(q,r,theta)(X)` using the threshold
`theta A_(q,r)(X)/M_(q,r)(X)`. The same mass split inside that class, followed
by the global energy upper bound, yields

\[
 H_{q,r,\theta}(X)\gg_{q,r,\theta,\epsilon}X^{1-\epsilon}.
\]

By `(25)`, its threshold is asymptotic to
`theta sigma_q(r)V X^(1/24)`. This is a cutoff-scale statement: it neither
gives one fixed target high multiplicity along all larger cutoffs nor controls
the pointwise lower tail below a fixed fraction of the mean.

There is also a uniform anti-condensation form. Fix `delta>0`, and for each
`X` let `T_X` be any subset of `[1,X]` with `|T_X|<=X^(1-delta)`. Use `(2)`
with loss `delta/2`. Then

\[
 \sum_{n\in T_X}a(n)
 \le |T_X|^{1/2}\left(\sum_{n\le X}a(n)^2\right)^{1/2}
 \ll_\delta X^{25/24-\delta/4}=o(A(X)).               \tag{26b}
\]

The estimate is uniform over the choice of `T_X`. The identical calculation
inside any fixed class `r mod q`, still bounding its square energy by the
global one, is `o(A_(q,r)(X))`. Thus no fixed-power-thin carrier can hold a
positive fraction of the representation mass. This remains compatible with
a hole set of density one, because it constrains tuple mass rather than the
number of zero fibres.

## 5. Classification firewall and exact leastness status

The conclusion is a **support-exponent and average-scale-rich-support**
theorem. It does not prove positive density, density one, `o(X)` holes,
eventual representation, finiteness of the hole set, or a pointwise asymptotic
for `a(n)`. In particular it does not classify all counterexamples.

THM-4027 supplies a separate useful restriction. Every residue modulo every
fixed modulus is locally represented. By adding arbitrarily large multiples
of a binomial input period to a chosen coordinate, one obtains arbitrarily
large exactly represented integers in that same target residue class.
Consequently no arithmetic progression has an exceptional tail. If the hole
set is infinite, it therefore cannot be eventually periodic or a nonempty
finite union of residue classes. This rules out a congruence-only
classification; it does not decide whether the hole set is finite or
infinite.

The exact leastness ledger is:

- **FINITE-EXACT:** THM-4026 proves
  `N=896315812331399` is a hole in the convention `(1)`.
- **CITED, not imported as a repository prefix certificate:** the published
  search checked the full prefix through `2*10^12`.
- **OPEN:** the interval

  ```text
  2,000,000,000,001 <= n <= 896,315,812,331,398
  ```

  contains `894,315,812,331,398` integers and has not been exhaustively
  discharged here.
- **OPEN:** whether `N` is least, and the finiteness, infinitude, density, or
  structural list of all holes.

Sparse CRT classes and local-density factors from THM-4027 can rank a search,
but universal modular solubility means they are not necessary-and-sufficient
hole tests. Searching selected CRT slices cannot certify the full leastness
interval.

## 6. A boundary lemma for any least hole

There is nevertheless an exact carry/height restriction at a hole's
neighbours. Let `n` be a hole in the canonical convention `(1)`.

1. If `n-1` is represented, then **every** representation of `n-1` has
   `x>=4`, `y>=6`, and `z>=8`. Indeed, if one of the higher coordinates were
   its canonical zero index `3,5,7`, respectively, replacing it by the first
   positive index `4,6,8` would add exactly one and represent `n`.
2. If `n+1` is represented, then **every** representation of `n+1` avoids
   `x=4`, `y=6`, and `z=8`. Lowering such a first-positive coordinate to its
   canonical zero index would subtract exactly one and represent `n`.

For a least hole the first hypothesis is automatic, so the first statement
is a genuine necessary boundary condition. At THM-4026's target, the exact
neighbour controls are

```text
a(N-1)=89,       a(N)=0,       a(N+1)=67.              (27)
```

Hence all 156 neighbour representations satisfy their corresponding shield.
This invariant retains a coordinate's height relative to the zero-to-one
carry. A residue-only projection destroys precisely that information.

## 7. Reproduction and scope

From the repository root, run

```text
python -B 04-computation/sun_2468_energy_support_thm4036.py
python -B 04-computation/sun_2468_energy_support_independent_audit_thm4036.py
```

The primary companion checks the four falling-factorial expansions, all 24
peeling orders, the support exponent arithmetic, and the exact leastness gap.
The independent companion checks endpoint strictness and the divisor bound,
the exact diagonal/off-diagonal decomposition at three finite cutoffs, and
classwise Cauchy--Schwarz modulo `q<=8`. The proof of `(2)--(4)` is the
argument above; the scripts are exact hostile controls. **QED.**
