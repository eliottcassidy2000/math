---
id: THM-2781
title: "Terminal-tail perfect-power rigidity and response count"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Over a characteristic-zero field, let f(0)=1, deg(f)<=d, write
  alpha=a/b in lowest positive terms, and suppose N=d alpha is integral.
  Then the d-1 coefficients of f^alpha in degrees N+1,...,N+d-1 vanish if
  and only if f is a b-th power; equivalently every coefficient above N
  vanishes.  Cubic and quartic hostiles show that d-2 zeros do not suffice
  uniformly, while an unreduced exponent shows lowest terms are essential.
  The theorem unifies the THM-2110 cubic two-response and THM-2778 quartic
  three-response mechanisms but does not construct a response bank, derive a
  Keller chart, or prove JC(2)/DC(2).
source: root/terminal-tail-perfect-power-rigidity-2026-07-28
depends_on: []
related:
  - THM-2110-cubic-fiber-degree-thirteen-faber-tail-gate
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2778-all-degree-complete-chosen-sheet-split-exact-prefix-closure
script: 04-computation/terminal_tail_perfect_power_rigidity_thm2781.py
output: 05-knowledge/results/terminal_tail_perfect_power_rigidity_thm2781.out
script_sha256: f911eabd8ac9b235180cb413bfbfe43b6b38de60d0d2740a1467733e8f2340c9
output_sha256: 325d695e53f860f239479d0ee6f9fbd46efd2dd068e06fe7d50c4bd822e3ddc7
hash_basis: LF-normalized bytes
---

# THM-2781 -- terminal-tail perfect-power rigidity

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The two missing coefficients in the depressed-cubic Faber gate and the three
missing coefficients in the depressed-quartic gate are instances of one
formal recurrence theorem.  The number of consecutive coefficients is
controlled by the degree bound, while the perfect-power conclusion is
controlled by the reduced denominator of the exponent.

## 1. The terminal-tail theorem

Let `K` be a characteristic-zero field.  Fix positive integers `d,a,b` with

```text
gcd(a,b)=1,             N=da/b in Z.                  (1)
```

The integrality in `(1)` implies `b|d`.  Let

```text
f(z)=1+p_1z+...+p_d z^d in K[z],                     (2)
```

where the actual degree may be less than `d`, and define its constant-one
formal power by

```text
f(z)^(a/b)=sum_(n>=0)c_n z^n in K[[z]].              (3)
```

Then the following are equivalent:

1. `c_(N+1)=...=c_(N+d-1)=0`;
2. `c_n=0` for every `n>N`;
3. `f=h^b` for a unique `h in K[z]` with `h(0)=1`.

Thus `d-1` consecutive terminal zeros are a complete universal certificate
for the denominator-`b` perfect-power locus, even when `p_d=0`.

## 2. The recurrence closes the entire tail

Put `alpha=a/b`.  The formal identity

```text
f(f^alpha)'=alpha f' f^alpha
```

gives, for every `n>=1`,

```text
n c_n=sum_(r=1)^d p_r ((alpha+1)r-n)c_(n-r),         (4)
```

with negative-index coefficients understood as zero.  Assume condition 1.
At `n=N+d`, every term with `r<d` uses one of

```text
c_(N+1),...,c_(N+d-1),                               (5)
```

while the only remaining predecessor is `c_N`, whose multiplier is

```text
(alpha+1)d-(N+d)=d alpha-N=0.                        (6)
```

Hence `c_(N+d)=0`.  For every later `n`, all predecessors `c_(n-r)` have
index at least `N+1` and are already zero.  Induction in `(4)` proves
condition 2.  This step uses the declared degree bound, not nonvanishing of
the top coefficient.

## 3. UFD multiplicities detect the denominator

Under condition 2, the series in `(3)` is a polynomial `g` of degree at most
`N`, and

```text
g^b=f^a.                                               (7)
```

For every irreducible `pi in K[z]`, unique factorization gives

```text
b v_pi(g)=a v_pi(f).                                  (8)
```

Since `gcd(a,b)=1`, equation `(8)` forces `b|v_pi(f)` for every `pi`.
Therefore `f=uH^b`.  Because `f(0)=1`, normalizing `h=H/H(0)` absorbs the
unit and gives `f=h^b` with `h(0)=1`.  The constant-one formal `b`th root is
unique, proving condition 3.

Conversely, if `f=h^b`, then the constant-one branch in `(3)` is `h^a` and

```text
deg(h^a)=a deg(f)/b <= ad/b=N.                        (9)
```

Thus condition 3 implies condition 2, and condition 2 trivially implies
condition 1.  This proves the equivalence.

## 4. The universal count cannot be shortened uniformly

The theorem does not claim that every specialized family needs all `d-1`
zeros.  It says that no smaller count works uniformly over the admissible
degree/exponent data.

For `d=3`, `a/b=1/3`, and `N=1`, take

```text
f=1+3z+3z^2,
f^(1/3)=1+z+0z^2-(1/3)z^3+... .                     (10)
```

The first of the required two tail coefficients vanishes, but `f` is not a
cube.  For `d=4`, `a/b=3/2`, and `N=6`, take

```text
f=1+z^3,
[z^7]f^(3/2)=[z^8]f^(3/2)=0,
[z^9]f^(3/2)=-1/16.                                  (11)
```

Here two of the required three coefficients vanish, but `f` is not a
square.  This is exactly the response-active root-zero hostile in reduced
quartic degree six.

Lowest terms in `(1)` are also load-bearing.  Display `alpha=2/4` and take

```text
f=(1+z^2)^2.                                          (12)
```

Its formal `alpha` power is `1+z^2`, so the entire displayed terminal tail
vanishes, but `f` is not a fourth power.  The correct conclusion uses the
reduced denominator two.

## 5. Cubic and quartic response dictionaries

For the depressed cubic, let

```text
S=1+pz^2+qz^3,        alpha=n/3,       3 does not divide n.
```

Here `(d,b,N)=(3,3,n)`.  Vanishing of the two consecutive coefficients in
degrees `n+1,n+2` makes `S` a cube.  Its constant-one cube root has degree at
most one, and the missing linear term forces that root to be `1`; hence
`p=q=0`.  This is the all-degree adjacent-tail mechanism in THM-2110.

For the depressed quartic used in THM-2778, let

```text
f=1+2dz^2+qz^3+(d^2-s)z^4,
alpha=M/4,             M=4k-2.                       (13)
```

Now `(d_degree,b,N)=(4,2,M)`.  Vanishing in degrees `M+1,M+2,M+3` makes
`f` a square.  Write its constant-one square root as `1+uz+vz^2`.  The
linear coefficient gives `u=0`; comparison in `(13)` gives

```text
v=d,                 q=0,                 s=0.        (14)
```

Thus the common triple-response support is the single weighted-projective
point `P_infty`.  THM-2778 then adds exact-prefix localization, the
root-zero section obstruction, slope four, and the vertical contradiction;
none of those global steps follows from the present formal theorem alone.

## 6. Transfer rule and scope

For a prospective degree-`d` source-fibre chart, the reusable design is:

```text
d-1 consecutive Faber/Gauss--Manin observables
        -> terminal recurrence
        -> reduced-denominator perfect-power locus
        -> use missing low coefficients to classify that locus.         (15)
```

This is a response-count rule, not a promise that the required observables
exist or are constant along a Keller pair.  It does not derive split or
nonsplit chart entry, control exact-prefix residues, classify finite poles,
or prove `JC(2)` or `DC(2)`.

## 7. Exact hostile controls

Run

```bash
python3 04-computation/terminal_tail_perfect_power_rigidity_thm2781.py
python3 -O 04-computation/terminal_tail_perfect_power_rigidity_thm2781.py
```

The assertion-free rational companion checks `163` admissible parameter
rows with `1<=d<=10`, `326` constructed perfect powers, `83` nonpowers,
the terminal multiplier and complete predecessor window, integral-exponent
and vanishing-top-coefficient boundaries, and hostiles `(10)--(12)`.  Its
`2103` finite gates verify the algebraic interfaces; the all-degree
quantifier comes from Sections 2--3.

```text
script_sha256 = f911eabd8ac9b235180cb413bfbfe43b6b38de60d0d2740a1467733e8f2340c9
output_sha256 = 325d695e53f860f239479d0ee6f9fbd46efd2dd068e06fe7d50c4bd822e3ddc7
hash_basis    = LF-normalized bytes
```
