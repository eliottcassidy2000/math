# THM-2781 independent hostile audit

## Verdict

**THEOREM/PROOF: PASS.  CURRENT EXACT COMPANION: FAIL ONE ADVERTISED
HOSTILE GATE; REPAIR REQUIRED BEFORE PROVED-canon promotion.**

The terminal recurrence and the UFD equivalence are correct with all stated
quantifiers, including `d=1`, `b=1`, and `deg(f)<d`.  The cubic and quartic
specializations are correct.  Lowest terms, constant-term normalization, and
characteristic zero are genuinely load-bearing.

The sole defect found is in the canonical script's unreduced-exponent
control.  It checks only

```python
len(trim(unreduced)) - 1 == 4
```

as evidence that `(1+z^2)^2` is not a fourth power.  Degree four is compatible
with being the fourth power of a linear polynomial, so this gate verifies an
intermediate statistic rather than the advertised consequence.  The example
itself is valid and the theorem text is correct; the exact companion must
test the unique linear fourth-root candidate or the coefficient
contradiction.

## 1. Independent recurrence derivation

Write

```text
f=sum_(r=0)^d p_r z^r,       p_0=1,
g=f^(a/b)=sum_(n>=0)c_nz^n.
```

From

```text
f g'=(a/b)f'g
```

the coefficient of `z^(n-1)` is

```text
n c_n
 =sum_(r=1)^d p_r(((a/b)+1)r-n)c_(n-r).                (1)
```

Put `N=da/b`.  If

```text
c_(N+1)=...=c_(N+d-1)=0,                               (2)
```

then at `n=N+d` every `r<d` predecessor lies in `(2)`, while the `r=d`
multiplier is

```text
((a/b)+1)d-(N+d)=0.                                    (3)
```

Hence `c_(N+d)=0`.  For every later `n`, all `d` predecessors have index
strictly above `N`, so induction closes the whole tail.  A vanishing top
coefficient `p_d` causes no problem; the proof uses a declared degree bound,
not exact degree.

Characteristic zero enters twice: every positive integer `n` in `(1)` is
invertible, and the constant-one rational power is uniquely defined.

## 2. Independent UFD derivation

Once the tail closes, `g` is a polynomial and formal exponentiation gives

```text
g^b=f^a.                                                (4)
```

For every irreducible `pi in K[z]`,

```text
b v_pi(g)=a v_pi(f).                                   (5)
```

Since `gcd(a,b)=1`, `(5)` forces `b|v_pi(f)`.  Thus
`f=uH^b`.  Because `f(0)=1`, `H(0)!=0` and

```text
u=H(0)^(-b),        h=H/H(0)
```

give `f=h^b`, `h(0)=1`.  Conversely, `f=h^b` gives the constant-one branch
`g=h^a` and

```text
deg(g)=a deg(f)/b <= ad/b=N.
```

The root is unique because any two polynomial `b`th roots differ by a
constant `b`th root of unity, and their constant terms are both one.

## 3. Boundary probes

### `d=1`

The response block has length zero.  This is not an exception: `b|d` and
lowest terms force `b=1`, so every `f` is a first power and `f^a` has degree
at most `a=N`.  The equivalence is correctly vacuous.

### `b=1`

Every `f` is a first power and `f^a` is a polynomial of degree at most
`ad=N`; all terminal coefficients vanish automatically.

### Actual degree below `d`

For example, with declared `(d,a,b,N)=(6,2,3,4)`,

```text
f=(1+2z)^3
```

has actual degree three and `f^(2/3)=(1+2z)^2`.  The full declared
five-coefficient terminal block vanishes.  Independent exhaustive tests
also include every top-zero coefficient pattern through declared degree
five.

### Uniform insufficient-tail hostile

The cubic example in canon is the first member of a uniform family.  For
every `d>=2`, take

```text
f_d=(1+z)^d-z^d,          alpha=1/d,         N=1.       (6)
```

Modulo `z^d`, its constant-one `d`th root equals `1+z`, while comparison at
degree `d` gives

```text
[z^2]f_d^(1/d)=...=[z^(d-1)]f_d^(1/d)=0,
[z^d]f_d^(1/d)=-1/d.                                   (7)
```

The polynomial `f_d` has degree `d-1` and is not a `d`th power.  Thus
`d-2` zeros fail for every `d>=2`, not only somewhere in the cubic/quartic
range.  The audit checks `(6)--(7)` exactly for `2<=d<=12`.

The installed cubic hostile is exactly `d=3` in `(6)`.  The quartic hostile
`f=1+z^3`, `alpha=3/2`, has coefficients

```text
[z^7]=[z^8]=0,             [z^9]=-1/16,
```

and is not a square because its degree is odd.  Both are correct.

### Unreduced exponent

For displayed `alpha=2/4`,

```text
f=(1+z^2)^2
```

has polynomial formal power `f^(2/4)=1+z^2`.  It is a square but not a
fourth power.  The repaired conclusion uses reduced denominator two.

An exact non-fourth-power check is:

```text
deg(f)=4, so a constant-one fourth root would be 1+uz;
[z]f=0 forces 4u=0, hence u=0;
but [z^2]f=2, contradiction.                            (8)
```

### Constant term

The condition `f(0)=1` provides a `K`-rational, uniquely selected formal
branch.  Over `Q`, even the constant polynomial `f=2` has no square root in
`Q[[z]]`.  If a nonzero constant term is supplied together with a chosen
`b`th root in `K`, rescaling reduces to the normalized theorem; without
that sidecar the expression is not defined over `K`.

### Characteristic zero

The scope cannot simply be weakened to "`char(K)` does not divide `b`".
Over `F_2`, take

```text
d=3, a=2, b=3, N=2,       f=1+z^3.                    (9)
```

The constant-one cube root of

```text
f^2=1+z^6
```

is unique because `3=1` in `F_2`, and begins

```text
f^(2/3)=1+z^6+O(z^7).                                  (10)
```

Thus the required coefficients `c_3,c_4` vanish but `c_6=1`; `f` is not a
cube.  The induction fails precisely when its left multiplier reaches a
multiple of the characteristic.

## 4. Independent exact path

The scratch companion does not compute rational powers from recurrence
`(1)`.  It forms `f^a` and solves

```text
g^b=f^a
```

coefficient by coefficient; at degree `n`, the new coefficient occurs with
linear multiplier `b`.  It then independently checks `(1)`.

Its exact universe contains:

```text
42 admissible (d,a,b) rows with 1<=d<=5, 1<=a<=5;
all 3,408 polynomials with coefficients in {-1,0,1};
1,832 positive window-zero implications;
5,888 exact gates;
11 uniform insufficient-tail hostiles, 2<=d<=12.
```

It separately checks `d=1`, `b=1`, actual degree below the declared bound,
the unreduced repair, constant-term scope, the `F_2` hostile, and the cubic
and quartic specializations.  Normal and optimized runs byte-match the
stored scratch output; there are zero Python `assert` nodes.

Scratch artifacts:

```text
.scratch/thm2781_hostile_audit/thm2781_independent_hostile_audit.py
.scratch/thm2781_hostile_audit/thm2781_independent_hostile_audit.out
```

LF-normalized hashes:

```text
script e3e2c4535c25e1c8271ae953f31b1154c2aa957e23b855398b0d99e69cf93799
output 2e17edec5f29da09049246139f6930c4af9816cbd30207b9b81661d78198003b
```

## 5. Canonical replay and required repair

Current canonical artifacts replay identically in normal and optimized
modes and byte-match stored output.  They contain no Python `assert` nodes.

```text
theorem ad24caeff22b0aade1203c7c05a8c186c9ed30ec759176342ff482919fa20729
script  f911eabd8ac9b235180cb413bfbfe43b6b38de60d0d2740a1467733e8f2340c9
output  325d695e53f860f239479d0ee6f9fbd46efd2dd068e06fe7d50c4bd822e3ddc7
```

Required script repair:

```python
candidate_fourth_root = [Fraction(1), unreduced[1] / 4]
gate(
    power(candidate_fourth_root, 4) != trim(unreduced),
    "unreduced hostile is square but not fourth power",
)
```

The existing degree-four gate may be retained separately, but it cannot
carry the non-fourth-power claim.  This replacement leaves the transcript
and gate count unchanged; refresh the script hash in theorem metadata and
the results index.

Recommended documentation strengthenings, not correctness repairs:

1. state explicitly that the `d=1` response block is empty and forces
   `b=1`;
2. replace or supplement the two isolated sharpness examples with the
   uniform family `(6)`;
3. record the characteristic-two hostile `(9)--(10)` as the exact reason
   for characteristic-zero scope;
4. state that divisor/response counts use the declared bound even when the
   top coefficient vanishes.

After the exact-companion repair, the theorem is suitable for
`PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED` status.  It remains
a formal response theorem only; none of the chart-entry or `JC(2)/DC(2)`
scope boundaries change.
