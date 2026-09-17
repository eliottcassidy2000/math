# Arithmetic braids II: floor reciprocity distinguishes squarefree from prime

**Status: PROVED elementary all-modulus identities + FINITE-EXACT controls.**
The user confirmed that all three sums contain floor functions. No new
canon ID or literature-priority claim is made. These identities do not prove
Collatz, twin primes, Goldbach, or LRC(14).

## Inheritance and working board

The closest proved repository mechanism is the retained boundary in
[THM-2422, swap-fixed operation fibres](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md):
averaging or identifying paired objects requires keeping their exceptional
fibre. Its canonical hostile is `4`, whose repeated-factor diagonal cannot
be deleted in a primality test. The previous session's
[divisor balance note](arithmetic_braids_20260917_divisors.md) retains the
whole prime exponent profile; that profile will again distinguish support
from multiplicity. The corrected near miss here is the assumption that a
cube map modulo a prime must permute the nonzero residues. Modulo 7 its
image is only `{1,6}`; the relevant symmetry is antipodal, not bijective.
The least-used sidecar is the exact count of **lattice points on the boundary**.

The routed [THM-4068, squarefree Stern-packet balance](../../01-canon/theorems/THM-4068-squarefree-stern-packet-and-tournament-apex-balance.md)
also uses CRT, while explicitly warning that a chosen representative's
parity need not factor. No Stern imbalance estimate is imported here.
[THM-2474, squarefree first-collision saturation](../../01-canon/theorems/THM-2474-squarefree-first-collision-primitive-character-saturation.md)
has a different root-packet carrier; sharing the word squarefree does not
transfer its Fourier conclusion to these sums.

| Lane | Object and representation | Question | Hostile / retained coordinate |
|---|---|---|---|
| Anchor | Cubic and inverse cubic lattice sums | Exact scope of the three formulas | `n=4,6`; retain boundary points |
| Niche | Residue ring `Z/nZ` | Nilpotents versus zero divisors | Squarefree composite `6` |
| Wildcard | All odd powers and product arities | Which mechanism generalizes? | Degree 1, even degree 2, nonpermuting cubes mod 7 |
| Bridge | Prime exponent profile | Which defect sees repeated primes? | Compare `p`, `p^3`, `p^2 q r` |

## 1. The exact repair and the strongest domains

For an integer `n>=2`, define

```text
F3(n) = sum_(k=1)^(n-1) floor(k^3/n),
H3(n) = sum_(k=1)^((n-1)(n-2)) floor((nk)^(1/3)),
P2(n) = sum_(i,j=1)^(n-1) floor(ij/n).
```

An empty sum is zero. Let

```text
R3(n) = #{0<=k<n : k^3=0 mod n},
Z2(n) = #{1<=i,j<n : ij=0 mod n}.
```

The complete formulas are

```text
F3(n) = (n-2)(n-1)(n+1)/4 + (R3(n)-1)/2,               (1)
H3(n) = (3n-5)(n-2)(n-1)/4 + (R3(n)-1)/2,              (2)
P2(n) = (n-2)(n-1)^2/4 + Z2(n)/2.                      (3)
```

Consequently the first two proposed identities hold **if and only if `n`
is squarefree**, whereas the third holds **if and only if `n` is prime**.
Each is true for primes; the first two are stronger than prime identities.
Even moduli are included. The lower bound `n>=2` matters: at `n=1` all
three empty-sum formulas hold but 1 is not prime.

The boundary counts have closed multiplicative formulas. If
`n=product_(p^a||n) p^a`, then

```text
R3(n) = product_(p^a||n) p^(a-ceil(a/3)),                (4)
Z2(n) = sum_(i=1)^(n-1)(gcd(i,n)-1)
      = product_(p^a||n) p^(a-1)((a+1)p-a)-2n+1.        (5)
```

The corrections need not be integers by themselves. At `n=4` the cube
and product baselines are half-integers and the correction is `1/2`.
Dropping this even-modulus boundary would hide the first counterexample.

| `n` | `F3(n)` | `H3(n)` | `P2(n)` | Cube/root defect | Product defect |
|---:|---:|---:|---:|---:|---:|
|2|0|0|0|0|0|
|3|2|2|1|0|0|
|4|8|11|5|1/2|1/2|
|6|35|65|27|0|2|
|8|96|201|76|3/2|5/2|
|9|141|309|114|1|2|
|12|358|853|311|1/2|17/2|
|27|4554|12354|4239|4|14|
|30|6293|17255|5925|0|38|

For the original unfloored expressions, `p=3` already gives `3` instead
of `2` in the cube sum and `3` instead of `1` in the double product sum;
the two cube roots sum to strictly more than 2. Thus inserting the floors
is mathematically essential, not a cosmetic change of notation.

## 2. The common mechanism: antipodal residue pairing

Here is the reusable elementary lemma. Let `X` be a finite set with an
involution `iota`, and let `h:X->Z` satisfy

```text
h(iota(x)) = -h(x) mod n.
```

If `Z` is the number of elements with `n|h(x)`, then

```text
sum_(x in X) floor(h(x)/n)
 = (sum_(x in X) h(x))/n - (|X|-Z)/2.                   (6)
```

**Proof.** Write `[h(x)]_n` for the representative in `{0,...,n-1}`.
The paired residues have sum `n` unless both are zero, when their sum is
zero. Summing over the involution gives
`2 sum [h(x)]_n=n(|X|-Z)`. Subtract this from `sum h(x)` and divide by
`n`. Fixed points cause no problem: a fixed nonzero residue equals `n/2`
and has already been counted correctly. This covers even `n`. QED.

For cubes, use `X={1,...,n-1}`, `iota(k)=n-k`, and `h(k)=k^3`.
The exceptional count is `R3(n)-1`. Since
`sum k^3=n^2(n-1)^2/4`, (6) gives (1).

For products, use `X={1,...,n-1}^2`,
`iota(i,j)=(n-i,j)`, and `h(i,j)=ij`. The raw sum divided by `n` is
`n(n-1)^2/4`. Substitution gives (3).

Thus the common reason is an involution and its zero-residue set. There is
no need for powers to permute residues, no probabilistic independence, and
no tournament orientation. The product-zero relation is symmetric and
includes diagonal pairs when present; forcing it into a tournament would
discard the boundary being measured.

## 3. The cube-root identity is a rectangle reciprocity law

Put `M=(n-1)(n-2)`. The exact endpoint identity is

```text
(n-1)^3 = n M +(n-1).                                  (7)
```

Consider the rectangle of integer lattice points
`1<=j<=n-1`, `1<=k<=M`. The sum `F3(n)` counts its points with
`nk<=j^3`; the sum `H3(n)` counts its points with `j^3<=nk`. These two
regions cover the rectangle and meet exactly on `nk=j^3`.

The intersection contains `R3(n)-1` points. Indeed each nonzero cube-zero
residue `j` supplies the unique integer `k=j^3/n`, and the endpoint `j=n-1`
is never a zero residue. Therefore

```text
F3(n)+H3(n)=(n-1)M+R3(n)-1.                            (8)
```

Inserting (1) proves (2). For `n=2`, `M=0`, the rectangle is empty and
the same proof and formulas remain valid.

The equal correction in (1) and (2) is thus forced. The inverse graph does
not create a third unrelated prime identity: it exchanges the two sides
of the same lattice boundary. A strict-versus-weak inequality convention
would move the boundary correction from one side to the other.

## 4. Composite arithmetic: reduced rings versus fields

By the Chinese remainder theorem, `k^3=0 mod n` is equivalent to
`v_p(k)>=ceil(a/3)` for each `p^a||n`. There are
`p^(a-ceil(a/3))` possibilities modulo `p^a`, which proves (4).

This is 1 exactly when all prime exponents are 1. Equivalently,
`Z/nZ` has no nonzero nilpotents exactly when `n` is squarefree. There is
even a direct hostile construction: if `p^2|n`, then `x=n/p` is nonzero
modulo `n` and satisfies `x^2=0`, hence `x^3=0`, modulo `n`.

At fixed nonzero `i`, the congruence `ij=0 mod n` has `gcd(i,n)`
solutions modulo `n`, one of them `j=0`. Summing proves the first
identity in (5). If `n` is prime every term is zero. If `n` is composite,
any proper factorization `n=ab` gives a nonzero pair `(a,b)`, so `Z2>0`.
This proves the primality equivalence without analytic number theory.

For completeness, including zero in both coordinates gives

```text
#{(i,j) mod n : ij=0} = sum_(i=1)^n gcd(i,n).
```

This count is multiplicative by CRT. At `p^a`, the valuation strata
`v_p(i)=s<a` each contribute `(p-1)p^(a-1)` to the count, and `i=0`
contributes `p^a`; their sum is `p^(a-1)((a+1)p-a)`. Removing the
`2n-1` pairs with at least one zero coordinate proves the last identity
in (5). The equivalent divisor convolution is

```text
sum_(i=1)^n gcd(i,n) = sum_(d|n) d*phi(n/d).
```

The separator `n=6` is intrinsic: `Z/6Z` is a product of two fields, so
it is reduced, but `2*3=0` supplies zero divisors. Accordingly, the cube
and cube-root identities hold while the product identity has positive
defect. The three sums distinguish **reduced ring** from **field**.

## 5. Every odd power and every inverse lattice rectangle

For an integer degree `d>=1`, define

```text
Rd(n) = #{0<=k<n : k^d=0 mod n}
      = product_(p^a||n) p^(a-ceil(a/d)),
Fd(n) = sum_(k=1)^(n-1) floor(k^d/n),
Sd(n) = sum_(k=1)^(n-1) k^d,
Md(n) = floor((n-1)^d/n),
Hd(n) = sum_(k=1)^Md(n) floor((nk)^(1/d)).
```

The rectangle argument proves, for **every** integer `d>=1`,

```text
Fd(n)+Hd(n)=(n-1)Md(n)+Rd(n)-1.                         (9)
```

For odd `d`, antipodal pairing additionally proves

```text
Bd(n) = Sd(n)/n-(n-1)/2,
Fd(n) = Bd(n)+(Rd(n)-1)/2,
Hd(n) = (n-1)Md(n)-Bd(n)+(Rd(n)-1)/2,
Md(n) = ((n-1)^d-(n-1))/n.                              (10)
```

For **each fixed odd `d>=3`**, either zero-defect identity in (10) is
equivalent to squarefreeness of `n`. The degree-one exception matters:
`R1(n)=1` and `F1(n)=H1(n)=0` for every modulus. Even-degree antipodal
pairing is false; at `n=3,d=2`, the actual power-floor sum is 1 whereas
`S2(3)/3-(3-1)/2=2/3`. Reciprocity (9) still survives at even degree.

The proof also generalizes to odd integer polynomials in (6), but their
zero set need not test squarefreeness. For example `x^3-x` already has
the distinct roots `0,1,-1` at an odd prime. The power map's particular
zero fibre is indispensable to the ring interpretation.

For product arity `r>=2`, the same involution gives

```text
sum_(1<=i1,...,ir<n) floor(i1*...*ir/n)
 = (n-1)^r*(n^(r-1)-2^(r-1))/2^r + Zr(n)/2,             (11)
```

where `Zr` counts nonzero-coordinate tuples with zero product. It vanishes
exactly for primes: a composite factor pair can be padded by `r-2` copies
of 1. Thus this product extension continues to detect prime moduli, not
merely squarefree ones.

## 6. What the defects retain about the radical

The nilradical of `Z/nZ` is the set of multiples of `rad(n)`, of size
`n/rad(n)`. The finite-degree count `Rd` is generally only part of it:

```text
Rd(n) <= n/rad(n),
equality iff d >= max_(p^a||n) a.                        (12)
```

For example at `n=p^4`, `R3=p^2` but the whole nilradical has size
`p^3`. This hostile prevents calling the cube defect the entire nilradical
size without an exponent bound.

For the three exponent patterns previously classified in the divisor
balance note, all exponents are at most 3, so the cubic count does saturate:

| Balance pattern | `R3(n)=n/rad(n)` | Defect in each cube/root sum |
|---|---:|---:|
| `p` | 1 | 0 |
| `p^3` | `p^2` | `(p^2-1)/2` |
| `p^2 q r`, distinct primes | `p` | `(p-1)/2` |

This is a real link through exponent profiles. It does not identify
`F=S+U` with a zero-defect floor identity: two of its three pattern types
have strictly positive floor defects.

The scalar defect still forgets substantial information. Every squarefree
modulus, including all primes and arbitrarily complicated products of
distinct primes, has the same zero cube/root defect. Consequently it
cannot determine the modulus, factor support, arithmetic height, or
Collatz itinerary. A squarefree-density theorem describes the frequency of
zero defects; it does not by itself supply a decreasing dynamical height.

There are inexpensive dynamical hostiles as well. For the accelerated odd
map `T(n)=(3n+1)/2^v2(3n+1)`,

```text
33 -> 25:  cube/root defect 0 -> 2,
 9 ->  7:  cube/root defect 1 -> 0.                     (13)
```

An exact scan of all odd starts `3<=n<=500` confirms these are respectively
the first loss and first gain of squarefreeness. Thus the scalar defect is
not even monotone under one step of the relevant operation.

There is an all-length obstruction to extracting descent from the first
two identities. The companion
[squarefree symmetry note, Section 3](arithmetic_braids2_20260917_squarefree_symmetry.md)
proves that for every fixed `L>=1` there are infinitely many positive
integers `q` for which all

```text
n_j=3^j*2^(L+1-j)*q-1,       0<=j<=L,                  (14)
```

are squarefree. They form `L` consecutive strictly growing Collatz steps
with halving exponent one. Thus the cube/root identities can hold at
every vertex of an arbitrarily long growing trajectory segment. This is
a joint realization theorem, stronger than unrelated squarefree densities
at the individual positions. It still makes no assertion that one
infinite trajectory is squarefree, or that a squarefree realization can
also be chosen among separately constructed terminating completions.

## 7. Verification, transfer contracts, and stopping boundary

| Source -> target | Map | Preserved predicate | Lost data / required sidecar |
|---|---|---|---|
| Cubic graph -> inverse root graph | Exchange lattice axes in the fixed rectangle | Counts off the boundary | Equality points must be counted separately |
| Odd powers -> floor defect | Antipodal residue involution | Number of zero residues | Entire nonzero residue distribution |
| Prime exponent profile -> cube defect | `(4)` | Presence of repeated primes | Exact exponents may exceed cubic resolution |
| Multiplication table -> product defect | Zero-product pair count | Primality at zero defect | Pair locations; count alone is not a factorization |

Reproduction:

```text
python 04-computation/experiments/arithmetic_braids2_20260917_floor_reciprocity.py
python -O 04-computation/experiments/arithmetic_braids2_20260917_floor_reciprocity.py
```

The script uses exact integers and rational fractions; cube roots and other
roots are evaluated individually by certified integer binary search.
Checks remain active under optimized Python. It writes
`04-computation/experiments/arithmetic_braids2_20260917_floor_reciprocity.json`,
including the source hash, explicit universes, examples, and hostile controls.

The finite controls cover 2,994 odd-power parameter pairs and 751,494
residues, 160 inverse-rectangle pairs with 743,638 individually evaluated
root terms, 1,113,775 multiplication-table entries, and 129,974 higher-arity
tuples. Degree 1, even degree 2, `n=2`, `n=4`, squarefree composite `n=6`,
and nonpermuting cubes modulo 7 are retained as controls.

The new reusable mechanism is the exact boundary defect. The next useful
question is not whether a zero cube defect proves termination; it plainly
forgets the trajectory. It is whether a retained **joint** observable of
valuation depth, residue address, and boundary incidence has a proven law
under the particular dynamical operation being studied. No such global
Collatz law is claimed here.
