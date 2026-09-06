# All near-minimizers of the sharp signed-root stability quotient

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent audit](continuing3_20260906_stability_near_minimizers_audit.md)
accepts the full proof; producer 477 and referee 72 exact gates pass normally
and with optimization. The sharp constant is inherited from THM-4454.

Concurrent credit: incoming commit 28f41c846 proves the core classification
as [THM-4455 / three-atom minimizing-sequence rigidity](../../01-canon/theorems/THM-4455-three-atom-minimizing-sequence-rigidity.md).
Our complete classification, singular-boundary exclusions and signed-dust
correction are an independent reproof. The additions here are the moment
residual and moment-pair equivalences, explicit residual/distance estimates,
and the locally uniform entire-product equivalence. The parameterized dust
family realizes arbitrary finite separate first-moment limits as well as
divergence. Incoming THM-4455 also supplies a stronger negative-count bound
than the total-dust count consequence used below; it remains a useful sidecar.
Its d2 and Delta3 denote our d2 squared and d3 squared, respectively.

## 1. Inheritance, objects, and scope

The primary supplier is **THM-4454**,
`01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md`.
Its two regional envelopes and their strictness, not their finite control
tables, are used below. Its exact optimum is already PROVED. The recovered
predecessor `05-knowledge/results/signed_uniform_stability_empty_core_next_sep06.md`
classifies near-minimizers of the different quotient `J` and supplies the
one-atom denominator warning, mixed-sign dust, and entire-product closure
method. The quantitative predecessor
`05-knowledge/results/quantitative_stability_empty_core_morning_sep06.md`
already exhibits distinct two-atom and equal-three-atom limiting constants.
Its old open-optimum wording is superseded by THM-4454.

Anchor: classify the present global quotient. Niche: recover the actual
locally uniform polynomial-product limit. Wildcard: separate signed net
dust from its positive and negative first masses. The board retains top
positive roots, both quotient denominators, regional envelope zeros, the
square-weighted root measure, the first-moment defect, and the exponential
dust factor. The corrected near miss is that an unnormalized objective
going to zero does not identify the limit of its quotient. The least-used
sidecar is the signed first moment after square-mass compactification.

The map `r -> sum r_i^2 delta_(r_i)` preserves total square mass and the
third/fourth moments. It forgets the first moment carried by many tiny roots.
Integer root multiplicity restores three macroscopic atoms, while the
original normalization restores the signed exponential dust factor. No
mathematical implication for LRC(14) or general actual Laurent two-rung
separation is asserted.

## 2. Complete classification

For each index `n`, let `r^(n)` be a finite real list, of arbitrary varying
length, with

```
p1=sum r_i=1, p2=sum r_i^2=1, E=(1-p4)/2>0.
D=(5-8p3+3p4)/6, J=D/E,
c_*=(13-8sqrt2)/3,
K3=4sqrt3/[3(1+sqrt2)(1+sqrt3)].
```

Zero padding and both signs are permitted. Let `a>=b>=c>=0` denote its
three largest positive entries, inserting zeros if necessary, and set

```
d2^2=2-sqrt2(a+b), R=(J-c_*)/d2^2,
z=1/sqrt3,
d3^2=inf_permutations ||r-(z,z,z,0,...)||_2^2
    =2-2z(a+b+c),
M=sum_i r_i^2(r_i-z)^2=p4-2z p3+1/3.
```

The denominator `d2^2` is positive for every eligible finite list: zero
would force exactly two positive roots `1/sqrt2`, contradicting `p1=1`.

**Theorem.** The following are equivalent along the whole sequence:

1. `R -> K3`.
2. `a,b,c -> z` and the square mass outside these three roots tends to zero;
   equivalently `d3 -> 0`.
3. `M -> 0`.
4. `p3 -> 1/sqrt3` and `p4 -> 1/3`.
5. The actual products `G_n(s)=product_i(1+r_i^(n)s)` converge locally
   uniformly on the complex plane to
   `(1+s/sqrt3)^3 exp((1-sqrt3)s)`.

No assumption bounding `E` or `d2` away from zero is needed. In fact they
necessarily approach `1/3` and `sqrt(2-2sqrt(2/3))`, respectively.
For every `epsilon>0`, there is a dimension-independent `delta>0` such that
`R<K3+delta` implies `d3<epsilon`. This is a qualitative modulus, not a
new claimed sharp quantitative stability constant.

If `dust` means the entries outside the three selected positive roots,
then every near-minimizing sequence satisfies exactly

```
sum dust=1-a-b-c -> 1-sqrt3,
sum dust^2 -> 0.
```

Separate positive and negative dust sums are not determined. In particular,
the negative absolute dust sum tends to `sqrt3-1` **if and only if** the
positive dust sum tends to zero. Section 7 gives exact counterexamples,
including unbounded cancelling first masses.

## 3. The possible zeros of the inherited regional envelopes

Use `u=sqrt2`, `v=sqrt3`, `h=1/u`, and the constants of THM-4454:

```
A=2-u, B=u-1, gamma=3K3/(4u), C0=A+u gamma,
C=C0-gamma(a+b),
F_actual=1-C-p3+C p4=(3E/4)d2^2(R-K3).
```

Since `E<=1/2` and `d2^2<=2`, condition 1 implies `F_actual->0`.
This fact alone is not yet a quotient classification.

In the region `b<=z`, THM-4454 bounds `F_actual` below by the continuous
secant envelope on the compact domain `a>=b>=0`, `a^2+b^2<=1`.
Its only zeros are `(a,b)=(z,z)` and `(1,0)`, as explicitly proved there.

In the region `b>=z`, put `c0=sqrt(1-a^2-b^2)`, which is total remaining
square mass rather than the actual third root. The inherited nonnegative
three-root envelope has exactly two zeros:

```
(a,b,c0)=(z,z,z), or (h,h,0).
```

Here is the precise zero-set deduction from the inherited strict bounds.
For `0<c0<z`, its upper interval endpoint `a=b=x` has value
`c0^2(x-z)` times braces strictly greater than `1/2`; all factors are
positive. Its lower endpoint `b=z` is bounded below by
`gamma(1-a)(a-z)P(a)>0`. Strict concavity in `t=a+b` makes every interior
value positive as well. At `c0=0`, the lower endpoint is still positive;
the upper endpoint `(h,h,0)` is zero, and the chord bound is positive
everywhere else. At `c0=z` the interval collapses to `(z,z,z)`.

Thus any subsequence on which the actual objective tends to zero has a
further subsequence whose first two positive roots approach `(z,z)`,
`(1,0)`, or `(h,h)`. Only finite-dimensional compactness of the envelope
coordinates has been used; no first-moment compactness of the dust is assumed.

## 4. Both degenerate quotient alternatives are excluded

### One positive atom

If `a->1`, put `q=1-a^2=sum_(i>=2)r_i^2`. Then
`|sum tail^3|<=q^(3/2)` and `sum tail^4<=q^2`, independently of length
and signs. The exact one-atom quotient calculation recovered from the
uniform-stability predecessor gives `J->1`. For example dividing numerator
and denominator by `q` in

```
J=[5-8a^3+3a^4-8 sum tail^3+3 sum tail^4]
  /[3(1-a^4-sum tail^4)]
```

gives the limit `6/6`; the scalar numerator factors by `1-a`.
Consequently

```
R -> (1-c_*)/(2-sqrt2)=sqrt2-2/3 > K3.
```

This rules out the zero-energy alternative for a minimizing sequence.

### Two equal positive atoms

Suppose `a,b->h`. Put `q=1-a^2-b^2` and `w=a-b`. All other roots together
have square mass `q`, signed third moment `O(q^(3/2))`, and fourth moment
`O(q^2)`. Write `v0=w^2` and `t=a+b=sqrt(2-2q-v0)`. The exact top moments are

```
a^3+b^3=t(1-q+v0)/2,
a^4+b^4=[(1-q)^2+(2-2q-v0)v0]/2.
```

Taylor expansion at `(q,v0)=(0,0)` is uniform on `q,v0>=0`. For
`g=B-p3+A p4`, its linear part is

```
g=(7sqrt2/4-2)q+(2-11sqrt2/8)v0+o(q+v0).
```

The tail bounds make the error dimension-independent:
`q^(3/2)/(q+v0)<=sqrt(q)`, while the smooth top-moment remainder is
`O((q+v0)^2)`. Also `E->1/4` and
`d2^2=q+v0/2+O((q+v0)^2)`. Therefore

```
R=K_two+(K_dust-K_two) q/(q+w^2/2)+o(1),
K_two=(64-44sqrt2)/3,
K_dust=(28sqrt2-32)/3.
```

The mixing fraction lies in `[0,1]`, and `q+w^2>0` on every eligible
finite row. The constants satisfy `K_dust>K_two>1/2>K3`; the first
inequality reduces to `sqrt2>4/3`, the second to `sqrt2<125/88`, and
THM-4454 already gives `K3<4/9`. Hence `liminf R>=K_two>K3`.
This excludes the vanishing-distance alternative without assuming any
relation between imbalance and dust scale. It also records the full
dimension-independent local limiting profile at that boundary, rather
than extrapolating one particular two-atom family.

## 5. Recover the actual third atom

Sections 3 and 4 force `a,b->z` along the entire minimizing sequence.
The third-root conclusion still requires recovery of the information
discarded by the envelope. Put `f_C(x)=x-Cx^2`. Near `(a,b)=(z,z)`,
`C->(3-sqrt3)/2`, `f_C(b)->(sqrt3-1)/2>0`, and the derivative of `f_C`
on `[0,b]` is eventually at least `(2-sqrt3)/2>0`.

The exact secant loss after keeping the largest positive root is

```
S=a^2 f_C(a)+(1-a^2)f_C(b)-(p3-Cp4).
```

It is the difference of `F_actual` and the secant envelope. Both tend
to zero, so `S->0`, even if the secant envelope is slightly negative on
the `b>z` side. Its individual terms are nonnegative in this neighborhood.
Each negative root contributes at least `f_C(b)r_i^2`; each positive
root other than `a` contributes at least
`[(2-sqrt3)/2]r_i^2(b-r_i)`.

It follows that the total negative square mass tends to zero and that
`sum_(positive i other than a)r_i^2(b-r_i)->0`.
The positive square mass after deleting `a,b` tends to `1/3`. If `c` is
the largest remaining positive root, its contribution to the last bound
is at least `(b-c)` times that entire remaining positive square mass.
Thus `b-c->0`, and `c->z`. The remaining square mass is now
`1-a^2-b^2-c^2->0`, proving condition 2.

Conversely condition 2 gives `p3->3z^3=z` and `p4->3z^4=1/3`, because
the remaining moments of orders three and four are bounded by its square
mass to powers `3/2` and `2`. The limiting denominators are positive,
and direct substitution gives `R->K3`. This proves the central equivalence
without an unresolved boundary case.

## 6. Moment-measure and entire-product equivalents

Condition 4 implies condition 3 immediately. Conversely `M->0` localizes
the square-weighted root measure near `z`, with its integer multiplicities
retained. There is a useful explicit elementary bound. If `M<5/6144`,
consider the band `[7z/8,9z/8]`. Its outside square mass is at most `192M`.
It contains at most three roots since `4(7z/8)^2=49/48>1`, and
at least three since two roots can carry at most `2(9z/8)^2=27/32` while
outside mass is less than `5/32`. Thus it contains exactly three roots.
Their total squared displacement from `z` is at most `192M/49`, giving

```
d3^2 <= (9600/49)M,                     M<5/6144.
```

In the other direction the global elementary estimate
`M<=(1+z)^2 d3^2` follows by matching the three target entries and using
`|r_i|<=1` on the remainder. Hence conditions 2--4 are equivalent.

For condition 5, fix a complex disk `|s|<=R0`. Under condition 2, the
largest dust entry tends to zero. Eventually every dust factor has
`|r_i s|<=1/2`, so the logarithm near one satisfies

```
|sum_dust log(1+r_i s)-s sum_dust r_i|
    <= R0^2 sum_dust r_i^2 ->0.
```

The signed sum is exactly `1-a-b-c`; no bound on the separate positive
and negative first masses is used. Multiplication by the three macroscopic
factors proves the locally uniform limit in condition 5. Conversely local
uniform convergence gives the first four coefficients by Cauchy's formula;
Newton identities give condition 4. This recovers the actual product
observable, not only its first few moments.

The qualitative dimension-independent modulus follows by contradiction
from the sequential equivalence: otherwise rows with `R-K3->0` and
`d3>=epsilon` would violate condition 2.

## 7. Exact mixed dust: the missing first-moment coordinate

For `L>=2` and any real `c_L>=0`, set `T_L=3+c_L` and

```
Delta_L=3L^2+L(T_L^2-3+c_L^2)-c_L^2,
d_L=[L(T_L^2-3)-c_L^2]/[L T_L+sqrt(Delta_L)]
   =[L T_L-sqrt(Delta_L)]/(L-1),
S_L=T_L-d_L>0.
```

Take three raw roots equal to one, `L` roots `c_L/L`, and `L` roots
`-d_L/L`, and divide every root by `S_L`. The exact quadratic gives

```
S_L^2=3+(c_L^2+d_L^2)/L,
```

so the normalized list has `p1=p2=1` and positive energy. Positivity of
`S_L` follows from `Delta_L-T_L^2=(L-1)(3L+T_L^2+c_L^2)>0`.
If `c_L=o(sqrt L)`, then `S_L->sqrt3`, the three main roots tend to `z`,
and all other square mass vanishes. These are therefore genuine
near-minimizers of the sharp quotient.
Indeed `Delta_L/L^2->3` and `T_L/L->0` in the displayed radical formula
for `S_L`, which proves the asserted normalization limit directly.

For the simplest hostile `c_L=1`,

```
d_L=(13L-1)/(4L+sqrt(3L^2+14L-1)).
```

The positive dust sum tends to `1/sqrt3`, while the negative absolute
dust sum tends to `4/sqrt3-1`, strictly larger than `sqrt3-1`.
The first failed implication is replacing vanishing square mass by
vanishing positive first mass. The strongest survivor is the signed net
dust identity. More generally a fixed `c_L=c` gives limiting positive
dust mass `c/sqrt3` and negative mass `c/sqrt3+sqrt3-1`; taking
`c_L=L^(1/4)` makes both masses unbounded while preserving near-optimality.

If the dust has `N` nonzero entries and square mass `q3`, Cauchy--Schwarz
gives `N q3 >= (sum dust)^2 ->(sqrt3-1)^2`. Thus the number of dust roots
must diverge. Nothing forces equal dust magnitudes, one dust sign, a unique
multiplicity schedule, or a finite limiting total variation.

## 8. Verification and stopping point

The proof rests on the exact regional bounds of THM-4454, the uniform local
expansions and secant-loss recovery above. A standalone verifier checks
the new identities, boundary constants, exact normalized positive/mixed-dust
families, and original product coefficients. Its finite controls are not
used to infer the universal compactness or limiting statements.

The finite universe comprises twelve literal original-product controls
(`L=2,3,5`, `c=0,1,2,3`), eighteen bounded mixed-dust controls, six
unbounded-dust controls with `L=n^4,c=n`, eight rational one-atom controls,
and eighteen two-atom controls covering imbalance, dust, and mixed scales.
Every numerical inequality uses a rational outward square-root enclosure
with 256 binary fractional bits. Symbolic polynomial identities separately
certify normalization and the original doubled-product fourth coefficient.
The verifier does not import an inherited computation or infer a universal
envelope sign from a finite sample. Its **477 always-active gates** pass
normally and with Python optimization, with byte-identical LF outputs.

Reproduction (Python 3 with SymPy):

```text
python3 -B 04-computation/continuing3_20260906_stability_near_minimizers.py
python3 -B -O 04-computation/continuing3_20260906_stability_near_minimizers.py
```

Source and frozen output use this report's stem, with `.py` and `.out`
extensions. The optimized replay has identical output. Raw-LF SHA-256 hashes:

- source: `b5a2ce85df71837e5b73b3609feee0d9fcfaa14f425e22b81b563a7a3574c598`;
- each output: `3a0d5f4c830b62c55c48a9458f315f3ea601d68aeec40e9c634980ad2e6f8337`.

The proof and exact companions are filed without allocating a duplicate
theorem ID. The manifest preserves all frozen source and output identities.
The optimum remains the inherited closed value `K3`; the surviving new
question is an effective sharp modulus for this three-atom rigidity, not
reopening the constant or the separate dust-mass implication.
