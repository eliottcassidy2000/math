# The period-fifteen CRT bridge: eleven sources on the formerly discarded axis

**Status: PROOF CANDIDATE; exact controls passed; independent audit pending.**
This bounded supplement leaves the frozen session74-clocks packet unchanged.
The bridge was proposed by the root session after the all-period golden-fibre
proof. It is an explicit map, not an inference from the equality `L_15=4*341`.
No Collatz convergence, integer-cycle classification, or historical novelty is
claimed.

## 1. Exact quotient and its four-element kernel

Let `O=Z[phi]`, `phi^2=phi+1`, and `R_n=O/(phi^n-1)`. Put

```
phi^15-1=376+610phi=2A,       A=188+305phi,
N(A)=188^2+188*305-305^2=-341.
```

Evaluation `phi -> 81 mod 341` is a ring homomorphism because
`81^2-81-1=19*341`. It kills `A`, since `188+305*81=73*341`, and is
surjective on constants. The principal ideal `(A)` has index `|N(A)|=341`,
so the induced surjection is an isomorphism `O/(A) -> Z/341`.
Modulo two, `A=phi`, a unit. Thus `(2)` and `(A)` are comaximal. Since
`t^2+t+1` is irreducible over `F_2`, the Chinese remainder theorem gives

```
R_15 = O/(2A) ~= F_4 x Z/341 ~= F_4 x F_11 x F_31,
phi -> (omega, 4, 19),       omega^2=omega+1.
```

The explicit map sends `a+bphi` to
`((a mod2,b mod2), a+4b mod11, a+19b mod31)`.
An optional principal-factor certificate is

```
phi^15-1 = -2 phi^9 (phi-4)(7-2phi),
N(phi-4)=11,                N(7-2phi)=31.
```

In particular, using `(phi-19)` as the norm-31 ideal would be wrong:
its norm is 341. The correct norm-31 factor is `(7-2phi)`.

The quotient

```
Theta:R_15 -> Z/341,        a+bphi -> a+81b mod341
```

forgets precisely the `F_4` coordinate. It is surjective and every fibre
has four elements. Its kernel is the four-element ideal `F_4 x {0} x {0}`,
not a four-to-one claim about an unrelated integer sequence.

## 2. Actual rational sources and integer clocks

The independently audited all-period golden-fibre theorem identifies all
`1364=L_15` rooted cyclic no-`11` words with `R_15`, hence with actual
rational 2-adic standard-Collatz points fixed by `C^15`. Here
`C(x)=3x+1` on odd inputs and `x/2` on even inputs; these are points of
period **dividing** fifteen. The exact CRT action has component orders
`3,5,15`, because `omega`, `4 mod11`, and `19 mod31` have those orders.

The kernel's four actual sources are exactly

```
{0,1,2,4}.
```

Proof: `R_3 ~= F_4`, and same-source repetition from three to fifteen has
multiplier `Q_(3,5)=sum_(j=0)^4 phi^(3j)`, whose CRT coordinates are
`(1,0,0)`. The period-three standard Collatz sources are zero and the
cycle `1 -> 4 -> 2 -> 1`. Thus forgetting the `F_4` component erases this
cycle's phase and zero as one fibre. The other 340 fibres each also contain
four rooted rational sources, without a claim that they are ordinary integers.

For the eleven period-five sources, same-source repetition gives

```
Q=Q_(5,3)=1+phi^5+phi^10=38+60phi,
Q -> (0,3,0),              Theta(Q)=124 mod341.
J:R_5=F_11 -> R_15,        x -> Qx,
Theta(J(x))=124x mod341.
```

Therefore their projected image is exactly `31(Z/341)`, the eleven-element
additive-order-11 axis, with the source-specific phase `x -> 3x` in its
`F_11` coordinate. All ten nonzero sources lie outside the unit sector;
the zero source does too. In particular the source `1/13`, whose base
label is three, projects to residue 31, not to a freely chosen axis label.

`J` is additive and `O`-linear but not multiplicative:
`Q^2=3Q`, hence `J(a)J(b)=3J(ab)`. The axis ideal has its own identity
`e_11=4Q` (CRT coordinate `(0,1,0)`), while the same-source lift is `Q`,
not `e_11`. Forgetting this distinction silently changes the phase.

**Clock guard.** Repetition preserves the underlying source and minimal
period, not free traversal cost: a literal length-five presentation is
traversed three times, with its odd/even counts and cost multiplied by
three; a length-three presentation is traversed five times. The chosen
standard clock in both cases is fifteen. This has no identification with
the Terras eleven-step negative cycle.

## 3. The real connection to the bracket-341 carrier

After projection to `Z/341`,

```
K=<4,81>=QR_11 x QR_31,     |K|=5*15=75.
```

For a direct generation proof, use generators `4` of `QR_11` and `19` of
`QR_31`. Their orders are five and fifteen, and `4=19^12 mod31`.
Thus the two generators have coordinates `(4,19^12)` and `(4,19)`.
An element `4^a81^b=1`, with `a mod5,b mod15`, forces
`a+b=0 mod5`, `12a+b=0 mod15`, hence `a=b=0`. The generated group
therefore has 75 elements and equals the product of square subgroups.

This is exactly the unit-square image used by the bracket-341 arithmetic
lane, now reached by an explicit golden quotient. It does **not** mean
that the eleven period-five sources are the unit subgroup: there are 300
golden sources over `K` (four per residue), whereas **none** of the eleven
period-five sources survive that unit filter. On their conductor-11 axis,
the `K` action factors through `QR_11` and has a 15-element kernel.

**Full-ring guard.** Scalar four is not a unit of `R_15`: its `F_4`
component is zero. Also the full-ring scalar 81 has `F_4` component one,
whereas `phi` has component `omega`. Thus `phi=81` and the group
`<4,81>` are statements **after this quotient**, not identities of scalar
or clock actions on the full `R_15`. On the eleven-cell axis, `phi` and
scalar four do agree, since only the `F_11` component remains.

## 4. The discarded axis has exactly twelve units of bracket correction

For the bracket statistic at `q=341`, the parameter change `z=2n+1`
permutes all residues. Its exact residue term is

```
T_341 = sum_(z mod341) (-z^2 mod341),
```

with least nonnegative residues, not just support. Restricting the uniform
source to the eleven lifted sources gives every `z=124x`, `x in F_11`,
exactly once. Squaring yields

```
residue:       0   31   93   124   155   279
multiplicity:  1    2    2     2     2     2.
```

The zero source contributes zero tax. Thus the positive-residue weighted
sum is 1364, the negative-square weighted sum is 2046, and this axis's
term in the bracket correction is

```
2*2046/341 = 12.
```

The separate additive-order-31 axis has negative-square sum 4092 and
correction 24. Dropping both axes loses exactly 36, agreeing with the
independently recovered bracket-341 hostile. This is a source-weighted
bridge to an actual previously discarded contribution, not a resemblance
of prime factors.

The uniform axis sum is insensitive to permuting all eleven source labels,
but individual sources and nonuniform weights are not. For base source
label `x=1`, the true phase gives `-(124x)^2 mod341=310`, while the naive
replacement by `31x` gives 62. Arbitrary source weights must follow
`x ->124x ->(124x)^2` and be summed on equal images; they cannot be replaced
by uniform mass on `K`, on all 341 residues, or on all 1364 golden sources.
No bijection between the bracket parameters and all golden words is claimed.

## 5. Exact controls and scope

`verify.py` is portable, standard-library-only, and uses exact integers and
Fractions with explicit exceptions. Run it normally and under `python -O`.
The retained receipts agree byte-for-byte, SHA256
`1c5cb4453362dba92754d9b0cfeda495fa4505eb55cfefef7c1a6d86c9aa5693`.
There are **25,996 checks**, including every one of the 1364 rooted words,
every actual prescribed rational parity and return, all 341 four-element
fibres, the kernel's actual sources, both same-source lifts, the 75-element
unit subgroup, the axis multiplicities and weighted sums, and hostile
controls for dropping the `F_4` component or the source-specific phase.
The structural CRT and source-lift statements are proved above for these
fixed moduli; there is no unbounded computational or Collatz claim.
