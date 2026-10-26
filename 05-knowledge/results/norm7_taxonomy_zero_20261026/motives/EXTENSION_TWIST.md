# The missing minus-one twist in the prime-power cone sums

**PROVED transport identities over every odd finite field; FINITE-EXACT
receipts on the declared fields.** This extends the six raw cone sums of
[HYP-9391](../../../hypotheses/HYP-9391-the-seven-vertex-paley-census-and-four-new-weight-five-motives.md).
It does not identify the four proposed pure motives, their ranks, or conductors.

The inherited evaluator is explicitly restricted to primes `q≡3 mod4`.
Its constant orientation sign is correct there. Extending that sign unchanged
to F_(p²) would compute a different character sum. The missing coordinate is
`ε_q=χ_q(−1)`, a genuine quadratic twist rather than a numerical convention.

## 1. Local cone identity

Put `L_q(λ)=Σ_x χ_q(x(x−1)(x−λ))`. For distinct a,b,d in F_q and a cubic
vertex with k of its three neighbors having larger vertex labels, its literal
oriented edge-product sum is

    ε_q^k χ_q(b−a) L_q((d−a)/(b−a)).                    (E1)

Indeed, each edge from the cone vertex to a larger-labelled neighbor has
factor χ(neighbor−x)=εχ(x−neighbor). Substitute x=a+(b−a)y in the resulting
cubic: its scalar is (b−a)³, with character χ(b−a). This proves E1 for every
odd field. No tournament is required to define this character sum.

The six pinned graphs H02,H04,H05,H06,H09,H18 have independent cone triples
and orientation exponents `1,0,1,0,0,0`. Their graph edge lists are inherited
data pinned in the receipt. Define C_i(q) by summing their literal oriented
edge products over distinct four-point bases fixed at 0,1,s,t, and over all
three cone coordinates independently. Thus C is a cone sum, not the fully
injective seven-vertex graph sum; merge corrections remain separate.

With `s,t∉{0,1}`, `s≠t`, direct substitution into E1 gives

    C02 = ε Σ χ(s)χ(t)χ(t−1) L(s)² L(t/s),
    C05 = ε Σ χ(s−1)χ(t)χ(t−1) L(s)² L(t/s),
    C09 =   Σ χ(s)χ(s−1)χ(t) L(s)L(t)L(t/s),
    C18 =   Σ χ(s−1)χ(t−1) L(s)L(t)L(t/s).             (E2)

On q≡3 mod4 these are exactly the inherited Y02,Y05,Y09,Y18. On q≡1 mod4,
the first two have the opposite sign from the displayed fixed-minus formulas
in HYP-9391. This is an extension boundary, not an error within that hypothesis's
stated prime class.

## 2. The paired reductions extend with an arithmetic sign

For **every odd prime power**, the correct relations are

    C04(q)=ε_q C02(q),       C06(q)=C05(q).             (E3)

Here is a proof retaining the field character. Set

    A(u)=Σ_(s≠0,1) χ(s)χ(u−s)L(s)²,
    B(u)=Σ_(s≠0,1) χ(s)χ(s−1)χ(u−s)L(s)².

Using t=su followed by s↦1/s in E2, and allowing zero character factors
to impose the excluded diagonals, gives

    C02=ε Σ_(u≠0,1) χ(u)L(u)A(u),
    C04=  Σ_(u≠0,1) χ(u)χ(u−1)L(u)A(u),
    C05=  Σ_(u≠0,1) χ(u)L(u)B(u),
    C06=  Σ_(u≠0,1) χ(u)χ(u−1)L(u)B(u).

The elementary Legendre substitutions give
`L(1/u)=χ(u)L(u)` and `L(u/(u−1))=εχ(u−1)L(u)`.
In A or B substitute s↦s/(s−1); both transform by the same factor
`εχ(u−1)`. Hence L(u)A(u) and L(u)B(u) are invariant under u↦u/(u−1).
On χ(u−1)=−1, this involution flips χ(u), preserves that subset, and has
no fixed point (its only allowed fixed point is u=2, where χ(u−1)=1).
The subset therefore contributes zero. Replacing χ(u−1) by 1 in the two
weighted sums proves E3. The inherited prime-class proof is recovered at ε=−1.

## 3. Exact extension data and the Euler-factor gate

For p≡3 mod4, `ε_(p^r)=(-1)^r`: a quadratic extension changes the sign.
For example the independently reconstructed finite-field calculations give:

| q | ε | C02 | C04 | C05 | C06 | C09 | C18 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | −1 | −128 | 128 | 128 | 128 | 0 | 0 |
| 49 | 1 | 2928 | 2928 | 11760 | 11760 | 1904 | −2320 |
| 27 | −1 | −8064 | 8064 | 12672 | 12672 | 11520 | 4992 |

The arithmetic signs must be fixed before Newton reconstruction. If a trace
sequence a_r of a specified representation were known, its quadratic Euler
coefficient would be `(a_1²−a_2)/2`. For the raw C02 sequence at p=7, using
the wrong extension sign would change that expression from 6728 to 9656.
Neither value is declared to be a coefficient of an identified pure motive:
the raw cone sum includes pieces that still need exact separation.

Moreover even a correct a_1,a_2 does not establish rank four, a symplectic
pairing, an integral model, or conductor exponents. Those are additional
hypotheses of the [certificate lane](../certificates/README.md), not outputs
of a fitted trace distribution. HYP-9391 and its MISTAKE-671 correction keep
the rank and identification OPEN.

## 4. Certified zero and undefined have different result types

The Paley orientation `x→y iff χ(y−x)=1` is a tournament precisely when
ε=−1. If ε=1 it is symmetric, so this tournament census is undefined.
The checker returns `count: null` together with that reason, rather than zero.
The raw cone sums remain defined in both cases.

An adjacency-based recursion independently certifies `TT6(P_27)=0`.
It counts each transitive subset by its unique source. In F9 the Paley
**tournament** input instead fails its domain gate. This is decidable domain
validation, not the computability-theoretic nowhere-defined function.

## 5. Reproduction and controls

`python3 -B check_extensions.py` prints the receipt; repeat with `-O`.
Universe: p=3,5,7,11,13,17,19, each prime field and its quadratic extension,
plus F27. Finite fields are independently implemented by polynomial quotients;
no inherited producer code is imported. Quadratic moduli use the least
nonsquare; F27 uses x³−x−1. Positive and hostile controls include:

- Character multiplicativity, Frobenius powers, and inverse tables.
- Literal cubic edge-product sums versus E1 for every declared field q≤27.
- Boundary fibres `L(0)=−ε`, `L(1)=−1` and smooth-fibre Hasse bounds.
- For λ∈F_p\{0,1}, independent elliptic recurrence
  `L_(p²)(λ)=2p−L_p(λ)²`.
- Both q classes, both extension parities, and distinct zero/undefined cases.

The general proofs are E1–E3. The receipt supports their implementation and
applications; it is not evidence that an undeclared Euler factor or conductor
has been computed.
