# Five exact overlap profiles for inert odd-3-unit LRC tails

**Status: PROVED relative to the named inherited overlap/absorption inputs,
with FINITE-EXACT exhaustive heads; INDEPENDENTLY AUDITED.** Research date
2026-09-06. No theorem ID or external novelty/priority claim is requested.
LRC(14) remains **OPEN**.

The [independent referee](creative_20260906_inert_pareto_audit.md) accepts
the proof and rebuilds both finite heads using physical half-lift wall cells
and a separate arithmetic generator. All 9,505 independent gates pass.

## Inheritance and connection

The closest mechanism is [THM-4453 / inert-sum disjointness](../../01-canon/theorems/THM-4453-lrc14-inert-sum-five-ray-disjointness-and-dyadic-entry-closure.md):
arithmetic excludes five high-overlap ratios, giving the sufficient body mass
`4/91`. The less-used sidecar is the *joint* mass and largest-component profile
from [THM-4153 / third-tier Haar transfer](../../01-canon/theorems/THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer.md).
The hostile `(1,9)` retains inert sum but loses the 3-unit condition. The
corrected near miss is treating measure-zero walls as invisible to connected
components; [THM-4451 / strict component widths](../../01-canon/theorems/THM-4451-lrc14-dyadic-three-tail-strict-component-width-caps.md)
requires literal open topology. We keep every excluded endpoint.

The connection is the intersection of two existing carriers: arithmetic
ratio classes -> literal opposite-parity interval intersections -> their
two-coordinate Pareto envelope. It preserves upper bounds on both overlap
measure and component width. It forgets location, endpoint owners and body
phase. The restoration sidecar for the consumer is the actual compact body
safe set, its mass `M`, its largest component length `L`, and the tails'
common scale `g`. The conclusion is sufficient entry, never an equivalence.

Exact-statement, fraction, and synonym searches in current canon/results
found no prior five-profile envelope. Its ingredients are inherited.

## 1. Carrier and theorem

Let `0<p<q`, `gcd(p,q)=1`, `gcd(pq,6)=1`, and assume every prime dividing
`p+q` is `2 mod 3` with exponent at most two. There is **no height bound**.
Let `C_(p,q)` consist of `y in R/Z` for which both phases `y/2` and
`(y+1)/2` have at least one of the two speeds `p,q` at distance **strictly
less** than `1/14` from an integer. Write

    m(p,q)=mu(C_(p,q)),   b(p,q)=largest open-component length,

with `b=0` for an empty carrier.

**Theorem.** Every pair `(m(p,q),b(p,q))` is coordinatewise bounded above by
one row of this table; its five rows are precisely the maximal profiles:

| Attaining `(p,q)` | Measure `m` | Largest component `b` |
|---|---:|---:|
| `(7,13)` | `2/49` | `1/49` |
| `(19,25)` | `138/3325` | `37/3325` |
| `(5,41)` | `12/287` | `2/287` |
| `(5,53)` | `78/1855` | `2/371` |
| `(1,67)` | `20/469` | `2/469` |

In particular `m<=20/469`, equality only at `(1,67)`, and `b<=1/49`,
equality only at `(7,13)`. Every row also belongs to the actual THM-3818
atlas `p+q<=356`. Sharpness here concerns the **carrier profile**, not the
minimum sufficient body mass or the sharp LRC threshold.

## 2. Proof: two complementary finite reductions

The inherited exact formula (THM-4153, equation (6)) gives

    m(p,q)=2/49 + 2(B2(u_-)-B2(u_+))/(pq),
    u_-={1/2+(q-p)/14}, u_+={1/2+(q+p)/14},
    B2(u)=u^2-u+1/6.

Since `osc(B2)=1/4`, `m<=2/49+1/(2pq)`. Now

    20/469 - 2/49 = 6/3283,
    1/(2*274) < 6/3283.

Thus all `pq>=274` are strictly below the claimed mass maximum. The
complete admissible head `pq<=273` consists of the following 21 pairs:

| `p` | All eligible `q` in the product head |
|---:|---|
| 1 | 19,43,49,67,91,109,115,163,169,187,211,229,235,241 |
| 5 | 17,29,41,53 |
| 7 | 13,37 |
| 11 | 23 |

Substitution gives the unique maximum `20/469` at `(1,67)`. The frozen
output lists every exact mass. Completeness follows by iterating
`1<=p<=floor(sqrt(273))`, `p<q<=floor(273/p)`, with precisely the stated
coprimality, 3-unit and prime-factor filters. No other filter is inherited.

For the width, work on `0<=y<=1`. Each component is a positive intersection
of two intervals

    ((7k-1)/(7p),(7k+1)/(7p)),
    ((7l-1)/(7q),(7l+1)/(7q)),       k and l of opposite parity.

Indeed `py` near even/odd integers means respectively that the first/second
half-lift is bad for `p`. A single odd speed cannot make both half-lifts bad,
so opposite parities are necessary and sufficient. Teeth within either
family are separated by positive gaps. These intersections are therefore
the actual components, and `b<=2/(7q)`. The endpoints `0,1` are not in the
carrier, so no circular merger occurs there. For `q>=15` this is strictly
below `1/49`. The only admissible pair with `q<=13` is `(7,13)`, whose two
components are

    (1/7,8/49), (41/49,6/7).

Both have length `1/49`, proving the width assertion.

Finally, every `q>67` has `b<2/469` and, by the mass result, `m<20/469`.
It is strictly dominated by `(1,67)`. The remaining universe is exactly the
46 eligible pairs `p<q<=67`. Literal interval intersection gives the table
as its complete Pareto frontier. All five rows are attained and mutually
incomparable. Exhaustive domination of all 46 rows is checked with rational
arithmetic in the standalone certificate. This proves the all-height claim;
the height cutoffs, rather than a extrapolated scan, handle the infinite tail.

## 3. An inclusive mass-and-width entry criterion

Let `H` be any nonempty finite set of positive integers and

    G_H={y: ||hy||>=1/14 for all h in H},
    M=mu(G_H), L=largest connected-component length of G_H.

Let the tails be `a=gp,b=gq`, for positive odd `g` and an admissible pair
above. If, for **each** row `(m_i,b_i)` of the table,

    M>=m_i  OR  gL>=b_i,                              (1)

then `2H union {a,b}` has a common phase of clearance at least `1/14`.
It is sufficient to use either scalar gate `M>=20/469` or `gL>=1/49`.
Intermediate uses of (1) can succeed while both scalar gates fail.

**Proof.** Strict failure would give `G_H subset m_g^(-1)(C_(p,q))`.
The right-hand side is proper open; the left-hand side is compact.
Whenever `G_H` is nonempty, containment forces `M<m(p,q)`: equality would
make their open difference empty, contradicting a nontrivial clopen subset
of the circle. Every closed component lies inside one open component of
the inverse image, giving `gL<b(p,q)` whenever `L>0`. An empty or
zero-length `G_H` cannot satisfy (1), whose thresholds are all positive.
Choose a table row dominating the actual profile. The corresponding
disjunction (1) contradicts one of these strict inequalities. QED.

The [recovered physical eleven-body](creative_20260906_physical_control.md)
`H={1,3,4,5,6,7,8,9,10,11,13}` has `M=5939/140140`, `L=11/728`.
At `g=1` it satisfies (1) while `M<20/469` and `L<1/49`: the first four
profiles are paid by mass, the fifth by width. Thus the joint gain is
physically realized. All odd common tail scales preserve this certificate.
Its full closed geometry is independently rebuilt, including six isolated
safe points. This is a new certificate application to a recovered core,
not a first safety claim. As with THM-4453, actual decoder applications must
use the labelled pair component; the displayed body/tail split is a physical
doubled-body example and has a crossing decoder edge.

## 4. Improved original-body gates

[THM-4450 / absorbed-label hierarchy](../../01-canon/theorems/THM-4450-lrc14-absorbed-label-overlap-hierarchy-and-component-address-decoder.md)
supplies the same absorption estimates used by THM-4453. For a body of ten distinct positive integers
`C` and odd 3-units `r,a,b` with admissible reduced `(a,b)`:

* For `4C union {2r,a,b}`, take `H=2C union {r}` and use
  `mu(G_H)>=mu(G_C)/2`. The sufficient body gate improves from `8/91` to
  **`40/469`**, inclusively.
* For `2C union {2r,a,b}`, take `H=C union {r}` and use
  `mu(G_H)>=mu(G_C)-8/63`. The sufficient body gate improves from
  `20/117` to **`716/4221`**, inclusively.

If `r` already belongs to `C` in the second case, absorption has no loss
and the displayed weaker estimate is immediate; otherwise the cited
absorption theorem applies with its stated new-label hypothesis.

The first scalar reduction is `2/67` of the old threshold, about 2.99%.
There is no universal lower bound on the actual body's mass or component
length reaching these criteria. Arbitrary arithmetic tail pairs, actual
body entry below these gates, and LRC(14) remain **OPEN**.

## 5. Hostiles, symmetries and reproducibility

The 3-unit restriction is essential: `(1,9)` has inert sum ten but mass
`4/63>20/469`. The exponent cap is essential: `(5,11)` has sum `2^4` and
mass `18/385>20/469`. Dropping inertness also admits `(1,11)` with mass
`4/77`. These are boundary controls for inherited assumptions, not new
corrections to canon. The empty carrier `(1,3)` checks zero conventions.
Interchanging tails and reflecting phase preserve the profile. Common odd
dilation preserves measure and divides component width by its scale;
forgetting that scale would invalidate (1).

[Source](../../04-computation/creative_20260906_inert_pareto.py) and
[frozen output](creative_20260906_inert_pareto.out) use exact integers and
fractions with always-active gates. Universes: the complete 21-pair product
head, complete 46-pair width head, all 582 primitive odd pairs `p<q<=75`
for Bernoulli-versus-literal-geometry checks, and all 548 odd-3-unit ratios
in the `p+q<=356` inert atlas. The latter is a consumer check, not the
universe of the theorem.

    python3 -B 04-computation/creative_20260906_inert_pareto.py
    python3 -B -O 04-computation/creative_20260906_inert_pareto.py

The method cards used are **Retain the local profile until its global
counting weights are known**, **Search the statement before the method**,
and **Preserve the selected side, not only the walls**. The new object
changes the anchor's entry gates. For the walk, Laurent-observer and Smith
lanes, its lesson is only to retain coupled coordinates; it provides no
map preserving their separate target predicates.
