# Near-square targets and consecutive counts inside the actual core grammar

**PROVED fixed-chart identities and an elementary dyadic-square classification;
FINITE-EXACT controls; all-even-square coverage and Chan–Pak remain OPEN.**
No prime-infinitude assertion or external novelty claim is made.

This continuation connects the user's Landau/Fermat input to actual
tournament constructors. Primality is a later predicate on the constructed
integer; it does not repair a missing Hamilton-count witness.

## 1. Recover the operations, and distinguish two functions called g

The closest existing source is **THM-5115**,
[synchronized-endpoint-gap-pairs-and-exact-repair-rigidity](../../../../01-canon/theorems/THM-5115-synchronized-endpoint-gap-pairs-and-exact-repair-rigidity.md).
It gives actual strong same-order pairs at every prescribed even difference,
including a pure-count pair with identical three rooted sidecars. Arbitrary
common O/Z tails preserve that difference, but their lower endpoint is still
constrained by an exact scalar equation. Its
[complete inverse](../../session244_even_gap_anchoring_20261110/selector/ROUND4.md)
already distinguishes exact anchoring from prescribed gaps and from modular
surjectivity. Universal even differences are inherited from Sessions199/221.

There is a potentially dangerous notation collision. Let f denote the
Session199 canonical-source function (called g in those papers):

    f(1)=f(2)=1,
    f(2m)=f(m),
    f(2m-1)=f(m)+m-(m mod2), m>=2.

Let g denote the current one-apex compiler's function. Then exactly

    f(n)=g(n-1), n>=1.

The even and odd recurrences prove this by induction. For example f(2)=1,
whereas g(2)=3; unshifted reuse would give incorrect count formulas.
The source is [Session199, section2](../../../../session199/hidden-structure/PROOF.md).

The old pure-count adjacent pair, with selected difference2 and arbitrary
binary tail, consequently has lower endpoint

    62M+9f(M)-14 = 9g(A)+62A+48, A=M-1>=0,

and upper endpoint2 more. This is a real controlled family, not arbitrary
anchoring. The old theorem excludes61 and173 from its full two-index family
at every order. It already records why even all fixed-modulus residues
would not overcome this absolute-value restriction.

Other recovered square sources have different types. The
[slack-spectrum square-minus-diagonal map](../../split_word_metagraph_spectra_20261021/session23-slack-spectrum/README.md)
has h=1+N^2 for an orbit-count h=H/|Aut| and N a power of two. With only
singleton atoms it recovers ordinary reversed-arc Fermat-type counts;
with triangle atoms the raw Hamilton count carries an extra automorphism
factor. This is not a constructor for all even n^2+1.
The [Gaussian square-norm transport](../../collatz_square_norm_transfer_20261007f.md)
is a marked Collatz guard coordinate, and
[quadratic Hamilton memory](../../session172_memory_dynamics_20261103/squares/NOTE.md)
is a one-vertex response theorem for undirected graphs. Neither supplies a
tournament Hamilton count simply by sharing a quadratic formula.

## 2. The first near-square failure, and two genuine repairs

Every Fermat number is in the one-apex image by g(2^m)=2^m+1. This does
not extend to every even square plus one. The first missing value is

    101=10^2+1,

followed by197=14^2+1 and785=28^2+1. The finite-target exclusion is complete
over arbitrary address lengths: repeatedly delete trailing binary ones
using g(2q+1)=g(q). For the resulting even address2q,
`g(2q)=g(q)+q+(q mod2)>=q+1`, so g(2q)<=X implies q<X.
Enumerating0<=q<X therefore finds every one-apex value through X.

For X=1,000,001 the [receipt](near_square_receipt.json) checks all500
targets n^2+1 with positive even n<=1000. Exactly231 are absent from this
one-apex image. This is a finite exact result about a restricted grammar,
not a spectrum omission. The small positive controls include5,17,37,65.

The current five-core bank repairs101 already with C3, out-cut O={0},
left address0 and right address16. The literal9-vertex output is

    (250,508,505,240,224,192,128,256,121).

An inherited different repair is the neutral C=6 source word
`11,01,00,00,10` in the [bounded-response atlas](../../ready_return_pentagon_20261103/ATLAS.md).
It has8 vertices, H101 and rooted profile
`(s,q,F-H,delta)=(17,17,6,6)`. The checker independently constructs and
recounts both witnesses. Its ready response differs from the C=2 five
families; this is an actual boundary enlargement.

There is a sharp all-size positive/negative statement for the first repair.
For C3 with O={0}, the positive coefficients are

    (h,d,L,i,o,s,e)=(5,3,12,6,3,2,1).

Consequently at left address0 and right address B,

    H=3g(B)+3B+2+(B mod2).

At B=2^m, m>=1, this becomes H=6*2^m+5. It is of the form n^2+1
for positive integer n **exactly** in the three cases

    (m,n,H)=(1,4,17), (4,10,101), (5,14,197).

Proof: n is even; write n=2t. The equation is
`t^2-1=3*2^(m-1)`. At m1 it gives t2; m2 gives the impossible square7.
For m>=3, t is odd and a=(t-1)/2 satisfies
`a(a+1)=3*2^(m-3)`. Since a and a+1 are coprime, their odd member is1
or3. The first possibility gives no factor3. The second forces the pair
{2,3} or {3,4}, giving m4 or m5. This proves the exhaustion.

Thus the same C3 dyadic chart repairs the first two failed even-square
targets101 and197, yet has only three near-square values on this dyadic
slice. It provides neither an infinite prime family nor all near squares.

## 3. Three consecutive odd counts from the two H9 cores

Let S5a=(2,4,9,19,7), S5b=(2,4,25,3,11), and k>=0. The following
three explicit one-apex charts give

    9g(k)+82k+62,
    9g(k)+82k+64,
    9g(k)+82k+66.

| Output constant | Core | Apex out-mask O | Left address A | Right address B |
|---:|---|---:|---|---|
|62|S5a|5|2k+1|0|
|64|S5b|5|2k+1|0|
|66|S5a|14|0|2k+1|

These are actual consecutive odd counts, on equal vertex orders. The
middle/lower channel is the odd-address specialization of the concurrent
alternative lane's exact minus-two switch; that lane retains literal
transport data. Here the additional third chart extends it to an interval
of three odd values. No vertex correspondence between all three charts is
inferred merely from equality of their orders.

The proof is direct substitution into the positive compiler. On the left
odd ray, g(2k+1)=g(k), so

    Phi(2k+1,0)=d*g(k)+2i*k+(h-d+i+s).

On the right odd ray it is

    Phi(0,2k+1)=d*g(k)+2o*k+(h-d+o+e).

The three rows have d9, the chosen i or o equal41, and constants62,64,66.
The exact seven-coordinate vectors are in the receipt. Direct child-DP
checks for every k0..15 provide48 literal controls, in addition to the
symbolic coefficient proof.

The first examples are71/73/75,153/155/157,253/255/257. The last endpoint
at k2 is the Fermat prime257, but this numerical instance is not the
mechanism: the three charts work at every k regardless of primality.

**Anchoring boundary.** The lower value is at least82k+71. Thus membership
of any target h in the union of these triples has a complete finite
inverse: test0<=k<=floor((h-71)/82) and the three equations. In particular
101 is not in this triple family although section2 constructs it. The
three formulas do not make k an unrestricted choice of lower endpoint.
At most3X/82+O(1) values through X can arise, so this family alone cannot
cover all sufficiently large odd numbers.

## 4. What these comparisons do and do not transfer

| Source | Target/map | Preserved predicate | Missing coordinate / decisive test |
|---|---|---|---|
| Old canonical-source f | Current address g | Exact counts under n->n-1 | Omitting the shift fails already at2. |
| One-apex Fermat construction | Even-square target family | Actual H at each explicit address |101 refutes universal near-square coverage in that image. |
| C3 dyadic right flank | H=6*2^m+5 | Literal adjacency and count | Factoring t^2-1 proves only three near squares on this slice. |
| H9 chart changes | Three consecutive odd values | Exact interval and equal order | The lower endpoint remains9g(k)+82k+62; use the complete inverse, not a gap-only argument. |
| Neutral C6 source | H101 with a named root | Full profile and defect6 | H agreement with a different101 constructor does not transfer readiness. |

The finite-state selector obstruction in [SELECTOR](../cover/SELECTOR.md)
also excludes an ordinary-binary rational witness selector covering all
even-square-plus-one outputs: those have order sqrt(X) distinct values
through X, exceeding every fixed polylogarithm. This does not exclude an
unbounded-register construction of them or five-core membership itself.

The useful next operation is an anchor-changing switch whose complete
endpoint data permit repeated interval extension. A fixed2-gap, even a
fixed three-odd interval, leaves that obligation open. The rooted response
and higher path-cover tensors identify what must survive the switch.

## Reproduction

```text
python -B 05-knowledge/results/chanpak_five_core_bridges_20261110/history/near_square_bridge.py
```

[near_square_receipt.json](near_square_receipt.json) records the complete
one-apex target universe, all231 near-square misses, the two101 witnesses,
4,096 f/g-shift checks, exact chart coefficients and48 literal consecutive
count controls. The proof of the three-case dyadic-square classification is
the factorization above, not a finite search over m.
