# Guarded inverse paths: edge content and primitive parameter-five completion

**Status: PROVED elementary identities, guarded constructions and obstructions;
FINITE-EXACT controls.** Collatz convergence and completeness of the known
odd-parameter cycle lists remain **OPEN**. No novelty claim is made.
The main extension here is the two-edge completion theorem in §5, which
retains arbitrary powers of five as well as the inherited binary/ternary
coordinates. It is an existence theorem for ancestors, not a cover of all
fixed starting integers.

## 1. Inheritance, input corrections, and the live board

The closest proved mechanisms are the inverse fibre `R_b(x)=4x+b` and
parameter content in
[signed cycles, §§2–4](arithmetic_braids2_20260917_signed_cycles.md), and
the one-edge finite-word completion in
[inverse completion, §§1–3](arithmetic_braids2_20260917_inverse_completion.md).
The nine parameter-minus-five cycle witnesses below are inherited from
the former note, not discovered or exhaustively classified here.
The canonical hostile is three disjoint `3n-1` basins, each dense in the
same binary/ternary residue quotients. The corrected near miss is
identifying that local support with global inverse reachability. The
least-used sidecar is content modulo the additive parameter.

Anchor: guarded reachability for `3n-5`. Niche: the exact gcd of one edge
and its Pythagorean image. Wildcard: completion in the parameter's own
prime-power coordinate, where content is a genuine obstruction.

| Live concept | Preserved predicate | Lost coordinate / decisive hostile |
|---|---|---|
| Odd acceleration | Exact valuation and odd endpoint | One division gives `5->8`, still even |
| Inverse braid | One fixed target and all ternary source residues | Later basin of that target |
| Edge gcd | Primitive endpoint pair after division | The three-part of next edge content can grow |
| Pythagorean triple | Edge magnitudes and squared content | Orientation, signs, exponent and parameter |
| Finite word | One binary source cylinder | Fixed integer height and infinite continuation |
| Mixed completion | Prescribed word and source residues modulo powers of 3 and 5 | Disjoint basins can share every such cylinder |

The pasted input's one-division shortcut is not an odd-only map:
`(3*5+1)/2=8`. Under actual odd acceleration, an output of `U_1` is never
divisible by three, but the rows one and five modulo six do not alternate:
`11->17` stays in row five. Negation interchanges rows one and five and
preserves row three; it conjugates parameters `b` and `-b` in forward
time, without turning escape into a sink. The three known negative `U_1`
cycles and positive `U_(-1)` cycles are the same reflected list, not two
independent lists, and their completeness is not proved by the input.

## 2. The exact inverse guard, including the parameter-three exception

For odd nonzero parameter b, let

```text
U_b(x)=(3x+b)/2^v2(3x+b),     x odd, 3x+b != 0.
```

A proposed odd target y has predecessor x with exponent k exactly when

```text
k>=1,  2^k y=b mod3,  x=(2^k y-b)/3.                    (I1)
```

Once this congruence holds, x is automatically odd; since y is odd, k is
the actual valuation. Positivity is an additional inequality, not part
of the congruence. If `3` does not divide b, y must be prime to three and
k has one prescribed parity. If `3|b`, y must instead be divisible by
three, and every positive k is allowed. This last distinction is essential
when generalizing the familiar parameter-one inverse graph.

In every admissible case,

```text
3(4x+b)+b=4(3x+b),
U_b(4x+b)=U_b(x),              k -> k+2.                 (I2)
```

For `3` not dividing b, this is the whole inverse fibre after choosing
its initial exponent. The inherited identity
`v3(R_b^j(x)-x)=v3(j)` makes it a full cycle on odd residues modulo
`2*3^s`. It still lies over one fixed target. For `3|b`, the exponent
spacing in the full fibre is one and the same full-ternary conclusion
does not apply to (I2).

## 3. One edge has a sharp gcd law; triangle content is its square

**PROVED for every defined odd edge `x->y`:**

```text
gcd(x,y)=gcd(x,b)=g,
gcd(y,b)=gcd(3x,b).                                     (I3)
```

The first equality follows from `2^k y=3x+b`, since x is odd. The second
follows by removing only powers of two from a number whose gcd with odd
b is unchanged. Thus the first formula has no parameter-three exception.
The exception concerns the *next* edge content:

```text
v3(gcd(y,b))=min(v3(x)+1,v3(b)).                         (I4)
```

All other prime parts are conserved. The three-part increases one level
per step until it reaches the full three-part of b. The smallest hostile
to unrestricted content conservation is `U_3(1)=3`: the edge content is
one, and the next state's content with b is three. On a cycle the three
part must already be saturated, recovering the inherited cycle-content
theorem for every odd parameter.

Put `X=x/g`, `Y=y/g`. They are coprime odd integers and satisfy the same
edge relation at parameter `b/g`. For unequal absolute values their
Pythagorean image is primitive:

```text
( |XY|, abs(X^2-Y^2)/2, (X^2+Y^2)/2 ).                  (I5)
```

The identity is immediate by squaring. An odd prime dividing its first
and second entries would divide both X and Y; two cannot divide its odd
first entry. Before normalization, the triple therefore has gcd exactly
`g^2`. Its hypotenuse plus and minus the even leg are precisely the odd squares
`X^2,Y^2`. The boundary `abs(X)=abs(Y)=1` has a zero leg and is not a
nondegenerate right triangle. In particular the portal `1->-1` at `b=-5`
must not be called a primitive positive triangle.

At `b=+/-1`, every distinct positive edge gives a primitive triple. At
`b=-5`, content-five edges give exactly 25 times a primitive triple;
content-one edges give primitive triples, subject to the same zero-leg
exception. This is an exact representation, but does not select a
descending orbit: the Pythagorean identity holds for every legal edge.
The [edge/angle companion](odd_square_20260921_edges.md) retains the
direction and scale coordinates needed for its geometric comparisons.

## 4. The exact parameter-minus-five decomposition

For `b=-5`, content is conserved and equals one or five. Odd dilation gives

```text
U_(-5)(5z)=5 U_(-1)(z),
Anc_(-5)(5v)=5 Anc_(-1)(v),                            (I6)
```

where Anc includes the target itself and uses all signed odd states. The
reverse inclusion in the basin identity uses conserved content: any
ancestor of `5v` is a multiple of five. This is an exact all-height
reduction of one stratum, not a classification of the other stratum.

The nine inherited witnesses divide as follows; the replay checks every
arrow, exponent, ordered carry, and content.

| Content | Primitive parameter | Known cycles, labelled by smallest absolute node | Least periods |
|---|---|---|---|
| 5 | -1 | `5`, `25`, `85`, `-5` | `1,2,7,1` |
| 1 | -5 | `-1`, `-19`, `-23`, `-187`, `-347` | `1,3,3,17,17` |

The content-five list is the fivefold copy of the four known signed
`U_(-1)` cycles. The five content-one cycles have minimal supporting
parameter magnitude five; they do not reduce to integer parameter-one
cycles by odd dilation. This table is a list of verified witnesses, not
an exhaustive all-period census.

The only positive-to-negative portal is `1->-1`. Consequently none of
the other negative cycles can be reached from a positive start. Positive
multiples of five stay positive by (I6), while a positive content-one
orbit can enter the negative half-line only through one. Whether every
positive content-one start reaches that portal remains open here.

## 5. Two extra edges complete any word in every compatible 3/5 cylinder

For a finite positive exponent word `w=(k_1,...,k_L)`, including the empty
word, set `K_i=sum_(j<=i)k_j`, `K=K_L`, and

```text
B=sum_(i=0)^(L-1)3^(L-1-i)2^K_i;
K=B=0 for the empty word.
```

The inherited binary cylinder for this exact prefix at parameter minus
five is

```text
n=(2^K+5B)3^(-L) mod 2^(K+1).                         (I7)
```

All intermediate guards, not just formal endpoint integrality, are
equivalent to this cylinder. For example, write the full carry as a
prefix carry plus `2^K_j` times an odd remaining carry. Reducing (I7)
modulo `2^(K_j+1)` forces exact valuation `K_j` of the j-step numerator.
Thus this is one class of relative density `2^(-K)` among the odd integers.

**PROVED mixed completion theorem.** Fix an odd target u prime to 15,
`s>=0`, `t>=1`, any residue r modulo `3^s`, and any unit residue a modulo
`5^t`. There are infinitely many odd n, with the same sign as u, whose
entire path to u has that sign, such that:

* n realizes w followed by two exponents `(k,ell)`, with `k in {1,2}`;
* `n=r mod3^s`, `n=a mod5^t`, and n retains the binary cylinder (I7);
* every step is an exact `U_(-5)` step, and all nodes are coprime to five.

Proof. Put `N=L+2` and, for either `k=1,2`, let

```text
C_k=2^(K+k)+3*2^K+9B,
E=K+k+ell,
n=(2^E u+5C_k)/3^N.                                   (I8)
```

These are exactly the composed affine equations for the proposed path.
The prescribed source residues are equivalent to

```text
2^E=(3^N r-5C_k)/u mod3^(N+s),
2^E=(3^N a-5C_k)/u mod5^t.                             (I9)
```

Both right sides are units. Two generates the full unit groups modulo
`3^m` and `5^t`, of orders `2*3^(m-1)` and `4*5^(t-1)`: modulo the primes
its orders are two and four, and the elementary lifting identities
`v3(4^j-1)=1+v3(j)` and `v5(16^j-1)=1+v5(j)` give the claimed orders.

Let the two discrete-log conditions in (I9) be `E=e_3 modP_3` and
`E=e_5 modP_5`. Their periods have gcd two. The parity of e_3 is
`K+k+epsilon`, where `u=(-1)^epsilon mod3`. The parity of e_5 is
independent of k: the right side modulo five is `3^N a/u`, because the
carry term contains five. Choosing one of `k=1,2` therefore makes the
parities agree. The ordinary noncoprime CRT now supplies E, in an
arithmetic progression of period

```text
P=lcm(P_3,P_5)=4*3^(N+s-1)*5^(t-1).                   (I10)
```

Take E sufficiently large that `ell=E-K-k>=1`. The first congruence
makes (I8) integral. Reading the affine relations backwards, the
divisibility by `3^N` successively implies every intervening division by
three is integral. Each resulting numerator is even times an odd integer
plus five, so each intermediate node is odd; thus all exponents are exact.
As E grows along (I10), every intermediate node has leading term of the
same sign as u. Discarding finitely many E therefore gives the asserted
same-sign paths. Conserved content makes every node prime to five.

This is an explicit construction; no orbit-convergence hypothesis is
used. The program uses finite discrete-log tables and CRT, then checks
the path independently by backward guards and actual forward valuations.

**Two additional edges are sharp for universal residue prescription.**
For the empty prefix and target `u=1`, every one-step predecessor has
even exponent and is `(2^ell+5)/3`, hence is two or three modulo five.
The zero-step target is one modulo five. Neither realizes source class
four. Two steps do: `459 -> 343 -> 1` with exponents `(2,10)`.
No claim of smallest source height is attached to this witness.

## 6. What this proves about inverse reachability, and what it cannot prove

For every `H>=1`, `s>=0`, `t>=1`, every odd class modulo `2^H3^s5^t`
that is prime to five contains infinitely many same-sign ancestors of
every target u prime to 15. To see this, read a finite word from any
integer in the desired binary class until its total exponent is at least
`H-1`; (I7) then fixes that binary class. Apply the theorem for the other
two prime-power coordinates. No integer zero numerator occurs at b=-5.
If a prefix was already prescribed, extend that prefix using a representative
in its compatible binary class before applying the same construction.

In particular each of the **five known primitive negative cycle basins**
meets every such class among negative integers. These basins are disjoint,
since their verified terminal cycles differ. Even prescribing a fixed
halving prefix, together with arbitrary fixed powers of two, three, and
five compatible with that prefix, cannot separate these basins.
For positive endpoint u=1 the construction gives positive ancestors of
the portal, which then reach `-1`. It gives no positive ancestors of the
other four primitive negative cycles; the sign obstruction in §4 remains.

The content obstruction is equally real: inverse generation from five
never reaches one, and inverse generation from one never reaches five.
Within the parameter-one positive map, saying that every integer is an
ancestor of one is exactly the original convergence claim; a full inverse
residue braid does not prove it.

Finally, for each fixed prefix, target, and residue prescription, (I8)
grows geometrically as E advances by P. This particular witness family
has only `O(log X)` sources of absolute size at most X, hence natural
density zero. Full finite-residue support is compatible with sparsity,
large heights, and several disjoint basins. A proof of universal backward
coverage must supply a height-controlled certificate for each fixed
integer, not just another ancestor somewhere in its residue class.

## Reproduction and finite scope

[Script](../../04-computation/experiments/odd_square_20260921_inverse.py)
and [JSON](../../04-computation/experiments/odd_square_20260921_inverse.json):

```text
python 04-computation/experiments/odd_square_20260921_inverse.py
python -O 04-computation/experiments/odd_square_20260921_inverse.py
```

The exact mixed universe contains all 13 prefixes of lengths zero through
two with exponents in `{1,2,3}`, targets `1,-1,-19,-23,-187,-347`, all
nine residues modulo nine, and all 20 units modulo 25: 14,040 completions.
Another 14,040 controls take a larger member of the same infinite family.
Every division guard and actual forward valuation is checked. Separate
controls cover odd parameters through absolute value 21, edge content and
triangle content, the parameter-three exception, inverse parity guards,
all nine inherited cycle witnesses, and the sharp two-edge example.
All checks use explicit exceptions and remain active under optimization.

Independent proof audits: **PASS** from the parent and a separate peer on
the two-exponent CRT, parity compatibility, backward guards, sign boundary,
binary refinements, and sharpness witness. The peer also replayed `main()`
read-only and obtained exactly the saved decoded JSON. Normal and optimized
script runs produce byte-identical JSON.
