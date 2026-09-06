# Full three-node three-jet dyadic Smith form from eight attained costs

**PROVED BY EXACT SYMBOLIC IDENTITIES / INDEPENDENTLY AUDITED.** The statement
has an exhaustive finite symbolic-minor proof, not an inference from integer
node samples. The independent polynomial, full-Smith and boundary audits
passed. It is a continuation of
[THM-4443, arbitrary-jet precision and dyadic unit boundary](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md),
not a new theorem namespace.

## 1. Inheritance and precise statement

THM-4443 gives the largest dyadic Smith exponent L for three nodes with
three Hasse jets each. Its hostile has identical metric but losses18/19.
The least-used sidecar is now the valuation of **every intermediate minor**,
retaining the global dilation weight separately from the close-pair depth.
The corrected near miss is to guess an unbounded affine formula from its
first few scales: the seventh determinantal exponent switches slope after
e=3d+1. The live concepts are complete inverse; intermediate ideals;
symbolic monomials; two independent depths; attained lower envelope;
intrinsic unit class.

For any three distinct integer nodes, write their dyadic pair depths as
(e,e,e+d), e>=0,d>=1. Let L be the exact THM-4443 value:

```text
d>=2: L=8e+5d-1;
d=1:  L=max(7e+4,8e+Gamma(t)),
```

where Gamma is its four-valued intrinsic unit class modulo16. Let D_j
be the valuation of the gcd of all j-by-j observer minors, with D_0=0.
The full formula is

```text
D_0=D_1=D_2=D_3=0,
D_4=e,
D_5=min(5e,4e+1,3e+d+1),
D_6=min(9e,7e+d+1),
D_7=min(13e+d,12e+4d+1),
D_8=27e+9d-L,
D_9=27e+9d.                                         (1)
```

Thus all nine Smith exponents are D_j-D_(j-1). All unit dependence is in
the last two factors, through L. In particular **metric plus worst loss
determines the full dyadic partition in this complete three-node three-jet
problem**, although metric alone fails. This is not true in the incoming
four-node two-jet examples and is not asserted at other primes/node counts.

## 2. Normalize over the local integral ring, without discarding unit data

Translate a closest-pair endpoint to zero. Divide the variable by the odd
unit part of the outsider's difference; this is an invertible integral
coefficient and Hasse-jet change over Z_2. The normalized nodes are

```text
0, 2^e, 2^e a,       a in Z_2, v_2(a)=d>=1.          (2)
```

This normalization preserves the full dyadic Smith form. The parameter a
need not be an ordinary integer. All polynomial valuations below apply to
every a in2Z_2 of the stated depth, so no rational-unit lift is lost.

Clear the three identity rows at zero. The residual6x6 matrix has columns
j=3,...,8, and rows indexed by(node,r), node in{1,a}, r=0,1,2:

```text
R_((node,r),j)=2^((j-r)e) binom(j,r) node^(j-r).       (3)
```

Write rows0,1,2 for the three jets at1 and rows3,4,5 for those at a.
For row subset I and degree subset J, every permutation term has common
dilation weight W=sum(J)-sum_(i in I)r_i. Its determinant is

```text
2^(We) P_(I,J)(a),        P_(I,J) in Z[a].           (4)
```

Thus a polynomial term c a^b gives the lower valuation bound
We+bd+v_2(c). Cancellation can only increase the valuation. The three
cleared unit factors mean rank-r residual determinantal valuations are
exactly D_(r+3).

## 3. Complete symbolic certificate and attained envelopes

The standalone
[certificate](../../04-computation/three_node_threejet_smith_overnight_hexagon_sep05.py)
enumerates **every** minor of residual ranks1,2,3,4:36,225,400,225 minors,
886 total. It expands each determinant by every permutation and retains
all nonzero integer monomials. This is a finite polynomial identity check
for arbitrary a,e, not a sampled valuation bank.

Set d=1+d', where d'>=0. Each nonzero term has nonnegative-domain affine
cost represented by(W,b,v_2(c)+b). The complete coordinatewise-undominated
sets of these triples are respectively

```text
r=1: {(1,0,0)},
r=2: {(5,0,0),(4,0,1),(3,1,2)},
r=3: {(9,0,0),(7,1,2)},
r=4: {(13,1,1),(12,4,5)}.                            (5)
```

For every discarded monomial the certificate exhibits an element of(5)
whose three coordinates are no larger. Therefore its value at every
e,d'>=0 is no smaller. Taking the minima of(5) proves the lower bounds
in(1) for D_4,...,D_7 simultaneously at every integral depth. No limiting
argument or finite cutoff on either depth is used.

The bounds are attained without cancellation by these eight exact minors:

| Rows I | Column degrees J | W | P_(I,J)(a) | exact valuation |
|---|---|---:|---|---|
| 2 | 3 | 1 | 3 | e |
| 0,2 | 3,4 | 5 | 3 | 5e |
| 1,2 | 3,4 | 4 | 6 | 4e+1 |
| 2,5 | 3,4 | 3 | 18a(a-1) | 3e+d+1 |
| 0,1,2 | 3,4,5 | 9 | 1 | 9e |
| 1,2,5 | 3,4,5 | 7 | 30a(a-1)(2a-1) | 7e+d+1 |
| 0,1,2,5 | 3,4,5,6 | 13 | 3a(a-1)(5a^2-5a+1) | 13e+d |
| 1,2,4,5 | 3,4,5,6 | 12 | 90a^4(a-1)^4 | 12e+4d+1 |

Every displayed parenthetical factor is a dyadic unit because a is even;
the constants have the indicated exact valuations. The certificate expands
each displayed factorization independently and checks it against the full
determinant. Thus every affine cost in(5) is an actual minor value, not
only a monomial lower bound. This proves equality for D_4,...,D_7.

The confluent determinant gives
D_9=3^2 sum_(i<j)v_2(x_i-x_j)=27e+9d. By definition the last Smith
exponent is D_9-D_8; the already-proved exact inverse law in THM-4443
therefore gives D_8=D_9-L. This completes the proof of(1).

## 4. Boundary controls and what the new operation changes

At e=0, (1) gives six unit factors followed by d,3d+1,5d-1, exactly the
independent two-node three-jet block plus the separated singleton.
At d=1 the final two factors retain Gamma even though all earlier
determinantal valuations forget it. For the earlier metric hostile at e=2,
the final pairs are(18,18) and(17,19), as required.

An initially tempting fit D_7=13e+d holds through e<=3d+1 but fails
afterwards. At e=5,d=1 the actual D_7 is65, not66. The responsible minor
uses the four derivative rows and has cost12e+4d+1. This is the first
failed implication: a shallow dilation bank did not see the second slope.
The repaired form retains every symbolic monomial weight and both scales.
The full candidate was not canonized before this test.

The primary checker also compares all886 symbolic minors with independently
evaluated literal determinants at a=2,4,6, and compares the complete Smith
partition on180 literal9x9 matrices: e=0,...,8, d=1,...,5 and t=1,3,5,9 at
the nodes(0,2^e t,2^(e+d)). These are positive/hostile controls, not the
unbounded proof. The independent
[referee](../../04-computation/three_node_threejet_smith_referee_overnight_hexagon_sep05.py)
reconstructs all886 minors using literal Hasse derivatives and a different
symbolic determinant algorithm, all eight attained costs, the intrinsic
unit class, and175 full Smith matrices with signed/unit-affine controls.
It passes1,110 gates. Primary and referee both agree under normal and
optimized Python; the primary has4,943 gates.

```bash
python3 04-computation/three_node_threejet_smith_overnight_hexagon_sep05.py
python3 -O 04-computation/three_node_threejet_smith_overnight_hexagon_sep05.py
python3 -B 04-computation/three_node_threejet_smith_referee_overnight_hexagon_sep05.py
python3 -B -O 04-computation/three_node_threejet_smith_referee_overnight_hexagon_sep05.py
```

Raw LF SHA256 manifest: primary source
`6359f824bd1301d3407c43bc9d6645616896a21842dee0b7573b2ae1a72934f9`;
[primary output](three_node_threejet_smith_overnight_hexagon_sep05.out)
`02f6d66898122c4ebcf857b9ffcfd2fd80c1f88f6864de0654e57a9210afb889`;
referee source
`c476bf118cc3f35468e2ea6b0895de83f682f8f02ae3a8bfd0e3a221c8a210e2`;
[referee output](three_node_threejet_smith_referee_overnight_hexagon_sep05.out)
`762c7d10e5c556113585f812b33ed90462117a8d1bab7c2e387ed08cf7b57edf`.
