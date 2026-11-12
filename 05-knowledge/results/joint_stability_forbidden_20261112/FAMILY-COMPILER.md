# A complete repair compiler for the binary birth-clock library

**PROVED for the declared finite library of infinite families;
FINITE-EXACT input and interface searches. Chan–Pak OPEN.**

All540 distinct normalized binary coefficient laws from the complete
stable-birth bank through ten now have every lower head repaired. An
additional320 distinct higher-mode laws have odd slopes at every length,
so their original tails already start at the first exponent. These are
actual tournament constructions. The remaining478 even-slope higher-mode
axes are not claimed to have uniformly completed heads.

## B1 A one-sided bound for the slope, separately from the defect

Fix a marked tournament, an unselected vertex p and its fixed exterior.
Write d+ and d− for the exterior out- and in-degrees of p. In a cover
counted by U, the second path starts at a selected vertex. Thus at most
one run of p can start a path, while at most two runs can end paths.
Distinct other runs require distinct exterior incoming or outgoing
neighbors. Therefore the complete slope has the finite expansion

```
U(a)=sum_(r=1)^R u_r r!S(a,r),
R=min(d−+1,d++2),       u_r>=0.                         (B1)
```

Hamilton paths starting selected, the other part of U, satisfy a stronger
bound and are included. The exterior vertices are single fixed vertices;
the claim does not make this bound uniform when unrelated exterior
modules also grow. Values at a=1,...,R recover all coefficients by the
triangular Stirling matrix. This is a direct application of the inherited
path-run substitution mechanism, keeping the asymmetric selected-start
condition that distinguishes U from delta.

[affine_births.py](affine_births.py) applies(B1) to every1,379 single birth
clock in the1,150-birth bank. [audit_affine_births.py](audit_affine_births.py)
independently recomputes every full signed and slope kernel by the separate
split-cover skeleton algorithm. All1,379 entire tensors agree, including
the absence of coefficients beyond the proposed bound.

| Exact mode class | Birth axes |
|---|---:|
| Constant odd part |20|
| Nonconstant binary affine slope |561|
| Higher modes, odd slope |320|
| Higher modes, even slope |478|

Different axes can induce the same coefficient law; these are not counts
of disjoint value sets. The20 constant-odd-part axes yield the14 constants
1,3,15,25,39,45,65,75,117,149,169,189,195,507. Their full rays are already
settled by [FORBIDDEN.md F3](FORBIDDEN.md), with7 and21 still excluded when
they are the requested output values.

## B2 Normalize a binary clock without losing its valuation

If all u_r for r>=3 vanish, then U(a)=q*2^a+b, where q=u_2 and
b=u_1−2u_2. The nonconstant-odd-part cases have q>0 and b!=0. Put

```
j=v2(b), v=v2(q), a0=max(1,j-v+1),
A=q*2^(a0-1)/2^j, B=b/2^j.
```

Here A is a positive integer, B is odd, and for every k>=1 at port length
a=a0+k−1 we have

```
U(a)=2^j*c_k,   c_k=A*2^k+B odd,   c_(k+1)=2*c_k−B.     (B2)
```

The source order is k+n+a0−2. The tail supplies all exponents s>=j+1.
This derivation treats the finitely many earlier port lengths separately;
it does not assume their valuations follow the eventual formula.
The561 axes give540 distinct normalized(A,B) laws. Where several axes
give the same law, the compiler retains the smallest available j.

## B3 Compile each unpaid head into a physical interface obligation

For s=1,...,j, the missing head is

```
1+2^s*c_k = D+C*2^k,
D=1+2^s*B, C=A*2^s.                                   (B3)
```

Choose an already realized positive factor G dividing gcd(|D|,A).
Write A/G=a*2^v with a odd. An actual marked graph with

```
delta=D/G, U=a*2^ell
```

gives(B3) after a positive tail of length t=k+s+v−ell and an ordered
join with the G witness. It therefore covers all k>=max(1,ell−s−v+1).
Any earlier k are a finite remainder. The literal graph, its selection,
ell, factor witness and tail length are all retained; a matching pair of
integers without an actual graph is rejected.

The minimum-j540-law bank has786 head obligations. The baseline search
uses every ordinary tournament through9 and all strong order10 sources
admitting a unit selection, but allows **all** selections when looking
for the new auxiliary interfaces. It retains4,475 interfaces for994
queries and closes784 head obligations, leaving exactly two. Its128
finite remainder entries are all in the inherited complete order-eleven
spectrum, with maximum value68,671. No larger-prefix coverage theorem is
silently used.

The sources of this compiler are actual endpoint interfaces; its targets
are actual Hamilton-count families; its map is positive-tail construction
followed by ordered join. The preserved predicate is the exact count
identity. Projecting to(D,oddpart(U)) loses ell, source order, selection
and graph identity, so all four are required sidecars.

## B4 The two remaining obligations expose the wrong source filter

For c_k=595*2^k−3, the unpaid s=2 family is2380*2^k−11.
For c_k=11007*2^k−173, it is44028*2^k−691. No matching physical interface
occurred in the declared baseline bank. This was a failure of that bank,
not evidence that any of these targets was impossible.

An independent C++ search of the first541,148 ordinary order10 classes
finds both missing interfaces:

| Rows and selection | H | delta | U |
|---|---:|---:|---:|
|[64,1,35,23,199,539,430,303,63,479], O={2,3,6}|1527|−11|9520|
|[512,37,385,71,143,93,279,363,571,254], O={5,6,8}|7237|−691|44028|

Neither underlying tournament has **any** unit selection. Restricting
auxiliary sources to unit-bearing tournaments was the precise obstruction.
Both complete endpoint tables are replayed by two independent Python
engines. The search stops at two positive witnesses and makes no complete
order10 nonunit-interface census claim.

The first interface uses tail length k−2 for k>=3. Its two finite
remainders4749 and9509 have these explicit ten-vertex witnesses:

```
4749: [512,129,11,147,103,79,271,373,575,254]
9509: [384,97,43,531,519,345,157,62,222,487].
```

The second interface uses tail length k for every k>=1. Consequently
every786 head obligations is closed. Combined with(B2), this proves

> For every one of the540 stored(A,B) laws, every k>=1 and every s>=1,
> the integer1+(A*2^k+B)*2^s is an actual tournament Hamilton-path count.

The complete library is [completed_bulk_heads.json.gz](completed_bulk_heads.json.gz).
The final independent audit recounts786 chosen interfaces and1,571 whole
constructed tournaments, with1,572 exact affine checks. All-length
correctness comes from the proved constructor identities, not extrapolation
from these finite controls. Finite remainder witnesses inherited from the
order-eleven spectrum remain explicitly distinguished from new literal
adjacency witnesses.

## B5 Another320 complete families from parity alone

For every r>=2, r!S(a,r) is even. Therefore(B1) immediately gives

```
U(a) = u_1 (mod2) for every a>=1.                       (B4)
```

Exactly320 higher-mode birth clocks have odd u_1. For each, set c_a=U(a).
The original tail realizes every1+c_a*2^s with a,s>=1: there is no head
debt. Their full Stirling coefficient lists are distinct, so these are
320 different parameterized coefficient sequences. They are not claimed
to have disjoint images or to be inequivalent under every reindexing.
The remaining478 higher-mode axes have even slopes. Their valuations
can vary with a, and(B4) alone does not close their lower heads.

## B6 A complete coefficient-overlap metagraph

The540 binary families also admit a full all-length equality classification.
For two laws A*2^k+B and C*2^l+D with B!=D, equality implies

```
min(k,l) <= v2(B-D).                                   (B5)
```

Enumerate that finite range in either coordinate and test whether the
remaining ratio is a positive power of2. This decides every collision,
including arbitrarily large parameters, without a sampled length bound.
For B=D, two laws intersect exactly when A/C is a power of2, in which
case they are the same eventual dyadic ray. No such pair occurs among
the540 laws here.

The resulting coefficient-overlap graph has540 vertices,257 edges,
312 connected components and cycle rank29. There are211 collision values,
the largest557,585, and maximum value multiplicity4. There are258 isolated
laws. Component sizes are completely recorded in the receipt. Its26
triangles all have a common coefficient in this bank; that is a finite
fact, not a general Helly theorem. Abstract dyadic rays2^k+101,
2^k+99,2^k+103 already give pairwise intersections103,105,107 with no
common value, so a general pairwise-to-joint inference is false.

Let M(c) be the number of these laws containing a collision value c. The
exact coefficient image count at any positive height X is

```
sum_(A,B) max(0, floor(log2((X-B)/A)))
   - sum_(collision c<=X) (M(c)-1),                    (B6)
```

where a nonpositive ratio contributes zero. Integer division and bit length
implement this without floating point. The finite collision correction is
essential; numbers of parameter words and distinct values differ.

In particular the number at dyadic height2^b is **exactly540*b−6078 for
every b>=21**. This threshold is proved by monotone inequalities for each
ray, plus the complete collision list, not inferred from a linear-looking
sample. At2^20 the count is4,723. In general(B6) is540 log_2 X+O(1).
[coefficient_intersections.py](coefficient_intersections.py) verifies the
closed formula against literal image unions through2^20 and supplies the
complete finite graph and value-labelled collision fibers.

This metagraph is different from the tournament contraction graph. Equality
of c does not identify marked tournaments or preserve their endpoint
interfaces. A valid graph-to-number projection retains the source witness
and valuation on every fiber. Switching constructors at a common number
does not generate coefficients outside this finite union by itself.

## B7 The remaining universal obstruction

These860 nonconstant laws, the finite coefficient band and the full birth
classification are substantial constructive libraries. Each still comes
from finitely many fixed templates. For the1,152 stable cylinders the
inherited path-height bound gives at most

```
967*L^2+145*L^3+38*L^4+2*L^5,   L=floor(log2(X-1)),
```

literal unit-tail outputs below X. The new finitely many head interfaces
add only another finite-template contribution. Their asymptotic sparsity
is proved and prevents this library alone from implying Chan–Pak.

The useful next hypothesis is a **recursive grammar that creates new
birth cores and has a finite repair-interface description**, rather than
a fixed finite set of cores. The present experiment proves complete head
closure for the binary library, and pinpoints why its initial source
filter failed. It does not prove that the same finite repair bank works
for every later birth or for the478 even-slope higher-mode axes.

The exact first gap of the completed coefficient sublibrary remains53,559.
This checks the full binary laws and all positive lengths of the320
higher-mode odd laws, not just a length sample: their nonnegative Stirling
coefficients make U nondecreasing and provide a certified stopping point.
It excludes none of the other possible seeds or the478 even-mode axes.
Moreover1+2*53,559=107,119 is explicitly realized by an eighteen-vertex
tournament from the inherited five-core attachment compiler. Its complete
adjacency rows and literal Hamilton count appear in
[library_boundary.json.gz](library_boundary.json.gz). A missing seed address
in this library is therefore not a globally forbidden Hamilton count.
