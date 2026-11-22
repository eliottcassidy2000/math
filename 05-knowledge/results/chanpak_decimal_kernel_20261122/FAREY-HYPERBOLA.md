# A target-selected lattice search for two Farey flanks

**FINITE-EXACT.** In the canonical tagged tower with left tag mask 12, no pair of binary flank addresses at spine lengths **12, 16, 20, 24, 28 or 32** realizes the lower count \(N=2^{60}(2^{61}-1)-1\). The other orientation always adds two. Each exclusion covers independent arbitrary words on both sides, including either or both empty. All lengths divisible by four from 36 onward were already excluded. Therefore only spine lengths **4 and 8** remain open in the canonical tower. [Inherited tower definition and large-length bound](FERMAT-RECEIVER-SEARCH.md)

**Remaining computation.** On lengths 8 and 4 the current geometric schedule estimates roughly 270 billion and 3.1 trillion tangent blocks, respectively, for an unrestricted sweep, in addition to the near-axis work. These are costs of this implementation. The searches below close stated sectors at those lengths. The [one-run and unrestricted-flank inverse](REVERSED-PHASE.md) and [fixed-sum inverse](CUBE-ROOT-PHASE.md#inverting-an-entire-antidiagonal-at-exponent-61) give complementary exclusions.

## Removing a congruence before choosing either flank

The lower chart from the seven slot readings is
\[
 N=k+t x+c y+qxy+r(x\bmod2)+e(y\bmod2)+h\bigl(f(x)+f(y)\bigr).
\]
Write \(g=\gcd(h,c,q,e)\). Regardless of the second flank, \(N\equiv k+tx+r(x\bmod2)\pmod g\). For this target every admissible left word at length 12 is \(141\) or \(170\) modulo 290. For every other tower length from 4 through 32 it is \(0\) or \(1\) modulo 10. Here \(g=145\) at length 12 and 5 otherwise.

Thus replace \(x\) by \(2gX+b\) for each allowed \(b\). More generally write the two addresses as \(s_x X+o_x,s_y Y+o_y\), with one scale 1 and zero offset. Expansion leaves positive coefficients \((Q,C,T,K)\) on \(XY,Y,X,1\). An even scale absorbs a fixed parity correction into \(K\). The unknown nonlinear readings still have coefficient \(h\). Choose the smaller of \(X,Y\) first. Its bound is \(\lfloor\sqrt{(N-K)/Q}\rfloor\), and the second orientation excludes the diagonal to count a solution only once. After normalization these bounds at lengths 20, 24, 28 and 32 are approximately 61 billion, 6.8 billion, 753 million and 84 million for each surviving residue. Both raw words are reconstructed before testing the chart.

## Why thin hyperbola strips can be decoded together

Put \(D=(N-K)/Q,\ \alpha=C/Q,\ \beta=T/Q\). Both offsets lie between 0 and 1 in the normalized charts. Dropping the nonlinear terms and the remaining parity corrections gives the exact real upper boundary
\[
 Y_0(X)=\frac{D+\alpha\beta}{X+\alpha}-\beta.
\]
For \(X\ge L\ge1\), \(Y\le D/X\) and \(0\le f(z)\le z\). Consequently a solution lies in \(Y_0(X)-G\le Y\le Y_0(X)\), where
\[
 G=\frac{h}{Q}\left(s_x+s_yD/L^2\right)+
       \frac{h(o_x+o_y)+r+e}{QL}.
\]
Parity coefficients already absorbed in \(K\) can be dropped from this bound.

On a block \(X_0,\ldots,X_0+W-1\), put \(\lambda=-Y_0'(X_0)\). Pick a reduced rational \(a/b\) near \(\lambda\), with \(b\le W\). For any integer anchor \(y_*\) near \(Y_0(X_0)\), set \(F=Y_0(X_0)-y_*\). Every solution has an *integer* line index
\[
 z=b(Y-y_*)+a d,\qquad d=X-X_0,
\]
in the interval
\[
 b(F-G)+\min(0,(a-b\lambda)W)
 \le z \le b(F+B)+\max(0,(a-b\lambda)W),\qquad
 B=(D+1)W^2/L^3 .
\]
Indeed \(Y_0\) is convex with second derivative at most \(2(D+1)/L^3\). On any surviving line, \(ad\equiv z\bmod b\) fixes the progression of offsets and the line determines \(Y\) exactly. The code takes \(W\) of the order \(L/D^{1/3}\) and sweeps adjacent fractions of the Farey sequence of order \(W\) as the slope decreases. Its neighboring fractions also give \(a^{-1}\bmod b\) without another Euclidean algorithm. The outer receiver uses a determinant of tagged path counts. The second Farey choice uses curvature at the prescribed target. Any reduced rational satisfies the interval proof.

Near the axes the code uses congruences first. After fixing the allowed residue, divide the coefficients of the normalized equation and the translated target by their common content. At length 12 this removes 145 and exposes the full modulus \(4097=17\cdot241\) in the direct row sieve. A divisor of \(h\) removes the two nonlinear readings. At fixed address residue and both parities this gives one arithmetic progression for the opposite address or eliminates that branch. For each surviving outer word the interval above rejects an empty progression near \(Y_0\). For longer opposite words the exact integer stage splits off \(b\) low bits, squeezes the remaining quotient by \(0\le f(z)\le z\), and intersects it with the necessary slope congruence modulo \(h\). It uses the [low-word concatenation identity](../chanpak_e61_fibonacci_target_20261122/SUFFIX-INVERSE.md). The fixed-address step always tests the literal exact height of a candidate.

The floating-point operations enclose the geometric strips. All reported heights use unsigned 128-bit integer arithmetic. For the vectorized row test the absolute rounding allowance is \(2^{-47}(D/(mL)+2)\) for a progression of spacing \(m\). The omitted shift \(\alpha\beta/(X+\alpha)\) gets an additional \(2/L\) in unscaled coordinates. For the slanted lines the anchor residual numerator \(N-K-TX_0-(QX_0+C)y_*\) is formed as a signed integer before division. With IEEE binary64 unit rounding error \(2^{-53}\), converting and dividing the positive chart entries costs at most six units of relative error. The slope \((y_*+F+\beta)/(X_0+\alpha)\) uses another two additions and a division. The possible cancellation in \(a-b\lambda\) is charged in **absolute** error proportional to \(a+b\lambda\). The implementation uses twice the positive bounds on curvature and thickness, and expands each line interval by \(2^{-46}[1+b(1+|F|+B'+G')+W(a+b\lambda)]\) using the enlarged bounds \(B',G'\). Thus neither rounding nor an inaccurate initial choice of the nearest rational can discard a solution. The searched normalized small coordinates are below \(2^{52}\). The tangent stage caps \(D/L\) at \(3\cdot10^{15}\).

At length 12, the upper bound for the smaller normalized coordinate is 910,755,735,559 for both residues of the physical left word. Both forward charts and both reverse charts above physical right address 128 returned no hit. The [unrestricted small-flank inverse](REVERSED-PHASE.md) already covers physical right addresses 0 through 128 at that spine. Running the command below from zero performs the combined test without using this overlap.

## Reproduction

From this directory:

~~~sh
python3 -B hyperbola_controls.py > /tmp/farey-e61-charts
c++ -O3 -std=c++17 farey_hyperbola.cpp -o /tmp/farey-e61-hyperbola
python3 -B hyperbola_controls.py --check /tmp/farey-e61-hyperbola
for n in 12 16 20 24 28 32; do /tmp/farey-e61-hyperbola "$n" < /tmp/farey-e61-charts; done
~~~

The controls compare the complete inverse with direct pair enumeration on small planted targets and test long nonconstant flanks in both orientations, including a pair in the tangent regime. For another lower height pass its decimal integer after the spine length. Extra arguments select direction (0 or 1), inclusive normalized small-coordinate range, and the allowed physical left residue. Each returned pair is checked against the integer chart. [Exact inverse](farey_hyperbola.cpp) · [Generator and controls](hyperbola_controls.py)

## A target inverse along every increasing pair of words

**PROVED.** Put \(P(x,y)=k+qxy+tx+cy\) on the physical, unnormalized coordinates. Every lower count belongs to the envelope
\[
 [P(x,y),\,P(x,y)+h(x+y)+r+e].
\]
In each of the eight charts \(q>h\) and \(q+t+c>r+e\). Increasing both coordinates by at least one increases the lower envelope bound by at least \(q(x+y+1)+t+c\). This exceeds the entire previous interval's width. Thus a prescribed absolute height has **at most one candidate state on any chain with both coordinates strictly increasing**. The last state with \(P\le N\) is the only one to test. A synchronized bit-output controller starting from two positive words has this property on every emitted bit pair. The statement concerns the words emitted by one controller and places no bound on its computation time between emissions.

For example, choose a positive rational direction in lowest terms \((a,b)\) and grow both words by a common integer scale \(j\). Its unique possible positive scale is
\[
 \left\lfloor
 \frac{\sqrt{(ta+cb)^2+4qab(N-k)}-(ta+cb)}{2qab}
 \right\rfloor .
\]
At each of the lengths 4, 8, 12 and 16 the search over all **9,185,685** primitive directions of product \(ab\le10^6\) is empty. At the unique scale none passes both the residual interval and divisibility by \(2h\). In particular any positive pair at lengths 4 and 8 must respectively have \(\gcd(x,y)\le1{,}076{,}565{,}685{,}641\) and \(\gcd(x,y)\le135{,}481{,}007{,}739\), since otherwise \(qxy\le N-k\) would force the primitive product back into the cleared range. Run the following from this directory. [Exact ray inverse](farey_directions.py)

~~~sh
python3 -B farey_directions.py 1000000 4 8 12 16
~~~

### A least set of supersingular primes after fixing absolute scale

Before reducing modulo \(2h\), require only the real envelope at the chosen scale and the actual parity terms:
\[
 0\le N-P(x,y)-r[x\ {\rm odd}]-e[y\ {\rm odd}]\le h(x+y).
\]
Among the same \(9{,}185{,}685\) primitive directions per chart this leaves exactly **1,493, 356, 80 and 7** pairs at the four successive spine lengths. Each of the 15 absolute supersingular primes has a separate witness in the resulting pool of **1,936** chart-address pairs. Yet simultaneous reduction at **17 and 47** rejects every one of them. The least possible number of primes for this finite pool is two, and nine distinct two-prime choices suffice. The marked upper count is the lower plus two, so a lower match also matches the upper coordinate over the same field. For the separate pools at lengths 12 and 16 one prime suffices, respectively 71 and 7. The test enumerates subsets of primes in increasing cardinality, reducing the actual chart differences on the same absolute candidates. This is a restriction on this bounded set of directions after the real inverse, compatible with the all-modulus local lifting theorem on unrestricted words. [Full set of local word inverses](LOCAL-WORDS.md)

## Every elementary cellular rule on small periodic tapes

**FINITE-EXACT.** For each of the four small charts, run every radius-one binary cellular automaton with one fixed local rule on a periodic ring of width \(1\le w\le14\). Begin with any labelled initial state, select two sites, and on each tick append the observed bits to two addresses initially equal to 1. Every retained pair consequently has equal binary lengths. Rotation lets the left site be site zero when all initial states are allowed. The right site can have any displacement, including zero. At **every possible time** these controllers miss exponent 61 in the four charts. The conclusion is specific to this family of synchronized observations.

Indeed at length \(\ell\) both coordinates are in \([u,2u-1]\), \(u=2^{\ell-1}\). Their lower count lies between \(P(u,u)\) and \(P(2u-1,2u-1)+2h(2u-1)+r+e\). Comparing these endpoints to \(N\) permits precisely lengths **50, 47, 44 and 41** at the successive spine lengths \(4,8,12,16\). Iterate the rule only to the permitted times. For each chart there are \(256\sum_{w=1}^{14}w2^w=109{,}052{,}416\) rule/state/displacement presentations. None even reaches the parity-corrected real envelope at the required time. Rules, observations and initial states can give the same word pair, so this count is of presentations. An independent planted control with rule 110, width 5, state mask 19, displacement 2 and nine emitted pairs gives words 832 and 576 and lower height 1,100,614,525,773 at spine length 4. The census at this target recovers ten presentations through width 5.

~~~sh
c++ -std=c++17 -O3 farey_elementary.cpp -o /tmp/farey-cellular
/tmp/farey-cellular 5 1100614525773 < /tmp/farey-e61-charts
/tmp/farey-cellular 14 < /tmp/farey-e61-charts
~~~

[Spacetime census](farey_elementary.cpp) · [Finite word controllers and all-prime local inverses](LOCAL-WORDS.md)

## Two different balanced sectors at lengths 4 and 8

**FINITE-EXACT.** At length 8 there is no hit with \(\min(x,y)\ge0.99\max(x,y)\). At length 4 there is no hit with \(\min(x,y)\ge0.999\max(x,y)\). Both physical addresses may be arbitrary words. To reproduce these regions put \(S=(N-k)/q\). The chart satisfies \(qxy\le N-k\le q(x+1)(y+1)\) at any solution. Thus \(\min(x,y)\ge\alpha\max(x,y)\) forces
\[
 \sqrt{\alpha S}-1\le x\le\sqrt{S/\alpha}.
\]
Since \(x=10X+i,\ i=0,1\), the large addresses in this range have \(X\le y\); the forward chart covers the whole sector. The first two commands use outward-rounded endpoints for \(X\). With the composite-divisor prefilter they check 73,195,203 and 266,868 final physical address pairs, respectively.

The normalized balanced end of the search has different physical slopes. With \(U=\min(X,y)\), the final 34,359,738,368 possible values of \(U\) at length 4 and final 68,719,476,736 at length 8 are excluded. Here both address orientations are run. The lower part of the full hyperbola at these two spines remains open. The [height census by physical residue](FAREY-PHASE-CENSUS.md) bounds the total number of targets across all such sectors.

~~~sh
/tmp/farey-e61-hyperbola 4 E61 0 107602780617411 107710491108520 < /tmp/farey-e61-charts
/tmp/farey-e61-hyperbola 8 E61 0 13480196806844 13616360410954 < /tmp/farey-e61-charts
/tmp/farey-e61-hyperbola 4 E61 -1 340405772222361 340440131960728 < /tmp/farey-e61-charts
/tmp/farey-e61-hyperbola 8 E61 -1 42774158359912 42842877836647 < /tmp/farey-e61-charts
~~~

## Prescribed differences along a Berggren direction

**FINITE-EXACT.** At lengths \(n=4,8\), choose a primitive direction \((a,b)\) from the following nine. Take \((1,1)\), take the six ordered choices of different entries from \(3,4,5\), or take the two orders of the actual Berggren parameters \((2^n+1,(3(2^n+1)-1)/5)\). For **every** such direction an exact E61 solution would require at least three binary one-runs in \(\lvert bx-ay\rvert\). There is no restriction on the binary complexity of the addresses \(x,y\) themselves. The last two directions use the attained endpoint counts of the relevant tagged receiver. [Physical Berggren parameters](BERGGREN-RECEIVERS.md)

For \(v\ge0\), the line \(bx-ay=v\) has a canonical initial nonnegative point. Choose \(0\le y_0<b\) with \(a y_0=-v\bmod b\), then \(x_0=(v+a y_0)/b\). Every other nonnegative solution is \((x_0+a j,y_0+b j)\) with \(j\ge0\). Interchange the coordinates to treat negative offsets. On that line let
\[
 W=N-P(x_0,y_0),\qquad u=t a+c b+q(a y_0+b x_0).
\]
If \(W\ge0\), the chain theorem makes the unique candidate
\[
 j=\left\lfloor\frac{\sqrt{u^2+4qabW}-u}{2qab}\right\rfloor.
\]
Thus choosing a slope and difference supplies **both** flank addresses from the prescribed absolute height. For \(v\ge0\), the upper bound \(v\le b\lfloor(N-k)/t\rfloor\) follows from the physical left slope. Each one-run offset is \(2^j-2^i\); each two-run offset chooses four strictly increasing endpoints. Zero is tested on the nonnegative side only.

At lengths 4 and 8 there are respectively 81,877,008 and 66,616,292 lines across the nine directions and both signs. At the selected scale 9,527 and 1,445 pairs survive the parity-corrected real envelope; respectively four and zero survive the necessary divisibility by \(2h\). All four miss the literal integer reading. Different slopes can describe the same pair. For a control at each direction, the program plants both signs of a unit determinant on a noninitial point of the line. [Exact determinant inversion](farey_sparse_determinants.py)

~~~sh
python3 -B farey_sparse_determinants.py 4
python3 -B farey_sparse_determinants.py 8
~~~

**The ordinary address gap needs at least four runs.** At \((a,b)=(1,1)\), testing exactly three one-runs extends the preceding exclusions. Each signed difference now chooses six strictly increasing binary run endpoints inside the same target bound. There are 2,625,088,745 additional signed lines at length 4, with 886,763 envelope survivors and 412 passing divisibility by \(2h\). At length 8 the counts are 1,804,806,008 lines, 143,111 envelope survivors and exactly one passing that divisibility. Every final integer read misses. The native replay locates the quadratic root by a floating seed, then increments or decrements until the two defining *integer* inequalities hold. Planted addresses \(32,11\) and their transposes exercise both signs at each spine. [Three-run determinant inverse](farey_gap_runs.cpp)

~~~sh
c++ -std=c++17 -O3 farey_gap_runs.cpp -o /tmp/farey-three-run-gaps
for n in 4 8; do
  for s in 0 1; do /tmp/farey-three-run-gaps "$n" "$s" < /tmp/farey-e61-charts; done
done
~~~
