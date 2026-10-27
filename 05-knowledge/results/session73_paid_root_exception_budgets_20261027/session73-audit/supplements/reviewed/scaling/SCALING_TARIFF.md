# Addendum: a forest-to-forest complement with linearly growing rank-cost change

**PROVED elementary; finite optimality controls independently use subset DP.**
This is an additive hostile to over-reading the small-n tariff census, not a
change to the frozen main packet or a new ROOT claim.

For n≥5, let S_n have natural rank0,…,n−1 and a single reversed pair
{1,e}, e=largest even number below n. Let P list the evens increasingly and
then the odds increasingly. The parity-chain steps are forward; the only
rank descent e→1 is exactly the supplied reverse edge. Thus P is Hamilton
in S_n. The source has f(S_n)=1: its lone exception gives the upper bound,
and 1→2→e→1 is a directed triangle (e≥4).

Write U_n for the tournament whose natural-rank exceptions are all pairs
{i,i+2}, i=0,…,n−3. This is a disjoint union of two paths, hence a forest.
Path edges of P are precisely those length-two pairs together with {1,e}.
Consequently

    tau_P(S_n) = converse(U_n).

Both source and target therefore have minimum cycle debt0, but

    f(U_n)=n−2−floor(n/4),                           (A1)
    f(tau_P S_n)−f(S_n)=n−3−floor(n/4).              (A2)

**Proof of(A1).** The directed triangles of U_n are exactly

    {i,i+1,i+2}, 0≤i≤n−3,
    {i,i+2,i+4}, 0≤i≤n−5.

Indeed write the two consecutive gaps in an arbitrary sorted triple as a,b.
A triangle occurs exactly at(a,b)=(1,1) or(2,2); checking the rule “only
length-two pairs reverse” proves both directions.
Assign weight1/2 to each first-kind triangle and1/4 to each second-kind
triangle. Each edge belongs to triangles of total weight at most1:
a length-one edge meets at most two first-kind triangles; a length-two
edge meets one first-kind and at most two second-kind; a length-four edge
meets at most one second-kind; all others meet none. Every order has at
least one backward edge in every directed triangle, so its number of
backward edges is at least the total triangle weight

    (n−2)/2+(n−4)/4 = 3n/4−2.

Integrality yields f(U_n)≥ceil(3n/4−2)=n−2−floor(n/4).
For the reverse bound split the natural order into consecutive blocks of
four and a remaining block of size<4. Inside each full block[a,a+1,a+2,a+3]
use order[a+2,a,a+3,a+1]. Its cost is1 instead of the natural2. Keep the
order of all blocks and the remainder fixed; every cross-block relative pair
is unchanged. The natural total cost was n−2, and there are floor(n/4)
independent savings. QED.

**First finite counterexample to the tempting half-size extrapolation.**
At n10, S has only exception{1,8}, P=(0,2,4,6,8,1,3,5,7,9), source f1,
target f6, so |Δf|=5>floor((10−1)/2)=4. The main packet's exhaustive n≤7
maxima0,1,1,2,2,3 do not justify that stronger all-order bound.
The main theorem's n−3 tariff is unaffected and not claimed sharp.

This family proves a linear-in-size complement charge is genuinely needed
if the charged quantity is feedback-arc count. It simultaneously warns that
edge debt and cycle debt are different observables: beta remains0 here.
All paths and the single source exception are literal named objects; no
relabelled-state averaging or source-weight substitution is used.

Reproduce with `python3 verify_scaling.py > receipt.json` and
`python3 -O verify_scaling.py > receipt.optimized.json`, compare byte-for-byte,
and run `sha256sum -c SHA256SUMS`. Full subset DP checks n5…18; literal
triangle classification/packing and explicit rank certificates check n5…40.
