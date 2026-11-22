# Two-ended suffix inversion at the top order

**PROVED.** The cyclic triangle receiver cannot straddle the perfect-shaped center at **any odd exponent**, at any address or total order, with the inherited six-vertex marked core. At exponent 61 the 123-vertex problem therefore has only two receiver types left. Their arbitrary suffixes remain **OPEN**. In each, both the number of bit changes and the number of equal consecutive bit pairs must be at least 22. These tests refer to the construction with a fixed marked core. [Three original charts and the order bound](SELECTOR.md)

## Remove a receiver before choosing its suffix

Write the prescribed lower count as \(N\), the address as \(y\), its parity as \(j\), and the chart as \(K_j+cy+hf(y)\). Every contribution to \(f\) is even. On writing \(y=2u+j\), a necessary condition is
\[
                 2\gcd(c,h)\mid N-K_j-cj.
\]
This is the exact affine congruence after relaxing \(f(y)/2\) to an arbitrary integer. For the triangle the chart is \((c,h,K_0,K_1)=(411,75,367,392)\), so its two possible residues are 1 and 2 modulo three. For every odd \(p\) the lower target \(2^{p-1}(2^p-1)-1\) is divisible by three. This proves the unrestricted obstruction for the triangle.

At exponent 61 this test rejects both parity branches for 214 of the 610 normalized strong receivers through order nine at **all** widths. Another 24 survive only on the even branch, and 372 survive on both. Of the 214 full rejections, 99 already hold modulo three for every odd exponent. These counts include the triangle. Run `python3 -B suffix_filter.py` from this directory.

## Squeezing both ends while reading just the low bits

For any fixed low word write \(y=Bz+r\), \(B=2^b\), \(b\ge1\), and put \(v=\operatorname{rev}_b(r)\), retaining all leading zeros. The earlier [binary concatenation formula](../chanpak_gap_star_ports_20261122/README.md) gives
\[
  T := N-K_j-cr-hf(r)
   =Dz+h\big(f(z)+[r<B/2](z\bmod2)\big),\qquad D=Bc+h(B-1-v).
\]
Indeed the zero bits in the suffix contribute \((B-1-v)z\), and the correction at the branch between the two words contributes \([r<B/2](z\bmod2)\). The recurrence gives \(0\le f(z)\le z-[z\text{ odd}]\) by induction. This puts the quotient in the **closed** interval
\[
 \left\lceil T/(D+h)\right\rceil\le z\le\left\lfloor T/D\right\rfloor,\qquad Bcz\equiv T\pmod h .
\]
Intersect with the required bit width and move both interval endpoints inward to the quotient residue class. Except at a terminal case with \(b=k\), \(z\) has exactly \(k-b\) bits. High bits common to these new endpoints are then forced. Every bit change already present in the low word and every bit change already present inside the common high prefix is permanent. A restriction on the total number of changes can therefore prune before the middle bits are known. For a prescribed fixed reference word take the exclusive-or with the reference before counting changes. At small remaining interval size directly test all quotients in the residue class using the original chart.

The two extensions of a low word need no fresh evaluation of \(f\). Appending a new most significant bit \(d\) to that word replaces \((r,v,f(r))\) by \((r+dB,2v+d,f(r)+d(B-1-v+[r<B/2]))\). These disjoint branches cover all allowed words. The replay compares with exhaustive uncompressed solutions at widths 4, 7, 10 and 13, including attainable targets.

**FINITE-EXACT.** At top order the unrestricted triangle is already impossible. For the singleton chart \((81,25,43,68)\) and four-vertex chart \((1021,125,915,940)\), both of the following searches return zero exact flips. The run cap is 22 including the first run.

| Width | Reference | Prefixes examined | Final quotient candidates |
|---:|---|---:|---:|
| 115 | all zeros | 642,985,660 | 544,270 |
| 112 | all zeros | 403,902,382 | 6,005,087 |
| 115 | alternating bits | 2,274,355,982 | 5,974,215 |
| 112 | alternating bits | 52,518 | 0 |

Exclusive-or with an alternating word exchanges changes and equal consecutive pairs. For width \(k\) these two disjoint exclusions cover exactly \(2\sum_{i=0}^{21}\binom{k-1}{i}\) distinct words starting with one, including words that already violate an earlier absolute bound. This is 110,334,209,188,095,564,170,048 words at width 115 and 60,086,979,220,266,461,252,416 at width 112. As a separate test, no address has at most 18 runs after exclusive-or with the periodic word 0011 read from the low end. The search has not exhausted the middle transition counts. [Exact two-ended replay](suffix_runs.cpp)

From this directory run:

~~~sh
c++ -O3 -std=c++17 suffix_runs.cpp -o /tmp/chanpak_suffix_runs
/tmp/chanpak_suffix_runs 22
/tmp/chanpak_suffix_runs 22 01
/tmp/chanpak_suffix_runs 18 0011
~~~

## How many potential addresses share the correct decimal suffix?

Let \(k\) be the prescribed width. Take every fixed-width address in the necessary interval
\[
 \max\{2^{k-1},\lceil(N-\max K_j)/(c+h)\rceil\}
 \le y\le
 \min\{2^k-1,\lfloor(N-\min K_j)/c\rfloor\}.
\]
For the four-vertex chart require also the five forced terminal ones. In this universe the exact numbers of addresses whose lower count ends in **2175**, the suffix of \(E_{61}-1\), are:

| Receiver | Matching addresses |
|---|---:|
| Singleton, width 115 | 1,548,134,166,998,502,827,632,877,103,933 |
| Four vertices, width 112 | 47,676,574,974,289,207,102,343,161 |

These are counts of decimal matches prior to the run tests. For the calculation track \(y\bmod 10000\) and \(f(y)/2\bmod(10000/\gcd(2h,10000))\) while appending bits at the right. A one keeps the second state and a zero adds \(\lceil y/2\rceil\) to it. Its needed residue is determined by \(y\bmod10000\) here since twice the second modulus divides 10000 in both charts. A digit DP below each endpoint counts the admitted words once. The control computes direct small-width counts in addition to the full outputs. Run `c++ -O3 -std=c++17 suffix_decimal.cpp -o /tmp/chanpak_suffix_decimal && /tmp/chanpak_suffix_decimal`. [Decimal replay](suffix_decimal.cpp)

## The same two top receivers at every odd exponent

**PROVED.** For every odd \(p\ge11\), every normalized exact construction with the six-vertex core would have order at most \(2p+1\). At that order the only possible receivers are the singleton and the rooted four-vertex tournament \((2,12,9,1)\). The four-vertex receiver still requires eight initial zeros after the leading one and five terminal ones. Primality of \(2^p-1\) is unnecessary.

Use the inherited fixed-order bound \(n\le4(N+10)/(28\cdot2^{m-7})\) for a nonsingleton receiver of order \(n\) in a construction of total order \(m\). Here \(N+10<2^{2p-1}\). At \(m\ge2p+2\) this gives \(n<16/7\), impossible for a strong nonsingleton receiver. The singleton address at any such order has \(y\ge2^{2p-7}\); its slope 81 alone exceeds the target. At \(m=2p+1\) the same inequality gives \(n<32/7\). The triangle has just been excluded. The three rooted order-four profiles have slopes 1,021, 1,046 and 1,071. For these the highest bit is \(Y=2^{2p-11}\), so any slope of at least 1,024 overshoots. For slope 1,021, writing \(y=Y+u\), positivity of the intercept gives \(u/Y<3/1021<2^{-8}\) and \(f(y)<3Y/125\). A zero among the lowest five bits would contribute at least \(Y/32\) by the binary recurrence. This proves the claimed digit conditions for the whole exponent family.

## Exact Boolean formulation of the unrestricted remaining words

The [unrestricted top-order census](TOP-ORDER-CENSUS.md) closes these two charts for every odd exponent from 11 through 37. It also gives exact strong edges at neighboring centers, including an edge centered at \(E_{13}+2\).

Writing \(b_i\) for the bit at position \(i\), introduce an integer \(A_d=\sum_j b_{j+d}(1-b_j)\) for every distance \(d\ge1\). The counting recurrence gives \(f(y)=\sum_{d\ge1}(2^{d-1}+[d=1])A_d\). Introduce a Boolean for each conjunction and sum the integer equation \(K_0+(K_1-K_0)b_0+cy+hf(y)=N\) in base 128 using seven-bit limbs and explicit nonnegative carries. The [OR-Tools CP-SAT model](suffix_cp.py) also imposes the proven terminal congruence and high/low prefix interval. It leaves every admissible suffix available. Any returned address is checked against the original recurrence. The planted 18-bit control returned 131,615 for each chart.

With 70 seconds per unrestricted chart, CP-SAT returned **UNKNOWN** for both. An independent Boolean carry-save encoding in CaDiCaL likewise returned unknown with a 20,000-conflict budget per chart. These runs settle neither existence nor nonexistence. To reproduce the CP-SAT attempt use `python3 -B suffix_cp.py four --seconds 70` and substitute `singleton` for the other chart. The `--export` option writes a solver model. With the time option omitted the solve has no set time limit.

The supplied paper by Pasten and Vidaux studies Diophantine definability and gives a conditional negative answer to Hilbert's tenth problem over quadratic number fields assuming \(abc\) over that field. For these particular finite receivers there is an unconditional exact integer encoding with a known bound on every variable. Write \(G(b)\) for the quadratic height expression above minus the target. The one quartic equation
\[
 G(b)^2+(b_{n-1}-1)^2+\sum_{i=0}^{n-1} b_i^2(b_i-1)^2=0
\]
has exactly one integral solution per admissible full-width address satisfying the chart. All summands are nonnegative on integers, so every solution lies in the Boolean cube and has leading bit one. This identifies exactly what an integer decision procedure would need to settle at widths 115 and 112. The complete smaller-width tests exploit the quotient interval. [Pasten and Vidaux, original paper](https://arxiv.org/abs/2611.24748) · [Exact smaller-width census](TOP-ORDER-CENSUS.md)
