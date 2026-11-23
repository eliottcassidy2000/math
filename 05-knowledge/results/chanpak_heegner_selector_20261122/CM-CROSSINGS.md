# One residue and one square root recover two independent binary runs

**PROVED with finite-exact source data.** Let \(J_d=744-j((1+i\sqrt d)/2)\) for the five discriminants below. At \(d=19,43,67,163\), \(J_d-1\) and \(J_d+1\) occur on strong tournaments differing by one arc. Each crossing belongs to a family with two independently variable positive widths. Given an odd target, one residue and at most one integer square root decide membership on the union of these four families and recover the exact graph. The scope of the inverse is the displayed union. [CM values and the near-integer](https://people.maths.ox.ac.uk/greenbj/papers/ramanujanconstant.pdf) · [Proofs for the other discriminants](README.md)

| \(d\) | \((-j)^{1/3}\) | Exact midpoint | Right address \(y\) | Strong order |
|---:|---:|---:|---:|---:|
|19|96|885,480|736|19|
|43|960|884,736,744|2,438,701|30|
|67|5280|147,197,952,744|650|31|
|163|640,320|262,537,412,640,768,744|1,086,757,889,205,407|58|

The last midpoint is the integer nearest Ramanujan's constant. At discriminant 11 the CM expression and rounding give 33,512 and 33,506 respectively. The construction of a tournament below takes the specified integer as input. The CM calculation identifies the targets.

## The literal reversals

Number the core vertices from zero, with bit \(j\) in row \(i\) denoting \(i\to j\). An apex beats the vertices in the tag mask. Put a transitive chain of length \(\operatorname{bitlength}(y)\) after the core, with all core-to-chain arcs forward. The apex contacts on the chain are the reverse complement of the ordinary binary expansion of \(y\), with a one denoting an outgoing apex arc. Reverse just the indicated core pair to get the upper graph.

|\(d\)|Lower outgoing rows|Tag mask|Pair|\(h,c,K\)|
|---:|---|---:|---|---|
|19|6,68,24,19,35,143,61,95|40|1,2|233,975,1051|
|43|4,5,56,3,11,27,63|80|2,3|23,356,337|
|67|72,37,97,22,71,153,170,31 plus receiver below|5 on core, 1 on receiver|0,2|604183,225699127,224054575|
|163|6,32,10,35,15,85,31|11|0,1|59,240,197|

On the selected parity of \(y\) (even, odd, even, odd respectively), the lower count equals \(K+cy+h f(y)\). The upper count has the same coefficients and constant \(K+2\). Here \(f(0)=0, f(2v+1)=f(v), f(2v)=f(v)+v+[v\text{ odd}]\). The four sets of coefficients come from the tagged subset-path recurrence. A separate verification reconstructs the actual graphs, checks strong connectivity and the arc, splits each Hamilton path at its apex and counts directly on the ordered blocks. For 163 a common source in the original eight-vertex core could be removed because it starts every path. No claim of minimum order is made for these four pairs. [Binary star construction](../chanpak_gap_star_ports_20261122/README.md) · [Independent block counter](../chanpak_decimal_kernel_20261122/e15_pair.py) · [Replay and target decoder](cm_crossings.py)

### How the Fibonacci receiver lands at 67

Insert the twelve-vertex receiver with rows \(1126,4092,4088,4081,4065,4032,3968,3841,3585,3073,2048,1\) immediately after the eight-vertex core and before the chain. All interblock arcs point forward. Its sole tag is vertex zero. It has 1,187 full Hamilton paths, of which 1,186 start outside the tag. On the eight-vertex core the heights with an apex are 1,415 and 1,417. The ordinary count 509 and the unsigned tagged start read 1,322 agree, while the signed start read changes by two in the convention of the receiver theorem. Consequently the joined gap is \(2(1187-1186)=2\) for *every* right address. At the even address 650 the lower equation is
\[
 224054575+225699127\cdot650+604183\,f(650)=147197952743 .
\]
This illustrates what is gained by selecting a receiver with one tagged start. The remaining 1,186 paths contribute to the slope that sets the absolute count. It also explains how the witness was found. Generate strong rooted Fibonacci receivers, keep small arc pairs with this signed change, and invert the slope equation at the requested integer. The script stores the selected receiver explicitly. [Receiver theorem and grammar](../chanpak_e31_physical_20261122/README.md) · [Recursive receiver generator](../chanpak_e61_fibonacci_target_20261122/check.py)

**Choosing the receiver from the target itself.** In fact the twelve-vertex receiver here has a transitive root deletion. Let \(v\) be the binary apex contacts on that ordered chain, leading one first, and set \(H=1+f(v),\ U=v+1,\ E=[v\text{ even}]\), \(m=2^{|v|}-\operatorname{rev}_{|v|}(v)\). The two-slot read \(Q\) is computed by four states which track whether a before-apex last contact is valid and an after-apex part has started. If the original core has reads \(h,c,e,k\), the joined lower chart becomes
\[
 K_v+(cQ+hU)y+hHf(y)+hE[y\text{ odd}],\qquad
 K_v=Hk+h(U-E)+c(Q-H-m)+e(H-1).
\]
This is the ordinal join applied to the [word receiver's explicit profile](../chanpak_e31_physical_20261122/README.md). The other orientation adds two for every \(v,y\), since the leading one forces exactly one tagged full start. For midpoint \(M\), only words \(1\le v\le(M-1)/h\) can work: the other summands of \(K_v\) are nonnegative and \(U-E\ge v\). For each word and suffix parity the linear residue modulo \(2\gcd(cQ+hU,hH)\) eliminates incompatible targets, and the [exact suffix inverse](../chanpak_gap_star_ports_20261122/replay.py) decides the survivors. This tests the whole one-root-over-a-chain subclass on this fixed core, stopping at its first solution. On input **147197952744 alone** it returns the first receiver \(v=1634\) and suffix \(y=650\), with no supplied receiver, suffix, or width bound. [Runnable target-to-receiver search](receiver_from_target.py)

## Two freely variable widths through four CM integers

Read the binary address as a fixed prefix of value \(w\), followed by \(a\) ones, \(b\) zeros and \(\ell\) final ones, where \(a,b\ge1\). Write
\[
 P=w+1,\quad S=2^\ell c+h,\quad
 C=K+1+c(2^\ell-1)+h(f(w)+1).
\]
The exact midpoint is
\[
                    M(a,b)=C+(S2^b-h)(P2^a-1).
\]
Indeed the first run gives the odd value \(u=P2^a-1\) with unchanged digital read. Appending \(b\) zeros increases that read by \((2^b-1)u+1\). The last ones preserve it. Substitution proves the formula and the same physical edge works at every pair of widths.

|\(d\)|\(w\)|\(\ell\)|\(h\)|\(S\)|\(P\)|\(C\)|Widths at \(J_d\)|
|---:|---:|---:|---:|---:|---:|---:|---|
|19|2|0|233|1208|3|1751|3,5|
|43|152418|1|23|735|152419|2537019|2,1|
|67|162|0|604183|226303310|163|297160719|1,1|
|163|4245148004708|5|59|7739|4245148004709|212736837811695|1,2|

**Modular selection.** Take \(m=26,896,077=27\cdot13\cdot19\cdot37\cdot109\), a divisor of \(2^{36}-1\). Compute \(M(a,b)+\epsilon\) for the four rows, \(a,b\bmod36\) and \(\epsilon=\pm1\). All \(4\cdot36^2\cdot2=10,368\) residues are distinct. The certificate enumerates those entries and can test every smaller divisor of \(2^{36}-1\) to show this is the least separating divisor. This extends the earlier [five one-width rays](../chanpak_ternary_buckets_20261122/PHASE-DECODER.md) to two independent widths. Modulus 37 has a separate role in the [Heegner-prime tournament](HEEGNER.md). Here its role is verified by the residue table for these particular source coefficients.

**Absolute recovery.** For an odd input \(T\), an occupied residue specifies \(d,\epsilon,a\bmod36,b\bmod36\). Put \(Z=T-\epsilon-C\). Each row satisfies \((2S-h)(2P-1)>2SP\), so
\[
                SP2^{a+b-1}<Z<SP2^{a+b}.
\]
Consequently \(L=a+b\) is forced by the bit length of \(\lfloor Z/(SP)\rfloor\). With \(X=2^a\), expansion gives
\[
       hP X^2-(SP2^L+h-Z)X+S2^L=0.
\]
Demand a nonnegative square discriminant, an integral power-of-two root, positive complementary width and the two specified phases. Exact substitution gives both necessity and sufficiency. Within a row at most one root can work, because two power-of-two roots would force their product \(S2^L/(hP)\) to be a power of two. The ratios \(1208/699,735/(23\cdot152419),226303310/(604183\cdot163),7739/(59\cdot4245148004709)\) all have a surviving odd factor in numerator or denominator. The four endpoints from the two [isolated order-11 edges](INVERSE.md) also have four distinct residues outside this table and can be tested by exact equality. Across rows and signs uniqueness follows from the modular table. Thus this grammar gives each of its odd outputs a unique address and sign. Other constructions and other isomorphism classes may realize the same height.

As a sharp test of the domain, neither neighbor of \(2^{p-1}(2^p-1)\) belongs to these four families for **any** positive integer \(p\). Reducing \(p\) modulo 36 suffices since 2 is a unit, and none of the 72 tested residues is occupied. In particular increasing the two run lengths in these Ramanujan constructions cannot hit exponent 61. Membership in the full Chan–Pak spectrum stays open outside established results.

## How many different integer heights and ordinary classes this supplies

Every pair of positive widths gives two distinct strong isomorphism classes. To see strongness at every width, the nonempty tags reach the entire deletion core by core arcs in each of the eight certified orientations. This already follows from strongness at the four displayed widths, since no suffix vertex sends a forward arc back into the core. Every core vertex reaches the last chain vertex, which beats the apex because the leading bit of the address is one. All chain vertices are reachable from the core and reach the apex through the last vertex. Uniqueness of the decoded height now distinguishes different width pairs, endpoints and sources even after all markings are forgotten.

In the four rows the total order is \(n=n_0+a+b\) with \(n_0=11,27,29,55\), respectively. Including the four exhibited order-11 endpoints, the size series of these particular ordinary classes is exactly
\[
 4z^{11}+\frac{2(z^{13}+z^{29}+z^{31}+z^{57})}{(1-z)^2}.
\]
The number of odd endpoints through height \(X\) on this grammar is \(4(\log_2 X)^2+O(\log X)\): the absolute height window bounds \(a+b\) by \(\log_2 X+O(1)\) in each row, and there are \(L-1\) positive ordered pairs of sum \(L\). Thus the exact decoder's domain has density zero. For any one of its endpoints at order \(n\), the full isomorphism-class series at that height dominates \(z^n/(1-z)^2\) coefficientwise by independently adjoining source and sink strings to the strong core. The distinct converse witnesses at the two order-11 edges double this latter lower bound there. [General isomorphism-fiber formula](../chanpak_arc_fiber_grammar_20261122/README.md)

### The complete class series in each of the 10,368 residues

For an occupied endpoint residue \(r\) modulo 26,896,077, let \(a_0,b_0\in\{1,\ldots,36\}\) be its two positive phase representatives and let \(n_0\) be the order offset of its selected core. Inside this grammar, the entire size series at residue \(r\) is
\[
                   G_r(z)=\frac{z^{n_0+a_0+b_0}}{(1-z^{36})^2}.
\]
Indeed \(a=a_0+36i, b=b_0+36j\) with independent nonnegative \(i,j\), and distinct choices give distinct heights and classes. Consequently at order \(n=n_0+a_0+b_0+36k\) there are exactly \(k+1\) such classes. At every other order the coefficient is zero. Each of the four isolated residues instead has the single term \(z^{11}\).

This yields an exact uniform sampler of the **marked** constructions with a specified endpoint residue and order. Draw \(i\) uniformly from \(0,\ldots,k\) and put \(j=k-i\). If labels are wanted, apply a uniform permutation of the vertices. Every marked pair has trivial automorphism group: its automorphisms act freely on the Hamilton paths of each completion, so their order divides \(\gcd(H,H+2)=1\). Thus each chosen pair has exactly \(n!\) labelings. To sample on a union of compatible residues, first choose a residue with weight its computed coefficient at order \(n\). This addresses the same source of bias as sampling Latin squares by construction histories: different partial choices have different numbers of completions. Here the exact coefficient gives the correction explicitly. [Marked-edge rigidity](../chanpak_arc_fiber_grammar_20261122/README.md) · [Mira Cult's exact Latin-square sampling proposal](https://github.com/Mira-Cult/math/tree/30ddaaad27b48df2002225c3e3f40cffd02c5253/notes/combinatorics/exact-latin-square-sampling)

From this directory, run the following for the finite certificate, literal midpoint counts, minimum separating divisor, and a target-only reconstruction.

~~~sh
python3 -B receiver_from_target.py 147197952744
python3 -B cm_crossings.py --minimal
python3 -B cm_crossings.py 262537412640768743
~~~
