# Which single-arc constructions can follow a perfect-number center?

**PROVED.** Write \(A=2^{p-1}\) and \(E_p=A(2A-1)\), with \(p\ge3\) an integer. In the one-apex spike construction, the consecutive unit-gap centers on opposite sides of \(E_p\) are exactly
\[
                        E_p-A/2,\qquad E_p+5A/2 .
\]
The physical star seam at exponent 31 escapes this interval. Its own digit addresses also have an obstruction to extrapolation: a fixed finite list of binary words, repeated affine numbers of times, cannot carry that same marked-core edge through \(E_p\) for unbounded exponents on an arithmetic progression. The proofs give a necessary numerical test on any one-sided analogue.

## 1. The exact gap for a bare transitive spine

Use the established single-site formula. Orient the fixed spine forward, and let the apex contact word be \(L b R\). Let \(u\) encode \(L\) in binary and \(w\) encode the reverse complement of \(R\). Put \(f(0)=0,\ f(2v+1)=f(v),\ f(2v)=f(v)+v+[v\text{ odd}]\). In particular \(0\le f(v)\le v\), by induction. The counts for \(b=0,1\) are
\[
 2uw+f(u)+f(w)+1+u+[u\text{ odd}],\quad
 2uw+f(u)+f(w)+1+w+[w\text{ odd}].
\]
The unit gap means the absolute difference is two. Exchanging the arms if necessary, take the first count as the lower; now \(w+[w\text{ odd}]=u+[u\text{ odd}]+2\) and \(u+1\le w\le u+3\). The center \(C\) exceeds the first count by one. [Single-site formula and physical words](../chanpak_apex_flip_local_global_20261122/README.md)

If \(u\le A-3\), then \(w\le A\) and
\[
 C\le 2(A-3)A+(A-3)+A+(A-2)+2 < E_p-A/2.
\]
For \(u=A-2\) the only possible values of \(w\) are \(A-1,A\). Use \(f(A-2)=A/2,\ f(A-1)=0,\ f(A)=A\) to obtain the larger center \(E_p-A/2\) at \(w=A\). For \(u=A-1\), the two possibilities are \(w=A+1,A+2\). Since \(f(A+1)=A/2\), the first gives \(C=E_p+5A/2\). The second gives \(C\ge2(A-1)(A+2)+A+2>E_p+5A/2\). Finally, \(u\ge A\) forces \(w\ge A+1\), so \(C\ge 2A(A+1)+A+2>E_p+5A/2\). Both extremizing address pairs are positive and therefore construct strong tournaments. This proves the consecutive-center claim for **every** exponent, regardless of primality.

There is also an exact boundary if the two tournaments admit the same apex and its deletion is transitive on both sides. A flipped spine arc then swaps adjacent spine vertices, since any intermediate vertex would make a directed triangle after reversal. Equal apex contacts give no change. Otherwise regard the swap as \(L01R\to L10R\). The one-apex count is one plus the sum of \(g(d)=2^{d-1}+[d=1]\) over inversions at distance \(d\). The new inversion contributes 2. Earlier ones and later zeros contribute distance increments \(g(d+1)-g(d)=2^{d-1}-[d=1]\). Thus the increase is \(2+u-[u\text{ odd}]+w-[w\text{ odd}]\). A unit gap requires \(u,w\in\{0,1\}\), leaving at most four spine vertices after discarding initial zeros and terminal ones. These swaps produce bounded centers (at most 12 by the inversion formula), below the perfect-number intervals for \(p\ge3\).

## 2. Every nonconstant periodic word contributes positive drift

Consider a binary address, possibly padded by leading zeros. A template consists of a finite concatenation of fixed words and powers \(V^{an+b}\) of fixed nonempty words with integer \(a>0\) and eventually nonnegative exponents; different blocks may have different words, slopes and shifts. Restricting \(n\) to a progression permits rational affine exponents whenever they are integral on that progression. An exponent \(p=p_0+kn\) of the perfect-number family may advance by any integer \(k>0\).

For a word \(V\) of length \(l\), let \(K=2^l\), \(d=\mathrm{bin}(V)\), \(z=\mathrm{bin}(\mathrm{revcomp}(V))\), and let \(e\) indicate that its first bit is zero. Appending one copy to an integer address \(w\) gives
\[
          w'=Kw+d,\qquad f(w')=f(w)+zw+f(d)+e[w\text{ odd}].
\]
For the empty suffix the corresponding update is the identity. This is the binary split identity for nonempty words. After the first copy the incoming parity is fixed by the last bit of \(V\). Writing \(\epsilon=[d\text{ odd}]\), the drift per extra copy in \(f(w)-z w/(K-1)\) is
\[
              \beta(V)=f(d)+e\epsilon-\frac{dz}{K-1}.
\]
[Exact digit splitting](../chanpak_gap_star_ports_20261122/README.md)

**Lemma.** The drift is zero for a constant word and positive for every nonconstant word. In fact, if \(N_1\) is the number of cyclic clockwise \(1,0\) transitions, then \(\beta(V)\ge l N_1/(2^l-1)\).

**Proof.** Let \(N_j\) count ordered pairs of a 1 followed by a 0 at cyclic displacement \(j\). The extra term \(e\epsilon\) completes the adjacent inversion count of \(V\) to \(N_1\). Expanding \(dz\) by pairs of binary digits, and subtracting the ordinary within-word exponential inversion sum in \(f(d)\), leaves weight \(2^{j-1}/(2^l-1)\) for every such cyclic pair. For a nonwrapping pair this uses \(2^l/(2^l-1)-1=1/(2^l-1)\); for a wrapping pair only the product term is present. Consequently
\[
 \beta(V)=N_1-\frac{\sum_{j=1}^{l-1}N_j 2^{j-1}}{2^l-1}.
\]
Cyclic binary symmetry gives \(N_j=N_{l-j}\) since the two transition directions across displacement \(j\) occur equally often. A pair counted by \(N_j\) crosses at least one of the \(N_1\) adjacent transitions along its \(j\) steps, giving \(N_j\le jN_1\). Therefore the numerator is at most \(N_1\sum_{j=1}^{l-1}j2^{l-j-1}=N_1(2^l-l-1)\). This proves the inequality and its asserted sign. ∎

It follows that the address of a template is a finite rational linear combination of \(2^{sn}\) (allowing \(s=0\)), while its \(f\)-value is such a combination plus
\[
                             n\sum_V a_V\beta(V).
\]
To see this, process each repeating block in order. After its first copy, the change in \(f-zw/(K-1)\) per copy is exactly \(\beta(V)\), and the new address after all copies is an affine function of the old address with factors \(K^{an+b}\). Incoming parity is eventually fixed by the last nonempty word. A fixed intervening word retains the incoming \(f\) with coefficient one. This proves the formula by induction, including repeated occurrences of the same word in different positions.

In a fixed two-flank marked-core chart the output is a sum of a constant, linear address terms, \(Qxy\), eventually constant parity terms and \(h(f(x)+f(y))\), where \(h>0\). Thus **no such pair of templates containing a growing nonconstant period can realize \(E_{p_0+kn}+C\) infinitely often for fixed \(C\)**. The sum of the two drifts stays strictly positive. A finite sum of distinct real exponentials minus the proposed target either has a largest nonzero exponential term, or reduces to a constant; in neither case can adding a positive linear term give infinitely many zeros. This includes templates whose two flanks grow at different rates.

## 3. Rigidity of the exponent-31 edge

For the first 62-vertex exponent-31 witness the lower and upper apex counts are 1021 and 1023. Both deletion counts are \(h=125\), both two-slot counts are \(Q=9728\), both start reads are \((523,548)\), and the end reads are \((3769,3884)\) and \((3784,3899)\). The signed end differences are equal. Substitution in the two-flank formula shows
\[
                 H_{\rm upper}(x,y)-H_{\rm lower}(x,y)=2+30x .
\]
In particular this flip gives a gap of two precisely when \(x=0\). The freely varying right address then has the adjusted-intercept form \(K_0+1071y+125 f(y)\) for a fixed parity. [Two-flank formula](../chanpak_gap_star_ports_20261122/README.md)

**Power-of-two interval test.** For any fixed integers \(c,h>0\), a one-sided output \(K_0+cy+h f(y)\) from a template as defined above can have a pure dyadic leading coefficient (in base \(2^n\)) only if \([c,c+h]\) contains a power of two, provided all growing blocks are monochromatic. Redundant leading zeros may be discarded; assume \(y\to\infty\).

Indeed a final fixed suffix of length \(l\) follows the last growing run. Immediately after that run the ratio \(f(w)/w\) tends to \(\delta\), where \(\delta=1\) for a growing zero-run and \(0\) for a growing one-run. A zero-run of length \(j\) replaces the address \(v\) by \(2^jv\) and adds \((2^j-1)v+[v\text{ odd}]\) to \(f\). A one-run leaves \(f\) unchanged and replaces the address by \(2^j(v+1)-1\). These facts and \(f(v)\le v\) give the limit for arbitrary incoming addresses, after removing redundant zero-runs. The trailing fixed suffix makes \(f(y)/y\to(\delta+z)/2^l\), by the split formula above. For an empty suffix take \(l=z=0\). The leading coefficient of \(y\) is dyadic rational, since each fixed block and each monochromatic affine run updates the address by multiplication by a dyadic monomial and addition of a dyadic polynomial. Thus the leading coefficient of \(cy+h f(y)\) is this positive dyadic number times \((c2^l+h(\delta+z))/2^l\). A pure dyadic power forces the latter numerator to be a power of two; dividing by \(2^l\) gives a power of two in \([c,c+h]\). If \(h\) is odd this sharpens to **one of the two endpoints**. Indeed that integer power \(2^t\) must satisfy \(2^l(2^t-c)=h(\delta+z)\). Therefore \(h\) divides \(2^t-c\), which lies between 0 and \(h\).

For the actual exponent-31 unit edge, \([1071,1196]\) lies strictly between 1024 and 2048. Together with the drift lemma this proves that **no fixed finite template of repeated right words** extrapolates this physical edge into exact perfect-number centers on an unbounded arithmetic progression of exponents, even if one requests hits just on an infinite subset of that progression. At each new exponent the digit pattern or the marked core must escape this prescription. This result leaves open infinitely many solutions by arbitrary right addresses of the same core.

The same conclusion holds for the other two witnesses at exponent 31. In the second 62-vertex example the signed start and end data again agree, and the slope of the difference in \(x\) is 30. The 60-vertex example has slope 66; its end reads pass from \((168696,170181)\) to \((168729,170214)\). In all three the two-slot counts agree, and the signed gap at \(x=0\) is two. Their complete one-sided tests are

| Witness order | \(c\) | \(c+h\) | Gap across its arc |
|---:|---:|---:|---|
| 62 | 1071 | 1196 | \(2+30x\) |
| 62 | 2267 | 2392 | \(2+30x\) |
| 60 | 53334 | 55083 | \(2+66x\) |

Each physical flip thus forces the external left address to vanish at a unit seam, and every one of its one-sided charts fails the test. The twenty earlier perfect-band witnesses at exponents 17 and 19 also all fail the endpoint test when their own left addresses are held fixed. Here the one-sided linear coefficient is \(c=Qx+e+o\) and the additive intercept may depend on \(x\) and the right parity. The check obtains all twenty values from the stored core and the actual left address. These assertions concern those specified templates. The constructive mechanism can change the internal receiver as the exponent grows. [All three physical exponent-31 pairs](../chanpak_e31_physical_20261122/README.md) · [Twenty prior band certificates](../chanpak_gap_star_ports_20261122/PERFECT-BAND.md)

The companion checker enumerates cyclic words, verifies the two exact one-site formulas against literal smaller tournaments, and checks the local exponent-31 reads. Run `python3 check.py` in this directory.

## 4. Two homogeneous growing runs already fail on a single apex

**Theorem.** No fixed one-apex word with at most two variable homogeneous runs of affine lengths can have count \(E_p+C\) for infinitely many exponents along a nonconstant progression, for any fixed additive constant \(C\).

The inherited two-run expansion has form \(C_0+UX+VY+WXY\), with \(X,Y\) powers of two from the two widths. A zero \(W\) leaves at most one active run. A single active run cannot supply both exponentials of \(E_p\). Otherwise comparing the two leading target rates forces both widths to have the rate of \(p\). Indeed the slopes of the widths add to the leading slope \(2p\), and an active secondary term must have the slope of \(p\). Write the widths as \(p-a,p-b\). The shifts may be negative and can absorb any fixed integral changes of width. [Two-run expansion and coefficient table](../chanpak_moonshine_prime_blocks_20261118/ONE_APEX_PRIME_SLOPE.md)

Let \(M\) be the fixed intervening word, \(K=2^l\) with \(l=|M|\), \(d=\mathrm{bin}(M)\), and \(z=K-1-\mathrm{rev}_l(d)\). The leading coefficient \(W\) must be a power of two. The outer address factors of \(W\) in the inherited table are thus also powers of two. In case 01 the two linear coefficients are positive. In case 10 the outer addresses plus one are powers of two, so both outer strings merge into the first one-run and final zero-run after removing redundant leading zeros and trailing ones. Now \(W=K\), forcing \(a+b=l+1\), and the magnitude of the negative linear coefficient at \(2^p\) is
\[
            (1+\mathrm{rev}_l(d))/2^a+(K-d)/2^b .
\]
For \(d=0\) the product of the two numerators is \(K\). Otherwise, if \(d\) has \(t\) trailing zeros, the reversed integer is at least \(2^{l-t-1}\) and \(K-d\ge2^t\), so the product is at least \(K/2\). AM–GM gives a sum at least one, exceeding the required \(1/2\). The same inequality applies for empty \(M\) with \(K=1\).

In case 00, the left outer address is a power of two and the right reverse-complement address plus one is a power of two. Absorb zeros to normalize the word to \(1\,0^{p-a} M 0^{p-b}\); initial zeros and terminal ones do not change its count. Unless both variable runs merge, trim \(M\) to begin and end in 1 and absorb the trimmed zeros into the widths. The linear coefficients of the two runs are then \(-\mathrm{rev}_l(d),d\), and the leading coefficient is \(K\). Thus
\[
        2^b\mathrm{rev}_l(d)-2^a d=2^l,\qquad a+b=l+1 .
\]
Both read integers are odd. For unequal shifts the 2-adic valuation on the left is \(\min(a,b)\), which cannot be \(l\) when \(a+b=l+1\), \(l\ge1\), and \(a\ne b\). For equal shifts \(l\) is odd and the equation requires \(\mathrm{rev}_l(d)-d=2^{(l-1)/2}\). Reversal of an odd-length word preserves its base-two value modulo 3, which excludes this equality. Converse followed by reversing the spine exchanges 11 and 00. This exhausts the cases.

Three growing runs can suffice: the prior face places marks at exponents \(2p-1,3p/2-1,p/2-1\) for every even \(p\ge4\), giving exactly \(E_p+1\). Thus the number of growing runs in that scalar construction is minimal. [Three-contact solution](../chanpak_perfect_wieferich_seams_20261122/PROOF.md)
