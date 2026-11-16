# Collatz rows, signed descent, and the induced tournament barriers

**Results.** The row picture has exact equations. Its simplest signed descent covers all odd target addresses arithmetically, but fails as a rule for induced tournament growth, even on legal prime counts. Within the one-apex class there is a sharp growth floor, a classification of its lowest carries, and an exact list of all possible binary fork parents.

## 1. Two-column tower and the mod-12 wheel

Number the rows starting with 1. Row \(2j+1\) has \(3j+1,3j+2\) on its left and right; row \(2j+2\) has the single point \(3j+3\). Write \(\rho(n)\) for the row of \(n\). The disjoint mod-12 classes are \(\{0,3,6,9\}\) (multiples of 3), \(\{2,4,8,10\}\) (remaining evens), and \(\{1,5,7,11\}\) (the prime-candidate odds). The last class is not a primality criterion.

In the shortcut Collatz map, all **forward odd** edges land on the right, since \((3n+1)/2\equiv2\pmod3\). The destination row is exactly the old integer \(n\), and the three cases are

| Old odd integer | Start row | Destination row | Row difference |
|---|---:|---:|---:|
| \(6j+1\), left | \(4j+1\) | \(6j+1\) | \(2j\) |
| \(6j+3\), short | \(4j+2\) | \(6j+3\) | \(2j+1\) |
| \(6j+5\), right | \(4j+3\) | \(6j+5\) | \(2j+2\) |

This gives exactly the three families starting every fourth row. In contrast, the ascending **reverse even** arrow is doubling. A short row \(r\) doubles to short row \(2r\). The independent chains begin at \(r\equiv2\pmod4\); the first starts at row 2 and has successive gaps \(2,4,8,\ldots\). On a long row \(2j+1\), its left point doubles to the right of row \(4j+1\), and its right point doubles to the left of row \(4j+3\). Inclusive lengths (end minus start plus one) are \(2j+1\) and \(2j+3\). The full length list has repeats \(1,3,3,5,5,7,7,9,9,\ldots\).

There is an exact shoelace selector for the examples \(4,5,7,11\): let \(n_i=3+2^i\). Then \(n_{i+1}=2n_i-3\), and \(\rho(x+3)=\rho(x)+2\) makes consecutive strands \([\rho(n_i),\rho(2n_i)]\) overlap in three inclusive rows. Their lengths are
\[
3+2\left(\left\lfloor2^{i+1}/3\right\rfloor-\left\lfloor2^i/3\right\rfloor\right)
=3,5,5,9,13,25,\ldots .
\]
There is a second length-9 strand \(13\to26\) in the full picture, but not as the next term of this particular selector.

For an odd modulus \(m\), the signed odd arrow \(n\mapsto(3n+\sigma)/2\) translates by \(z=n+\sigma\) to multiplication by \(3/2\) modulo \(m\). It is a permutation iff \(3\nmid m\), and a genuine permutation of the nonzero classes for primes \(m\ne2,3\). The even arrow under the same coordinate has \(z\mapsto(z+\sigma)/2\), not multiplication by the same unit. An odd residue class doesn't tell which arrow is legal at an integer lift. That is the exact limit of using the multiplication table as the dynamics.

The earlier [weighted interpretation](../collatz_order_geometry_20261116/README.md) keeps a different piece of information. For the pruned shortcut inverse tree in the first-hit basin of 1, live fork phases \(2,20;11;8;17,26\bmod27\) have length pairs \((2,2),(2,4),(3,4),(1,4)\) respectively. The unweighted skeleton forgets all four distinctions; the lengths recover the suppressed stretches and the sterile rays by that construction. This concerns the actual basin, not an assumption that every positive integer reaches it.

## 2. Arithmetic descent that does not depend on the Collatz conjecture

There are two distinct possible parent maps.

* **Signed shortcut inverse.** For odd \(p>3\) with \(3\nmid p\), the unique integral \(q=(2p-\sigma)/3\), \(\sigma\in\{-1,1\}\), is odd and smaller. In \(t=(H-1)/2\) coordinates the forward child is \(t\mapsto\lceil3t/2\rceil\). Indeed \(t=2j,2j+1\) gives \(3j,3j+2\); the ordinal join with a triangle gives \(3j+1\) via \(t\mapsto3t+1\). Thus these are exact complementary output lanes. If the signed moves could be lifted for every admissible parent/target (and a witness could be chosen each time), induction would suffice. Parent holes require only the already-realized restart counts 11, 31 and 63.
* **Greedy signed odd part.** For odd \(p>1\), choose \(\sigma=1\) at \(p\equiv1\pmod4\) and \(-1\) otherwise, then \(q=\operatorname{oddpart}(3p+\sigma)<p\). No such parent is divisible by 3. The only forbidden parent it can meet is 7, beginning at 19 among primes. A neighboring dyadic parent can bypass that specific hole (for example \(19=2\cdot9+1\)). This is arithmetic only, not a claim of universal physical lifting.

The barriers are not just failure to find a suitable orientation. Exhaust the fixed-Hamilton-path frames at orders up to seven and test every induced subset using the same subset DP. For the **first** rule, no \(H=13,17,19,23,29\) tournament contains respectively \(H=9,11,13,15,19\). For the **second** rule, no 17 contains 13, no 23 contains 17, and no 31 contains 23. These are all-order statements: a prime target has only one nontrivial strong component, the [Busch bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) makes it order at most seven for targets through 31, and ordinal singleton padding changes no induced count outside the core. Every tournament has a Hamiltonian path so isomorphisms into a fixed frame lose no cases. Counts in the receipt are frames, not isomorphism classes.

There is also a concrete mod-12 warning for a different proposal: any prime \(p>3\) has a unique \(p=6h+c\) with odd \(h\) and \(c=-5,-1,1,5\) in residue order \(1,5,7,11\bmod12\). Doubling a triangle-scaled source and adding only \(\pm1,\pm3\) misses both \(\pm1\bmod12\) prime lanes. For that design, the \(\pm5\) carries are not optional.

## 3. A general growth floor and the smallest physical carries

Use the [one-apex seam formula](../one_apex_all_residues_20261116/followup/README.md): at \(U|V\), set \(X=X(U)\) (binary value), \(Z=Z(V)\) (reversed complemented binary value), and \(K=H(U)+H(V)-1+[10]\). If \(H\) is the parent and \(H'\) is its zero-insertion child, then
\[
H=K+XZ,\qquad D=2H'-3H=X(Z+2)-K+2[11].
\]
One insertion is enough to analyze growth. A strong old spine starts 1 and ends 0. At an interior gap \(X,Z\ge1\) and \(K\le X+Z+2\), since \(f(n)=H(\operatorname{bin}(n))\le n+1\). For \(X\ge3\) this gives \(D\ge(X-1)Z+X-2\ge3\). For \(X=2\), the prefix is 10, so \(D=2Z+2-f(Z)\ge Z+1\ge2\). For \(X=1\), let \(q\) be the first bit of \(V\); then \(D=Z+1-f(Z)+3q\). If \(q=1\) this is at least 3. If \(q=0\), then \(Z\) is odd and \(f(Z)=f((Z-1)/2)\le(Z+1)/2\); \(D=1\) is possible only at \(Z=1\), namely \(H=3\to5\). The only active exterior zero insertion is at the right, giving \(D=2X-H\ge X-1\); except for 3 this is again at least 3 by odd parity. Reverse-complement covers insertion of 1.

Consequently **any count-increasing induced extension within the one-apex class from a parent with \(H\ge5\) has \(2H'-3H\ge3\)**, even if arbitrarily many vertices are added. To obtain a chain of strong extensions, strip neutral outer padding and choose the final eligible apex: it must belong to the old nontransitive core and restrict to an eligible apex there. In its final strong spine insert any new extreme 1 on the left or 0 on the right first when needed, then the remaining interior bits; this grows by strong words at each step. Any neutral ordinal singleton padding can be ignored. If the output is prime, \(D\ge5\) for \(H'\equiv1\bmod3\), and \(D\ge7\) for \(H'\equiv2\bmod3\). This is an obstruction at **absolute counts**, not residues.

One can classify every \(D\le7\) move once \(H\ge23\). Two active steps cannot fit, since they raise \(H\) by at least a factor \((3/2)^2\). For one step, the exterior bound forces \(H\le9\); \(X=2\) forces \(Z\le6\) and \(H\le21\); the \(X\ge3\) inequality gives \(H\le17\) on its allowed pairs. For \(X=1,q=0\) it gives \(Z\le13\) and \(H\le21\). Thus only \(X=1,q=1\) remains. Write \(g(n)=n+1-f(n)\). The exact append recurrence
\[
f(2n+1)=f(n),\qquad f(2n)=f(n)+n+(n\bmod2),\qquad f(0)=1
\]
gives \(g(2n)=g(n)-(n\bmod2)\) and \(g(2n+1)=g(n)+n+1\).
For even \(Z=(2v+1)2^r\), \(r\ge1\), it follows that \(g(Z)=v+g(v)\). The bound \(D=3+g(Z)\le7\) leaves \(v=0,1,2,4\). The exact rays are therefore
\[
(H,H')=(2Z+1-g,\ 3Z+3-g),\qquad
(Z,g)=(2^r,0),(3\cdot2^r,2),(5\cdot2^r,2),(9\cdot2^r,4).
\]
Conversely take \(V\) to be the reverse-complement of \(\operatorname{bin}(Z)\); insert 0 between the first 1 and \(V\). The first ray only outputs multiples of 3; the other three give the prime candidates \(9\cdot2^r+1,\ 15\cdot2^r+1,\ 27\cdot2^r-1\), with no infinitude claim about primes. For example \(23\to37\) and \(33\to53\) are exact. A fixed replacement of the signed carry does not reach arbitrary primes in this class.

More quantitatively, the number of possible output counts \(H'\le B\) with \(0<D\le C\) is \(O(C\log B)\) in this class. For \(X\ge2\) the bounds give \(XZ=O(C)\); there are \(O(C\log C)\) pairs, and a strong prefix and reversed-complemented suffix have no leading zeroes, so the pair fixes the words. Exterior cases have \(X=O(C)\). For \(X=1\), \(g(Z)=O(C)\), giving \(O(C)\) odd cores (because \(g(2v+1)\ge v+1\)) with \(O(\log B)\) trailing-zero lengths. The latter bound also follows directly from the strong floor \(H'\ge1+2^{n-2}\) on a spine of length \(n\ge3\). If two or more active insertions occur, \(H'\ge(3/2)^2 H\) implies \(H'\le 3C/2\). Use \(C\le2B\) without loss. The count bound is on the small-deviation method, not on the full spectrum.

## 4. Full classification of one-apex binary fork parents

Surprisingly, allowing different gaps does not enlarge the fork set. For \(H>5\), two zero runs cannot both admit the two signs: on crossing a 1 the zero deviation grows by at least \(X_{\rm before}+Z_{\rm after}\) as proved in the parent note. A jump of at most 2 at an interior crossing would force both to be 1, i.e. the word 110 with \(H=5\). A leading pad isn't a legal child. The one-bit case is the reverse-complement.

For opposite inserted bits at two separate cuts \(U|W|V\), write \(x=X(U),z=Z(V),A=X(W),B=Z(W),P=2^{|W|}\) and \(E=H(U)+H(W)+H(V)-2+[U|W=10]+[W|V=10]\). If 0 is at the left cut and 1 at the right, the sum of deviations from \(2H\) is \(x(1-B)+z(1-A)+\epsilon_0+\epsilon_1-2E\). A pure middle word lets one cut slide to the other. A mixed one has \(A,B\ge1\); for the sum to be zero requires \(E=1\) and both corrections 1, hence \(W\) starts 1 and ends 0, contradicting \(H(W)\ge3\).

In the other order slide the two cuts inward as far as possible, so if they remain distinct then \(W\) starts 0 and ends 1. The sum is now \(A+B+(P-B)x+(P-A)z+\epsilon_1+\epsilon_0-2E\). Here \(x,z\ge1,\ A,B\le P/2-1,\ E\le x+z+H(W)\) by \(f(j)+(j\bmod2)\le j+1\), and \(H(W)\le\min(A,B)+1\). The sum is at least \(2(P-\min(A,B)-3)\ge P-4>0\) for \(|W|\ge3\); for \(W=01\) the sharper bound is \(x+z>0\). This is incompatible with opposite signs. An unpointed two-apex parent has reverse-complement pointed words by the archived classification, so their sets of possible output counts coincide. At the tiny \(H=3,5\) cases direct checking suffices.

It remains to solve the same-gap arithmetic. The two deviations are \(X+[11]-K,Z+[00]-K\), as in the parent note.
* At boundary 01 both \(X,Z\) are even, so they would be \(2n,2n+2\) and \(K=2n+1\). The append recurrence instead gives \(K=2n+1+f(n)+f(n+1)>2n+1\).
* At 10 they are \(2n+1,2n+3\); \(K=2n+2\) requires \(f(n)+f(n+1)=2n+2\). Only \(n=0,1\) work: for even \(n\ge2\) the sum is at most \(n+1+n/2+1\), and for odd \(n\ge3\) it is at most \(2n+1\). The counts are 5 and 19.
* At 00 the nonsquare case has \(X=2n,Z=2n-3\), requiring \(f(n)+f(n-2)=n-(n\bmod2)\). Even \(n=2v\) is impossible by the recurrence (the sum exceeds \(2v\)). Odd \(n=2v+1\) reduces to \(f(v)+f(v-1)=2v\) and gives \(v=1,2\) from the 10 case. The counts are 23 and 79.
* At 00 the square case has \(X=2n,Z=2n+1\), requiring \(2f(n)=n+2-(n\bmod2)\). Even \(n=2v>0\) requires \(f(v)+(v\bmod2)=1\), impossible. Odd \(n=2v+1\) needs \(f(v)=v+1\), attained only at \(v=0\) and powers of 2 at least 2, since \(g(2v)=g(v)-(v\bmod2)\) and \(g(2v+1)=g(v)+v+1\). Boundary 11 is its reverse-complement counterpart.

Thus the **exact** fork-parent counts are
\[
\{5,9,19,23,79\}\ \cup\ \{(2^{j+2}+3)^2:j\ge1\}.
\]
All have literal witnesses in the verifier or the parent construction. In particular **79 is the last prime parent** admitting both induced binary children within one apex, at any gaps. The earlier next-step nonrenewal theorem now applies regardless of whether the initial gaps were chosen separately. This rules out a complete induced one-apex branching proof, not all possible ways to realize the primes.

## 5. Why a single local edit is also too coarse

Keep the apex contacts but flip the spine arc between two adjacent
vertices with opposite contacts, then read the new transitive spine in
order. This changes \(U01V\) to \(U10V\). If both graphs are strong,
the gap is interior. Direct use of the concatenation formula gives
\[
H(U10V)-H(U01V)=X(U)+Z(V)+1-p+q,
\]
where \(p\) is the last bit of \(U\) and \(q\) the first bit of \(V\).
The \(X,Z\) terms come from the pair's ordinary and reflected binary
values, both changing from 1 to 2. Internal \(H\) rises by 2; the
boundary adjacencies subtract \(p+1-q\).
For \(i=|U|,j=|V|\ge1\), the exact minimum is
\(2+[i\ge2]2^{i-1}+[j\ge2]2^{j-1}\), by the leading 1 in each
binary address. Thus at large spine length \(n=i+j+2\) the smallest
such edit at any position is of order \(2^{n/2}\); since strong
one-apex counts are of order \(2^n\), these edits cannot supply
small absolute corrections by themselves. This is a different
restriction from the induced-growth floor. A next constructive
test is to look for a second physical seam whose large response
cancels the first, or for a larger feedback core that permits a
genuinely independent carry. The [mirror follow-up](MIRROR.md) finds two such edits with exact net +2, including the twin-prime counts 10007 and 10009.

## Replay

Use python3 verify.py for the row, growth and fork controls, and compile descend.cpp with a C++17 compiler for the fixed-frame checks. Run the latter with no argument for greedy odd-part ancestry, or with any argument for signed-shortcut inverse ancestry. The finite output is retained in replay.txt. The row/growth/fork arguments and the bounded-order reduction explain exactly which all-order statements these checks support.
