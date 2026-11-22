# From the 3–4–5 triangle to a physical receiver

**PROVED.** The two boundary proportions that produce the exponent-15 flip are \(3/5\) and \(376/627\). The determinant equation selecting the second has a complete arithmetic parametrization. One branch starts with the Euclidean parameters of the 3–4–5 triangle. Infinitely many of its subsequent points admit an explicit tagged tournament receiver. A congruence explains exactly why changing this receiver alone cannot cross exponent 61 with the original left tags.

## The branch in the ternary tree

For the fixed left proportion \(3/5\), the solutions of \(5r-3h=-1\) with \(h,r>0\) are
\[
                    (h,r)=(5k+2,3k+1),\qquad k\ge0.
\]
They are coprime and have opposite parity. Consequently \((h^2-r^2,2hr,h^2+r^2)\) is a primitive Pythagorean triple for every \(k\). At \(k=0\) it is \((3,4,5)\). At \(k=125\) it is \((251753,471504,534505)\); the two Euclidean parameters at this point are the actual path counts \(627,376\) in the [exponent-15 certificate](E15-FAREY.md).

In Euclidean parameter coordinates two of the three Berggren children are \(L(u,v)=(2u-v,u)\) and \(R(u,v)=(u+2v,v)\). Their triple coordinates follow by substituting these parameters in Euclid's formula. For every \(k\ge1\), apply \(L\) exactly \(k-1\) times to \((2,1)\), then \(R\), then \(L\). The three successive states are \((k+1,k)\), \((3k+1,k)\), \((5k+2,3k+1)\). This gives the precise ternary-tree address of every member of the arithmetic branch. Only odd \(k\) can have \(h\) equal to a tournament height, and \(k=1\) gives the excluded height 7. [Earlier Berggren and tournament growth operations](../../../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md)

## Every numerator at a power-of-two-plus-one denominator

Let \(V_n=C_3[\mathrm{TT}_n,1,1]\), ordered cyclically as spine \(\to a\to b\to\) spine. Exactly one full path starts at \(a\), and one starts at \(b\). A path starting at spine rank \(i=0,\ldots,n-1\) must use a nonempty increasing initial run, then \(a,b\), then the remaining spine vertices in order. Each larger spine vertex can independently belong to either run. There are \(2^{n-1-i}\) such paths. Hence
\[
 H(V_n)=2^n+1,\qquad \{r_{V_n}\text{ over all tag sets}\}=\{0,1,\ldots,2^n+1\}.
\]
Here \(r\) counts full paths starting at a nontag. The binary digits of any \(r<2^n\) specify exactly which ranks to leave untagged when both singletons are tagged. The last two values use the singleton starts as well.

For any physical left flank with reduced critical proportion \(u/d=(h_A-e_A)/h_A\) in lowest terms **and with \(d=h_A\)**, \(0<u<d\), the determinant selects \(V_n\) exactly when \(u(2^n+1)\equiv\pm1\pmod d\). It then specifies an actual numerator \(r=(u(2^n+1)\mp1)/d\) and the rule above chooses its tags. If the unreduced counts have a common divisor larger than one, the unit determinant is impossible. For odd \(d>1\) all admissible \(n\) follow by one finite search modulo the multiplicative order of 2 at \(d\). The chosen sign determines the direction of the unit flip. This statement uses the middle edit with half-differences \((-1,1,1)\) of the exponent-15 proof.

For \(u/d=3/5\) and the negative determinant, the precise condition is \(4\mid n\). At positive \(n\) divisible by four, \(k=(2^n-1)/5\) gives the required arithmetic point. This realizes a physical receiver and a permanent unit edge for infinitely many points on the Pythagorean branch. The denominators include the Fermat primes \(17,257,65537\) at \(n=4,8,16\). The same receiver construction works at composite denominators. The modulus controls attainment of the determinant; the absolute target still depends on the left tags and exterior addresses.

## Changing the left tags changes which targets can be reached

Keep the actual four-vertex left tournament \(A\) and six-vertex middle edit \(D\) of the certificate, with middle tags 35. Let \(V\) now be **any** nonempty tagged tournament, and put \(h=H(V)\). Write \(j=x\bmod2\) for the exterior left-address parity. The seven-coordinate join gives the following entire modulo-five projection. The columns refer to the middle states in original and flipped orientation, for any exterior right address.

| Left tags | Original count modulo 5 | Flipped count modulo 5 |
|---:|---|---|
| 1 | \(h(2+j)\) | \(h(1+j)\) |
| 12 | \(h(1+2x+3j)\) | \(h(2x+3j)\) |

Indeed, in the coordinate order \((H,Q,c,t,e,r,k)\) the profiles of \(A\to D\) at tags 1 reduce to \((0,0,0,0,0,1,2)\) and \((0,0,0,0,0,1,1)\). At tags 12 the corresponding tuples are \((0,0,0,2,0,3,1)\) and \((0,0,0,2,0,3,0)\). Joining \(V\) multiplies each surviving entry by \(h\), and adding the exterior flanks gives the table. [Exact join law](E15-FAREY.md)

At tags 1, if \(5\mid h\) both heights vanish modulo five. Otherwise neither height vanishes. This excludes \(E_p-1,E_p+1\) for **every \(p\equiv1\bmod4\)** across the full class of replacement receivers at every order. The left retagging to 12 keeps \(h_A=5,e_A=2\), hence exactly the same Farey selector. For the branch \(5r-3h=-1\) it also gives \(h\equiv2\bmod5\), with the flipped middle the lower completion. The target residues \((0,2)\) are attained precisely when \(x\equiv0,1\bmod10\). In particular an empty left flank is permitted. This isolates a specific physical change that removes the all-receiver modulus-five obstruction while preserving the gap.
