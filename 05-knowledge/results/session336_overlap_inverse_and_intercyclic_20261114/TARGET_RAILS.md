# Two unbounded direct inverses grown from the target witnesses

The exact-overlap selector gives strong words of counts 385 and 573 which are not in the intercyclic class. Two homogeneous inflations turn each into an explicit target family. For $a,b\ge1$ put $A=2^a$, $T=2^b$, $R=3^b$.

| spine word | ordinary count | one locked frame |
|---|---|---|
| $1\,0^a30\,2^{b+1}1$ | $1+18T+12R+A(-6+18T+42R)$ | $(\alpha,\beta)=(a,b+1)$ |
| $2\,0^a30\,2^{b+1}12$ | $5+36T+24R+A(8+48T+36R)$ | $(\alpha,\beta)=(a,b+2)$ |

At $a=b=1$ these are the target-385 and target-573 witnesses respectively. Increasing either exponent simply inflates a vertex by a transitive block, preserving strongness. Each count lies in the tensor product of the bases $1,2^r,3^r$ for the two inflated lengths: fixed occupied-run sets give the three surjection polynomials. The exact subset DP on the two $3\times3$ grids verifies the displayed expressions, which suffices for all lengths.

For the first word the ordered four anchors are first 1, the 3, last interior 2, final 1. In the second they are first 2, the 3, penultimate 1, final 2. The first two allow a distinct outgoing contact for each hub, and the last two allow a distinct incoming contact. Between the first and second anchors are the $a$ doubly eligible zeros. Between the two inner anchors are respectively $b+1$ and $b+2$ triply eligible sites. The final anchors are adjacent. This is a locked frame guaranteed at the exact *ordinary* count, not a marked-fibre count substituted for it.

## Starting with the target

Given $N$, enumerate $b\ge1$ only up to $3^b<N$ (all counts above exceed $3^b$). Subtract the part independent of $A$ and divide by its positive coefficient. Test whether the quotient is an integer power of two at least 2. Each success returns a word and its overlap frame, with no search over spine words. If the application prescribes a minimum $\beta$, use the displayed $\beta$ to restrict $b$ before testing. This takes $O(\log N)$ big-integer tests per family and is complete **for these two families only**.

The first family's increments genuinely leave the earlier offset-one smooth rail. Indeed
$$
 (N-1)/6=3T+2R+A(-1+3T+7R).
$$
This has nonzero residue $-A$ modulo 3, so $v_3(N-1)=1$ throughout. If $b\ge2$, the first term is a multiple of 4, the second is 2 modulo 4, and the third is divisible by 4; thus $v_2(N-1)=2$. If $b=1$, the expression gives $N-1=12(6+13A)$, with $v_2=3$ for $a\ge2$. Consequently only $a=b=1$ gives a $2,3$-smooth increment. This is a verified mixed-count inverse, not a reiteration of the single smooth monomial rail.

Run **python replay_target_rails.py** for the two coefficient certificates and explicit target inversions. The broader search is still needed for targets outside the rails; no claim of complete Chan–Pak coverage follows.
