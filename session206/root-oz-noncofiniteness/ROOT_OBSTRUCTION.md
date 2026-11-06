# Root O/Z scalar noncofiniteness and the d=1 bank collapse

Status: independently audited proof. This does not exclude every finite physical bank.

Write Q(0)=0, Q(2m)=Q(m), Q(2m+1)=Q(m)+m+(m mod 2). The rooted-edge O/Z tail scalar is H(m)=1+2m-Q(m), for m>=0. Parameter m records tail letters from low to high bits, with O=0 and Z=1. Nonempty canonical tails end in Z; terminal O letters do not change the scalar.

## 1. An exact converse involution

For m>=1 let n be its bit length, m=2^(n-1)+a with 0<=a<2^(n-1), and let rev_k be reversal of exactly k binary digits (including leading zeros). Define

    T(m)=2^n-1-rev_(n-1)(a).

Then T is an involution of each dyadic block, and H(T(m))=H(m).

Proof: the physical rooted-edge construction consists of distinguished vertex c and a transitive chain v_0,...,v_n. Vertex c beats v_0. For i>=1, v_i beats c exactly when tail digit b_(i-1) is one. Thus the c-adjacency word in chain order is

    (0,b_0,...,b_(n-2),1).

Reverse every arc of this tournament, and list the transitive chain in reverse order. Its c-adjacency word is the reverse complement

    (0,1-b_(n-2),...,1-b_0,1).

This is exactly another rooted-edge O/Z construction, with parameter T(m). Reversing the order of a directed Hamiltonian path gives a bijection with Hamiltonian paths in the converse tournament. Hence the two scalars agree. Applying the operation twice proves the involution. In particular the arithmetic identity is

    Q(T(m))-Q(m)=2(T(m)-m).

The fixed-point equation says that the n-1 lower bits equal their reverse complement. For n-1 odd there are no fixed points; for n-1 even there are exactly 2^((n-1)/2).

## 2. Rigorous image-density obstruction

Let A(X)=#{m>=0:H(m)<=X} and I(X)=#{H(m):m>=0,H(m)<=X}. The parameter-capacity theorem frozen in `dependencies/PARAMETER_CAPACITY.md` gives

    A(X)=(log 2)X+O(sqrt X).

Because Q(m)<=m-1 for m>=1, H(m)>=m+2 (and the weaker H(m)>=m+1 also suffices). Consequently all parameters counted by A(X) have O(log X) digits. Summing the fixed-point counts in Section 1 gives O(sqrt X) such parameters. The finite set counted by A(X) is T-invariant, except m=0 which is treated as a fixed singleton. Each T-orbit maps to a single scalar. Therefore

    I(X) <= A(X)/2 + O(sqrt X)
         <= (log 2)X/2 + O(sqrt X).

All these scalars are odd. In particular the number of missing odd integers at most X is at least

    ((1-log 2)/2)X - O(sqrt X).

Thus the root O/Z image is not cofinite among odd integers. This is an all-order obstruction, not a conclusion from numerical density. It does not assert existence or the value of the image's asymptotic density, nor unbounded lengths of consecutive gaps. The numerical upper coefficient is approximately 0.346574 relative to all integers, or log 2 relative to odd integers.

## 3. Every physical seed with d=1 collapses into this same root image

Throughout this packet a physical seed is a nontrivial tournament V (at least two vertices) with distinguished c having one Hamiltonian path starting at c, and d=H(V-c), with the established R>d>0. The singleton is not a source in this statement. If d=1, then V-c is transitive. Indeed a nontransitive tournament has a directed triangle. Each of the triangle's three directed Hamiltonian paths extends to a Hamiltonian path of the full tournament by inserting the remaining vertices one at a time. The insertion always exists: insert the new vertex before the first path vertex it beats, or at the end if it beats none. The three final paths remain distinct because their relative triangle orders differ. Therefore a tournament with exactly one Hamiltonian path has no directed triangle and is transitive.

Write V-c as the transitive chain v_0,...,v_k. The unique possible Hamiltonian path starting at c is c,v_0,...,v_k, so c beats v_0. Every later v_i either beats c or loses to c, while losing to all earlier chain vertices. Hence V itself has a rooted-edge O/Z construction. Appending any further O/Z tail is still a rooted-edge O/Z construction. Therefore the entire scalar tail image of every physical d=1 seed is a subset of the root image.

Consequently every finite bank whose physical seeds all have d=1 has the same noncofiniteness obstruction, regardless of how many seeds it contains. The result even applies to an arbitrary infinite collection of such seeds.

## 4. A sufficient obstruction for banks that may contain the root

For a finite physical bank B, let B_> be the subset with d_s>1. If at least one member has d=1, replace the union of all such members by the root image as an upper bound. Combining Sections 2-3 with the frozen parameter-capacity bound gives

    #{represented scalar values <=X}
      <= X[(log 2)/2 + sum_(s in B_>) log(R_s/(R_s-d_s))/d_s]
         + O_B(sqrt X).

If the bracket is strictly less than 1/2, the bank misses a positive proportion of odd integers. This improves the elementary reciprocal-width obstruction and can exclude some root-containing banks. It remains a sufficient test only: a bracket at least 1/2 neither proves nor disproves coverage.

For banks with no d=1 member the unmodified parameter-capacity criterion applies. One can also omit any bank member whose tail image is already known to lie in another member's image; doing so only improves the sufficient test.

### A nontrivial root-containing example

From the rooted edge (h,R,e,d)=(1,2,0,1), two J operations followed by O give

    J:(3,3,1,1), J:(5,6,2,3), O:(5,12,0,3).

Here J sends (h,R,e,d) to (h-e+R,h+R,h-e,h). Thus the two-member physical bank consisting of the rooted edge and this JJO seed has image-count upper coefficient

    (log 2)/2 + (log(4/3))/3 < 1/2.

The strict inequality needs no numerical logarithm approximation: convexity gives log 2<3/4 by the trapezoidal bound for the integral of 1/x on [1,2], while log(4/3)<1/3. Thus the coefficient is less than 3/8+1/9=35/72<1/2.

The added seed genuinely enlarges the root image: its tail parameter m=2 gives 29, while the root image does not contain 29. For the latter exact finite membership check the inverse bound restricts root parameters to m<=27; direct use of the Q recurrence excludes all of them. This finite check only establishes the stated distinction of the two images; the noncofiniteness conclusion is the all-order inequality above.
