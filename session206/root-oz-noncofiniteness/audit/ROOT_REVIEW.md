# Independent root-image obstruction review

Verdict: PASS for the four theorems of `session206-tail-obstruction/ROOT_OBSTRUCTION.md`, under its inherited nondegenerate physical-source convention R>d>0 and nonempty V-c.

The converse tournament reverses the transitive-chain order and complements the distinguished-vertex adjacency word. This gives exactly the stated bit involution, and the Hamiltonian-path reversal bijection preserves the scalar. Its fixed points are anti-palindromic lower bits, with the stated counts. The admissible-input set is involution-invariant because the scalar is preserved. Taking orbit counts gives I(X) <= (log 2)X/2+O(sqrt X), and a positive lower density of missing odds.

Independent arithmetic check: write m=2^k+x, 0<=x<2^k, and y=2^k-1-rev_k(x). For k>=1, Q is invariant under reversal of exactly k bits because every pair weight depends only on distance. Complementing k bits gives

Q(complement(z))-Q(z)=2^k-2-z-rev_k(z)+lowbit(z)+topbit(z).

The identity Q(2^k+x)=Q(x)+rev_k(x)+topbit(x) then yields Q(2^k+y)-Q(2^k+x)=2(y-x). The cases m=0,1 are fixed separately.

The d=1 collapse is valid: a directed triangle supplies three distinct Hamiltonian paths that remain distinct after successive vertex insertion. Thus H(V-c)=1 forces V-c transitive. Unit start then forces c to beat the chain's first vertex, so the entire seed is a rooted-edge O/Z construction. All its further tails stay in that root image. The mixed-bank union bound consequently follows. The singleton boundary must be excluded or treated separately rather than silently using a nonexistent first chain vertex.

No exact root image density, long consecutive gaps, or general finite-bank impossibility follows.

Revision review: the author has now explicitly excluded singleton sources, resolving the scope note above. The added physical bank example is correct: JJO sends (1,2,0,1) to (5,12,0,3); m=2 gives 29. The root inverse bound H(m)>=m+2 confines any preimage of 29 to m<=27, and exact recurrence evaluation excludes these. The strict coefficient inequality can be certified without decimal arithmetic: log(2)<3/4 by the strict trapezoidal bound for the convex integrand 1/x on [1,2], and log(4/3)<1/3; hence its coefficient is less than 3/8+1/9=35/72<1/2.
