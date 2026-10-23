# The canonical diagonal grid, exact counts, and a typed phase quotient

Status: PROVED under the canonical three-unit convention in the Fibonacci packet. The current prefix does not settle the earlier conflicting color at35. No claim here resolves that historical choice.

## 1. Three unit colors form two macrotiles

The canonical substitution is sigma(C)=CA, sigma(A)=B, sigma(B)=CA. Every C is immediately followed by A and every A is preceded by C. Group CA into a length-two tile U and B into a length-one tile V. Then

sigma(U)=UV, sigma(V)=U.

Thus the three-color word is exactly the ordinary Fibonacci tile word with its long tile split into two differently owned unit positions. There are not three interchangeable color states: the two long-tile halves each have density phi^-2, while the short tile has density phi^-3.

This is compatible with the user's hidden-unit interpretation under the explicitly selected proper expansions, but it is not proof that any arbitrary nonconsecutive decomposition yields that word.

## 2. Exact positions, not a fitted prefix

Put alpha=phi^-2, beta=phi²=3−alpha. The kth occurrences of the three colors are

C_k=floor(k beta)−1,
A_k=floor(k beta),
B_k=2 floor(k phi)+k.

Here is a self-contained Beatty proof. The two sets floor(k phi) and floor(k beta) partition the positive integers: their cumulative counts through n are floor((n+1)/phi) and floor((n+1)/beta), whose sum is n because the two irrational arguments sum to n+1. The increment at each n is exactly one, proving disjointness and coverage. Use these as U and V positions in a binary word. Under U→UV,V→U, the kth newly created V is at position floor(k phi)+k=floor(k beta), so the word is fixed and begins U. It is therefore the Fibonacci fixed word.

Expanding its long tiles gives the kth C position floor(k phi)+k−1, followed by A. The kth V tile is at tile-position j=floor(k beta); among the earlier tiles exactly j−k are long. Its expanded position is j+(j−k)=2floor(k phi)+k. This proves every displayed formula.

The exact prefix counts at length H are consequently

c(H)=floor((H+2)alpha),
a(H)=floor((H+1)alpha),
b(H)=H−c(H)−a(H).

They imply bounded discrepancy directly, rather than only limiting frequencies. Combined with the Fibonacci fibre theorem, the total numbers of representations for1≤N≤H are

sum |D(N)|=3H−floor((H+1)alpha),
sum |R(N)|=4H−1−floor((H+1)alpha)−floor(H alpha).

These exact expressions sharpen the already proved asymptotics phi²H+O(1) and2phi H+O(1).

## 3. The triangular grid and a logarithmic exact algorithm

Number a row's unit cells from the left, j=1,...,N. Row coherence gives the color at(N,j) as the(N−j+1)st letter of the word. Every diagonal N−j=constant therefore has one fixed color. This is an exact diagonal slicing of the triangular grid, not a metaphor for its appearance.

Define S(H)=sum_(m=1)^H floor(m alpha), with S(0)=0. The color populations in all rows1,...,H are

T_C(H)=S(H+2),
T_A(H)=S(H+1),
T_B(H)=H(H+1)/2−T_C(H)−T_A(H).

In particular T_C(H)−T_A(H)=floor((H+2)alpha)=c(H): equal asymptotic color densities leave a precisely measured boundary excess.

For K=floor(H alpha), a lattice-point transpose gives

S(H)=K(H+1)−3K(K+1)/2+S(K).

Indeed S(H)=sum_(j=1)^K(H+1−ceil(j/alpha)), and ceil(j beta)=3j−floor(j alpha). The recurrence reduces the argument by the factor alpha<1 and therefore takes O(log H) iterations. All arithmetic can be exact: for m>0,

floor(m alpha)=(3m−isqrt(5m²)−1)//2,

where isqrt is the integer square root. The subtract-one is essential because m sqrt5 is irrational. Treat m=0 separately. No floating-point threshold decisions are used.

## 4. Zero-slot elimination is a relation, not reconnection

The boundary factor graph before deleting zero is A—Z—B—C, with the zero slot Z forced unselected. Correct existential elimination of Z leaves all four binary assignments(A,B) allowed: both constraints adjacent to Z are already satisfied. Reconnecting A directly to B instead imposes a NEW hard-core constraint and removes(A,B)=(1,1).

Thus ignoring a zero length is harmless to lengths but need not be harmless to admissibility. The removed A+B representation sheet is the exact missing joint constraint coordinate. This is a precise relational-join analogy to shared-port attachment choices, not a bijection between Fibonacci decompositions and tournament ports.

## 5. A genuine order-three quotient, and why it is not the original colors

In coordinates(C,A,B), the substitution incidence matrix is

M=[[1,0,1],[1,0,1],[0,1,0]].

It kills e_C−e_B. Its image consists of vectors(x,x,y), and on that image it acts by

T(x,y)=(x+y,x),    T=[[1,1],[1,0]].

Equivalently q(c,a,b)=(c+b,a), j(x,y)=(x,x,y), M=j q and q j=T. Modulo2, T permutes the three nonzero vectors in a cycle. This is a real C3 action on a quotient charge space. It is NOT a cyclic permutation of the original letters: C and B have already been identified, the letter substitution is not invertible, and the three letter frequencies are not equal.

The same distinction guards against importing a Collatz or {7,21} conclusion. A common three-state count is not a map preserving arithmetic trajectories or forbidden Hamiltonian counts. The exact retained connection here is a matrix quotient and a shared-incidence elimination rule, not a proof about either open problem.

## Reproduction

Run python3 check_colored_grid.py and python3 -O check_colored_grid.py. The checker validates exact word positions and diagonal counts through20000, compares floor sums directly, checks the quotient matrix and zero-slot relation, and uses inputs through501 digits for the shrinking recurrence. Independent audit uses a separately generated word and a different floor-sum control.
