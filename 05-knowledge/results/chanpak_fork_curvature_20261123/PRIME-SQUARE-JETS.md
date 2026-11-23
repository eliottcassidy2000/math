# All nonzero-base lifts at prime squares for the even physical controller

Let $F(t)$ denote the integer slot polynomial of a nonempty even-order tournament. The operations are two ordinal singleton sources and the cyclic sandwich, respectively
$$PF(t)=t^2F(t),\qquad WF(t)=t^2F(t)+t(F(t+1)-F(t-1)).$$
They, and all steps described here, are physical filler extensions inside the same uniform factor-seven exterior.

## Exact finite state and square reachability

For odd prime $p$ and $d=(p-1)/2$, the square state consists of
$$a_i=F(i)\pmod {p^2},\qquad b_i=F'(i)\pmod p,\qquad 1\le i\le d.$$
Its boundary has $a_0=b_0=0$, $a_{d+1}=a_d-pb_d$ modulo $p^2$, and $b_{d+1}=-b_d$ modulo $p$. Indeed $d+1=p-d$ and even parity gives $F(p-d)=F(d-p)$; Taylor expansion proves the two identities. For any integer input, even parity and Taylor expansion recover its value modulo $p^2$ from this same state (multiples of $p$ evaluate to zero). Direct differentiation gives the finite update:
$$
\begin{aligned}
a_i(P)&=i^2a_i,& b_i(P)&=i^2b_i+2ia_i,\\
a_i(W)&=i^2a_i+i(a_{i+1}-a_{i-1}),\qquad&
b_i(W)&=i^2b_i+2ia_i+(a_{i+1}-a_{i-1})+i(b_{i+1}-b_{i-1}).
\end{aligned}
$$
The $a$-equations are modulo $p^2$ and the $b$-equations modulo $p$.

For $p=3,5,7$, the [finite queue replay](even-jets.cpp) starts at the actual even transitive filler $F=t^2$. It reaches 18, 15,000, and 40,235,958 states respectively, exactly $(p^d-1)p^{2d}$ in each case. Neither move can enter zero base from nonzero base since both base matrices are invertible at these primes. This proves every square lift over a nonzero base, not just existence of one lift over each base. The calculation remains valid at 3 and 7 where the odd-parity matrix is singular.

## A congruence-kernel certificate at larger primes

A fiber difference is $z=(\Delta a/p,\Delta b)\in\mathbb F_p^{2d}$ at fixed base $v=a\bmod p$. Denote the ordinary even and odd wrap matrices modulo $p$ by $W_e,W_o$, differing by the sign of the reflected boundary. On differences, $P$ acts as $B_P=\operatorname{diag}(P,P)$, and
$$B_W=\begin{pmatrix}W_e&-dE_{dd}\\0&W_o\end{pmatrix}.$$
For $p=5,11,13,17,19,29$ both are invertible (the receipt computes exact inverses); they and the invertible base matrices make $P,W$ permutations of the entire square state set. Inverses in the reasoning can therefore be replaced by finite positive powers of the physical operations.

The power $P^{p-1}$ acts on every fiber as translation by $Cv$, with top entry $C_i v=((i^{2(p-1)}-1)/p)v_i$ modulo $p$ and bottom entry $C_{d+i}v=-2i^{-1}v_i$. The first is the Fermat quotient of the physical scaling and the second is its derivative. If a word $g$ acts by base matrix $A_g$ and on fiber differences by $B_g$, then $g^{-1}P^{p-1}g$ translates the fiber of $v$ by $B_g^{-1}CA_gv$. The constant term of an affine fiber map cancels in this conjugacy. Translations of multiples are obtained by repetition, and translations of sums by composition.

The [exact word-rank replay](even-jet-kernel.py), with fixed literal words, computes full ranks $4,10,12,16,18,28$ at $p=5,11,13,17,19,29$, respectively, over the fiber of $e_1$. Thus every lift over $e_1$ can be reached from any one lift there. Wherever the base nonzero orbit is transitive, transport to $e_1$, use this whole translation space, and transport back; the transports are bijections on fibers. This proves the entire nonzero-base square orbit at 11 and 13 using the separately certified base orbit, and applies with the independently certified base orbit results at 17,19,29 as well. Without that base certificate the kernel assertion alone concerns the base orbit of $e_1$, not every nonzero vector.

## No observational redundancy in the derivative state

From the first $d$ successive ordinary outputs under $W^0,\ldots,W^{d-1}$ one recovers all $a_i$ modulo $p^2$: at each first appearance the new coordinate has unit diagonal $(i-1)!$, and no derivative reaches coordinate 1 yet. For two states with those same $a_i$, look next at $W^d,\ldots,W^{2d-1}$. A difference from $b_j$ first reaches ordinary coordinate 1 at time $2d-j$. Its unique shortest route goes from derivative index $j$ up to $d$ along nonzero subdiagonal entries of $W_o$, crosses at index $d$ with coefficient $-pd$, then goes left to ordinary index 1 along the nonzero superdiagonal entries of $W_e$. Its coefficient is $p$ times a unit modulo $p^2$. Hence determine $b_d,b_{d-1},\ldots,b_1$ triangularly. The genuine fork gap separates the same outputs since the difference of its two exits is $4Q_1$.

Consequently each distinct attainable jet requires a distinct exact predictor state, regardless of the encoding. The full nonzero-base orbit supplies the lower bound $(p^d-1)p^{2d}$; zero-base completion at 3 and 5 is given below. The square results do not extend by themselves to higher prime powers, where additional Taylor coefficients enter.

Run `g++ -O3 even-jets.cpp -o even-jets; ./even-jets 3; ./even-jets 5; ./even-jets 7` and `python3 even-jet-kernel.py`.

### Zero-base completion at 3 and 5

There is also full zero-base synthesis at these two primes. At zero base write the state as $(u,b)=(a/p,b)$. Ordinal multiplication by any even filler of base vector $v$ scales both $u_i$ and $b_i$ coordinate by coordinate by $v_i$: its derivative term multiplies the zero base and disappears. Every nonzero $v$ is available by the foregoing proofs and zero $v$ by the cyclic filler from the base-prime packet. The zero state can also be obtained by multiplying any zero-base filler by itself.

At 3, the five explicit even seeds in `S3` of `even-zero-jets.py` have $(u,b)$ equal to $(0,0),(0,1),(2,0),(2,1),(2,2)$. Ordinal scaling by 2 (the even tournament $W(T_2)$ has $H=5$) makes these all nine states, without needing any further wraps. At 5, the three even seeds in `S5` have $(u_1,u_2,b_1,b_2)$ equal to $(3,3,4,3),(3,3,3,0),(0,4,2,3)$. Start with their coordinatewise base-scalings and apply the same $P,W$ closure, using the fiber matrices with zero base. The cumulative exact closure sizes are 577, 601 and 625, so all $5^4$ zero-base states are attained. The companion independently derives the seed jets by subset DP and falling-factorial path-cover inversion, and replays the closure. Together with the triangular separation this gives exact minimal memory cardinalities of **27 for modulus 9** and **15,625 for modulus 25**.
