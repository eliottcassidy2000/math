# Modulo 290: local pair wheel, an oriented filter, and the 145-state quadratic map

Scope: proved elementary results and separately identified finite computations. This is not a proof of infinitude of twin primes or of minimality of a tournament carrier. Prior work is the repository file `05-knowledge/results/jacobian_eth_prime_wait_20261102/PRIME-MRB.md` at commit `4c6f0ca`; it already establishes the least-prime clock and MRB filter, so those are inherited results, not novel claims. It also records the `U(290)/{±1} ≅ C2×C28` calculation and explains why a 56-element Cayley graph is not the nonregular 290-edge order-six metagraph.

## Pair wheel and triangular coordinates

Let V be the units modulo 290 and let the directed +2 edges be (r,r+2) whenever both endpoints lie in V. There are 112 vertices and 81 edges. More generally, for even m=2^a product p_i^e_i, the number of directed +2 unit edges is

    2^(a-1) product p_i^(e_i-1)(p_i-2).

Proof: independently modulo each odd p_i, exactly two residues, 0 and -2, are forbidden, and parity admits exactly the odd residues modulo 2^a. CRT multiplies the counts. In particular, 81=3·27. A vertex or edge represents a local residue condition, not primality. For p,p+2 both >29, every twin pair lands in the wheel; exceptional actual twin pairs involving a prime factor 5 or 29 include (3,5), (5,7), (29,31). The actual twin pair (3,5) is a hostile example to an unconditional assertion that all pairs are wheel edges; (7,9) is a hostile example to the converse.

In the +2 ordering of all 145 odd residues there are 33 deleted sites, but two adjoining deletions mean 31 surviving runs. The graph has 25 paths of four vertices, two of three, two of two, and two isolated vertices 27 and 263. Verification: each deleted multiple of 5 produces 29 equally spaced cuts on the 145-cycle and surviving runs of length four; four additional odd multiples of 29 are inside these runs. They split two runs into lengths 1+2 and two into 3+0, yielding the stated counts. Equivalently, the exact enumeration in `check.py` reconstructs every path independently. The wheel graph has the single wrap edge (289,1), representing possible pair residues (-1,+1), but no assertion that each residue edge yields a twin pair is justified.

For a pair of odd positive integers write p=c-h, q=c+h, 1≤h<c. The triangular region of (c,h) has local admissibility gcd(c-h,290)=gcd(c+h,290)=1. For each fixed h modulo 290 the number of admissible c modulo 290 is

    (4 if 5|h else 3) × (28 if 29|h else 27),

namely 112,108,84,81 in the four cases. Proof: c has parity opposite h; modulo an odd prime p, forbid ±h, which are one residue when p|h and two otherwise. CRT completes the count. This is an exact jagged periodic mask on the triangular array, not a primality detector. A finite-state cellular automaton or FRACTRAN/Turing-machine construction could implement this mask, but no assertion here identifies such an implementation with Conway's PRIMEGAME; finite modular masks cannot themselves decide primality.

## First ±1 primes and orientation

For each even m>2, both target classes are coprime to m and are distinct, so Dirichlet's theorem gives a prime in each class. At m=2 the classes coincide. Reproduction via elementary trial division shows first +1 prime 1451 and first -1 prime 2029, whereas the first twin pair crossing from -1 to +1 is (6089,6091). Thus first hits do not generally pair, and their existence does not give any twin-prime existence theorem. The published least union hit 1451 has eligible-prime rank 227, versus 56=φ(290)/2 orientation-merged classes (prior work). Neither a geometric waiting law nor a recurrence for twin primes follows from CRT or Dirichlet.

Let f(n)=n^(1/n)-1. Then f(n)=log(n)/n+O(log²(n)/n²), and f is eventually decreasing to zero. A periodic coefficient sequence a_n of mean μ obeys

    Σ_(n≤N) a_n f(n) = (μ/2)(log N)² + C_a + o(1).

Indeed the zero-mean part has bounded partial sums and converges by summation by parts; Σ_(n≤N) log(n)/n − 1/2 log² N converges, and the error series converges absolutely. For the alternating original sign (-1)^n, every unit is odd and the mean is -112/290=-56/145. Thus the inherited result strengthens to -28/145 log² N plus a constant and o(1). Filtering instead to the unoriented union ±1 gives -1/290 log² N plus constant and o(1). Assigning opposite weights to the two orientations has mean zero, hence converges in ordinary order; it does not converge absolutely since the absolute weights have density 2/290. Consecutive block differences f(290k+289)-f(290k+291) are absolutely summable, but this does not make the original termwise series absolutely convergent.

## A precise bridge to -7/4, -29/16 and 145

For the rational quadratic map f(x)=x²−29/16, the known exact orbit is

    -7/4 → 5/4 → -1/4 → -7/4.

The final arithmetic step is (1−29)/16=−28/16. Write x=a/4 with a odd and define F(a)=(a²−29)/4. It is an odd integer because odd squares are 1 mod 8. Moreover

    F(a+290)−F(a)=145(a+145) ≡ 0 (mod 290),

so F descends to the 145 odd residues modulo 290. The known 3-cycle becomes exactly 283→5→289→283; the intermediate residue 5 is not a unit, so this is not a cycle in the prime wheel.

The complete 145-state graph has eight cycles, all period three. Seven basins have 20 vertices apiece with preperiod histogram 3,7,10 at distances 0,1,2, and one basin has five vertices with histogram 3,2,0. The small cycle is (29,203,145); all eight explicitly appear in the checked receipt. The rational cycle lies in the 20-state basin.

Proof: odd residues modulo 290 identify with Z/145 ≅ F5×F29. In F5, F(a)=4(a²+1), all five points enter 0→4→3→0 after at most one step, and exactly 3 points are on the cycle. In F29, putting y=a/4 conjugates F to y². Nonzero elements form a cyclic group of order 28=4·7. Squaring kills the 4-part within two steps; the resulting seven odd-order elements consist of one fixed element and two 3-cycles, while zero is another fixed point. A point in the odd-order subgroup has four predecessors under the appropriate eventual projection; zero has one. Combining each period-d cycle with the mod-5 3-cycle gives gcd(3,d) product cycles. Thus the 29-zero factor gives one 5-state basin, the nonzero fixed factor gives one 20-state basin, and the two period-3 factors each give three 20-state basins. Also the max of the two factor preperiods gives 3,7,10 for each larger basin.

The finite modular graph says nothing by itself about boundedness of integer orbits. Here an independent real-growth argument completely classifies the odd-integer lifts: if |a|≥9 then F(a)≥13, and F(a)>a for a≥9, so the orbit increases without bound thereafter. Each of the eight odd values -7,-5,-3,-1,1,3,5,7 enters the known three-cycle, as can be checked directly. Thus exactly these eight integer numerators give bounded rational orbits x=a/4; this is not a classification of arbitrary rational starting values. A concrete hostile lift is 9→13→35→299→…, which diverges although its residues form the modular cycle 9→13→35→9. The +2 wheel is not invariant under F: 9 and 11 are wheel-adjacent but map to 13 and 23, which are both units and are not adjacent. In particular, the 145 states and the 28 in the finite-field unit order arise from the explicit CRT construction; this does not identify the order-seven tournament carrier, prove minimality, or give a causal interpretation of the coincident MRB coefficient.

## Reproduction

    python3 check.py > normal.out
    python3 -O check.py > optimized.out
    cmp normal.out optimized.out

Uses standard-library exact arithmetic and trial division, and explicit failure checks that remain active with -O. Floating-point MRB partial sums are only a numerical sanity check, not a proof. The path census, first-prime search and functional graph search are finite exhaustive tests. The proofs above are independently checkable without trusting these outputs. No external mathematical sources are required for the elementary claims except the classical Dirichlet theorem and cyclicity of finite-field multiplicative groups; for the former see the English translation of Dirichlet's original 1837 work at https://arxiv.org/abs/0808.1408 .

An independent audit of the mod-290, filter, and quadratic-map claims is recorded in [AUDIT.md](../prime290_audit/AUDIT.md). The separate [general wheel theorem](WHEEL_GENERAL.md) has its own [independent audit](../prime290_audit/WHEEL_GENERAL_AUDIT.md).
