# Compiling an exact Hamilton-path count into FRACTRAN

For each binary word w beginning 1 and ending 0, the accompanying `fractran.py` constructs a finite ordered FRACTRAN program and an initial integer. The program halts with the exact Hamilton-path count of the corresponding strong one-apex tournament in the exponent of prime 5. The number of fractions is 6(|w|-2)+2. This is a family of word-specialized programs, not a fixed universal FRACTRAN interpreter; no claim is made that the example below has Conway's 14-fraction prime-game behavior.

## Tournament and count

For w=w_1...w_n, orient the base vertices by i→j for i<j and add apex 0 with 0→i exactly when w_i=1. For a word beginning 1 and ending 0, the count is

    H = 2 + Σ_{j<i, w_j=1, w_i=0} 2^(i−j−1) + #(adjacent 01 transitions).

Thus appending 0 adds the integer represented by the preceding prefix, while appending 1 adds 1 exactly after a 0. The [binary bank note](NOTE.md) proves this path-splitting identity and the strong-connectivity condition.

## Compiler

Let B=2, T=3, H=5 denote register primes; their exponents are the binary-prefix integer, scratch, and path-count accumulator, respectively. Every control state is a distinct prime p_i,q_i,r_i,s_i not among B,T,H. Set initial exponents B=1, T=0, H=2 and activate the first state. The initial 2 is the prefix accumulator, not the Hamilton count of the two-vertex prefix tournament.

For each bit b of w except the initial 1 and final 0, let delta=1 precisely when this bit is 1 and its predecessor is 0. The following six fractions are placed in this order, with p_next the next bit's entry state:

    (q_i T^2 H^(1-b)) / (p_i B),
    p_i / q_i,
    r_i / p_i,
    (s_i B) / (r_i T),
    r_i / s_i,
    (p_next B^b H^delta) / r_i.

At p_i, the first two fractions repeatedly move one unit of B to two units of T, adding one to H exactly when b=0. Once B is empty, the third fraction changes phase. The next two repeatedly restore T to B, then the last changes to the next bit, adds b to B, and adds delta to H. The first applicable fraction rule supplies both zero tests: the fallbacks are reached only when the earlier transfer is impossible. Thus a bit sends (B,H) to (2B+b, H+(1-b)B+delta).

For the final 0, it suffices to add B to H and halt, rather than preserving B. Use two further fractions, q H/(p B) and p/q, and no fallback from p. The machine halts exactly when B=0. The scratch register is zero, so the final integer is p times 5 raised to the exact count. The construction is total: each transfer consumes a nonnegative counter, and there are finitely many bit phases. The literal integer representation can grow exponentially in the size of the input word; no polynomial-resource claim is made.

## A 14-fraction instance

For w=1010, the compiler produces initial integer 350 and the ordered list

    495/14, 7/11, 13/7, 34/39, 13/17, 19/13, 207/38,
    19/23, 29/19, 62/87, 29/31, 370/29, 205/74, 37/41.

The ordinary FRACTRAN first-applicable rule halts in 32 steps at 72,265,625 = 37*5^9. Independently counting paths in the specified tournament gives 9. Its count is a fixed finite example and is not evidence for universality.

## Check

`python3 fractran.py` and `python3 -O fractran.py` both literally execute the generated integer FRACTRAN programs and compare their final 5-adic exponents to independent subset/endpoint Hamilton-path dynamic programming for all 127 strong words of lengths 2 through 8. `python3 fractran.py 1010` prints the explicit instance. The code also requires the other register exponents to be zero at halting. The symbolic transfer argument proves the construction for all allowed finite words; these tests are controls.

The supplied simulator caps execution at 1,000,000 steps and may reject longer valid instances; the termination proof concerns the mathematical programs, not unbounded execution by this script.
