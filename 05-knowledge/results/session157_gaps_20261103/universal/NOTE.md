# All eligible residue classes from strong one-apex tournaments

Status: PROVED. For any positive integer M, every residue modulo M occurs as a strong tournament's Hamilton-path count if M is odd, and every odd residue occurs if M is even. The tournaments can be chosen with feedback vertex number at most one. This is a congruence-coverage theorem, NOT an exact integer-spectrum or positive-density theorem; it neither proves nor disproves the proposed exceptions 7 and 21. No claim of novelty is made.

## Explicit construction

Write M=2^e m with m odd. Let r be arbitrary modulo M, odd if e>0. Choose an integer a>=1 such that

    h(a) = (10*4^a+6a-1)/9 = r (mod 2^e).

Such an a exists: for large a, h(a)=(6a-1)/9 mod 2^e, and varying a modulo 2^(e-1) gives every odd class. If e=0 set a=1.

Let L be the least common multiple of ord_p(2) over odd primes p dividing m, or L=1 when m=1. Choose k>=max(2,e+1) with k=1 mod L, and write q=2^k. Then q=2 mod p for every p dividing m; in particular q-1 is relatively prime to M. For n>=1 take a transitive base with one-apex orientation word

    1 (0^(k-1) 1)^n (10)^a,

where the apex beats exactly the 1-bits. This is strong: the apex beats the first base vertex, which reaches every other base vertex; every base vertex reaches the last, which beats the apex.

## Exact path count and 2-adic preservation

Delete the apex and split any path into the two transitive pieces on either side. The two exterior placements contribute 2 paths. A 0 following a 1-bit at an earlier position j contributes 2^(i-j-1) for that pair, while an adjacent 01 contributes one. In left-to-right processing, let B be the sum of 2^(i-j-1) over preceding 1-bits. Update B to 2B+1 at a 1 and to 2B at a 0; a 0 adds its current B to H, and an adjacent 01 adds one. The initial 1 leaves B=1, H=2, and previous bit 1. Each block 0^(k-1)1 sends

    B -> qB+1,
    H -> H + (q/2-1)B + 1,

ending with previous bit 1. The condition k>=e+1 makes q/2=0 mod 2^e. Thus (H,B,previous)=(2,1,1) is fixed modulo 2^e by any number of blocks. The final suffix (10)^a therefore gives the same residue as the short word 1(10)^a, namely h(a). This includes e=0 vacuously.

For reference, all-n counts have a simple exact formula. Put

    B_j = (q^(j+1)-1)/(q-1),
    C_a = 2(4^a-1)/3,
    D_a = (4*4^a+6a-13)/9.

C_a and D_a are integers. The suffix contributes C_a B_n + D_a, and hence

    H(n) = 2+n+(q/2-1)*sum_{j=0}^{n-1} B_j + C_a B_n+D_a
         = A q^n + [q/(2(q-1))] n + C,

where A and C have denominators dividing (q-1)^2, because C_a,D_a are integers. The middle coefficient is a p-adic unit for each odd p dividing m.

## Simultaneous odd-prime CRT

For every p|m put L_p=ord_p(q); it divides p-1, so every prime divisor of L_p is smaller than p. Process the prime factors of m and of the various L_p in increasing order, fixing n modulo enough powers of each. When processing p|m, the already fixed smaller-prime coordinates determine n mod L_p. Write n=r_p+L_p t. For 1<=j<=v_p(m), q^(L_p p^(j-1))=1 mod p^j, so increasing t by p^(j-1) changes H modulo p^j by

    [q L_p/(2(q-1))] p^(j-1),

whose coefficient is a p-unit. It follows inductively that one can select n modulo p^v_p(m) to obtain the prescribed r modulo that prime power. Extend this coordinate arbitrarily to the precision needed by the orders for larger primes. CRT then joins all the selected coordinates. Increasing n by a multiple of their product preserves every odd-prime condition, while the 2-adic condition holds for every n; thus n>=1 can be ensured.

The construction produces arbitrarily large n, and H(n) grows without bound, so in fact every eligible residue contains infinitely many such strong one-apex counts. It still does not produce every eligible integer.

## Verification

Run `python3 check.py` and `python3 -O check.py`. The script independently builds literal tournaments and runs Hamilton subset/endpoint DP for small parameters; checks the derived formula; and executes the constructive CRT for every eligible residue for M<=250 and selected mixed/prime-power moduli. The modular computations are finite hostile controls for the proof, not a substitute for the symbolic recurrence or CRT argument.
