# A strong one-apex family with prime-power modular coverage

Status: PROVED. This is modular coverage, not exact integer-spectrum coverage. No assertion of novelty or existence of 7 or 21 is made. Strong one-apex coverage is restricted to residues odd when the modulus is even and nonzero modulo 3 when the modulus is divisible by 3; the all-eligible-residue extension uses ordinal sums and may be reducible.

For n>=1, let F_n have vertices r,y,x_1,...,x_{2n}. The subtournament y,x_1,...,x_{2n} is transitive in that order. Let r beat y and x_i for i odd; let x_i beat r for i even. It is strong: every x_i reaches x_{2n}, which beats r; then r beats y and y beats all x_i. (When i=2n use the trivial length-zero reachability to x_{2n}.) Deleting r leaves a transitive tournament.

Its Hamilton-path count is exactly

    h(n) = (10*4^n + 6n - 1)/9.

Proof. In a Hamilton path, split the transitive vertices around r. The two exterior cases contribute two paths. For nonempty sides let i be the maximum index before r and j the minimum index after r, counting y as index 0. Necessarily i=2k and j=0 or j odd. If i<j, then all transitive vertices are assigned and j=i+1; this gives n-1 paths. If i>j, vertices strictly between j and i may independently be assigned to either side, giving 2^(i-j-1). For each i=2k this sums over j to

    2^(2k-1) + sum_{l=1}^k 2^(2k-2l) = (5*4^k-2)/6.

Adding gives 2+(n-1)+sum_{k=1}^n(5*4^k-2)/6 = h(n). This determines all paths and is not an experimental fit.

For every odd prime p>=5 and every e>=1, the values h(n) attain every residue modulo p^e, with n>=1 and F_n strong. Put L=ord_p(4), so p does not divide L. For every k>=1, 4^(L*p^(k-1))=1 mod p^k, by the binomial theorem or LTE. Thus for n=Lt,

    h(L(t+p^(k-1))) - h(Lt) = (2L/3)*p^(k-1) mod p^k.

Here 3 is invertible modulo p^k and 2L/3 is a unit. At k=1 the map t -> h(Lt) is therefore a bijection modulo p. Inductively, changing t by a*p^(k-1), a=0,...,p-1, gives every lift of the previously chosen residue modulo p^k. Choosing positive representatives preserves the residues, and n>=1.

For every e>=1, h(n) attains every odd residue modulo 2^e. Indeed, when 2n>=e, 4^n=0 mod 2^e and h(n)=(6n-1)/9 mod 2^e; as n ranges over any complete residue system mod 2^(e-1), these are all odd classes. Shift n by multiples of 2^(e-1) until 2n>=e and n>=1.

Reproduction: `python3 check.py` and `python3 -O check.py`. The independent endpoint/subset DP checks the explicit tournament and count for n=1,2,3,4; transfer and reachable-residue checks provide further controls. Finite checks are not the proof of the all-n identity or the p-adic lifting.

## Simultaneous modular coverage, including powers of 3

The family gives the following stronger conclusion. For each positive M and each residue r modulo M, if r is odd whenever M is even, there is a tournament with Hamilton-path count congruent to r modulo M. The tournament can be chosen strong and one-apex whenever r is a unit at 3, meaning 3 does not divide r if 3 divides M. The general construction is an ordinal sum of directed triangles and one of the strong F_n. This does not imply any particular exact count or density.

First the strong family realizes every eligible residue r whose 3-part is a unit. We need two observations about h modulo powers of 3. The binomial expansion of 4^n=(1+3)^n gives

    h(n) = 1 + 4n + 10 sum_{j>=2} 3^(j-2) binom(n,j),
    h(n) = 1+n+binom(n,2) (mod 3).

Hence h(n) is 2 mod 3 when n=1 mod 3, and 1 mod 3 when n=0 or 2 mod 3. For k>=1, let u_k=(4^(3^k)-1)/3^(k+1). The binomial theorem gives u_1=7 and u_{k+1}=u_k mod 9; so u_k=7 mod 9. Therefore

    [h(n+3^k)-h(n)] / 3^k
      = (10*4^n*u_k+2)/3 = n (mod 3).

The division is exact. Thus, to lift a desired unit residue modulo 3^e, choose n=1 mod 3 for target 2 mod 3, or n=2 mod 3 for target 1 mod 3. At each later step, adding one of 0,3^k,2*3^k to n gives the three possible next digits, because n is a 3-adic unit. The same formula also shows that h(n) mod 3^e depends only on n mod 3^e.

To synchronize all prime powers, put L_p=ord_p(4) for each prime p>=5 dividing M, and let L be their least common multiple (L=1 if there are none). Each prime divisor of L_p is strictly less than p. Process the primes q dividing 6ML in increasing order and fix the residue of n modulo q^a to enough precision, namely a>=v_q(L), a>=v_q(M) for odd q, and a>=v_2(M)-1 for q=2. Choices are as follows:

* At 2, if 2^e divides M, choose n modulo 2^(e-1) to hit the desired odd residue using h(n)=(6n-1)/9 mod 2^e for n large. Extend this coordinate arbitrarily to precision a.
* At 3, if 3^e divides M, choose n modulo 3^e by the unit-lifting procedure above. Extend arbitrarily to precision a.
* At each p>=5 dividing M, all prime-power coordinates of n needed to determine n modulo L_p have already been fixed, since all divisors of L_p are smaller than p. Let n=r_p+L_p t. For each k>=1, changing t by p^(k-1) changes h(n) modulo p^k by (2L_p/3)*p^(k-1), a unit step. Thus t modulo p^e, and hence n modulo p^e, can be chosen to hit the desired residue. Extend arbitrarily to precision a.
* For all other q, fix an arbitrary coordinate.

CRT joins these mutually coprime coordinates. Adding multiples of their product makes n positive and large enough for the 2-adic simplification, without altering any specified output residue. This proves the simultaneous strong-family assertion. The proof requires the increasing-prime ordering: fixing n=0 mod the full L at the outset could incorrectly constrain the p-coordinate.

Finally suppose 3^b is the exact power of 3 dividing M and the desired residue has 3-adic valuation a<b. Choose a 3-adic unit u with 3^a u equal to r modulo 3^b, and choose the strong F_n simultaneously to have h(n)=u modulo 3^(b-a) and h(n)=3^(-a)r modulo M/3^b. All other eligible local conditions, including parity, are preserved. Take the ordinal sum of a directed triangles and F_n; every Hamilton path traverses the blocks in order, so its count is 3^a h(n). If r=0 modulo 3^b, use a=b and impose only the condition modulo M/3^b. If 3 does not divide M, use a=0. The result is an actual tournament in all cases, although the ordinal sum is not strong when a>0.
