# A strong one-apex family with prime-power modular coverage

Status: PROVED. This is modular coverage, not exact integer-spectrum coverage. No assertion of simultaneous coverage for a general composite modulus, novelty, or existence of 7 or 21 is made.

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
