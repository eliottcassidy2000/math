# Moving the length supplies every 19-adic residue without local switch choices

The [fixed-length physical bank](ADDITIVE_PORTS.md) only supplies a sparse set of 19-adic switch addresses. Its other degree of freedom, the length, gives the opposite result. Two constant-cut words, with **no choices within a word**, suffice to realize every unit modulo every power of 19. This is a modular result for the same block construction, not a new general modular-completeness claim.

Take the order-six block from the physical bank with $h=19,q=256$. Use the same apex out-cut $c$ on each of $n$ ordinal copies. With the response coordinates of the [port formula](ADDITIVE_PORTS.md), the split across distinct overlapping modules contributes
$$
s^*e^*\sum_{t=0}^{n-2}(t+1)q^{n-2-t}19^t .
$$
Here $t$ counts blocks strictly outside the visited overlapping interval; there are $t+1$ ways to distribute them on its two sides. All the other terms are multiples of $19^{n-2}$. For $n\ge k+2$ the tail of the formal geometric derivative vanishes modulo $19^k$, giving
$$
H_n(c)\equiv s^*(c)e^*(c)q^{n-2}(1-19/q)^{-2}\pmod{19^k}.
$$
In mask convention, cut 17 has $(s^*,e^*)=(94,187)$ and cut 4 has $(29,185)$, with products respectively $3,7\bmod19$. The former is a nonsquare and the latter a square. In both cases the literal tournament is strong: every module is strong, the apex can enter the first module and be reached from the last, and the ordinal arcs carry reachability forward.

The residue of $q$ has order nine modulo 19 and $q^9=1+12\cdot19\bmod 19^2$. Binomial lifting gives order $9\cdot19^{k-1}$ modulo $19^k$, so its powers exhaust the square subgroup of the units. The inverse-square factor is a unit square. Choose the appropriate cut from the target's square class; varying $n$ in its exponent class, with $n\ge k+2$, reaches the target residue. To realize nonunits, append an ordinal power of the 19-path block, multiplying by the needed power of 19 (these appended outputs need not be strong). Zero is reached by appending at least $k$ copies.

This sharply separates **exact integer coverage** from residue coverage. Every resulting positive count grows at least as $19^{n+t}$ with $t$ appended blocks. Thus the entire two-cut family with appendages supplies only $O((\log X)^2)$ integers up to $X$, despite meeting every residue modulo every $19^k$. A length adjustment can bypass the local selector for all modular targets and still miss almost every absolute target.
