# The clean guard has a complete dyadic clock but a sparse integer image

**Status:** proved for a fixed odd-order host; no unrestricted spectrum claim.

For a tournament (B), let (P_B(t)) count its ordered path covers by (t)
slots with empties allowed. The odd-cycle packing formula makes this an
integer-coefficient polynomial, with parity equal to the order of (B);
(P_B(1)=H(B)) is odd by Rédei. Add two guards (a,b) with
(b\to a\to B\to b), writing (G(B)). The
[earlier guard proof](../chanpak_guard_ternary_20261115/README.md)
gives (P_{G(B)}=L P_B), where
[
 (Lf)(t)=t^2f(t)+t(f(t+1)-f(t-1)),qquad
 H(G(B))=P_B(1)+P_B(2).
]
Let (h_k=(L^kP_B)(1)). The theorem below also applies to the evaluation
sequence from any odd (f\in\mathbb Z[t]) with (f(1)) odd.

**Theorem.** For (r\ge4), every consecutive block of (2^{r-2}) depths
starting at (k\ge r) maps bijectively to the residues
(h_4+4\mathbb Z\pmod{2^r}). More precisely, for distinct (k,k+d) with
(d>0) and (k\ge\max(4,3+v_2(d))),
[
v_2(h_{k+d}-h_k)=2+v_2(d).
]
There is a compatible bijective dyadic interpolation
(Phi_B:\mathbb Z_2\to h_4+4\mathbb Z_2) with
(v_2(Phi_B(x)-Phi_B(y))=2+v_2(x-y)) and
(h_k\equivPhi_B(k)\pmod{2^k}) for (k\ge4).
This interpolates the stable residue system; it does not assert equality
of the initial integers (h_k) with (Phi_B(k)).

## Stable residue dynamics

Read (L) on the (2^r)-component module of evaluations, with argument
indices modulo (2^r). Modulo 2, the even rows are zero and the odd rows
are the identity on odd coordinates plus a difference of even
coordinates. Thus (L\bmod 2) is idempotent (the off-diagonal block need
not vanish).

Over (\mathbb Z/2^r), split the module (M) into a nilpotent summand
and an invertible summand using the eventually stable image and kernel
of powers of (L). To see the decomposition, take a common large
exponent (j) at which both stabilize: (L^j) is invertible on its
image, so each vector differs from an element of the image by an element
of the kernel, and the intersection is trivial. Reducing this direct
sum mod 2 shows that the nilpotent restriction is zero and the
invertible restriction is the identity there, as both are idempotent.
In particular (L^r) kills the nilpotent summand, and
(D=L^2-I) sends the invertible summand into four times itself.

The only seed check is modulo 16. Every odd power of (t) has the same
evaluation vector as one of (t,t^3,t^5,t^7): above degree 3 even
coordinates vanish, and odd coordinates have fourth power 1.
The following table uses the displayed recurrence and records every
coordinate of the last column, not just (t=1). Put
(u(t)=8[t\text{ odd}]).
[
\begin{array}{c|rr|c}
f & (L^4f)(1)&(L^5f)(1)&(L^6-L^4)f\\\hline
t&7&3&u\\
t^3&5&9&u\\
t^5&13&1&u\\
t^7&13&1&u
\end{array}
]
Since (Lu=u\pmod{16}), summing with the coefficients of (f) gives,
for every (k\ge4), ((L^{k+2}-L^k)f\equiv u\pmod{16}). In
particular the vector is divisible by 8, its value at 1 has valuation 3,
and the successive scalar increments are 4 modulo 8. The latter
increments persist because the two-step increments vanish modulo 8.

Now take (d=2m) and work modulo (2^{4+v_2(m)}) if determining the
leading bit, equivalently modulo (2^{3+v_2(d)}). At depths at least this
precision the vector (g=L^kf) is in the invertible summand, so
(Dg\in8M_{\rm inv}) (project the divisibility in the full module
onto that summand). Expand
[
 L^{2m}g-g=mDg+\sum_{j\ge2}\binom mj D^j g.
]
For (j\ge2) the valuation is at least
(v_2(m)-v_2(j)+3+2(j-1)>v_2(m)+3).
The leading term at 1 therefore has exactly valuation (v_2(m)+3).
An odd (d) gives valuation 2 immediately modulo 8. This proves the
distance law. It implies the modular bijection by comparing two
depths inside one block; increments by the block length preserve the
residue using the same distance law. Deep representatives of each
depth residue class give the compatible inverse limit and its metric law.

## A target-first decoder, and why the clock itself cannot close Chan–Pak

For a host of order (n) and (k\ge1),
[
 h_k\ge1+2^{n+2k-2}.
]
Indeed (H(G(B'))=P_{B'}(1)+P_{B'}(2)), and the empty collection of
odd cycles contributes 1 and (2^{|B'|}) respectively. Thus if a target
(N>1) exists in this family, the only possible depths are 0 through
(K=\lfloor(\lfloor\log_2(N-1)\rfloor-n+2)/2\rfloor), including 0
regardless of the bound. Choose (r\ge4) with (2^{r-2}>K).
After checking depths below (r) (and 0), the modular inverse returns
**at most one** candidate in ([r,K]). One exact evaluation decides
it. The scalar algorithm does not presume a residue match is a witness.
Precision is (O(\log\log N)) bits; a modular vector table uses
(O(\log N)) entries for the least such (r). Exact evaluation can
propagate the displayed recurrence only in the triangular dependence
region between the initial arguments and (t=1), rather than count
Hamiltonian paths of a growing literal tournament by a subset census.

A sharper separation includes ternary information. For the singleton
seed, the mod-9 vectors at depths 3 and 4 are
((0,6,6,0,0,0,0,3,3)) and ((0,3,3,0,0,0,0,6,6));
(L) interchanges them. Modulo 8 the odd-depth scalar is 3 and the
even-depth scalar is 7 from depth 4 onward. Hence the singleton guard
chain has **all** compatible joint residues
(39,51\pmod{72}) at every further dyadic precision:
for each (r\ge3), every residue (N\pmod{9\cdot2^r}) reducing
to one of those two classes is attained arbitrarily deep. There are
only (O(\log X)) actual integer counts at most (X). This provides
a physical, one-parameter example where arbitrarily fine joint filters
leave a positive-density candidate set but settle none of its generic
absolute placement.

The distinction suggests a useful proof test: an exact target decoder
must bound the *physical* depth/size at which its local inverse is used.
A new guard word or a changing host must increase the number of
absolute baselines, not only refine the stable residue dynamics.

## Programming independent ternary and dyadic digits with only six towers

There is a compact upgrade which removes the compatibility condition
between the mod-9 digit and the dyadic depth parity. Start from an
odd-order core (B) with (h=H(B)) coprime to 3. Before wrapping, pad it
by an ordinal sum with (4j) singleton components, for (j=1,2,3).
The three input polynomials (f_j=t^{4j}P_B(t)) are identical as value
vectors mod 16: each odd coordinate has fourth power 1 and each even
coordinate contributes zero after multiplying by (P_B(t)).
The clock therefore has the same common mod-4 class for all three
towers. This class is (h\pmod4), since every monomial is now of degree
at least 5 and the seed table has value 1 mod 4 in those rows.

The prior [mod-9 guard calculation](../chanpak_guard_ternary_20261115/README.md)
is exactly
[
 (L^k f)(1)/3\equiv(-1)^{k-3}(f'(0)+f'(1))\pmod3\qquad(k\ge3).
]
For the padded cores (f_j'(0)=0) and (f_j'(1)=P_B'(1)+4jh);
these three derivatives exhaust the ternary values. For a prescribed
residue modulo (9\cdot2^r), first read the required depth parity
from the common mod-8 clock, then choose (j) from the prescribed
mod-9 residue, then choose the dyadic depth in that tower. Consequently
these three towers jointly hit every residue that is a multiple of
3 mod 9 and equals (h\pmod4), at every precision (r\ge3).

For (h\equiv1\pmod4) take the singleton; for (h\equiv3\pmod4)
take the five-vertex tournament with outgoing bit rows
((24,1,3,6,14)), whose (H=11) is checked by literal subset DP.
Thus **six literal towers cover every odd multiple of 3 modulo
(9\cdot2^r), for every (r)**, with each residue attained arbitrarily
deep. Those six towers still supply only (O(\log X)) absolute
counts at most (X). Selecting a new host changed the local source
of the ternary digit, but did not repair the integer placement.
