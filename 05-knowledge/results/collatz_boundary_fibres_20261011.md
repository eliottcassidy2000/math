# Real projective boundaries, arithmetic fibres, and a fixed-point obstruction

2026-10-11. **PROVED:** the elementary closure and obstruction statements
below. **FINITE-EXACT:** the companion controls. **OPEN:** a source-specific
completion rank for the unresolved Collatz children. Compactifying the state
does not establish convergence of a positive integer orbit.

## Inheritance and purpose

The [marked-ray geometry](simplex_sphere_cube_geometry_20261009.md) keeps
current value, original source, and a constant coordinate through simplex,
sphere, and cube norm sections. The
[hedgehog comparison](collatz_hedgehog_family_20260927.md), Sections 1–2,
already separates escape from cycles and interprets growth families as
shadows of negative rational cycles. We derive the particular Mersenne
identities here rather than importing an unspecified valuation-word clock.
[THM-4604, reducible Collatz words and the carry extension cocycle](../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md)
already says why a trace forgets carry. The
[debt observer](collatz_debt_observer_integration_20261009.md) keeps the exact
native word cells that such a quotient loses.

The new connection is an explicit boundary map: setting the constant
coordinate to zero removes the carry, while retaining a 2-adic driver keeps
the branch rule. A fixed point on that enlarged boundary gives a concrete
obstruction to one uniform continuous rank ansatz. The least-used coordinate
is the simultaneous real height / 2-adic address, not another scalar slope.

Board: marked ray, arithmetic address, affine carry, valuation guard,
forward shadow, source-specific rank. The cheapest hostile is two words
with the same multiplier and different native cells.

## 1. One real point hides every initial odd valuation

For a positive odd source set
\[
s(n)=\frac{(n,n,1)}{2n+1},\qquad b=(1/2,1/2,0).
\]
Then
\[
\|s(n)-b\|_\infty=\frac1{2n+1}.
\]
All large sources approach the same real boundary point. But for every
integer \(a\ge1\), the exact initial valuation is \(v_2(3n+1)=a\) on
the unbounded progression
\[
n\equiv(2^a-1)3^{-1}\pmod {2^{a+1}}.                 \tag{1}
\]
Indeed this is exactly \(3n+1\equiv2^a\pmod {2^{a+1}}\).
Every such residue is odd. Thus even the two valuation classes 1 and 2
preclude a continuous discrete-valued valuation extension to \(b\).
Positive constant coordinate retains the integer exactly; taking the real
boundary forgets its arithmetic branch.

There is a lossless boundary repair. Regard each integer as a 2-adic integer
and take the closure of
\[
\{(s(n),n): n\ge1\text{ odd}\}\subset\Delta_2\times\mathbb Z_2.
\]
Its part on the zero-marker boundary is exactly
\[
\{b\}\times\mathbb Z_2^\times.                      \tag{2}
\]
One inclusion follows because a vanishing marker forces \(n\to\infty\),
and a limit of odd 2-adic integers stays odd. Conversely, for any odd
\(z\in\mathbb Z_2\), let \(r_k\) be its canonical residue modulo \(2^k\)
and take \(n_k=r_k+2^{2k}\). Then \(n_k\to\infty\) in the real sense
and \(n_k\to z\) 2-adically. This proves (2), including all quantifiers;
finite residue samples are controls rather than its proof.

The other norm sections transport this statement by their continuous radial
homeomorphisms. This is a product with an arithmetic fibre, not a claim that
the real simplex alone can distinguish the branches.

## 2. A continuous coupled dynamics

Use the Terras convention
\[
T(x)=\begin{cases}x/2&x\text{ even},\\(3x+1)/2&x\text{ odd}.\end{cases}
\]
The same formula defines a continuous map on \(\mathbb Z_2\): parity is
clopen, and each numerator is divisible by 2 on its declared branch.
An output residue modulo \(2^k\) is determined by an input residue modulo
\(2^{k+1}\). A finite residue register consequently consumes precision;
it is not an autonomous substitute for the full 2-adic coordinate.

Let \(p=[x:n:c]\) be a nonnegative real ray, with \(n\) the original
source and \(c\) the marker, and retain two arithmetic coordinates
\((\xi,\eta)\in\mathbb Z_2^2\) for current and original values. If
\(a=\xi\bmod2\), define
\[
A_a=\begin{pmatrix}1+2a&0&a\\0&2&0\\0&0&2\end{pmatrix},\qquad
\mathcal T(p,\xi,\eta)=(\operatorname{section}(A_ap),T(\xi),\eta).
                                                               \tag{3}
\]
This is continuous on a compact norm section times \(\mathbb Z_2^2\).
Neither matrix annihilates a nonzero nonnegative ray, and the two parity
pieces are clopen. On an actual marked integer state
\((\operatorname{section}(x,n,1),x,n)\), it gives exactly
\((\operatorname{section}(T(x),n,1),T(x),n)\).

For precision, take the actual trajectory-monitor set
\[
\mathcal S=\{(\operatorname{section}(T^j(n),n,1),T^j(n),n):
                 n\ge1,\ j\ge0\}.
\]
Its compact closure is invariant under (3), by continuity and
\(\mathcal T(\mathcal S)\subset\mathcal S\). Points of the full ambient
product need not have consistent real/arithmetic coordinates; no such
consistency is inferred merely by writing (3).

## 3. The boundary is exactly where carry disappears

For an odd valuation word \(u\), write its affine value as
\((P x+B)/Q\), where \(P=3^{|u|}\), \(Q=2^{\sum u}\), and \(B\)
is its ordered carry. On the marked monitor the homogeneous matrix is
\[
M_u=\begin{pmatrix}P&0&B\\0&Q&0\\0&0&Q\end{pmatrix}.
\]
On the boundary \(c=0\), it acts as
\([x:n:0]\mapsto[Px:Qn:0]\). The carry disappears exactly, and where
\(n>0\), the ratio \(x/n\) is multiplied by \(\lambda=P/Q\).

For example, \(u=(1,3)\) and \(v=(3,1)\) have the same multiplier
\(9/16\), but carries 5 and 11 respectively. They agree on the entire
real boundary and differ on every positive-marker ray. Their exact native
source classes are disjoint:
\[
u:\ n\equiv19\pmod {32},\qquad
v:\ n\equiv29\pmod {32}.                                \tag{4}
\]
The congruences follow by imposing each valuation successively. They are
also checked by exact supplied-word replay over two full periods.

Consequently repeating a word gives the valid multiplier identities
\[
\lambda\mapsto\lambda^2,\qquad
J=\lambda+\lambda^{-1}\mapsto J^2-2.
\]
This is the precise surviving quadratic-dynamics connection. Neither
identity restores \(B\), authenticates a native word, or proves payment.
Here \(J\) has the displayed definition; it is not being confused with
the trace of a determinant-normalized affine matrix.

## 4. A fixed boundary point from actual Mersenne monitors

Let \(M_E=2^E-1\). Direct induction proves
\[
T^j(M_E)=3^j2^{E-j}-1\quad(0\le j\le E).                 \tag{5}
\]
For \(j<E\) the current value is odd; at \(j=E\) it is even.
For every \(1\le j\le E\), it exceeds \(M_E\). Thus the first \(E\)
Terras parity bits are 1, and no fixed universal bound on forward descent
time is possible. This is distinct from an accelerated odd valuation word,
whose leading run of valuation 1 here has length \(E-1\).

Now take an actual middle monitor with \(E=2k\):
\[
n_k=4^k-1,\qquad x_k=T^k(n_k)=6^k-1.
\]
We have \(x_k/n_k>(3/2)^k\), so its real marked ray tends to the current
axis \(e_x=[1:0:0]\). Also \(x_k\to-1\) and \(n_k\to-1\) in
\(\mathbb Z_2\), since their distances to \(-1\) are divisible by
\(2^k\) and \(2^{2k}\), respectively. Therefore
\[
z_*=(e_x,-1,-1)\in\overline{\mathcal S}.
\]
It is fixed by (3): \(T(-1)=-1\) and \(A_1e_x=3e_x\), which is the
same normalized ray.

**Uniform continuous-rank obstruction.** There is no continuous real
function \(h\) on \(\overline{\mathcal S}\) and constant \(\epsilon>0\)
such that
\[
h(z)-h(\mathcal Tz)\ge\epsilon                         \tag{6}
\]
at every actual positive nonterminal trajectory-monitor state.
The states just constructed have \(x_k>1\). Taking their limit in (6)
would give \(h(z_*)-h(z_*)\ge\epsilon\), a contradiction. The same proof
works if any fixed finite set of current labels is exempted from (6).

This rules out a particular forward proof ansatz. It does not exclude
strict but nonuniform drift tending to zero, a rank singular at the
boundary, an integer/source-dependent rank without such a continuous
extension, or a controller using authenticated inverse-child receipts.
It is not a positive integer counterexample: both arithmetic coordinates
of \(z_*\) equal the negative 2-adic integer \(-1\), and its real marker
is zero. Every actual state used to approach it is finite and marked.

## 5. Consequence for the next proof attempt

The [adaptive receipt rank](collatz_receipt_flow_geometry_20261009.md)
has a valid source-specific refinement theorem. To extend it, supply a
coarse rank, actual terminal receipts, and inclusion of every refined row
in the old certified row polytope. This note explains why a global smooth
picture alone cannot supply those premises: the boundary both loses carry
and contains an exact nonterminal fixed point for the forward extension.

A concrete **OPEN** ansatz is to use different ranks on arithmetic guard
cells, allow explicit blow-up near their negative-cycle boundaries, and
prove that paid child receipts transfer between cells before that rank
budget is exhausted. The required transfer inequalities must include the
ordered carry, original source, and actual terminal. A change of norm
section preserves them but cannot manufacture them.

## Reproduction and scope

The script checks odd residues through depth 9 with three unbounded lifts
each; initial valuations 1–16 on 16 lifts; currents 1–44 and original
sources 1–23 in all three sections; depth 1–10 precision controls; carry
5/11 controls and two native periods; Mersenne exponents 2–64; and midpoint
monitors \(k=1,\ldots,64\). All arithmetic is exact, and checks remain
active under optimization. The infinite results are proved above.

    python -B 04-computation/experiments/collatz_boundary_fibres_20261011.py
    python -B -O 04-computation/experiments/collatz_boundary_fibres_20261011.py

Both runs reproduce the saved output: **14,369 exact checks**. Independent
geometry review checked the boundary fibre, native residues, both clocks,
and the fixed-point obstruction. No ROOT certificate is generated here.
LF-normalized output SHA256:
`ac4ae5da62687049a0dd457f20c7b1289dd8802dcd9f7032ce1057c4c575b610`.
