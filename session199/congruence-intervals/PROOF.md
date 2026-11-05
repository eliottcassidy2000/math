# Fixed-length residue forcing and physical lower approximation

Session 199. **PROVED:** every compatible joint one-apex residue occurs in
every sufficiently long aligned binary interval, with an explicit uniform
length. Consequently, for each fixed modulus m, every sufficiently large odd
N has a physical balanced C1 count H at most N with H congruent to N modulo m
and N-H < C(m) sqrt(N). The construction has only one root change, including
its final balanced closure. **Not exact coverage:** this same source family
omits a positive fraction of odd integers.

## 1. Physical input and the unexpected connection

The frozen source compiler has seed t=(1,0,2,1) and chronological moves

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

Here D is the Hamilton-path count, E the root-ending count, U the reserved
cover count, and d the root-deletion count. These are actual rooted graphs,
not arbitrary integer quadruples. Every U>=2 has a zero-J source

    s(U)=(g(U), U mod 2, U, 1).

It is read from the binary expansion of n=U-1: discard the leading 1, map
bit 0 to Z and bit 1 to O, and apply that word to t. The source has
bit_length(U-1)+1 vertices. Put f(n)=g(n+1). The exact recurrence is

    f(0)=1,
    f(2x+1)=f(x),
    f(2x)=f(x)+x+(x mod 2).                       (1)

In particular f is odd and 1<=g(U)<=U. The inequality follows either from
the physical source cone or directly by induction in (1). Appending OJ to
any source gives a balanced strong tournament with exactly one Hamilton
cycle and H=D+2U. These physical facts are frozen in dependencies/FOUR_BLOCKS.md.

The least-used sidecar in the earlier modular proof is a zero-bit self-loop
at the residue n=0. It allows words of different lengths to be padded to a
common length without changing either residue. This adds a real interval
constraint to the modular construction; mere profinite density does not.

A seemingly unrelated prior route is THM-4469, the Mahler bridge for adjacent
Collatz blocks (dependencies/MAHLER_BRIDGE.md). Its useful mechanism is to
retain a word cylinder together with the arithmetic congruence, rather than
projecting to either one. Here a binary prefix is the real interval cylinder,
and a residue-correcting suffix is a finite exact sidecar. No Collatz or Mahler
statement is imported. The map sends a binary cylinder and target residues
to one physical source inside that cylinder. It preserves the prescribed
joint residues and the prefix interval, loses exact scalar target equality,
and needs the complete ordered suffix. The cheapest hostile test is an
arbitrary nonzero prefix, where the old from-zero loop formula alone does
not give the desired residues.

## 2. Every binary prefix admits a common-length correcting suffix

Let M>=2 be even. Write M=2^k q with k>=1 and q odd, and b=ceil(log2 M).
Choose a>=2 with 2^a=1 modulo q; for q=1 use a=2. Such a is effective:
for example 2 phi(q) works when q>1. Set

    ell=a+k,
    T=(M/2-1) ell,
    L=2b+T,   B=2^L.                              (2)

**Theorem 1.** For every positive binary prefix P, every residue u modulo M,
and every odd residue z modulo M, there is a binary suffix v of exactly L
bits such that the integer n with binary expansion P followed by v satisfies

    PB <= n < (P+1)B,
    n+1 = u (mod M),   f(n)=z (mod M).            (3)

Thus each aligned interval [PB+1,(P+1)B] of possible U values, P>=1,
contains a physical source with the specified joint (U,g(U)) residues.
The suffix is constructed explicitly, not found by an unbounded search.

**Proof.** Reading one binary bit updates (n,f) by (1). Modulo an even M,
the residue pair determines both transitions, including the parity term.
Let c be the residue -2^b P modulo M, chosen in [0,M-1]. Append its b-bit
word. It fits because 2^b>=M. The resulting state has n=0 modulo M; write
its second coordinate F0 modulo M. F0 is odd.

At n=0 modulo M, a zero bit leaves both residues unchanged. The word

    W=1^a 0^k

returns n to zero and changes f by

    Delta=(2^a-1)(2^k-1)+1.                       (4)

Indeed the a ones leave f unchanged and send n to y=2^a-1. The following
zeros add y(1+2+...+2^(k-1))+1. The final n is 2^k y=0 modulo M.
Moreover gcd(Delta,M)=2: modulo q the increment is 1; if k>=2 it is 2
modulo 4, and if k=1 it is even and coprime to q.

Let r be u-1 modulo M in [0,M-1]. Appending the b-bit word for r from
n=0 changes f by f(r)-1. Hence exactly one j in [0,M/2-1] satisfies

    F0+j Delta+f(r)-1 = z (mod M).                (5)

For M=2 this means j=0. Append T-j ell zero bits, then j copies of W,
then the b-bit word for r. All padding is at n=0, so it changes neither
residue. The suffix has b+(T-j ell)+j ell+b=L bits and has the required
endpoint. The original positive prefix remains unchanged, proving the
interval assertion and physical realizability. QED.

This is stronger than the earlier arbitrarily-large residue witnesses:
it controls source height inside every aligned interval, uniformly over the
entire existing binary prefix and both target residues. It does not assert
that every residue appears at every shorter length.

## 3. A shorter dyadic suffix

**Theorem 2.** For M=2^k, Theorem 1 holds with

    L_1=2,
    L_k=(5k^2+k-6)/2 for k>=2.                   (6)

This is a bound, not a minimal suffix length.

**Proof.** Reset and finish with k bits each as above. The loop 1^a 0^k
has increment 2-2^a modulo 2^k. To obtain -2j, 0<=j<2^(k-1), read j in
binary. Its bit i=0 uses the a=2 loop. Each bit i>=1 uses one a=i+1 loop
and one a=2 loop, since their combined increment is -2^(i+1).
All loops return n to zero. For k>=2 their total length is at most

    T_k=(k+2)+sum_(i=1)^(k-2)(2k+i+3)
       =(5k^2-3k-6)/2.

Pad the unused length with zero bits at n=0. Adding the k-bit reset and
k-bit finish gives (6). For k=1 there is no f-residue adjustment, so L=2.
The target adjustment in (5) is even, so it always equals -2j modulo M.
QED.

One can use either bound, whichever is shorter. For arbitrary even M the
simple choice a=2 phi(q) gives L=O(M^2); for dyadic M the sharper L is
O((log M)^2). No polynomial bound for B=2^L is claimed.

## 4. Congruence-exact lower approximation to every large odd target

Fix m>=1, put M=2m, choose any valid L from Theorem 1 or 2, and let B=2^L.

**Theorem 3.** For every odd integer N>=9B there is an explicitly constructed
balanced strong tournament with exactly one Hamilton cycle and count H such
that

    H<=N,   H=N (mod m),   0<=N-H<3 sqrt(BN).     (7)

It has at most floor(log2 N)+3 vertices. Its input to the final OJ closure
uses only O and Z; the whole displayed source word uses exactly one J.
In particular the available balanced C1 counts have gaps O(sqrt(N)) even
within any prescribed compatible congruence class, with a constant depending
on the fixed modulus.

**Proof.** Let Q be the largest power of two with BQ^2<=N. Then Q>=2 and,
putting x=sqrt(N/B),

    Q<=x<2Q.

Set

    P=floor(N/(B(Q+1)))-1.

Because x>=3 and Q<=x, x^2>=2(x+1)>=2(Q+1), so P>=1. Use Theorem 1 in
the aligned U interval [PB+1,(P+1)B], with target residues

    U=1 (mod M),   g(U)=N-Q (mod M).

The second residue is odd. Thus U is odd, g(U)+QU=N modulo M, and

    N/(Q+1)-2B < U <= N/(Q+1).                  (8)

Take s(U), append r=log2 Q-1 copies of O, and then append OJ. Every step
is a physical source operation. The resulting balanced C1 count is

    H=g(U)+QU.

Since 1<=g(U)<=U, (8) gives H<=N and

    N-H < N/(Q+1)+2BQ <= N/Q+2BQ
        = sqrt(BN)(x/Q+2Q/x) <= 3 sqrt(BN).

The last inequality holds for x/Q in [1,2]. The congruence modulo M is
stronger than the claimed congruence modulo m.

Finally, n=U-1<N/Q. The total number of graph vertices is

    bit_length(n)+1+r+2 = bit_length(n)+log2 Q+2.

Since nQ<N and Q is a power of two, bit_length(n)+log2 Q=bit_length(nQ)
<=bit_length(N)=floor(log2 N)+1. This proves the vertex bound. QED.

The proof can be read as a constructive real-and-modular interpolation:
a suffix chooses an actual source in a nearby U interval with the right
residue, and a physical dyadic tail converts the interval error and the
uncontrolled g(U) term into errors of opposite scale. Balancing them gives
the square-root bound. For fixed m, H/N tends to 1 uniformly in odd N.

## 5. Scope, hostile, and exact missing step

There is no implication from (7) to exact coverage. The zero-J readout
family {g(V)+2V:V>=2} is precisely the unit-deletion family, and the frozen
four-block proof gives its positive-density omission bound

    |{g(V)+2V} intersect [1,N]| <= cN+O(log N),
    c=sum_(j>=0)1/(2^(j+2)+1)=0.431166447...<1/2.

Our H is in that same family with V=QU/2. Thus the upper bound on gaps and
a positive density of missing odd targets are simultaneously true.
Taking m=N does not repair the issue: B(2m) depends on m, and the required
N>=9B(2m) is not supplied when the modulus grows with the target. Nor does
an even deficit N-H authorize adding that many Hamilton paths by a graph
operation. A useful next exact lemma would be a physical, height-controlled
repair of such deficits that preserves the appropriate surrounding state.
No such repair is asserted here. The four-block all-order inverse and the
{7,21} theorem remain open in this packet.

## 6. Verification and dependencies

The verifier has explicit failing checks that remain active under python -O.
It tests the finite-state endpoint construction on arbitrary modular states,
recomputes actual large integers via a separate halving recurrence, checks
every one-apex source in a small range by the O/Z compiler, verifies many
approximation outputs and the exact squared error inequality, and counts
small literal balanced graphs independently by subset path/cycle DP.
Normal and optimized receipts must agree. Finite tests validate code and
boundary conventions; the displayed inductions prove the all-order claims.

Dependencies are frozen local copies of the physical four-block proof and
the earlier JOINT_MODULAR_GADGET.md. MAHLER_BRIDGE.md records provenance of
the cylinder-plus-congruence idea only; it is not a theorem dependency.
