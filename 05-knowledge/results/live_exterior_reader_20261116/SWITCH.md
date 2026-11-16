# Choose the gap first, and leave an entire exterior free

The signed address-reversal identity is exact because J is invariant under
reversal: at Z=82 and g=9, H(b)-H(rev(b))=100(b-rev(b)). The change need
not be small. A difference of two such contrasts is still a *virtual*
signed sum. The following construction instead compares two actual
tournaments on the same vertices with at most two changed arcs.

## The physical gap theorem

For u>=2 let V=2^u. Every even d in [-V,V], positive, zero or negative,
has at least V/2 different *lower-state anchors* H0 of order Theta(V^2)
in Hamiltonian count, with H1-H0=d and at most two actual arc reversals.
Here lower-state means central bit zero, not necessarily the smaller
count. All outputs are strong, and their vertex order is O(u). For any
fixed u,d the two selected readers and the bit to flip are fixed,
independently of all V/2 exterior choices. On the explicit even strand
below the anchor spacing is Theta(V).

Use the first two literal seeds and their common synchronizing zero from
[variable readers](../variable_reader_coverage_20261116/PROOF.md).
They share g=9 and Z=82 and have I=-11,-9; the seeds differ in one arc.
Apply m additional zeroes, giving Z=82*2^m. Use an n-bit address with u
bits U above one variable bit, and i=u+m+3 bits L below it. Thus
n=i+u+1, W=2^i=8*2^m V, D=2^n=2VW. Write X=rev_i(L) and
e(t)=t+(t mod 2). Flipping the variable bit 0 to 1, and optionally
swapping the first seed for the second (s=0 or 1), gives the exact gap

    d=62*2^m V - 9(e(X)+e(U)) + 2s.                 (1)

Indeed the address increase is W, the reversed-address increase is V,
and the interactions from below and above the bit are respectively
e(X) and e(U), by the distance kernel. No far arcs change.

Choose m in {0,...,5} and s in {0,1} so that the numerator
62*2^m V-d+2s is divisible by 18. After dividing the congruence by
two, the first term 31*2^m V runs over the six units modulo 9;
their translates by 1 together with the units cover all residues.
Write K=(62*2^m V-d+2s)/9. It is even and V<K<W. Now take **every**
even U in [0,V), put X=K-U, and reverse X to obtain L. This proves
actual selection with an entire exterior choice still free. Each
interior even value e(X), e(U) also permits the adjacent odd preimage,
which gives more phases without changing (1).

For the bit-zero anchor apply the exact cross interaction 2XU (there
is one separating zero):

    H0=I+2ZD+4gWU-2ZVX-2gXU + B,
    B=g(2L-J(L))-Z rev_u(U)-gJ(U).

On the even strand X=K-U, the variable smooth term is
(452*2^m V-18K)U+18U^2. Its slope is at least 308*2^m V,
using K<W, and is at most 488*2^m V. Because 0<=J(t)<=t, all B
lie in an interval of length at most 235*2^m V. Consecutive even U
therefore give strictly increasing anchors with spacing between
381*2^m V and 1211*2^m V (loose integer bounds). This proves
distinctness and gives a binary search for bracketing a proposed baseline
inside the strand's range. The absolute scale follows directly from
D=2VW and Z<=82*32: H0>=ZD+I and H0<=2ZD+18D+I.
Strongness is the same final-zero argument as for the inherited seeds.

## Why this is not an exact-target theorem

Equation (1) fixes the *difference*, not the anchor. At fixed gap its
smooth terms move together, giving a monotone quadratic rather than a
tangent direction. Within this one strand the anchor gaps actually grow
like V. The dreamer's four-score contrast instead allows cancellation
between independent mirror-pair scales: W_(j+1)-2 W_j is a power of
two times 3. These are different mechanisms, and the theorem neither
identifies such a four-score contrast with one edit at its original
anchor nor proves every desired anchor lies on a strand.

The right next experiment is to impose exact N as a **second** constraint,
and vary the contact cut or synchronized core, not to fix one exterior.
After choosing d each strand has a full address list, and the remaining
equation is a monotone quadratic plus an explicit bit-reversal remainder.
