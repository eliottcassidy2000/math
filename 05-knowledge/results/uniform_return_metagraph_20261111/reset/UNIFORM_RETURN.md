# A uniform seven-retained-vertex analogue of the141→425 reset

**PROVED in this note, with exact literal controls.** The new operation maps
every member of the obstructed C4 family to its3H+2 target. It returns first
to an older ready family at a smaller count, then lands in a new neutral
defect24 family. Recursive ternary closure and Chan–Pak remain OPEN.

## 1. Inheritance and the changed operation

The source is T_a(k), k>=1, from
[the preceding C4/source note](../../chanpak_five_core_bridges_20261110/alternatives/NOTE.md):
its residual chain has word11^k,01,11,10; its ordered ports are a→b;
both ports beat the root r, which beats the chain. Write D=2^k. Then

    H(T_a(k))=34D+5,       s_r=q_r=9,       F_r-H=4.

The old141 reset deleted port b and used two insertions. Its uniform
formula31D+17 equals34D+5 only at D4. The present operation deletes the
MIDDLE tail vertex11 instead. It preserves the entire original11 prefix,
both ports, the root, and the two outside tail vertices01,10.

This deletion gives exactly the old ready-pentagon family BA_k, with

    H(BA_k)=14D+3,         s_r=q_r=5,       F_r-H=2.

This is the actual deletion map already observed in the preceding NOTE,
using the complete older [ready classification](../../two_entry_funnels_20261103/READY.md).
The new continuation differs from the standard two-birth ready reset:

1. Adjoin a global sink x.
2. Adjoin a second global sink y, so every old vertex beats x→y.
3. Adjoin z beating exactly x and the old01 tail vertex; every other
   vertex, including y, beats z.

The final three vertices contain the cycle x→y→z→x. The full count trace is

    34D+5 → 14D+3 → 14D+3 → 14D+3 → 102D+17
                                              = 3(34D+5)+2.          (1)

At k2 this is141→59→59→59→425. Seven of the eight original vertices and
every arc among them survive. At general k, k+5 of the k+6 old vertices
survive. The original root and BOTH ordered ports survive as well.
The two neutral intermediate sources are weak; the returned BA source
and the final target are strong. Their weak status is part of the program.

For T_b(k), first use the preceding count-neutral middle-tail/port edge
switch to T_a(k), then perform the same operation. Alternatively delete
one of its two01 tails to obtain the identical BA_k source.

## 2. Exact fixed core and the all-k proof

A canonical final graph at k1 has outgoing rows

    U=(4,1,10,51,71,151,175,287,127),

with root7, ordered ports a=6,b=5, and free prefix vertex8. Replace vertex8
by a transitive k-block, every copy having its old outside comparisons.
The six old BA1 vertices in order(prefix,B,A,a,b,r) map to

    (8,3,4,6,5,7).

The newborns(x,y,z) map to(1,0,2). This proves all preserved-incidence
claims for every k, without extrapolating finite counts. In append order
after deleting the middle tail, the outmasks are0,0,33·2^k.

Let W=U[{0,...,6}] and J={5,6}. Exact seven-vertex path-bin counts give

    H(W)=41,        E_J,W(2)=102,        C_J,W(2)=24.                 (2)

Here E counts ordered two-bin spanning path covers, the first nonempty
and ending in J; the second may be empty. C counts two nonempty bins,
the first ending at one port and the second starting at the other.
The implementation checks (2) directly from subset endpoint tables.

Prefix W by TT_k, and add r losing exactly to J and beating TT_k and
W\J. Deleting r from a Hamilton path gives the usual two-bin cover;
prefix vertices have two independent bins in E and are forced before
the terminal port in C. Thus the inherited polarized-cover argument gives

    H = 2^k E_J,W(2)+H(W)-C_J,W(2)=102D+17,
    s_r=q_r=H(W)=41,      F_r-H=C_J,W(2)=24.                        (3)

No transitivity of W is needed for this argument. Independently, the
positive free-module response is A1=221,A2=102,A_j=0 for j>=3; since
2!S(k,2)=2^k−2, it gives exactly the same formula.

Strongness follows directly: r reaches the prefix, the prefix reaches
every W vertex, and every W vertex reaches a port. The latter fixed
seven-vertex condition is explicit in U. Hence the finite graph check
and the structural argument agree at every k.

## 3. The new family belongs to the five-core hierarchy

Remove r and its two in-ports. Besides the transitive prefix, the residual
five-vertex rows are

    (4,1,10,19,7),

with H9 and isomorphism(1,0,2,3,4) to S5a=(2,4,9,19,7).
Thus the new neutral defect24 family has exact residual TT_k→S5a.
It is the next kind of structure the bounded-response atlas must retain:
the residual has become a nontransitive member of the original five-core
bank. Its gap is constant24, not a growing experimental statistic.

The inherited preparation operation that prepends11 preserves this
family and changes its count by H_next=2H−17. Therefore its order series
as distinct rooted classes is x^9/(1−x), and every root has the full
known profile(102D+17,41,41,102D+41). Distinct k have different orders;
no unrooted multiplicity or automorphism assertion is inferred here.

There is also a different, simple multiplicative realization of (1):
delete the first01 and middle11 tails of T_a(k), leaving the ready A_k
family of count6D+1, and take its ordinal join with any fixed H17 graph.
The product is17(6D+1)=102D+17. This explains a fixed divisor of the
target family. It does not replace the stronger one-deletion construction,
and it gives no information about the coprime ray102D+19.

## 4. What has and has not returned

The source-return is a genuine all-k program with a temporary count
decrease. It does NOT restore34D+5 before taking the branch; the older
single-module obstruction explains why that stronger architecture fails.
The source with H425 here is a ten-vertex graph, distinct from the preceding
eleven-vertex425 witness. Neither retains the deleted middle vertex.

The final root is neutral24-ready, not1/2-ready. This first reusable
growth family does not yet supply a3H+2 or3H+4 return from each of its
own outputs. Applying the standard ready reset to its retained root
would give3H+48. Changing only the final cut in the specified two-sink
schedule exhausts256 base cuts: the uniform3H+2 cut is unique, and no
uniform3H+4 cut exists in that schedule. The broader+4 and recursive
return questions remain open.

Reproduce with `python -B uniform_return.py`. The receipt retains actual
rows, all old-vertex maps, full traces, fixed-core counts, root responses,
and literal checks for k1..6. `classify_modules.py` and `scan9.py` record
the independent finite classifications that led to this construction.
