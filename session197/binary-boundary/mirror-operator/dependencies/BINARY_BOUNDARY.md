# The physical zero-clone boundary and a quadratic self-converse family

Session 197. **PROVED exact physical deletion identity.** The positive-length rank-one binary mirror formula does not automatically survive deletion of the whole marked clone module. Its exact defect is the square of a physical endpoint count. A literal single-edge flip proves this directly. The construction also gives strong self-converse tournaments with one Hamilton cycle and counts 2(2^k+1)^2-1, including 49, 161 and 2177=7*311.

## 1. The graphs and four side responses

Let S be a nonempty tournament of order t with a marked vertex a. Put

    h=H(S), e=H_start_at_a(S), f=H_end_at_a(S),
    q=sum_(nonempty directed paths P ending at a) H(S-P),

with H(empty)=1. Every path is an ordered vertex word.

Take a copy O of S and a copy I of its converse, with matched labels. Orient all cross arcs I->O except the single arc O_a->I_a. Between them put a transitive module Z of m vertices, with I->Z->O. Denote this physical tournament by F_m, allowing m=0: F_0 is simply the graph on I and O after the one cross-edge reversal. Swapping matched I/O labels and reversing the order of Z is an explicit anti-involution. Thus every F_m is self-converse, including F_0.

The independently reviewed positive-length rank-one identity is

    H(F_m)=h^2+q^2*2^m,  m>=1.                         (1)

The result proved here is

    H(F_0)=h^2+q^2-e^2.                               (2)

Equivalently, for all m>=0,

    H(F_m)=h^2+q^2*2^m-e^2*indicator(m=0).              (3)

The indicator is explicit; there is no ambiguous 0^0 convention. The ordinary generating function is h^2/(1-z)+q^2/(1-2z)-e^2. In particular deletion from m=1 lowers the path count by exactly q^2+e^2.

## 2. A direct single-edge flip proof

Before reversing I_a->O_a, the graph is the ordered join I->O and has h^2 Hamilton paths. A path lost by the reversal must use I_a->O_a as its unique crossing from I to O. It consists of a Hamilton path in I ending at I_a and one in O starting at O_a. Converse identifies the first endpoint count with e, so exactly e^2 paths are lost.

Every gained path uses the new arc O_a->I_a. It has the unique form

    I_prefix, O_prefix, O_a, I_a, I_suffix, O_suffix,

where the full O segment before the reversed edge is a nonempty directed path ending at O_a, and O_suffix is a Hamilton path on its unused O vertices. I_prefix may be empty. The full I segment after the edge is a nonempty path starting at I_a, and I_prefix is a Hamilton path on its unused I vertices. O_suffix may be empty. Converse turns the latter choices into the same q arrangements as on O.

Any two choices give a valid gained path: all unmarked cross joins point I->O, and the only missing such arc I_a->O_a cannot be needed at another join because both marked endpoints already occur at the central reversed edge. Conversely the positions relative to that unique O-to-I edge recover both choices. There are exactly q^2 gained paths. Subtracting e^2 and adding q^2 to h^2 proves (2), also when t=1 and both unused complements are empty. This proof does not extrapolate (1) to zero or invoke an odd-parity filter at a nonexistent two-cycle.

## 3. The alternating response and elementary consequences

The independently proved endpoint reciprocity identity is

    d=sum_(P ending at a) (-1)^|P| H(S-P)=-e.

It explains why an odd-cycle parity calculation of the deletion defect gives -d^2: the defect is exactly -e^2. This is an interpretation of a now independently established physical formula, not a positivity assumption.

If q_even and q_odd denote the same positive response restricted to even and odd path lengths, reciprocity gives

    q=q_even+q_odd, e=q_odd-q_even,
    H(F_0)=h^2+4*q_even*q_odd.

Consequently q>=e, q and e have the same parity, H(F_0)>=h^2, and H(F_0)=1 modulo4. Equality H(F_0)=h^2 holds exactly when a is a source vertex of S. Indeed an in-neighbour supplies a two-vertex path ending at a with positive complement count, whereas a source has no nontrivial path ending there. These assertions include weak side tournaments and singleton sides.

For a transitive side with j predecessors of a, h=1, q=2^j and e=indicator(j=0). The exact boundary is therefore 1+4^j-indicator(j=0). The simplest failure of naive positive-length extrapolation is S=K1: all positive lengths give 1+2^m, while F_0 is a two-vertex transitive graph with H=1, not 2.

## 4. Hamilton cycles at the boundary

For t>=2, every Hamilton cycle of F_0 must use the sole O-to-I arc O_a->I_a. Anchoring the cycle at that arc leaves exactly one Hamilton path through I starting at I_a and one through O ending at O_a. Both have f choices, by converse. Their other cross join is valid because at t>=2 its two endpoints cannot be I_a,O_a. Thus

    C(F_0)=f^2, t>=2.                                (4)

For every t>=1 and m>=1, the clone module occupies one increasing run in a Hamilton cycle: two runs would require two disjoint O-to-I transitions, but only the marked edge exists. The complementary old path traverses all O vertices ending at O_a and then all I vertices starting at I_a. Hence C(F_m)=f^2 for every positive m. The sole zero-length exception is t=1, where F_0 has only two vertices and C=0, while f=1.

In particular whenever f>0 and t>=2, the boundary graph has a Hamilton cycle and is strong. This supplies strongness directly, without relying on deletion preserving strongness in general.

## 5. An explicit quadratic C=1 family

For k>=1 let S_k have vertices 0,...,k+1, every arc pointing forward except k+1->0, and mark a=0. Put t=k+2.

There is one Hamilton path avoiding the backward edge. A path using it is determined by the arbitrary subset of the k interior vertices placed, increasingly, before k+1; the rest lie increasingly after 0. Thus h=2^k+1. A Hamilton path starting at 0 or ending there is unique, so e=f=1. A path ending at 0 is either the singleton path or an arbitrary increasing subset of interior vertices followed by k+1,0. Its complement is transitive. Therefore q=2^k+1=h.

The resulting boundary graph has order 2k+4 and

    H(F_0)=2(2^k+1)^2-1, C(F_0)=1.                    (5)

It is strong by (4) and self-converse by the explicit paired anti-involution. A direct spanning cycle is O_0,I_0,I_(k+1),I_k,...,I_1,O_1,...,O_(k+1),O_0.

The first six counts are

    17,49,161,577,2177,8449.

At k=2,3,5 these are respectively 7*7, 7*23, 7*311, with prime quotients. The verifier supplies the literal 14-vertex k=5 adjacency rows and exact subset-DP count, not merely an identity-derived large number. The general divisibility condition is

    7 divides 2(2^k+1)^2-1 iff k=0 or2 modulo3.

This is not a primality condition: at k=6 the output is 8449=7*17*71.

For the same fixed side S_k, every positive clone length has H=(2^k+1)^2(1+2^m). None is prime or 7p: it has a nontrivial odd square factor and is at least three times that square, whereas a number 7p has no such square factor unless p=7, in which case its value is only49. Thus allowing the zero-clone boundary genuinely changes the arithmetic output class, despite preserving C=1 for these t>=3 graphs. This does not claim a returning inverse grammar or full prime/7p coverage. The family is sparse.

## 6. Sources, attribution and verification

The positive-length identity (1) is the already independently reviewed `../binary-rigidity/PROOF.md`, SHA256 `65102468958b198f8538022288238c40d4e3ca33684cf0d419bc86b0993681a9`. The endpoint reciprocity proof is `session197-odd-cycle-path-polynomials/ENDPOINT_RECIPROCITY.md`, SHA256 `4761f31465aabbe5e079c73131cabbece27dabca299a5900e5065baeb221a74c`, independently reviewed in `session197-bridges/endpoint-reciprocity-audit/REPORT.md`. The single-edge flip proof was independently identified during that comparison. No literature-priority claim is made, and the known broad C=1 structural language is not claimed to be newly classified here.

Compile `verify.cpp` with a standard C++17 compiler and run it. The standalone checker counts the boundary graph and paths using its reversed edge separately for every rooted labelled side through order five and 192 deterministic order-six rooted controls. It checks the exact gain/loss partition, endpoint reciprocity, cycle-square law, equality and parity boundaries, positive-length contrasts through side order four, and six literal members of (5) through order16. The k=5 witness and the prime-quotient checks are included in the deterministic receipt. Optimized and unoptimized assertion-enabled runs agree. All-order conclusions come from the physical bijections, not the finite controls.
