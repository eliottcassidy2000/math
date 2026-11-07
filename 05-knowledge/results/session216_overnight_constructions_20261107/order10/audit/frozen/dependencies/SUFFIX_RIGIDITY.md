# 3-adic suffix rigidity gives a complete finite search for transitive-tail universal pairs

Status: PROVED candidate; finite bank result independently replayable, audit pending.
This is an all-length theorem for a specified construction family, not an obstruction to arbitrary source choice or J-containing extensions. Chan–Pak completeness remains OPEN.

## 1. The suffix theorem

Use nontrivial rooted source states p=(h,s,e,L,R,d,C,b,T). Compare two sources with the same positive b. For any two O/Z words u,v of the same length k, suppose their resulting T coordinates agree. Put a=v_3(b). Then u,v have a common suffix of length at least max(0,k-a).

Proof. Encode O by0 and Z by1, chronologically z_0,...,z_(k-1). The inherited transfer b'=2b, T'=3T-zb gives exactly

T_k=3^k T_0-b P_k(z),
P_k(z)=sum_(i=0)^(k-1) z_i 2^i 3^(k-1-i).

Write b=3^a b_* with 3 not dividing b_*. If k>=a, equality implies P_k(u)-P_k(v) divisible by3^(k-a). Reduction modulo3 determines the final bit, because its coefficient2^(k-1) is a unit and the bit difference is in {-1,0,1}. After this difference vanishes, divide by3 and repeat. The last k-a bit differences vanish. For k<a the claim is empty. No positivity assumptions on T_0 beyond its integrality are needed. QED.

## 2. Full-profile collision cutoff

Suppose the two sources have identical initial b, and after arbitrary O/Z words their eight NONSTART coordinates agree. Their word lengths are equal, since the final b is the initial b times2^length. There is then a collision of their full eight coordinates after prefixes of the same length at most a+1. More precisely, any collision of length k>a+1 already holds after the first a+1 letters of the actual two words. All remaining letters form a common suffix.

Proof. Section1 gives a common suffix starting no later than position a. After the first letter of this common suffix, both endpoint coordinates e agree: O sets e=0, and Z sets e=d. Their d coordinates agree, since d is invariant and final full-profile equality includes d. Thus at prefix length a+1 the two e coordinates agree, whether or not other coordinates yet agree.

Both O and Z are injective on pairs of nonstart states having equal e before the operation. For O this follows from h,L,d unchanged and R,C,b,T multiplied by2,2,2,3. For Z recover d, b=b'/2, R=(R'+d)/2, T=(T'+b)/3, L=L'-b, C=(C'-2T+2L)/2, and h=h'-R+e. The input e must be supplied; equality of it is exactly the needed sidecar. Therefore peel the common suffix backward from final equality down to prefix length a+1. QED.

Conversely, any eight-coordinate collision persists under every common O/Z/J suffix, because none of those eight-coordinate transfer formulas involves s. If s_1=s_2+2, with their residue and width constraints, this is precisely the universal-pair criterion from Session203. Thus an unbounded search for universal pairs by asynchronous O/Z padding of two sources with common b reduces to a finite prefix search. A cutoff test at exactly a+1 suffices: append identical O moves to any earlier collision to reach that depth.

The sidecar e cannot simply be dropped. O erases e; Z erases it while retaining h-e. The extra +1 in the bound is needed for this loss and is not silently omitted.

## 3. Fixed deletion-core consequence

Fix any nonempty tournament Q. Adjoin a new designated root with arbitrary incidences to Q, and then adjoin an arbitrary transitive tail. All old nonroot vertices beat tail vertices; the root-to-tail incidences may be arbitrary and independently chosen on the two sides. Every such construction is exactly a root-cut source followed by an O/Z word. All initial sources have b=B_Q(2), so the complete test for a universal pair is: enumerate root cuts, then all O/Z words of length v_3(B_Q(2))+1, and match the eight nonstart coordinates plus the label invariant s+2(s mod3). Retain positive starts at most26 with nonzero residue. The label invariant prevents losing candidates by overwriting different starts in the same residue class.

This criterion is complete for TWO CUTS OF THE SAME FIXED Q, allowing arbitrary independent transitive tails. It does not classify different deletion cores, arbitrary J words, or changes of root. It says nothing about finite-word or O/Z-only carrier collisions, which are weaker targets than full-eight equality.

## 4. Exact finite bank certificate

The unchanged supplied bank BANK.txt contains7412 deletion tournaments through order8; we make no claim that it is class-complete. The C++ root compiler evaluates all1,823,702 root incidences. There are63,546 retained bounded-start sources. The cutoff test evaluates219,976 leaves, with maximum cutoff7, and finds ZERO universal pairs. It stores208,504 distinct full-profile-and-label keys. The global table compares more cases than required (including some different cores), but the all-length conclusion is asserted only separately within each fixed bank core.

Combining Sections1–3 with that finite exact computation: for every deletion tournament Q actually supplied in BANK.txt, no two root cuts with arbitrary independent transitive tails produce a universal O/Z/J paired source within the stated start budget. This is an all-tail theorem about this finite family of cores, not an unrestricted no-go theorem or a spectrum hole.

Reproduction:
- Compile root_orientation_tail_complete.cpp with g++ -O2 -std=c++17.
- Run the executable on BANK.txt; stdout equals empty ROOT_TAIL_COMPLETE.txt and stderr equals ROOT_TAIL_COMPLETE_COUNTS.txt.
- Run verify_suffix.py, which checks the integer digit rigidity and explicit inverse maps.
- Compile root_interface_control.cpp, run on BANK.txt, and compare ROOT_INTERFACE_CONTROLS.txt. Then run verify_root_interface.py. All212 emitted profiles are independently reconstructed as literal tournaments and counted with the inherited subset-DP engine.

## 5. The root-incidence compiler

For each subset A of Q, let H_A count its Hamilton paths, E_A(i) its paths ending at i, and S_A(i) its paths starting at i. H_empty=1; endpoint counts of empty sets are zero. Let x_i=1 mean root->i. For fixed Q define

EL_i=sum_A E_A(i) H_(Q\A),
SR_i=sum_A S_A(i) H_(Q\A),
ST_i=sum_A S_A(i) B_(Q\A)(2),
QH_ij=sum_(A disjoint B, A union B=Q) E_A(i) S_B(j),
QC_ij=sum_(A disjoint B) E_A(i) S_B(j) H_(Q\(A union B)),

where A,B are nonempty in QH,QC. Then

s=sum_i x_i S_Q(i), e=sum_i(1-x_i)E_Q(i),
L=H_Q+sum_i(1-x_i)EL_i,
R=H_Q+sum_i x_i SR_i,
T=B_Q(2)+sum_i x_i ST_i,
h=s+e+sum_ij(1-x_i)x_j QH_ij,
C=2[H_Q+sum_i(1-x_i)EL_i+sum_i x_i SR_i+sum_ij(1-x_i)x_j QC_ij],
d=H_Q, b=B_Q(2).

Proof: split the path containing the new root into its before-root and after-root portions. For h these partition Q; for C the remaining vertices form the other labeled bin. The root lies in either bin, giving factor2. For L,R,T the root occupies the endpoint prescribed by that coordinate, leaving one or two ordinary bins. These cases are disjoint and exhaustive. This supplies a direct positive integer compiler, not a fit to sampled profiles.
