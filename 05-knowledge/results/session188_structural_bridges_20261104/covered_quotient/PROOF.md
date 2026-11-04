# Unequal transitive blocks in a fully covered cyclic quotient

**PROVED:** Let T=C3[TT_a,TT_b,TT_c], where all three physical blocks are transitive and a,b,c>=2. Let U be obtained by reversing the quotient triangle while retaining each child's interior orientation and every old label. If the three sizes are not all equal, T and U have different full named two-port responses. An explicit battery uses singleton masks in one block. The proof gives exact inequality, not necessarily a residue-two separation.

This is the simultaneous reversal of the whole core and all three child modules, a fully covered laminar family. It does not settle arbitrary child tournaments, general quotients, or arbitrary fully covered families. Equal sizes are addressed by the separate equal-profile argument, not by a deletion shortcut here.

## 1. Exact start-count tomography

Attach x->y, with a_u=[u->x] and c_u=[u->y], and write F_T(a,c) for the full Hamilton count. Let H be the old Hamilton count, S_u and E_u count old Hamilton paths starting and ending at u, and C_uv count ordered spanning two-path covers whose first path ends at u and whose second starts at v. Both cover slots are nonempty. Then

    F_T(0,0)=H,
    F_T(e_u,0)=H+2E_u+2 sum_v C_uv,
    F_T(0,e_u)=H-S_u+E_u+sum_v C_uv.

For the second formula, y is either first or immediately follows x. In either position the nontrivial x-incidence is an old path ending at u, followed either by nothing or by a disjoint old path. The remaining case starts x,y and contributes H. In the third formula, x is a universal source, so deleting it leaves a one-port insertion. Relative to the universal-source y case, exactly S_u paths are lost and E_u+sum C_uv paths are gained. Hence

    2S_u = F_T(e_u,0)-2F_T(0,e_u)+F_T(0,0).       (1)

All terms are actual two-port queries. This exact extraction retains the complete response; no cubic or quartic term is discarded. Formula (1) is an elementary instance of the inherited endpoint tomography, not a new reconstruction claim.

## 2. A general cyclic-run identity

Let T=C3[A_0,A_1,A_2] with quotient 0->1->2->0. For any nonempty child A_i, let p_i(k) be its number of ordered spanning covers by k nonempty paths, and set p_i(0)=0 and p_i(k)=0 above its order. Let U=C3^op[A_0,A_1,A_2]. Write S_i(T) for the sum of the physical start counts over A_i. Indices below are cyclic, with j=i+1 and l=i+2. Then

    S_i(T)-S_i(U)
      = sum_(k>=1) p_i(k)
          [p_j(k)p_l(k-1)-p_j(k-1)p_l(k)].       (2)

To prove this, collapse each maximal child run in a Hamilton path. Because the quotient is a directed triangle, the resulting word is forced to cycle. For a path beginning in block i, its length is 3k, 3k-1 or 3k-2. The run-count vectors are respectively (k,k,k), (k,k,k-1), and (k,k-1,k-1), in order i,j,l. Each ordered cover of each child supplies its successive runs uniquely. Reversing the quotient exchanges j and l; only the middle case changes. This proves (2) for arbitrary child tournaments.

## 3. Monotone Stirling ratios

Let s(m,k) be the Stirling number of the second kind. Every nonempty vertex subset of a transitive tournament has exactly one Hamilton path, so

    p_m(k)=k! s(m,k).                              (3)

We need the elementary ratio property: whenever b>c>=1 and k>=1,

    p_b(k)p_c(k-1)-p_b(k-1)p_c(k) >=0.             (4)

For k=1 both sides of the minor vanish. For 2<=k<=c, divide by the positive denominators. The factorials cancel from the ratio p_b(k)/p_c(k), so it suffices that s(b,k)/s(c,k) increases with k.

Here is a self-contained proof. Every Stirling row is nonnegative, has consecutive positive support 1,...,m, and is log-concave. Induct using s(m+1,k)=k s(m,k)+s(m,k-1). If a is a nonnegative log-concave row without internal zeros and b_k=k a_k+a_(k-1), then

    b_k^2-b_(k-1)b_(k+1)
      = k^2(a_k^2-a_(k-1)a_(k+1))
        +a_(k-1)a_(k+1)
        +(k+1)(a_k a_(k-1)-a_(k-2)a_(k+1))
        +(a_(k-1)^2-a_(k-2)a_k) >=0.

The wider minor in the third term is nonnegative by the successive-ratio form of log-concavity. Extend rows by zeros at both ends; the boundary cases obey the same identity. The one-entry initial row establishes the induction.

Consequently

    s(m+1,k)/s(m,k)=k+s(m,k-1)/s(m,k)

is strictly increasing over 1<=k<=m. Multiplying these positive increasing adjacent-order ratios proves the asserted ratio property for b>c. At k=c+1 the minor in (4) is nonnegative because p_c(k)=0, and above that both products vanish. Thus (4) includes the support boundary. Finally its k=2 value is strictly positive:

    p_b(2)p_c(1)-p_b(1)p_c(2)=2^b-2^c >0.         (5)

## 4. The separating battery

Suppose a,b,c>=2 are not all equal. Choose a block i such that the other two block sizes differ. In (2), all summands have the same weak sign, by (4), and the k=2 summand is strictly signed because p_i(2)>0 and (5) holds. Therefore S_i(T)!=S_i(U).

Whole quotient converse preserves H: converse the whole tournament, then converse each child module to restore its interior, operations which preserve scalar Hamilton counts. Subtracting (1) between T and U and summing over physical u in block i gives

    2[S_i(T)-S_i(U)]
      = sum_(u in A_i) [F_T(e_u,0)-F_U(e_u,0)
                       -2F_T(0,e_u)+2F_U(0,e_u)].   (6)

The right side is nonzero. At least one of the 2|A_i| individual named queries therefore separates. This is a constructive battery, not a claim that all these queries differ or that their differences are nonzero modulo four.

## 5. Boundaries of the argument

Equal transitive sizes make every start-block difference in (2) zero; this method alone does not settle that case. A one-port response of a vertex deletion does not in general recover the deleted card's separate endpoint counts, so deleting one vertex and reapplying (2) is not a valid two-port proof.

Neither Stirling ratios nor their deletion analogue extend automatically to arbitrary child tournaments. For A with outgoing rows (8,1,3,6), deleting physical vertex 0 gives profiles

    q(A)=(5,22,36,24),  q(A-0)=(1,6,6),

at k=1,... . The k=2 minor q_2(A)q_1(A-0)-q_1(A)q_2(A-0)=22-30=-8. This exact small obstruction is retained in the checker.
