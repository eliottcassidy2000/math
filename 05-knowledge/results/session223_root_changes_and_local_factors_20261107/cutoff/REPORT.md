# Two-root-change cutoff: exact digital condition and a certificate boundary

Status: PROVED elementary reductions below; FINITE-EXACT regression checks. The cutoff conjecture remains OPEN. No new census or infinite family of holes is claimed.

Sources: `root_change_carmichael_chamber_20261107/PROOF.md` and its audited correction file; Session220 `05-knowledge/results/session220_binary_envelopes_and_selectors_20261107/selectors/PROOF.md` for the earlier fixed-bank selector obstruction, including its asynchronous extension. This report uses the physical moves and source-prefix convention of the former. The closing J is not counted. The report is a limited research outcome, not a substitute for an all-order covering theorem.

## 1. Exact digital middle-block condition

Write an O/Z block B=b_1...b_n with b_i=1 for Z and 0 for O, and define

- L=2^n;
- q=sum_i b_i 2^(n-i), the ordinary binary integer;
- a=sum_i b_i 2^(i-1)=rev_n(q);
- c=sum_(j<i) b_i b_j 2^(i-1-j)+sum_(i=2)^n b_i b_(i-1);
- i=b_1 and e=b_n=q mod 2.

For n>=1 its exact action is

    B(D,E,U,d)=(D-iE+aU-cd, ed, LU-qd, d).                 (1)

Indeed the U recurrence is U_j=2U_(j-1)-b_j d. At a Z birth D increases by U_(j-1)-E_(j-1), and E_(j-1)=b_(j-1)d after the first birth. Summing proves (1), including the adjacent-digit term in c. Omitting that term is unsound: B=ZZ has c=2, not 1.

Fix arbitrary outer O/Z words A,C and let s1=t.AJ=(D1,E1,U1,d1). For nonempty C, compute its parameters (i_C,a_C,c_C,q_C,L_C,e_C) above and set

    h=a_C+L_C, r=c_C+q_C+e_C,
    alpha=1-i_C+h-r, beta=i_C-1, gamma=1+h.

Then

    F(s.JC)=alpha D+beta E+gamma U.                       (2)

For C empty the same formula holds with (alpha,beta,gamma)=(1,0,2), since F(s.J)=D+2U.

Combining (1) and (2), the exact equation for w=A J B J C with B nonempty is

    N=alpha D1-alpha i E1
      +(alpha rev_n(q)+gamma 2^n)U1
      -(alpha c(q)+gamma q-beta(q mod 2))d1.              (3)

Every integer q in [0,2^n) is permitted, including leading zeros; it specifies one physical block. Conversely every two-J source word has this form, with B empty handled directly by (2). Thus existence of A,C,n,q satisfying (3), together with the zero/one-J layers, is necessary AND sufficient for bounded-depth membership. All outer words remain quantified, not replaced by a finite bank.

The useful coordinate is c(q), not just q and its reversal. It records weighted pairs and adjacency. The target condition is not a simple bilinear factorization in independent free integers: q, rev_n(q), c(q), and the leading/trailing bits are coupled. Consequently neither a residue-cover computation nor treating their extrema as an interval proves coverage. No useful positive descent closing this coupled condition was found in this pass.

## 2. All mixed O/Z periods have unavoidable positive drift

Let P be any nonempty mixed O/Z word of length k, b=2^k. Let A be any fixed physical prefix (it may contain J), and C any fixed suffix. Put d_A=d(t.A)>0. There are rational constants X,Y and a strictly positive rational delta_P, and a positive integer z_C, such that for every m>=1

    F(t.A P^m C)=X b^m + z_C d_A delta_P m + Y.           (4)

Here z_C=F(e_D.C), where e_D=(1,0,0,0) is used only as a linear perturbation, not asserted to be physical.

Proof: across one period U maps to bU-q d. Its virtual fixed point is U*=q d/(b-1). For a mixed period every cyclic rotation of U*/d lies strictly between 0 and 1. E is periodic from the end of the first copy onward. D has eigenvalue 1 and increments affinely in U, hence its solution is an exponential plus a linear drift. The drift divided by d is the sum of virtual increments at the Z positions, namely U*/d minus the preceding Z-bit. Partition the cyclic word into Z-runs. A run of length l, entered at virtual U/d=u after an O, contributes

    u + sum_(j=1)^(l-1) 2^j(u-1)
      =1+(2^l-1)(u-1).

Its ending virtual U/d is 1+2^l(u-1)>0, so its contribution is greater than 2^-l. Thus delta_P>0. The U and E coordinates have no linear term, d is constant, and the D-coordinate linear term is d_A delta_P. Linearity makes the coefficient after C equal to z_C d_A delta_P.

To show z_C>0, use coordinates x=D-E, E, W=U-d, d. The moves are

    O:(x+E,0,2W+d,d),
    Z:(x+W,d,2W,d),
    J:(W+d,x,W+d,x+E).

Starting from e_D=(x,E,W,d)=(1,0,0,0), the cone of nonnegative coordinates with E<=d is preserved. A nonzero vector in this cone remains nonzero under each move. Its F=x+W+d is positive. This proves the assertion for every suffix, including arbitrary J budgets.

For example P=OZ has b=4 and delta=2/3. The linear term cannot be discarded as a lower-order error when exact equality is the target.

## 3. Consequence for finite-state certificate attempts

Fix K. Let R be ANY rational relation accepted by a finite-state two-tape transducer whose tape labels are finite words over {O,Z,J} and finite binary words, including empty words. Every accepted pair (w,v) must satisfy #J(w)<=K and value_binary(v)=F(t.w). Leading zero output padding is permitted. The transducer can select an arbitrary subset of valid witnesses; neither uniqueness nor recognition of all valid pairs is assumed.

**Theorem.** There is a constant R0, depending only on the transducer, such that every selected source word has at most R0 O/Z runs. Consequently the number of distinct selected outputs in [1,X] is O((log(X+2))^R1) for some constant R1. In particular such a literal-word finite-state certificate cannot cover a positive-density set of odd integers.

Proof: expand labels into single-symbol or epsilon edges. Any loop on an accepting path contains no J, since pumping would violate K. If its source label P contains both O and Z, pumping gives source words A P^m C. The binary output words have the form u v^m z. If v has length l>0, their numerical values have the form T 2^(lm)+S. If l=0 they are constant. Equation (4) has a nonzero linear term in m and is incompatible with either form, for any l. Formally exponential-polynomial sequences with distinct characteristic roots are linearly independent; alternatively divide by the larger exponential and then subtract the surviving exponential to expose the unbounded linear term. Thus no productive loop has mixed source labels.

In the trimmed graph each strongly connected component therefore contains only O-labeled edges, only Z-labeled edges, or only source-epsilon edges. Otherwise strong connectivity creates a mixed closed walk. Components are visited in DAG order; only boundedly many intercomponent edges occur. This bounds the number of source runs. J letters occur only on intercomponent edges.

The inherited height inequality F(t.w)>=2^(#O+#Z)+2 gives source length <=log2(X)+K. With boundedly many runs and at most K J letters there are only polynomially many in log X such words. Multiple binary encodings of the same integer do not increase the number of distinct outputs.

This extends Session220's already asynchronous fixed-bank selector obstruction to the literal bounded-root-change word grammar. Asynchrony itself is not new here. It does NOT exclude arithmetic constructions with unbounded carries, output-only automata, different witness encodings, or cofinite output support. It is specifically a warning against interpreting a finite-state word-output selector as the needed cutoff certificate.

## 4. Regression checks and limits

Run `python3 verify_drift.py` and `python3 verify_blocks.py`. The first checks all 494 mixed periods through length 8, four physical prefixes and five suffixes, using exact fractions. It checks 19,760 second-difference identities and positivity. The second independently checks digital block formulas against literal transitions, and the two-J scalar formula.

No conjectural cutoff, no new largest hole, and no positive asymptotic coverage claim follow from this report. The useful positive reformulation is (3). Its unresolved obligation is finding physical outer blocks and a compatible middle digit string uniformly for every sufficiently large odd N. Finite-state pumping cannot remove its unbounded drift coordinate.
