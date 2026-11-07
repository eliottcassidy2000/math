# Independent audit: fixed-H5 register and persistent pairs

Verdict: PASS for the exact register spectrum, joint endpoint support, reset-created pair, and same-order pair stable under two-sided padding. The separate four-choice H5 menu remains exploratory and is not promoted by this audit.

Input REGISTER.md SHA-256:

    eeaa666670b9065cdfcd906f7f11574f67469851efc5aa04303d90ed5b69fb8c

## Analytic checks

The fifteen core-cut values match literal vertex-subset path/cover DP on C=(8,1,3,6). The sink update is the independently audited physical cover recurrence, not an inferred scalar operation on arbitrary Hamilton-path counts.

For alpha>=17, inverse descent is alpha/2 for even alpha and (alpha-5)/2 for odd alpha. These are smaller integers at least5. A predecessor8 can arise only at16 or21, both direct seed values; the low starting range5..16 is covered except8, with10 obtained from5. No seed or allowed predecessor produces8. This proves the all-length spectrum, not merely its finite computational check. Every forward step at least doubles alpha, so a realizing path from a seed>=5 has at most floor(log2(alpha/5)) append steps. The stated body-order upper bound follows.

After a nonempty word, the final e is0 or5. The last bit is forced by the parity of alpha and its predecessor is unique. Applying the exact alpha spectrum gives precisely the stated e=0 and e=5 sets. The e=1,2,3,4 values occur only at the fifteen initial core cuts, exactly as listed. Original e=5 at alpha27 is already included in its asserted infinite set.

The reset-created bodies have (h,alpha,tau,e)=(5,35,24,5) and(5,35,26,5). A common append word preserves equality of h,alpha,e, hence preserves the tau difference2. The arbitrary-module formula therefore gives count difference2H(U); transitive U has H(U)=1. This is a returning marked-body operation, not an induced +2 operation or a completeness proof.

The same-order bodies have the full seven states

    (5,4096,3807,3219,2984,5,5),
    (5,4096,3807,3219,2986,5,5).

The local precursor states agree in (h,q,alpha,gamma,tau,e)=(5,512,223,398,168,5) and differ only in s=0 versus2. One selected prepend resets s to5 and adds the old s difference to tau; subsequent common prepends preserve this difference. All other coordinates agree. Common transitive padding at either end preserves the remaining coordinates and tau difference2 by the exact prepend/append updates. General common ordinal factors multiply this difference by their Hamilton counts under the audited composition law.

The displayed original pairs have strong closures. Common suffix appends preserve this. Arbitrary common prefixes still preserve the arithmetic pair, but need the usual initial-cut condition if strongness is additionally desired. The proof does not claim every unselected-prefix extension is strong.

## Independent controls

The checker uses the previous independent three-core audit's literal vertex-subset endpoint DP, copied locally for portable reproduction. It recomputes all initial states rather than trusting a state table.

- Exact support and full joint(alpha,e) support checked through alpha100000 against all claimed membership classes.
- The reset-created closed graphs counted literally:51 and49.
- Both eleven-vertex bodies' full states recomputed by subset DP.
- Both thirteen-vertex closures counted literally:4635 and4633; reachability confirms strongness.
- Both local precursor states recomputed literally, including s=0 and2.
- All49 common two-sided transitive continuation choices of total length<=3 preserve the claimed full-state relation. The all-length conclusion follows from the algebra, not these controls.

The finite register test is a hostile control on the exact proof, not a substitute for induction. Omitting the endpoint coordinate would identify local precursor states that behave differently under a selected prepend; the explicit0-versus2 witness prevents that loss.

## Reproduction and domain

    g++ -O2 -std=c++17 audit.cpp -o /tmp/register-audit
    /tmp/register-audit

The program's fixed bounded checks fit signed64-bit integers. The unbounded theorem is mathematical, not a claim that this finite C++ program can evaluate arbitrary-size states without overflow. The separate numeration executable has its own signed64-intermediate qualification. No large-target formula receipts are used to certify these small literal pair counts.

No source edits, stage edits, or publishing were performed. The density and three-core portable releases were not changed.
