# A safety–entropy tradeoff for the five-core arithmetic compiler

Session253, 11 November 2026. **NEW PROVED-IN-NOTE bridge and FINITE-EXACT restricted-language probes.** Chan–Pak, LRC(14), and Mahler 3/2 remain OPEN. No canon theorem ID or literature-firstness claim.

## Useful result

A natural transfer from Mahler to Chan–Pak fails for a precise reason. Mahler-safe binary words give a useful linear height bound for the flank statistic g. However, imposing this safety on BOTH flanks of ANY fixed finite core bank makes its positive integer output set have density zero. The same failure already holds for the weaker Fibonacci language forbidding11. This is not a statement that the five cores fail with unrestricted flanks.

The productive next probe was asymmetric: retain a Mahler-safe flank and an unrestricted flank. It succeeds through4095 but FAILS first at52477. The obstruction has only four canonical certificates, and reveals exactly where the Fibonacci/no11 relaxation loses Mahler safety. A still weaker one-no11-flank restriction survives the tested range through65535.

## Inheritance, scope, and recovered sidecars

- Closest proved mechanisms: the Session251 fixed-core compiler; Session252 seven-coordinate algebra; THM-4654 safe slack automaton and entropy; THM-3848 renewal count.
- Hostile: arbitrarily long words of ones have g=1 but unbounded address. Scalar g alone does not control address height.
- Corrected near miss: the symmetric returning clock 1^L0^L rediscovered here is only a specialization of Session197's 1^a0^k modular loop. Joint modular universality was already proved there and recovered in Session251/gems. No novelty is assigned to it.
- Least-used sidecar: the growth exponent of the accepted address language, combined with the positive bilinear AB term, not just the number of automaton states.
- Current concurrent input: `chanpak_five_core_bridges_20261110/cover/SELECTOR.md` proves a polylogarithmic output bound for rational asynchronous binary witness relations. Our result has a DIFFERENT hypothesis: the relation can be completely arbitrary, but both input languages must be thin and coercive. In particular it applies to the nonregular strict finite Mahler-safe language, without assuming any finite-state selector.

Concept board: safe carry language / Fibonacci language / integer height / bilinear positive compiler / density of output support / asymmetric flank freedom.

## 1. A sharp elementary coercivity lemma

Use the inherited scalar recurrence

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod2).

Let F be the nonnegative integers whose ordinary MSB-first binary word has no adjacent ones; the word of0 is taken empty. Then

    3g(n) >= n                       (n in F),
    3g(n) >= 2n+1                    (n in F even).       (1)

Proof by simultaneous induction. The claim holds at0. For n=2q, if q is even, the even induction bound gives

    3g(n) >= 2q+1+3q >= 4q+1.

If q is odd, the odd bound and the parity bonus give

    3g(n) >= q+3q+3 >= 4q+1.

For odd n=2q+1 in F, q is even (including q=0), so

    3g(n)=3g(q)>=2q+1=n.

All prefixes used in the induction remain in F. This proves(1).

The coefficient1/3 in the unrestricted no11 bound cannot improve: set n_k to the binary word 1(01)^k. Then n_k=(4^(k+1)-1)/3 and 3g(n_k)-n_k=2k+2, by the two-digit recurrence. Hence g(n_k)/n_k tends1/3. This sequence is NOT asserted Mahler-safe for large k.

## 2. The target itself bounds the product of restricted addresses

The concurrent positive compiler, inherited from THM-016 and the Session251 compiler, is

    H_y(A,B)=h+d(g(A)-1+g(B)-1)+L AB
              +iA+oB+s(A mod2)+e(B mod2),              (2)

where h,d,L are positive integers and i,o,s,e are nonnegative integers. This holds for every nonempty finite tournament core and cut. The proof below needs only these signs, not a five-core classification.

Since g>=1 and h,d,L>=1,

    H_y(A,B)>=g(A)+g(B)-1+AB.

If both A,B are in F, (1) gives

    H_y(A,B)>=(A+B)/3+AB-1,
    (A+1)(B+1)<=3H_y(A,B)+4.                           (3)

The universal constant in(3) does not depend on the core. In contrast, the unrestricted axis A=0,B=2^k-1 can have constant g(B), so no analogous coercive statement follows without an input condition or address normalization.

## 3. Thin-language hyperbola theorem

Let S,T be subsets of F, containing0 if desired, with

    #{a in S:a+1<=Y}=O(Y^alpha),
    #{b in T:b+1<=Y}=O(Y^beta),

for positive alpha,beta. For a FIXED FINITE bank of core/cut rows, the number of distinct values H_y(A,B)<=X with A in S,B in T is

    O(X^max(alpha,beta))                  if alpha != beta,
    O(X^alpha log(X+2))                   if alpha = beta. (4)

Proof. Bound distinct outputs by input pairs and the finite bank size. Put Y=3X+4. In the dyadic range 2^j<=A+1<2^(j+1), (3) forces B+1<=Y/2^j. Its number of pairs is at most a constant times

    2^(alpha j) (Y/2^j)^beta.

Sum j from0 to floor(log2 Y). A geometric sum gives the unequal-exponent case; the equal-exponent case has O(log Y) equal-size terms. This argument retains the actual positive target height. It is not just modular reachability.

### Fibonacci consequence

The number of no11 words of length k is F_(k+2)=O(phi^k), with phi=(1+sqrt5)/2. Therefore S=T=F has exponent

    alpha_F=log_2(phi)=0.69424...<1.

Every fixed finite core bank restricted to no11 on both flanks has only

    O(X^alpha_F log(X+2))

outputs<=X. It cannot cover all sufficiently large odd integers, even after removing finitely many exceptions.

### Mahler consequence

Let S_M consist of integers whose finite MSB-first binary word c satisfies every strict suffix inequality

    sum_(j>=0) c_(i+j)(2/3)^(j+1)<1,

where the sum stops at the end of the word. By THM-4654, S_M is contained in F. THM-3848/4654 give the inherited word estimate a_k=O((3/2)^k). Hence

    alpha_M=log_2(3/2)=0.58496...,
    #{H_y(A,B)<=X:A,B in S_M}=O(X^alpha_M log(X+2)).      (5)

The exponent is sharper, but the key obstruction already comes from no11. Finite safety here is a word predicate; no infinite Mahler Z-number or ordinary-integer carry stabilization is assumed. We feed the words to the Chan–Pak binary compiler as addresses, not as a claimed conjugacy of its dynamics with Mahler.

## 4. Exact new controls and the asymmetric survivor

`probe.py` imports the inherited exact scalar compiler and all cuts of the literal five cores. It tests ONLY new word restrictions, not an enlarged inherited spectrum. There are82 distinct signatures. It checks(1) for all no11 addresses<=65535 and computes the restricted support at255,1023,4095.

At X=4095, the number of distinct odd outputs is:

- unrestricted:2046 (only7,21 absent);
- at least one no11 flank:2046;
- at least one Mahler-safe flank:2046;
- both no11 flanks:1957;
- both Mahler-safe flanks:1871.

Both two-sided restrictions already miss47 and63. The full first-hole lists are in `results.json`. These are finite support statements, not evidence that the density-zero exponent is optimal.

The exact inverse enumerator normalizes a degenerate g-axis by deleting trailing ones until the remaining address is even or zero. This preserves prefix-closed no11 and finite-safe predicates; hence it does not discard the possible restricted certificates. Other branches are bounded directly by the positive coefficients, as in the inherited finite inverse.

### A new minimal asymmetric hostile,52477

The apparently promising one-safe-flank conjecture above was tested immediately, rather than left as a recommendation. `asymmetric.py` performs the exact new restricted inverse through65535. The unrestricted bank has32766 outputs; the one-safe-flank restriction has32765. Its UNIQUE extra missing target is52477. Requiring only one no11 flank loses no target in this range.

`hostile.py` exhausts the exact inverse at52477. There are exactly four canonical certificates:

- S5a, cut6, (A,B)=(11,55), order16;
- S5a, cut19, (A,B)=(55,11), order16;
- S4, cut10, (A,B)=(38,42), order17;
- S4, cut10, (A,B)=(42,38), order17.

For S5a, the two words1011 and110111 both contain11, whose Mahler weight is10/9. For S4,38=100110 contains11;42=101010 avoids11 but contains10101, whose Mahler weight is266/243>1. Thus the first obstruction is not captured by the Fibonacci adjacency rule: the next genuinely Mahler forbidden pattern matters. Independent subset Hamilton-path DP on the literal16-vertex S5a construction gives52477.

The inverse list is complete even allowing arbitrarily large address words. Away from the degenerate axis positivity supplies the inherited finite bound. An axis witness can be reduced by removing trailing ones while preserving finite safety, and appears in the canonical inverse. Since no canonical safe witness exists, no omitted unbounded safe axis witness can exist.

The safe-side minimum among all other targets is at most5. More precisely the first extra holes when the safe side is restricted to values<=K are:

- K=0:2893;
- K=1:14681;
- K=2:19813;
- K=4:30769;
- K=5:52477, the sole remaining extra hole through65535.

Thus, on this range, all safe words beyond the tiny set{0,1,2,4,5} add no support. This is finite evidence only. It does NOT prove that a bounded flank bank has density zero, positive density, or eventual coverage. With one flank fixed, the unrestricted other address still has full binary entropy and the two-thin-languages theorem does not apply.

The strongest survivor is an **OPEN one-no11-flank sufficient conjecture** for the five-core bank. The first52477 hostile instead pinpoints a useful positive surgery target: transport the S4 witness with one10101 safety defect to a bounded new core or modified cut, retaining the actual count. A proof allowing only one finite-safe flank is already refuted for the original five bank; adding a new core would be a genuinely different statement.

## 5. Exact bridge contract and LRC boundary

Source: finite Mahler-safe carry word, or its coarser no11 language.
Map: interpret the word as the MSB-first integer address for one Chan–Pak flank, and evaluate the physical compiler.
Preserved: literal word, word length, positivity of the resulting construction, and the proved height estimate.
Not preserved: Mahler carry evolution, oriented real phase, or ordinary-integer terminality.
Needed sidecar: second flank's actual word, core/cut, and exact H.
Decisive test: the hyperbola output count, not residue coverage.

For LRC, THM-3848 already solves the mixed-power positive control at loneliness2/5 but loses oriented Mahler safety. The present result does not repair LRC entry, arbitrary rungs, or synchronization. The useful cross-lane lesson is now quantitative rather than rhetorical: a restrictive safety language can make witnesses height-controlled and still make universal target coverage impossible. LRC finite-clock proposals similarly require their own target-preserving entry argument; this counting theorem must NOT be exported to LRC without a map giving(2) and(3).

Reproduce from this directory: `python3 probe.py`. Its compiler import searches repository parents for the pinned Session251 `ternary` directory. No Git mutation was performed by this lane. `asymmetric.py` and `hostile.py` reproduce the new restricted-range and four-certificate tests; their JSON outputs are saved beside this note.
