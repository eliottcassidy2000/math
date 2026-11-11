# A finite core bank still needs unbounded arithmetic memory

2026-11-10. **PROVED elementary extension of an inherited obstruction.**
This concerns a proposed witness-selection architecture. It does not refute
the five-core covering conjecture, nor Chan–Pak. No literature priority is
claimed. The one-input mechanism is inherited from
`session220_binary_envelopes_and_selectors_20261107/selectors/PROOF.md`.

## Statement and exact representation

Fix any finite collection of nonempty finite tournaments, with any finite
selection of apex cuts. Allow arbitrary transitive left and right flanks.
The inherited Session251 compiler assigns to each bank row y and two
nonnegative binary integers A,B an actual tournament count H_y(A,B).

A rational binary witness selector is a finite-state, possibly asynchronous
three-tape relation on ordinary binary representations of (A,B,N), contained
in the union of the equalities N=H_y(A,B). It may be partial, nondeterministic,
unlabelled by y, and use leading zeroes and epsilon transitions. It need not
select every witness or a unique witness. Every successful path must encode
a valid equality for at least one bank row. Binary word boundaries are part
of the representation, not values in an unbounded arithmetic register.

**Theorem.** Such a selector produces at most O((log(X+2))^R) distinct
positive outputs N<=X, for a constant R depending on the finite selector
and bank. In particular it cannot certify every sufficiently large odd N,
or any output set of positive lower density.

Finite-state recognition of output support alone is not excluded. Neither
are an algorithm with unbounded registers, a non-rational witness encoding,
or a target-dependent number-theoretic argument. The theorem is about the
relation retaining the actual two ordinary-binary flank addresses.

## 1. The positive compiler and target bounds

Write d=H(C). The inherited alternating-path identity THM-016,
`01-canon/theorems/THM-016-claim-b-alternating-sum.md`, gives k_-=-d and
converts the signed endpoint terms to nonnegative Hamilton endpoint counts.
Consequently the general compiler has the form

    H_y(A,B) = h + d(g(A)-1 + g(B)-1) + L AB
               + i A + o B + s(A mod2) + e(B mod2),       (1)

where h=H(C+x), d,L>0, and i,o,s,e>=0. The function g is inherited:

    g(0)=g(1)=1;
    g(2q+1)=g(q);
    g(2q)=g(q)+q+(q mod2).

For mixed contact both i,o are positive, because singleton core paths have
positive residual weights. Thus N<=X implies A,B=O(X). For uniform contact
one of i,o is zero; the other is k_+>0. Away from the degenerate axis the
term LAB again bounds both inputs by O(X). On the degenerate axis the
count is exactly d*g(v). Repeatedly remove trailing binary ones from v,
using g(2q+1)=g(q), until v is even (or zero). Its value is unchanged and
the reduced v=2q satisfies g(v)>=q+1, so v<=2X/d. This normalization never
increases the number of binary runs. These observations apply uniformly
over the finite bank.

The same bounds can instead be obtained from Session251's signed formula
and its finite-target inverse. No unproved coverage assumption is used.

## 2. A mixed loop has positive drift on either input tape

For a fixed prefix u, repeated nonempty binary word w of length k, and
fixed suffix z, the integer v_m=value(u w^m z) has the form

    v_m = a 2^(km) + b.

The inherited mixed-period lemma proves

    g(v_m) = c 2^(km) + D_w m + e,                      (2)

with D_w>0 if w contains both binary digits. If w is monochromatic,
D_w=0. For an empty loop word v_m is constant and we also put D_w=0.
The fixed suffix does not change D_w: appending a fixed word adds to g
an affine function of the old prefix and a fixed parity term, leaving
the coefficient of g exactly one.

For completeness, positivity is not merely numerical. At the virtual
fixed prefix -r/(2^k-1), where r=value(w), every cyclic prefix lies strictly
between -1 and 0 when w is mixed. A cyclic zero-run of length l starts
after a one and contributes 1-(2^l-1)q, where 0<q<2^-l. This is positive.
Sum over all zero-runs. All-zero and all-one loops have zero drift directly
from the recurrence. This is the Session220 argument, now used on both
flank tapes.

## 3. The bilinear term does not absorb this drift

Pump any loop on a successful path of the three-tape selector m times.
Each input A_m,B_m and output N_m is an exponential in m plus a constant,
or is constant if its loop tape word is empty. The product A_m B_m is a
sum of exponentials and a constant; it has no linear-in-m term.
The parities are fixed for m>=1. Therefore for each fixed bank row y,

    H_y(A_m,B_m)-N_m
       = sum_j c_j lambda_j^m + d(D_A+D_B)m + c_0,       (3)

where the distinct nonconstant bases lambda_j are positive powers of 2.
If either input loop word is mixed, d(D_A+D_B)>0.

Expression (3) can then vanish for only finitely many integer m. Indeed,
the largest base having a nonzero coefficient eventually dominates every
smaller exponential and the affine term. If no exponential coefficient
survives, the remaining affine function has positive slope. There are
only finitely many bank rows, so even allowing the successful bank row to
change with m cannot provide a valid witness for every pump count.
This contradicts the pumping property. Hence no successful loop is mixed
on either input tape. Unlike a growth-rate-only argument, this proof also
allows the two input loop lengths to differ and allows one input to stay
constant. The essential coefficient is the positive d multiplying both
g-values; no cancellation can remove the sum of their drifts.

## 4. Bounded runs imply sparse outputs

Trim the finite automaton to successful states. Each strongly connected
component has monochromatic non-epsilon labels on each input tape:
otherwise internal edges of both digits, connected by paths in the SCC,
would form a mixed closed walk. An accepting path visits each SCC at most
once in the condensation DAG. Hence both input integers have uniformly
bounded numbers of binary runs, with a bound depending only on the number
of states. Epsilon edges and leading zeroes cause no problem.

For each accepted N<=X, choose one valid bank row and apply the axis
normalization in Section1 if necessary. This produces a pair A',B'=O(X)
with no more runs than the original pair and the same N. The number of
binary integers at most O(X) with at most r runs is polynomial in log X:
choose their length, first digit, and at most r-1 run boundaries. The
number of normalized pairs is therefore also polynomial in log X.
Multiplying by the finite bank size bounds the number of distinct outputs,
proving the theorem. The normalization need not itself be performed by
the selector; it is a counting map on its valid witnesses.

## Consequence for the five objects

The five bank labels cannot serve as all the memory of a positive-density
constructor when its flank witnesses are encoded in ordinary binary.
Adding any fixed finite number of further core labels does not repair
this architecture. This does not say that five cores fail to cover the
target set. A proof of their coverage can retain an unbounded drift/carry
coordinate, use a richer selector, or prove existence without a finite-state
witness relation. The arithmetic register is precisely the information a
bare finite label forgets.

`check_selector.py` gives independent rational-arithmetic controls of the
drift calculation and of the two-input exponential annihilator. These
finite controls are regressions for the proof, not its infinite justification.
