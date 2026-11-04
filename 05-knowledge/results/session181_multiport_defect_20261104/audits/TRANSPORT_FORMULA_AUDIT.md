# Independent transport formula review

2026-11-04. Formulas(1)–(6) in `../transport/PROOF.md` are mathematically correct. The deleted-root splitting proof, squarefree birth update, dominant-prefix recurrence and multiport affine criterion were checked independently. The actual candidate arithmetic in Section6 also checks:18583−3690=14893, the offset is8308−2076·2^s, and the next affine constant targets are44681 and44683.

One proof-convention repair was requested: in Section3, `A−S_O(U)=S` assumes U is nonempty. For singleton B, the formula remains correct, since B_t is transitive and h=f=A=1, S=0. The proof should explicitly handle that case first, then assume U nonempty. No mathematical formula needs changing.

The author's full direct-enumeration verifier was replayed and its output exactly matches the frozen `verify.out`. A separate implementation, `multiport_formula_audit.py`, independently checks formula(5) on all3,069 combinations of a labelled suffix of order1 through4, a nonempty physical port set and prefix length1 through3. Suffix F and J are calculated from permutations with one cut, while the full target count and root two-slot response are calculated by backward endpoint recursion. The test intentionally includes nonstrong suffixes, verifying the algebraic identity beyond the strong-construction scope.

These exact evaluator laws and finite controls do not prove any unbounded two-digit realization or native guard legality.
