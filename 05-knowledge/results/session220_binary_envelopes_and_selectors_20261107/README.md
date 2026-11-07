# Session220: exact binary envelopes and a constructive-selector obstruction

**Checkpoint 1: PROVED scoped theorems with independent audits. Chan–Pak completeness, the three-core inverse conjecture, and eventual F2 coverage remain OPEN.** No unrestricted coverage extension or literature-priority claim.

This proof-focused session connects the two strongest inherited arithmetic construction carriers. It obtains exact optimization and a limitation on how a coverage proof can encode witnesses. Positive density and ribbon-construction investigations continue separately.

## 1. One binary function controls both carriers

Define f(0)=0, f(2u)=f(u), and f(2u+1)=f(u)+u+(u mod 2). For the inherited one-apex function g and e(v)=1−(v mod 2),

    g(v)+f(v)=v+e(v).

Consequently every fixed-bank join expression is exactly

    a g(v)−b e(v)+c(v+1)
      =(a+c)v−a f(v)−(a−b)(v mod 2)+(a−b+c).

This is the same quadratic carrier as a fixed suffix of Session216's three-core construction. The identity preserves the exact count and the binary witness; it does not identify the two physical families or transfer their coverage statements.

## 2. Exact maxima and minima without enumerating words

For K,h>0 and S>=0, the maximum of Kp−hf(p)−S(p mod 2) on 0<=p<2^m occurs at p=0 or p=2^m−2^j. An internal zero followed by a one can be rounded upward while weakly decreasing f and strictly decreasing binary reversal. The candidate list reduces further to a constant number near an explicit square-root balance.

This gives exact, parity-periodic square-root deficits, including all ties. For the full two-word three-core formula at fixed padding lengths, both maximizing words have at most one run of ones. Thus the maximum needs only polynomially many evaluations, instead of exponentially many binary words.

Every quadratic pair coefficient in the same formula is nonpositive. Therefore its minimum on **any fixed-bit cylinder** is an exact minimum-cut problem with polynomially many vertices and edges. The verifier includes capacities larger than 10^100; augmentation uses a proved capacity bound, not a numeric sentinel.

[Proof, sharp formulas, scope and reproduction](extrema/PROOF.md).

These are exact interval hulls, not occupancy certificates. A concrete legal rectangle has support {517,561,635,669,763,797,871,895}: its hull contains519 but the rectangle does not. This does not exclude519 from the full grammar.

## 3. Finite-state witness selectors are necessarily sparse

Fix any finite bank of a,c>0 and b>=0. A finite-state transducer selecting valid ordinary binary pairs

    (v,N),  N=a g(v)−b e(v)+c(v+1),

can produce only polylogarithmically many distinct outputs N<=X. The transducer may be asynchronous, nondeterministic and partial, and need not report which bank row is used.

The proof is structural: every mixed repeated binary period contributes a strictly positive linear drift to g in addition to its exponential term. A pumped finite-state output cannot reproduce that drift. Thus every useful input loop is monochromatic, forcing a uniformly bounded number of binary runs.

This excludes a finite-state construction retaining the actual F2 witness as a route to all-target coverage. It does **not** exclude a proof recognizing only output support, a different non-rational witness encoding, or an algorithm with unbounded integer registers. Fixed radix changes do not evade the theorem. The known streaming compiler with registers evaluates witnesses but still does not select one for every target.

[Selector theorem and exact boundary](selectors/PROOF.md) · [Independent audit](audit/AUDIT.md).

## 4. Inheritance and remaining proof obligation

Source snapshot: acbad6c8636cf8cbb675b4c947fa19e33ac0df23. Session214's induced-growth and surgery barriers, Session215's fixed-cap composition, Session216's three-core formula and density theorem, Session217's exact inverse and run obstruction, and Session219's Fermat/composition work remain in their stated scopes. The curvature obstruction to uniform target halving is credited to Session216, not promoted as new.

The concept board is exact occupancy, binary curvature, repeated-word drift, nonuniform composition, and ribbon state spaces. The successful connection preserves the common nonlinear coordinate; removing it loses the very drift that blocks finite-state selection. Sharp extrema discard interior occupancy, which must be restored by a genuinely different invariant.

Neither a new physical all-odd constructor nor the needed target-dependent inverse has been found. Theorems here improve exact optimization and rule out a concrete proof architecture; they are not a solution of Chan–Pak.

## Reproduce

    python3 verify_all.py

The portable runner replays both producer verifiers and independent audits, checks exact receipts, and verifies SHA256SUMS. Assertions in normal Python runs are load-bearing; optimized runs are only smoke tests. The audit records the repaired zero-length cylinder boundary and independently checks arbitrary fixed-bit domains and very large capacities. No large coverage scan is claimed.
