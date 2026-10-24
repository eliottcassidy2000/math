# Base-nine Lucas modulo nine: extension packet

Status: PROVED as a consequence of inherited TheoremB at p3, rather than independently of that theorem. Its exceptional lift identity has been independently reproduced above.

Modulo9, Q=1+t and H_3=1+6t (d0=1,d1=15,d2=315). Since t^sigma=t³ modulo9, TheoremB gives G=s(1+6t)G(t³). Apply this twice:

G=s(t)s(t³)(1+6t+6t³)G(t9) modulo9.

Now (1+t)(1+t³)=(1+t)^4−3t(1+t)². The normalized square root is therefore
s(t)s(t³)=(1+t)²−(3/2)t=1+5t+t² modulo9;
the omitted quadratic perturbation is divisible by9, and all denominators are powers of2, hence units. Multiplication gives

G(t)=(1+2t+4t²+3t³+3t4+6t5)G(t9) modulo9.

Consequently c(n) modulo9 equals the product of digit weights [1,2,4,3,3,6,0,0,0] of the base9 expansion of n. In particular c(9n)=c(n) modulo9 for EVERY n. The coefficient test confirms this digit rule independently through8000. The simple Cartier swap itself does not lift to9, as the explicit d1 counterexample shows. No novelty claim: this is an immediate algebraic consequence of the inherited inert branch, newly extracted here.

For general inert p, c(p²n)=c(n) modulo p² remains only FINITE-EXACT in this packet. Tested ranges (p,n_max): (3,888),(5,320),(13,47),(17,27),(19,22),(31,8),(41,4),(47,3),(59,2),(61,2),(73,1),(83,1),(89,1). This is distinct from the inherited OPEN split-prime strengthening c(p²n)=c(pn) modulo p4. Do not substitute one for the other or infer an all-inert theorem from the finite data.

## Corrected companion and exact defect budget

From G=P9 G(t9), applying C3 gives C3G=(1+3t)G(t³) modulo9, because C3P9=1+3t and G(t9)=G((t³)³). In contrast H=(1+6t)G(t³) modulo9. Thus C3G-H=-3tG(t³) modulo9. Also C9G=G implies C27G=C3G modulo9. These formulas include n0.

The base9 digit rule yields an exact support criterion modulo9: c(n) is nonzero modulo9 iff no digit is6,7,8 and at most ONE digit lies in3,4,5. Digits0,1,2 have unit weights1,2,4. Digits3,4,5 each have weight divisible by3 but not9. Two such digits force a factor9. This is a precise one-defect budget, not an analogy.

The mod27 strengthening is false: c9=26337699740, c1=2, so c9-c1=18 modulo27. Neither C9G=G nor the base9 product identity extends verbatim modulo27. The naive C3G=H already fails modulo9, with difference6 at n1.

## Dependency and freeze ledger

This extension is separate from the prime-dichotomy packet. The initial packet was temporarily amended before the separation request arrived, then restored to its pre-addendum scope. Only the final manifest hashes transmitted to the auditors identify authoritative frozen versions.

The sole non-elementary proof input is TheoremB, inert p3 case, in inherited-proof.md (sections2–3), together with the exceptional lift identity in section6. The latter is independently checked by check_modular.py past its stated Sturm bound32, but the modularity and whole TheoremB remain inherited dependencies. The original proof source is repository commit7818149bf0a1eca778549fdfd1dde34f479421fd, path05-knowledge/results/cooper_s7_weight_one_root_frobenius_20261018.md. This is a PROVED CONSEQUENCE OF AN INHERITED THEOREM; the independent audit of that theorem is separate. No new literature originality asserted.

Reproduction: python3 replay.py (works from any working directory). It verifies the closed exact file manifest, executes each verifier under normal Python and -O, and compares outputs to both committed fixtures. No external data, installed mathematical package, network, or repository access is needed for this finite replay. The general inert-prime table above records a separate finite exploration; only the p3 theorem is promoted by this extension.
