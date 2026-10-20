# Dependencies and primary-source scope

20 October 2026.

1. The sole external classification input is the characteristic-sequence
   theorem for characteristic-zero plane curves with one place at infinity.
   The exact numerical use is: after elementary degree reduction to a
   nonprincipal pair, both initial pole orders are at most conductor+1;
   the characteristic gcd and primitive-generator rules then classify the
   possible initial pairs for semigroup <3,4>.

   A. Sathaye and J. Stenerson, *On Plane Polynomial Curves*, author text
   dated 10 September 2014, Section 2.2 (ground-field convention),
   Section 3.1(1)–(5), especially (4), printed pp.8–12:
   https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf.
   Directly reopened on 20 October 2026. The published chapter is *Plane
   Polynomial Curves* (1994), pp.121–142,
   https://doi.org/10.1007/978-1-4612-2628-4_6.
   No byte-identical identification of those two versions is asserted.

   S. S. Abhyankar and A. Sathaye, *Uniqueness of plane embeddings of special
   curves*, Proc. Amer. Math. Soc. 124 (1996), 1061–1069,
   https://doi.org/10.1090/S0002-9939-96-03254-6, printed p.1062 Section 2,
   identifies the initial coordinate pole orders and elementary reduction.
   The precise original paragraph was separately checked by producer and
   reviewer of the inherited fixed packet. A new DOI fetch in this session
   returned an internal retrieval error; no fresh successful PDF retrieval
   is claimed. We do not invoke a general equivalence-classification result.

2. Audited local dependency: ../low_delta_source_conductor_exclusion_20261020/
   LOW_DELTA_QUADRATIC_AND_RELATIVE_GLOBALIZATION.md (staging source
   session13-low-delta/). Its trace-dual conductor identity proves the exact
   conductor of each plane polynomial curve and the codimension/semigroup
   calculations. Its actual-source claims are used only in the final
   scoped consequence, not in the moving entire proof itself.

3. Audited local dependency: the separate fixed cubic/quartic continuation,
   FINITE_EMBEDDING_CLASSES_AND_FIXED_CYLINDERS.md and
   EXACT_EMBEDDING_COUNT_APPENDIX.md (staging source
   session13-cubic-quartic/fixed-continuation/). The concrete terminal-pair
   calculations and all fixed-cylinder exceptional cases are inherited.
   The descent of elementary reductions from an algebraic closure to K
   is supplied explicitly in the present Section 3, not attributed to a
   global relative-lifting theorem.

4. The moving quadratic-cylinder packet supplies the inheritance pattern,
   not an omitted theorem: the monic-conductor descent and moving-period
   calculation needed here are proved in full in Sections 4–5. Its fixed
   finite-period lemma underlies the already audited fixed-cylinder
   consequence in item 3.

No Mason–Stothers theorem, Riemann–Hurwitz theorem, generic-to-global
cubic claim, or full degree-six Keller exclusion is imported.
