# Independent audit: idoneal information channels

**PASS, with two scoped publication corrections.** The reviewed source note was
`5b7790cc4adf658ab14cf893941fdbe0bd7a2d2c04aed135141cc105c6c1c229`.
The published overlay changes only the Kani citation and the ambiguous `nonROOT`
sentence. The complete original producer packet is retained under
`../provenance/original/`; `../provenance/approved-corrections.patch` is the exact
diff. The original odd-map counterexample and all mathematical tables are unchanged.

## Independent proof review

For a finite abelian class group A, every quotient fibre of A -> A/A² has
|A²| elements. Uniform conditional entropy is therefore log2|A²|. The
identification with genera and the idoneal criterion are classical inputs,
not conclusions inferred from a finite sign table.

At discriminant -252, the complete reduced-form list has four classes, and
(8,-6,9) is not its inverse (8,6,9). Thus the group is C4 rather than C2².
Our independently implemented proper composition additionally gives its square
(7,0,9). Discriminant -252 = 6²(-7) has conductor 6; the maximal order at -7
has only (1,1,2). The conductor cannot be discarded. The displayed genus signs
also follow directly modulo 3. Their two-element fibres lose exactly one bit.

At discriminant -1092, independent b-first enumeration reproduces all eight
printed forms, including the boundary representatives. Canonical reduction of
each inverse returns the same form, proving exponent two. For every literal
printed representative n we found a primitive representing pair and independently
computed the three signs using n modulo 4, 3, and 7. The latter formulas are the
quadratic reciprocity evaluations of (-4/n), (-3/n), and (-7/n), respectively.
All eight triples occur. Classical genus theory identifies these characters as
homomorphisms; the explicit finite composition check below is an additional
control, not a replacement for that cited theorem.

The three selected edge pairs are exactly the three perfect matchings of labelled
K4. Enumerating all 64 graphs independently gives exactly the eight regular graphs
listed by this code: the empty and complete graphs, three matchings, and their
three complementary cycles. Their edgewise XOR is addition of the three bits.
The literal class eta=(6,6,47) has code 111; actual form multiplication verifies
that eta translation is graph complement. Inversion and order reflection both
act identically here. Every staircase tournament has the advertised marked path
and is fixed by converse followed by reversal of that path. No inference about
unmarked isomorphism classes, other paths, tournament regularity, or metagraph
connectivity follows. For rank r>0 the abstract extension uses an injective linear
map into a sufficiently large reversal-fixed edge space sending a chosen nonzero
eta to its all-one vector; extend eta to a basis. Rank zero is correctly excluded.

The joint mod15/mod35 register is (3,28) for both 63 and 273, but the actual
unstopped primitive transitions are (95, clock1) and (51, clock4). Direct odd-map
calculation gives U(273)=205 with valuation2. These are different typed clocks,
not stopped arrival clocks. The flag counterexample U(3)=5, U(7)=11 and
U(17)=13 is correct. The original phrase calling 3 and 7 `nonROOT` was ambiguous:
both lie in the inherited arrival frontier. It is replaced by the explicitly
unstopped odd-map statement. No frontier membership is used in this counterexample.

## A genuinely separate composition control

`independent_idoneal.py` imports no producer code. It enumerates reduced forms
with b as the outer variable, implements canonical Gauss reduction, and finds
proper determinant-one substitutions making leading coefficients odd, coprime
to the discriminant, and pairwise coprime. For equivalent forms (a,b,c) and
(A,B,C), both b and B are even. CRT selects h with h=b/2 mod a and h=B/2 mod A.
The Dirichlet product is the reduction of

    (aA, 2h, ((2h)^2-D)/(4aA)).

We check determinant one, discriminant, reduction back to each original class,
product integrality, closure, every inverse, all 16 C4 and all 64 C2³ products,
the XOR law for all 64 products, and actual eta translation. Complete product
tables are printed in the receipt. Other controls parse the literal Markdown
form table rather than duplicating it silently, verify primitive represented
integers and signs, enumerate all regular labelled K4 graphs, construct actual
T5 arc sets, and replay the distinct P and U clocks.

The independent normal and optimized runs are byte-identical and perform 11,339
explicit checks using exceptions rather than removable assertions. The producer
normal and optimized runs likewise agree at 3,135 checks. The bounded recovery
of 65 candidates through 1848 is explicitly not a global completeness theorem.

From this directory run:

    python3 independent_idoneal.py
    python3 -O independent_idoneal.py

Compare with `independent.normal.json` and `independent.optimized.json`. The
checker pins the staged proof bytes in its output; it also accepts `--proof`.
From the parent directory the producer commands are:

    python3 verify_idoneal_channels.py
    python3 -O verify_idoneal_channels.py

## Literature and scope audit

The originally supplied Kani author URL redirects to an author homepage. The
verified publisher PDF is [Kani, Idoneal Numbers and Some Generalizations,
Annales des sciences mathematiques du Quebec 35(2), 197–227](https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf).
Directly read printed pages 200–202 identify genera with square-subgroup cosets
and give Theorem6/Remark7's equivalent criteria. This is why the citation changes.
[Brandhorst–Sonel–Veniani, Remark1.5](https://arxiv.org/html/2003.08914v7)
distinguishes GRH completeness and possible exceptions; the all-orders bound
must not be replaced by the squarefree one-extra bound. The
[Elsenhans–Kluners–Nicolae abstract](https://arxiv.org/abs/1803.02056)
retains its stated zero-free hypothesis.

The public OpenAI TeX blob 2000faddbbebac5de0ecfe0b962534ea61a852d5 was directly
read and does state the 65-number consequence attributed to the repository's
adopted input. This audit does not re-prove its analytic theorem. None of the
explicit order computations, maps, clock counterexamples, or entropy argument
uses completeness. No new priority, analytic, idoneal-completeness, universal
Collatz, native annulus energy, or six-colour forcing claim is approved here.
