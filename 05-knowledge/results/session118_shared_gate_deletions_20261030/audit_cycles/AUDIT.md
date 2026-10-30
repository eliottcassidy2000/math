# Independent audit: repeatable tight GP/P2 prefixes

**Verdict: PASS, with the producer's stated scope restrictions.** This audit
checks `../worst_cases/TIGHT_CYCLE_HOSTILE.md` against the pinned elimination
rules. It does not import the producer's Python implementation and does not
certify the external candidate's full theorem.

## Structural checks and priority

The fair affine system has two equations A=x+a and C=x+b per block and an
identity auxiliary matrix. The complete charged incidence graph, not just a
row-count word, is retained. Every object has an output descendant, every gate
has two distinct nonconstant pins, and every AND with affine parents has two
independent affine forms. No unpaired parameter has a unique equation child.
Thus no D1-D4 or M step is forced before the chosen G=xy.

After GP sets G=0 on the three-point domain, its sole child K=G+z reroutes to z.
Exactly G and K disappear. The x parameter still has two distinct equation
children A,C; y still has the gate child B=yb. All other parameters retain an
occurrence. Therefore P0 and PC are absent, P2 applies to x, and PE does not
preempt it. In particular no D4 or newly movable affine-parent XOR intervenes:
the output XOR receiving z has nonlinear B as its other parent.

Safe x=0 deletes A,C and substitutes a,b into H. H becomes ab, a nondegenerate
AND, not a constant or unary gate. The surviving B,T,H products keep every
output combination nonaffine at the syntactic preparation stage. All surviving
parameters are used and none has a unique equation child. The next unprocessed
G remains the first nonlinear gate in the prescribed topological list. This
justifies the consecutive prefix, including every earlier-priority closure.

Each block transition has exact changes

    GP: (-2,-2,+1,-1), P2: (-2,+1,-1,0)

in (s,I,q,a). There is no double-counting of shared descendants. The common
child H of the two affine equations survives both substitutions.

## Function and loop eligibility

Independent expansion gives

    Q=x+xa+xb+ab+yb+ac+xy+z.

Its polar radical is precisely {0,11,32,43}, with coordinates ordered
(x,y,a,b,c,z) from least to most significant bit. Hence polar rank is four.
The z direction is radical and Q(z)=1. On k blocks, the polar matrix is block
diagonal, rank 4k, and the radical linear restriction remains nonzero.

The quadratic-flat lemma is valid: for a constant affine flat a+U, U is
isotropic and its radical intersection lies in the kernel of Q on the radical.
These give the upper bound N-r-1 when polar rank is 2r. Conversely a lift of a
quotient Lagrangian, enlarged by the radical kernel, has that dimension. The
restriction of Q is linear and annihilates its radical intersection, so it can
be canceled by a polar functional B(a,-), producing a constant coset. Thus the
bound is exact, not merely a necessary dimension estimate.

For this family d=4k-1. Before each chosen GP and before its following P2,
a>=d+2. The actual function therefore supplies the needed affine-disperser
hypothesis and loop guard. The counterexample does not rely on an impossible
choice of d for a simple affine output.

## Scope that must remain visible

The initial ten-object block has a seven-gate implementation

    (x+b)(x+a+y)+ac+z.

The six-gate residual block also has a five-gate implementation

    b(a+y)+ac+z.

Consequently this is a counterexample to a universal positive-surplus or
bounded-run assertion for *arbitrary syntactically reduced states with the
specified legal choices*. It does not refute a strategy that resynthesizes
semantically, chooses different eliminations, or establishes additional
minimum-circuit structure. Original-circuit minimality cannot simply be
assumed for every restricted intermediate state.

Also d is linear in n, not logarithmic. The family does not decide whether a
low-d hypothesis excludes long tight prefixes. These qualifications are
essential and correctly present in the producer report.

## Independent certificate

`check_independent.py` builds its own three stage descriptions and audits
incidences and preparation obstructions for block sizes 1 through 25. It checks
325 exact cycles / 975 stage graphs, both loop guards, the direct-sum polar
ranks, all 64 one-block truth assignments, both smaller resynthesis identities,
and the radical by direct polarization. It imports no producer code.

Normal Python and `python3 -O` produce byte-identical `normal.json` and
`optimized.json`. Explicit exception checks survive optimization. The producer's
six-variable affine-flat census is additional evidence, not reimplemented here;
the all-k dimension claim is audited analytically above.

## Safe resynthesis corollary

[FACTORING_REWRITE.md](FACTORING_REWRITE.md) supplies the missing constructive
step: each residual block's exclusive common-factor cone saves one gate, and
all k blocks together save k distinct gates without altering outside users.
The independent checker additionally validates 1,016 general factoring truth
assignments and 1,056 whole residual/replacement DAG assignments. This gives
surplus for the exhibited family, not automatic surplus in every legal trace.


## Audited affine-state continuation

After the safe rewrite, M moves U=y+a into the fair system and D4 legally
exchanges y with U and deletes the old auxiliary y/equation. The input map
must retain y=w+a. This saves a second charged object per block without
changing I,q,a. The final representation has 5k-1 gates on 5k essential new
parameters and is minimum in those coordinates. In original coordinates,
the exhibited circuit still needs k XORs to form w and has 6k-1 gates.
The checker verifies the affine map, its inverse, matrix invertibility and
semantic agreement independently. This is not a free input-coordinate change.
