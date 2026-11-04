# A second literal +2 return, with a complete two-edit suffix boundary

**FINITE-EXACT new strong realization; PROVED all-size transport and isolation.**
Source pin: eliottcassidy2000/math, commit `4ad839b20dd280cdcfbadf86cb62b2be3c9615ad`, `05-knowledge/results/triangle_return_170807_20261103/{interface.json,zero_ready.json,ZERO_READY.md}`.
No closed ternary grammar or global Hamilton-path spectrum is claimed.

## 1. The new marked suffix

The inherited 170807 tournament has physical root39, sole in-neighbor31,
prefix50→51, and strong suffix in the order

    [1,2,4,6,7,8,14,21,22,31,32,33,40,54].

Its marked suffix responses are `(h,F)=(2195,42153)`, where
`F=E_31(B,2)` counts ordered spanning two-bin path covers, first bin
nonempty and ending at31, second possibly empty.

Reverse internal suffix edges on unordered pairs `{4,31}` and `{6,22}`.
The resulting strong suffix B' has

    H(B')=5655, E_4(B',2)=31673.

Retain physical root39 but change its sole in-neighbor from31 to4 by
reversing the two attachment edges `{39,31}`, `{39,4}`. Append two new
vertices55,56 to the prefix, giving50→51→55→56. Every prefix vertex beats
every suffix vertex;39 beats the prefix and all suffix vertices except4.
The new full tournament is strong by the canonical zero-ready grammar,
retains all seventeen old named vertices, and has

    H =5655+16*31673=512423=3*170807+2.

`candidate.json` supplies its complete nineteen-vertex incidence relation.
`build.py` reconstructs it from the pinned source, recording all four old
edge changes and two births. Independent native compilation succeeds with eight actual complements,
four temporary helpers100..103 (each deleted), and the two permanent births.
Completed counts are170807→170845→304027→428445→132347→259039→512423;
every completed repair is strong and peak order is19. The named guard
paths are recorded in the independent native audit. Joint and phase
certificates are separate audits: old spare6 has changed an incident edge
and must not be assigned its old port response by identification.

## 2. All-size arbitrary-prefix transport

For every nonempty tournament A, take source Z(A,B,31), and construct
Z(A→TT2,B',4). Both are strong zero-ready tournaments with the same actual
root39. By the inherited identity

    H(Z(A,B,w))=H(A)H(B)+P_A(2)E_w(B,2),

and by `H(A→TT2)=H(A)`, `P_(A→TT2)(2)=4P_A(2)`, their difference is

    H_out−3H_in =233 P_A(2)−930 H(A).                 (N1)

The product formula for P follows by restricting each bin to A and TT2:
each of the two transitive vertices chooses its bin independently, and
all concatenations are legal across the directed cut.

For an n-vertex A, cutting a Hamilton path at each internal position
injects `(n−1)H(A)` objects into ordered nonempty two-bin covers. Adding
the two empty-bin possibilities yields

    P_A(2)≥(n+1)H(A).

Consequently for n≥4, (N1) is at least235H(A)>4. For the only tournament
isomorphism classes of orders≤3, `(H,P_A(2))` are

    TT1:(1,2), TT2:(1,4), TT3:(1,8), C3:(3,12).

Their offsets are respectively−464,2,934,6. Therefore this actual
all-size suffix/prefix transformation gives +2 precisely for A=TT2,
and never +4. The new count is a positive next step, not a reusable
return engine. On transitive A=TT_s the offset is233*2^s−930.

## 3. Complete bounded-suffix / unbounded-prefix classification

Start with the fixed fourteen-vertex inherited suffix B. Permit at most
two arbitrary internal edge reversals, retain any resulting strong
suffix, select any of its fourteen physical vertices as the new sole
in-port, and attach ANY nonempty tournament A as prefix. Root39 remains
outside, with its actual attachment edges set by that selected in-port.
No constraint on |A| or its internal structure is imposed.

`search.cpp` exhausts1+91+binom(91,2)=4187 suffix edits, of which3815 are
strong. For each strong suffix and physical mark it independently builds
subset endpoint path counts and computes all exact `(h,F)`. The necessary
integer equation is

    h*x+F*y ∈{512423,512425}, x=H(A), y=P_A(2).

Here x is positive odd, y positive even, and y≥2x. Testing all
`1≤x≤target/(h+2F)` yields precisely53 candidate records, all y≤198.
For every n-vertex tournament A,

    P_A(2)=sum_(S subset V) H(A[S])H(A[V\S])≥2^n,

with H(empty)=1, because every tournament has a Hamilton path. Hence
all candidate prefixes have n≤7. `prefix.cpp` exhausts all2,131,019
labelled tournaments of orders1 through7 and records their exact profiles.
The sole matched profile in the53 records is `(x,y)=(1,16)`, on the
suffix/mark of section1, for target512423. Such a prefix has H=1 and is
transitive, so y=2^n forces n=4. There is no512425 realization anywhere
in this precisely specified grammar. This is stronger than a finite
prefix-bank exclusion: the bound proves that every larger prefix is
irrelevant to these two exact targets.

The mathematical scope is still narrow: three internal suffix flips,
other suffix changes, additional vertices inside the suffix, native
operations leaving this grammar, or a different target remain open.

## 4. Immediate fixed-suffix obstruction

Even with no search, the unchanged original marked suffix cannot reach
either child through ANY prefix: the positive Diophantine equations
`2195x+42153y=512423` and `=512425` have no positive integer solutions.
This removes transitive-prefix restrictions from the inherited next-step
clock obstruction, but does not constrain a changed suffix.

## 5. Exact computation conventions

The search fields `flip_a`, `flip_b`, and `mark` are zero-based
indices, not physical names. Edge IDs enumerate lexicographic pairs
(i,j), i<j, in the displayed fourteen-name suffix order; mark indexes
that same order. Thus IDs31,40 mean physical pairs{4,31},{6,22},
and mark2 means physical vertex4. A−1 edge ID denotes no second flip.

Rows store outgoing-neighbor masks. Let D[S,v] count paths on S ending
at v. Initialize D[{v},v]=1 and use

    D[S,v]=sum_(u in S\{v}, u→v) D[S\{v},u].

Then `H[S]=sum_v D[S,v]`, `H[empty]=1`, and

    E_w(B,2)=sum_(S containing w) D[S,w] H[V\S].

Both subset DPs are exact integer calculations. The largest actual
sums in these menus are far below signed64-bit capacity; no modular
counting or statistical inference is used. Saved finite counts can be
recomputed with `reproduce.py`; the external independent prefix atlas
provides a second implementation rather than an assumption.
