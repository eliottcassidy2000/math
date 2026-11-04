# Retaining q=3: the old cap has no entry, and relocation overshoots

2026-11-03. **PROVED deleted-joint star formula and positive relocation
floor; FINITE-EXACT complete censuses of two fixed named carriers.**
This is a separate relaxation of the q1 repair target. It keeps q=3.

In the original C2 children, varying only inner0's incident edges while
retaining root99 end defect0 gives six H1409 states and three H1427 states.
None has a start defect1 or2 at any vertex. Moving physical200 after the
terminal triangle produces positive entries but cannot have either target
count: every strong state already has at least1536 Hamilton paths.

## 1. Exact inheritance and the two actual carriers

The input graphs are the actual C2 children in
[`entry_selector/receipt.json`](../../../c4_port_returns_20261103/entry_selector/receipt.json),
not a replacement by equal counts. Every physical name survives. Deleting
inner mark w=0 leaves

    D=TT(99,10,11,2,1,3,4,12,200) -> C3(201->401->400->201).

The alternative D' moves200 to the right:

    D'=TT(99,10,11,2,1,3,4,12) -> C3(201->401->400->201) ->200.

Only {200,201},{200,400},{200,401} change outside the w-star. There are
4096 star cuts for each carrier. Both deleted graphs have exactly three
Hamilton paths. The board is source names / positive path-cover floor /
star count compensation / exact root guard / terminal entry / deletion
resource. The corrected near miss is that retaining q3 is a necessary
escape from q1 obstruction but does not itself regenerate a positive entry.

## 2. A star formula over any fixed deleted graph

For a fixed tournament D and x_i=[w->i], let E_i,S_i be its Hamilton endpoint
counters, and J_ij its ordered spanning two-path-cover counter with first
path ending at i and second starting at j. Then

    H(D+w)=sum_i (1-x_i)E_i + sum_j x_j S_j
            +sum_(i,j) (1-x_i)x_j J_ij.                 (1)

The terms put w last, first, or between two nonempty paths. This is an
exact bijection of actual paths and ordered covers, not a sufficient
certificate or an equality of quotient statistics. A fresh deleted-graph
profile supplies every coefficient. Formula(1) then classifies all4096
cuts for each fixed D without enumerating isomorphism classes.

Strongness is separately checked by forward and reverse reachability.
For original D it is equivalent to w->99 and at least one cap vertex->w;
for D' it is equivalent to w->99 and200->w. The former has1792 strong
stars and the latter1024. In either strong family, root99 has sole incoming
vertex w. The cover concatenation used in the parent packet gives

    E_99(2)-H=S_w(T-99,1)=3[ w->10 ].                   (2)

Thus the root guard leaves896 and512 states, respectively. Every guarded
state is independently counted by full literal subset path DP. Formula(1),
strongness, root guard, every start/end defect and every ordered joint are
therefore checked by actual named graphs, not by a reused old phase gauge.

## 3. The positive-cover floor explains relocation failure

More generally let D'=TT_l -> C3 -> z, l>=1, with r the first prefix vertex.
Every strong one-apex extension has w->r and z->w. All ordered covers
counted by J_D'(z,r) concatenate to distinct Hamilton paths through w.
The other l-1 prefix vertices choose their slot independently, and the
triangle contributes its ordered two-path cover count P2(C3)=12. Hence

    H(T) >= J_D'(z,r)=2^(l-1)*12.                       (3)

This is a positive contribution, so no cancellation with the other terms
of(1) is possible. For the actual l=8 carrier the floor is1536, exceeding
both1409 and1427 before imposing any root-end or entry guard. The computed
coefficient J_D'(200,99)=1536 independently verifies the structural count.

Relocation does supply entries. Physical z=200 has sole out-neighbor w.
Cutting the distinguished start slot at its unique outgoing edge gives

    sigma_z=S_z(2)-H=E_w(T-z,1)
                       =#{cap vertices beating w}.     (4)

After removing z, a Hamilton path ending at w is the forced transitive
prefix followed by a triangle Hamilton path ending at an in-neighbor of w.
Each cap endpoint contributes exactly one. Therefore k=1 or2 cap incoming
ports give the requested positive entry, while k=0 or3 do not. This
entry mechanism coexists with the larger positive-cover lower bound.

## 4. Complete finite results and saved actual witnesses

| Deleted carrier | All cuts | Strong/root-end0 | States with some sigma1/2 | H1409 / H1427 hits |
|---|---:|---:|---:|---:|
| TT9 -> C3 | 4096 | 896 | 0 | 6 / 3 |
| TT8 -> C3 ->200 | 4096 | 512 | 384 | 0 / 0 |

The original carrier's no-entry statement here is a complete finite result
on these896 states; it does not rely on extrapolating the finite census.
The separately proved [structural packet](../../structure/PROOF.md) supplies
the all-order no-entry theorem for TT_m->C3 with m>=3 and its smaller-prefix
boundary controls.
The relocated carrier has minimum guarded H1683 and minimum positive-entry
H2075, which are census extrema, not claims of optimal all-order formulas.

`census.json` contains every star word with its physical bit order, actual
deleted carriers, coefficient profiles, strongness and typed entries.
`witnesses.json` preserves all nine original-carrier exact-count witnesses,
including the original children themselves. Each has actual rows, names,
incoming/outgoing star cuts, all changed physical pairs relative to its
own source, every port defect, full E1/S1/E2/S2/J, and fresh golden/F49/F343
observer values. These are source-relative graph descriptions; there is
no native repair program or intermediate-strongness claim.

All1408 guarded actual graphs were independently profiled; the two inherited
children are positive controls. No graph is accepted solely because its
Hamilton count agrees with the target. The requested combination is still
missing: exact original count, root guard, preserved names and positive
entry in one actual state.

Reproduce from repository root:

    python 05-knowledge/results/compensation_returns_20261103/star/q3/reproduce.py

Normal and `-O` runs must agree byte-for-byte. The parent fixed-q1 packet
and all older source packets are read only. This bounded experiment is the
last carrier search in this packet; the surviving structural question is
how to add an entry without paying the full positive cover mass in(3).
