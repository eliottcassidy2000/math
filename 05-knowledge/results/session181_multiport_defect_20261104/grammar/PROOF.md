# Named bad-join defect transport and the prototype continuation boundary

2026-11-04. **PROVED exact unbounded-state transport; FINITE-EXACT bounded extension census.** This is not an all-depth two-digit realization or a global spectrum theorem.

## Provenance and scope

Session180's named squarefree endpoint algebra and paired prototype are the starting point. Current main was also checked at `9f5bdfaea3a2debeaac8d076c0367d5006aca71a`, specifically `05-knowledge/results/zero_ready_strong_switch_20261103/README.md` and `DEPTH.md`. That newer packet already supplies two generations of both branches on the actual 170807 lineage. The small prototype census below is about a different, clock-uniform source family; it neither supersedes nor contradicts those stronger actual-source certificates.

## 1. The defect is an explicit joint-port observable

For a named tournament B on V, let G_ij be the squarefree generating polynomial of all simple paths from i to j, with vertex weight z_v and z_v²=0. Put e_w=sum_i G_iw, s_v=sum_j G_vj and P=sum_ij G_ij. The marked two-slot count is F_w=[V]e_w(1+P); the Hamiltonian count is h=[V]P.

Define Q_w=sum_{v→w}s_v. Then

    b_w := F_w−h = [V] e_w Q_w
                   = sum_{v→w} J(w,v),
    J(w,v) := [V] e_w s_v.                         (1)

Proof. An ordered spanning two-slot cover whose first slot ends at w is either joinable (the second slot is empty, or its first vertex v satisfies w→v) or has a bad join v→w. Joining the two slots in the first case is a bijection to Hamiltonian paths, with inverse splitting a path just after its occurrence of w. The second case is counted by (1). This proves nonnegativity without a numerical comparison, and displays exactly which physically named pair correlations make up the defect. It holds for arbitrary tournaments, not only strong ones.

## 2. All-order birth transfer for the changing defect

Birth a fresh name x, writing a_i=[i→x] and c_i=[x→i]. Every polynomial on the right is evaluated on the **old** named B. Define

    L=sum_i a_i e_i,        R=sum_i c_i s_i,
    R_w=sum_j c_j G_jw,
    K_w=sum_{v→w} sum_i a_i G_vi + c_w.

Here c_w is a constant polynomial, and 1 is the empty path slot. The exact transfer is

    H(B+x) = [V](1+L)(1+R),
    b_w(B+x) = [V]((1+L)R_w Q_w + e_w(1+R)K_w).     (2)

Proof. Every path using x splits uniquely at x, with its incoming and outgoing pieces possibly empty. For the defect use the already proved endpoint update to obtain

    e'_w=e_w+z_x(1+L)R_w,
    Q'_w=Q_w+z_x(1+R)K_w.

For the second identity, each old in-neighbor v contributes the update of s_v; x contributes s'_x=z_x(1+R) exactly when x→w. Multiply, use z_x²=0, and extract [V∪{x}]. This proves (2). The formula works after arbitrary previous births: keep the updated full matrix, the current physical incidence relation and the actual mark each time. Thus it is a genuinely unbounded exact state grammar, but the state can grow exponentially and the formula does not force either requested affine digit.

Edge reversal, deletion and re-marking can use Session180's full-matrix rules before extracting (1). Native Hamilton-path-complement guards still require separate literal witnesses; algebraic reversals are not automatically legal native programs.

For a clock-preserving ternary lift, the two independent constraints remain

    h'=3h+d,       b'=3b−d,       d∈{2,4}.

Checking (2) or merely preserving a finite-field codec is not a proof that either target constraint is realizable.

## 3. Complete two-birth census from the previous uniform prototype

The Session180 prototype outputs, with physical mark0, are

    C2=[56,33,51,54,2,16],  (h,F)=(17,21),
    C4=[56,1,19,22,2,30],  (h,F)=(19,21).

Enumerate all 64 incidence cuts for a fresh vertex6 and all128 cuts for fresh vertex7. All old edges and names are fixed. Each source therefore has8192 labelled two-birth extensions. Retain strong output suffixes, and require the original physical mark0 to satisfy(h',F')=(3h+d,3F).

- From C2 there are exactly six +2 extensions and no +4 extension.
- From C4 there are no extensions for either digit.

The six first-source cuts (old vertices beating the newborn) are

    (15,15), (15,79), (15,111), (47,15), (47,47), (47,111).

For example the first gives rows[248,225,243,246,2,16,48,112], with(h,F,b)=(53,63,10). Hence this is one more uniform +2 lift on every transitive-prefix clock, preserving all old names and the port.

Without the strongness test, C4 has two misleading +2 responses: cuts(57,127) and(63,57). Each has a global sink (vertex7 and vertex6 respectively). They cannot serve as a strong marked suffix in the canonical zero-ready grammar. This is an explicit check against accepting the response pair alone.

A larger independent menu (`search_allmarks.cpp`) allows **any of the eight output vertices as the mark**, including either newborn. It has12 +2 marked realizations from C2 and zero +4, and no realization of either digit from C4. The added six +2 realizations use marks2,6,7. Every surviving marked response is independently permutation-checked in `verify.py`. Thus re-marking alone does not restore the missing branch in this menu.

The negative result is only for exactly two births with the old induced graph retained. It does not exclude more births, edge rewrites, a different marked suffix with the same scalar response, or the newer actual-source branch certificates.

## 4. Reproducible verification

`search.cpp` exhaustively computes the above census by integer subset DP and independently checks directed reachability for strongness. At order8 all counts fit uint64_t: F is at most8·8!, since a permutation and a cut determine every two-slot cover. No larger-order arithmetic claim is made.

`verify.py` uses explicit permutation enumeration, not DP, to check every displayed survivor. Separately it exhausts every labelled tournament of orders1 through4, every existing marked vertex and every incidence cut for one new vertex. For each case it verifies both (1) and (2) against direct permutation counts; its squarefree matrix is itself generated by ordered-tuples enumeration. Checks use explicit exceptions and remain enabled under Python optimization. The receipt is `verification.json`.

Commands: compile `search.cpp` with the existing C++17 compiler and direct output to `search.out`, then run `python verify.py` and `python -O verify.py`. No installs or repository startup helper are used.
