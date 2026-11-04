# Complete C4-port acceptance and the first persistent gaps

2026-11-03. **PROVED:** a finite, target-dependent reduction and exact
coefficient inverse for the declared ordered atom language.
**FINITE-EXACT:** its complete fixed-marginal classifications at
f=21,63,189,126,252,504,1008. All sixteen C4 port masks fail to repair
the common-clock 24/44 targets, or their first four fixed output-clock
shifts. **OPEN:** larger atom languages, unbounded clock shifts,
regenerated entry, and the global tournament Hamilton-path spectrum.

Starting pin: `1b50a1d763`. The closest proved mechanism is the
[seven-coordinate ordered-component law](../../session182_physical_returns_20261104/induction/PROOF.md).
The inherited hostile is the missing 24/44 intercepts in the
[complete singleton/triangle inverse](../../multiport_defect_grammar_20261103/COMPONENT_LANGUAGE.md).
The corrected near miss is treating fresh target acceptance as a physical
conversion of the preceding source. The least-used sidecar here is the
nonnegative coefficient equation for a complete component skeleton.
The live board is C4 port masks, joint correction j, coefficient ideals,
terminal-port strongness, output clock, and actual source compatibility.

## 1. Declared language and full C4 port bank

Every word is an ordered sum of nonempty atoms; every earlier component
beats every later component. The atoms are a singleton T, a directed
triangle C, and the fixed strong four-vertex tournament Q with adjacency
rows [8,1,3,6]. A row's bit v records an edge to vertex v.

T has its two port subsets. C uses one representative of each cardinality
0,1,2,3, since triangle rotation takes any such subset to any other of the
same cardinality. Q uses all sixteen actual port masks. Its four physical
vertices are never collapsed into cardinality classes. The final atom has
nonempty ports; interior atoms may have no ports. This is precisely the
terminal-SCC condition making the usual rooted attachment with any
nonempty transitive prefix strong.

Use the inherited signature sigma=(h,p,f,g,j,e,s): Hamilton count, ordered
two-slot path-cover count, marked first-slot end response, marked
second-slot start response, their joint response, and full Hamilton
endpoint/start counts. The Q signatures all have (h,p)=(5,32):

| Port mask | f | g | j | e | s |
|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | 7 | 5 | 0 | 1 | 1 |
| 2 | 9 | 6 | 0 | 2 | 1 |
| 3 | 16 | 11 | 4 | 3 | 2 |
| 4 | 6 | 9 | 0 | 1 | 2 |
| 5 | 13 | 14 | 2 | 2 | 3 |
| 6 | 15 | 15 | 6 | 3 | 3 |
| 7 | 22 | 20 | 12 | 4 | 4 |
| 8 | 5 | 7 | 0 | 1 | 1 |
| 9 | 12 | 12 | 4 | 2 | 2 |
| 10 | 14 | 13 | 2 | 3 | 2 |
| 11 | 21 | 18 | 10 | 4 | 3 |
| 12 | 11 | 16 | 4 | 2 | 3 |
| 13 | 18 | 21 | 10 | 3 | 4 |
| 14 | 20 | 22 | 12 | 4 | 4 |
| 15 | 27 | 27 | 22 | 5 | 5 |

The singleton and triangle signatures are respectively

    T_b: (1,2,b,b,0,b,b),
    C_k: (3,12,3k,3k,k(k-1),k,k).

In the receipt C uses actual masks 0,1,3,7, rather than k. All 22 atom
signatures are counted both by literal subset path DP and, independently,
by listing all directed permutations on every subset and pairing the
ordered slots. The table retains the orientation gauge and complete
port mask; f alone is not substituted for the word's joint signature.

## 2. Finiteness without an arbitrary order cutoff

For two successive component signatures x=(h,p,f,g,j,e,s) and
y=(H,P,F,G,J,E,S), the proved ordered law is

    x*y = (hH,pP,fH+pF,gP+hG,
           jH+hJ+eS+gF,hE,sH).                          (1)

**Finite reduction theorem.** Fix a positive target marginal F_*. Every
admissible word with marginal F_* has at most
1+floor(log_2 F_*) atoms. More precisely, the product p of the complete
prefix before its final atom is at most F_*. Therefore enumeration of all
such prefixes, followed by every positive-port final atom, is complete.

Proof. If the final atom has signature (H,P,F,G,J,E,S), then F>=1 because
its port set is nonempty. For a nonempty preceding word, (1) gives

    F_* = f_prefix H + p_prefix F >= p_prefix.           (2)

For an empty prefix use p_prefix=1 and f_prefix=0 just in (2); its endpoint
coordinates need not be defined. Every atom multiplies p by at least two:
the multipliers are 2,12,32. Thus a prefix with l-1 atoms has
2^(l-1)<=p_prefix<=F_*. Every valid word is consequently visited by the
finite prefix enumeration. Checking its actual marginal and final port
set gives necessity as well as sufficiency. There is no externally chosen
vertex or word-length cutoff. QED.

An extension of a prefix to a later nonempty terminal port has strictly
larger f than that prefix: its concatenated suffix has F>0, H>=1, and the
prefix has p>=1. Hence the checker may safely stop extending a prefix
once f>=F_*, after first recording any currently valid terminal word.
This pruning is separate from, and unnecessary for, the finiteness proof.

For these particular atoms the smallest nonzero F values are 1,3,5, equal
to their H values. Thus the stronger last-atom guard is
p_prefix*min(F_last)<=F_*. Also

    p_full <= (32/5)F_*,       h_full <= F_*.

The first inequality uses max(P/F)=32/5 over positive-port atoms; the
second uses p_prefix>=h_prefix and F_last>=H_last. Finally each atom has
P>=2^(its vertex count), so an n-vertex suffix satisfies

    2^n <= p_full <= (32/5)F_*.                          (3)

These are derived completeness bounds, not empirical cutoffs.

## 3. Independent coefficient inverse and the exact ideal guard

Fix a component-kind skeleton A_1,...,A_l, forgetting its port masks but
retaining whether each atom is T, C, or Q. Its H_i and P_i are fixed. On
expanding (1), its full marginal is exactly

    f = sum_i w_i F_i,
    w_i = (product_(r<i) P_r)(product_(r>i) H_r).          (4)

Every weight and every available F_i is a nonnegative integer. The final
F_l must be positive. The independent inverse first lists all skeletons
allowed by (2), then solves this finite weighted equation using the full
port bank at each position. Each solution still carries its actual port
masks; it is subsequently evaluated with all seven coordinates in (1).

At a partial choice, let D be the gcd of all remaining weighted available
F values. A necessary condition is that the remaining target lie in D Z.
The checker also keeps the minimum and maximum remaining sums. These are
safe rejection guards. Neither gcd divisibility nor lying within the
interval is claimed sufficient: actual choices are enumerated and their
full nonnegative sum is checked. This keeps the coefficient ideal separate
from its bounded positive realization problem.

This skeleton inverse is independent of the direct prefix search: it does
not use the latter's intermediate f pruning. Their complete sets of words
agree at every target in the receipt. The simple mixed-radix inverse from
T/C alone is not used; Q has P/H=32/5 and generally nonintegral F/H.

## 4. Common-clock targets: C4 does not repair 24 or 44

For an actual suffix sigma and common transitive prefix length t>=1, the
full attachment has

    H_t = h-j+2^t f.                                    (5)

The initial TT3/all-port source has (f,h-j)=(7,-5). Following a digit d in
{2,4} requires f->3f and a->3a+d for intercept a=h-j. The classifications
are complete over this atom language at their prescribed f:

| Level | f | All candidate words | Accepted digit words | Missing digit words |
|---|---:|---:|---|---|
| 1 | 21 | 5 | 2,4 | none |
| 2 | 63 | 9 | 22,42 | 24,44 |
| 3 | 189 | 18 | 422 | 222,224,242,244,424,442,444 |

At f=63 the two extra C4-containing words are

    T_1 T_1 T_0 Q_{mask4}: (5,256,63,201,46,1,5), a=-41;
    T_1 T_1 Q_{mask9}:     (5,128,63,108,52,2,5),  a=-47.

The desired 24/44 intercepts are -35 and -29; neither is present. At f=189,
the five C4-containing words give intercepts -97,-105,-97,-123,-141 and
likewise repair none of the seven missing requested targets.

The new bank does have a positive uniform output at level one. The single
atom Q_{mask11} has sigma=(5,32,21,18,10,4,3), hence intercept -5. It gives

    -5+7*2^t -> -5+21*2^t = 3*(-5+7*2^t)+10.

Thus the extra atom is not inert: it supplies a different digit. No native
operation or reusable-entry statement is inferred from this identity.

## 5. Four fixed output-clock shifts also fail

Allow an output prefix clock t-k, for fixed k>=1 and t>=k+1. A target
intercept a remains unchanged, while its required suffix marginal becomes
f'=63*2^k. The exact finite inverse gives:

| k | f' | All words | Prefix states | 24/44 accepted? |
|---:|---:|---:|---:|---|
| 1 | 126 | 12 | 403 | neither |
| 2 | 252 | 20 | 1,059 | neither |
| 3 | 504 | 26 | 2,819 | neither |
| 4 | 1008 | 42 | 7,555 | neither |

Every row includes all atom word lengths permitted by (2), not just words
of the same order as the starting source. This table does not establish
nonacceptance for unbounded k. A physical realization of a successful
clock shift would also have to account for the old prefix names, for
example by absorbing them into the suffix. No such source-preserving
transformation is asserted or required by the present acceptance test.

## 6. Reproduction, positive controls, and precise limitations

Run from the repository root:

    python 05-knowledge/results/c4_port_returns_20261103/language/check_language.py
    python -O 05-knowledge/results/c4_port_returns_20261103/language/check_language.py

Both modes use explicit RuntimeError checks and write `receipt.json` only
in this directory. The receipt contains the complete 22-atom bank and all
132 accepted fixed-f words, each with literal adjacency rows, actual port
mask, full signature, intercept, and target acceptance indices. Checks are:

- 22 atom signatures by independent permutation enumeration and subset DP.
- Equality of direct prefix and coefficient-equation inverse word sets at
  all seven declared marginal targets.
- All 132 word signatures independently recounted from literal adjacency.
- All 32 unshifted words attached at clocks one and two, with independent
  strongness and Hamilton counts, for 64 full-attachment controls.
- Inherited positive 2,4,22,42,422 acceptance and explicit 24/44 hostility.

Normal and optimized stdout and JSON bytes agree. The inherited checker
is imported read-only and its normalized SHA256 is recorded. No earlier
receipt or proof is rewritten.

The exclusion is exactly about this bank and these fixed slopes and clocks.
The output is an actual target graph if accepted, not a source-compatible
native program. Conversely these missing typed targets do not imply any
missing Hamilton count in unrestricted tournaments. Regenerated entry,
all-order two-digit closure, and the global spectrum remain separate open
obligations.
