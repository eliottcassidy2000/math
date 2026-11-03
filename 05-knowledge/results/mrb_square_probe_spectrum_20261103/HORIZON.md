# A finite-horizon repair: sparse births observe every future count

**PROVED all-size, fixed-horizon synthesis; FINITE-EXACT controls.** This
extends the packet's one-step observation bank, not its attainable-count
range. The input T is always the same actual named source. No finite
observer independent of horizon, positive combination, or coverage theorem
is asserted.

The incoming [S162 star-filter theorem](../session162_mrb_square_20261103/mrb_star/NOTE.md)
already proves that reversing edges incident to one vertex gives a Boolean
Hamilton-count polynomial of degree at most two: a Hamilton path uses at
most two such edges. We combine that local fact with the recovered
square-scaled pair theorem. This is not a new proof of local star degree.

## H1. Tensor interpolation with actual named histories

Fix T on n>=1 named vertices and k>=1 future births x_1,...,x_k. Let S_j
be the old incoming set of x_j among the n+j-1 earlier names. Write
F_T(S_1,...,S_k) for the raw Hamilton count of the resulting tournament.
Then F is a Boolean polynomial of degree at most two in *each birth block*
separately, though its total degree may grow with k.

Indeed the orientation indicator of a candidate Hamilton word is a product
of edge literals. The variables of birth block j correspond to a subset of
edges incident to x_j, of which that word uses at most two. Summing word
indicators preserves the separate degree bound. The blocks partition the
new arcs; an edge between two newcomers belongs only to the later block.

For a target set S with t=|S| define weights on A subset S, |A|<=2:

    lambda_empty(S)=(t-1)(t-2)/2,
    lambda_{i}(S)=2-t,
    lambda_{i,j}(S)=1.                                  (H1)

Zero weights may be omitted. Ordinary quadratic interpolation in each
block gives

    F_T(S_1,...,S_k)
      = sum_(A_1,...,A_k) [product_j lambda_Aj(S_j)]
          F_T(A_1,...,A_k),                             (H2)

where each |A_j|<=2. One can prove this by expanding a multilinear monomial:
in each block the formula is exact on1, each variable, and each pair.
Tensoring preserves it on every allowed monomial. Every term on the right
is the count of an actual birth program on T, with the same vertex names
and exactly k births. Replacing an earlier cut does not relabel later ports.

The full sparse bank has

    product_(j=0)^(k-1) [1+(n+j)+binom(n+j,2)]             (H3)

programs, versus 2^(kn+k(k-1)/2) unrestricted programs. This is finite for
fixed n,k and polynomial in n for fixed k. It is not bounded independently
of k. The coefficients can be negative: (H2) is an observation identity,
not a positive assembly of one source from the other sources.

## H2. All sparse histories fit in one finite arithmetic host

Assume T is supplied with a faithful square-sum orientation chart as in
README P1; otherwise retain its old arcs separately. There is one finite
positive-integer host containing a selected copy of every sparse program
in (H3), with identical old source labels across all copies.

Start with the source label pool P_0. At stage j apply P1 to the *entire*
current pool, rescaling it once, and append its candidates for all incoming
subsets of size0,1,2. Store the map from each such subset of pool names to
its candidate's new pool name. Different candidates have no square edges
between them. The pool grows from m to

    m+1+m+binom(m,2)=binom(m+2,2).                       (H4)

An actual history chooses one candidate per stage, whose neighbor subset
is the supplied A_j among its actual previously selected names. P1's exact
nonedges ensure that no other selected old vertex becomes an unintended
incoming neighbor. Later global scales preserve all earlier incidences.
At the end every selected history shares the same scaled P_0. The relative
order is the source order followed by birth order; unused pool vertices do
not participate in its Hamilton paths. This proves the claim inductively.

For C3 and two future births the pool sizes are3,10,66. This one66-label
host contains all77 sparse two-birth programs as selected five-vertex
graphs. Their counts determine all128 unrestricted two-birth counts.
The code decodes all77 actual selected graphs and checks every selected
pair for arithmetic adjacency, not just their path counts.

This common-host construction grows rapidly, and its numeric scales do
not stabilize with increasing horizon. It therefore does not give an
infinite integer embedding or a constant-size recursively closed carrier.
It supplies a finite, source-preserving hierarchy that repairs the
one-step triangle blindness at the next level.

## H3. The free pair is a terminal privilege

README P2 saved a pair probe using equality of the two constant cuts.
That equality applies when the observed newborn is the *last* birth.
It need not survive a fixed future context.

For a minimal two-coordinate witness, take the old edge0->1. Adjoin x with
incoming set A, then adjoin y with incoming set {0} among0,1,x. The final
four count values are

    A:      empty   {0}   {1}   {0,1}
    F(A):      3     1     5       1.                  (H5)

In particular F(empty)!=F({0,1}). Applying the terminal constant-cut
constraint here would incorrectly predict3 for the last entry. The even
smaller one-coordinate witness starts from a singleton and has final
values3 and1. The quadratic degree bound survives; the missing constraint
does not. Accordingly (H3) retains every pair at every nonterminal stage.
For each fixed sparse prefix, the terminal stage can still omit one pair.

This is a direct refinement lesson: a datum recoverable before gluing may
become independent after gluing. The incoming source-return results and
the earlier square-clock shared-port failures are instances of the same
need to state the future operation language. Their objects are not
identified with one another.

## H4. Tests and the remaining task

`horizon.py` checks all128 two-birth continuations from each named triangle
orientation and all4096 three-birth continuations from C3. It evaluates the
right side of (H2) using only sparse-program counts and compares with direct
subset DP on the requested graph. Independent full-permutation enumeration
recomputes all4096 three-birth values and checks the frozen digest. It also
checks the66-label common host and both terminal-constraint hostiles.

```sh
python3 horizon.py > horizon.json
python3 -O horizon.py > /tmp/mrb-square-horizon-optimized.json
cmp horizon.json /tmp/mrb-square-horizon-optimized.json
python3 audit_horizon.py > horizon_audit.json
```

The target now separates into two obligations: finite future observations
can be made faithful by this hierarchy; finding a bounded family of actual
sources closed under all required target branches is still OPEN. The
incoming strong returns31->95->285 and31->97->291 prove particular returns,
while the65-atom exclusion of their six next exact source states shows
that simply replaying that alphabet is insufficient. A larger observation
bank diagnoses that failure; it does not supply the missing source.

The final incoming two-ear computation continues the zero branches to855
and873 and exhaustively excludes the positive offsets for those supplied
parents with at most two added vertices. See INCOMING item8: this sharpens
the source-selection boundary without changing the interpolation theorem.
