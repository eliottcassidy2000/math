# Delayed marking: exact four-coordinate block switching and its boundary

**PROVED all-order block formulas and positive finite-alphabet transfer;
FINITE-EXACT source-image omissions. Complete realization remains OPEN.**
The block-count formula itself is inherited; the marked interior profile
and its closed command language are the present specialization.

The closest mechanism is [Session150's rooted transitive block](../session150_hamiltonian_spectrum_20261103/gadgets/RESULT.md).
The hostile is [the equal-profile H37 pair](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md):
four coordinates do not close arbitrary older-root switching. The corrected
near miss is the balanced reset's factor obstruction. The least-used sidecar
is the position of the next marked vertex inside a newly created block.

## B1. The full marked block response

For an arbitrary nonempty rooted tournament `(T,v)`, retain
`z=(h,s,q,F)` with the same definitions as in the previous adaptive packet.
Put `a=h-s`, `d=F-q`. Add a transitive block X by repeated births at v, so
`v->X->T-v`. Choose w in X, with A vertices of X preceding w in its directed
order and B following it. Equivalently A are newer births and B older births;
`k=A+B+1`. Mark w for the next operation. Then

    h' = a+q+(2^k-1)d,
    s' = 2^B d                         if A>0,
         a+(2^B-1)d                    if A=0,
    q' = h                             if k=1,
         a+q+(2^(k-1)-1)d              if k>=2,
    F' = 2^A a + (2^(B+1)3^A-2^A)d + 2^(A+1)q
         + 1_(B=0)(s-q).                                (B1)

These identities retain a literal source and a specified mark; they do not
declare an arbitrary integral profile realizable. On strong sources every
birth remains strong, and the new root indegree is exactly A+1.

**Proof.** The inherited unmarked block formula gives h'. Deleting w gives
the same block construction with k-1 vertices, proving q', including its
separate k=1 boundary. For a path starting at w, its initial X-run can use
any subset of the B older vertices. If A>0 a nonempty post-v X-run is
necessary, and the old path before v must be nonempty; this contributes
`2^B(F-q)`. If A=0, taking all X initially contributes h-s, while every
proper initial subset contributes F-q, giving s'.

For F', write a two-slot cover as `(P,Q)`, with P ending at w. If P contains
old vertices it contains v, its last old vertex. Each newer X vertex has
three possible bins: before P's old path, after v in P, or before Q's old
path. Each older X vertex has only the first and third bins. There are
F-q old covers whose first path ending at v is nonsingleton, giving
`2^B 3^A d`. For the singleton old path v, the first bin must be empty,
giving `2^A q`.

Otherwise P is entirely in X. It contains w and a subset of the A newer
vertices. The remaining Q is a Hamilton path of T plus a block of size B+r,
where r newer vertices remain outside P. Therefore this case contributes

    sum_(r=0)^A binom(A,r) H(T_(B+r))
      = 2^A(a+q-d) + 2^B 3^A d + 1_(B=0)(s-q).

The correction is precisely the difference between H(T_0)=h and the
positive-block expression a+q. Summing the two cases proves F'. All bins
carry their forced transitive order, so there are no factorial weights.

## B2. A finite positive alphabet for fresh-block changes

Let R keep the old root after one birth. Let J_A add A+1 vertices at that
root, then mark the first-born one. Thus J_0=N. The general command(B1)
is exactly R repeated B times followed by J_A; B is not an additional
unbounded alphabet. With an indegree cap d, the finite alphabet is
`R,J_0,...,J_(d-1)`, acting on sources whose marked root initially meets
the cap. All individual births keep the old root's indegree unchanged;
the new mark has indegree A+1.

All these maps are nonnegative integer linear maps in physical coordinates

    (a,t,s,delta) = (h-s, q-s, s, F-h-q+s).

The inherited endpoint-defect theorem gives delta>=0. Explicitly,

    R = (2a+delta, 0, t+s, delta),
    N = (a+t+s+delta, s, a, 0),
    J1 = (3a+t+s+2delta, a+t+s, a+delta, a+t+2s+delta),
    J2 = (7a+t+s+6delta, 3a+t+s+2delta, a+delta,
          7a+5t+6s+5delta).

For general A>=1 put `D=2*3^A-3*2^A+1`. Then J_A has coordinates

    ( (2^(A+1)-1)a+t+s+(2^(A+1)-2)delta,
      (2^A-1)a+t+s+(2^A-2)delta,
      a+delta,
      D a+(2^(A+1)-3)t+(2^(A+1)-2)s+(D-2^A+2)delta ).

All coefficients are nonnegative. A=1 checks directly; for A>=2 use
`3^A>=2^(A+1)` for the last coefficient. Every command has rank at most3:
R lands on t=0, N on delta=0, and J_A for A>=1 lands on
`h'-q'=2^A s'`. Different commands land on different charts, so this does
not imply a single three-dimensional linear observer for all words.

Indeed prefixes `empty,R,N,J1` from rooted C3, tested against contexts
`empty,R,RR,NR` (read chronologically), have Hamilton-response matrix

    [[3, 5,  9,  9],
     [5, 9, 17, 15],
     [5, 9, 15, 17],
     [9,23, 45, 41]],    determinant -32.

Thus dimension4 is the minimum rational linear dimension for this language
when J1 is allowed. The cap1 C3 language remains rank3. This is an observer
minimum, not an assertion about nonlinear parametrizations or source images.

## B3. The hostile survives the exact compression

Start from all four rooted versions of the strong order-four tournament.
Complete bounded saturation under R and every allowed J_A still omits
39 and prime67. The receipt includes capped2, capped3 and **unbounded A**
versions through199. Strict H growth makes these permanent omissions in
the specified language: every derivation of67 stays below67, and the
block formula bounds every possible A at each state. There is no silent
vertex-order cutoff or restriction to a finite list of A values.

Both values occur in the broader adaptive language:

    H39: roots (0,3,0,0),
    H67: roots (0,1,3,4).

The first explicitly returns to the original0 after using3. The second
needs a previously retained vertex outside the most recent fresh block.
The missing coordinate is an older name eligible for return. Adding a
new marked block is now fully controlled, but it cannot restore a name
that the grammar has discarded. The two-live-port language in
[LIVE_PORTS.md](LIVE_PORTS.md) keeps precisely this extra option.

The procedural path was: test whether a fixed-root block restores a small
observer; derive all marked coordinates; expose its finite alphabet and
positive charts; then challenge actual coverage at the first old-root
targets. The fit suggested the formula; the bin decomposition proves it.
No unrestricted four-coordinate closure is asserted.

The producer checks850 rooted block profiles using85 rooted presentations
from all ordinary source classes of orders1..5 and block sizes1..4. The
independent audit uses exact subset/end-path recursion on a separate universe, including a regular
five-core and the adaptive H37 hostile. Reproduce `python3 -I -B
block_switch.py` and its `-O` version; stdout equals block_switch.json.
