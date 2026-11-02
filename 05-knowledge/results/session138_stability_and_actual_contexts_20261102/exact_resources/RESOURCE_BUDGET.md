# Bounded actual-source budgets fail even profinite universality

**Status: PROVED corollary of the inherited actual-source normal form.**
This is a resource obstruction, not an omitted residue for the unrestricted
spectrum and not a restriction on substitution width alone. The observation
was supplied by the parent Session138 lane and sharpened here to explicit
budget and modulus bounds.

## 1. Fixed finite affine-dyadic banks

Let A be any subset of the integers contained in J exceptional points and K
rays

    {a_i*2^k+d_i : k>=k_i},       1<=i<=K,

where a_i,d_i,k_i are fixed integers and k_i>=0. For every r>=1 put
`M_r=2^r-1`. Since 2^r=1 modulo M_r, each ray has at most r residues. Hence

    |A mod M_r| <= J+Kr.                                  (1)

There is an effective r with J+Kr<2^r-1. At that odd modulus A misses a
residue, and consequently A is not profinitely universal among odd integers.
Parity imposes no residue restriction modulo an odd modulus: every residue
has an odd integer representative. Negative offsets or ray coefficients cause
no issue. Formula(1) is a cardinality upper bound, not an assumption that all
r residues are distinct.

The inherited Session137 constructors `BOUNDED_INTERIOR.md` expresses the
H-values of every fixed finite set of actual source budgets `(h,L)` as exactly
such a finite union, with actual bounded-order bases. Therefore such a source
bank is **not even profinitely universal**, stronger than its O(log X)
integer density-zero consequence. In particular the Session137 width-three
all-modulus controller must escape every bounded `(h,L)` bank when viewed
through these actual-source coordinates. Bounded substitution width does not
bound these source budgets.

## 2. Explicit modulus from height and interior bounds

Fix H>=1 and I>=0, and allow every actual source with parent height h<=H and
interior coordinate L<=I. The proof of the inherited normal form retains all
nonsingleton SCCs. Replacing its coarse core-order bound h by the inherited
Session133 exact budget K(h) gives a uniform base-order bound

    N(H,0)=K(H)+1,
    N(H,I)=K(H)+3+floor(log_2 I),          I>=1,          (2)

where

    K(H)=max_(b=0,1,2;5^b<=H)
              [4b+3 floor(log_3(H/5^b))]
         <=3 log_3 H.

Enumerate all labelled parent tournaments of orders1..N and all their cuts.
The number of this deliberate overcount of possible actual bases is

    C_N=sum_(n=1)^N 2^(n(n+1)/2).                         (3)

Each base supplies at most one exceptional native output and one affine
positive-tail ray. Thus the entire bounded-budget source class satisfies

    |H-values mod (2^r-1)| <= C_N(r+1).                  (4)

Choose any integer r with

    C_N(r+1)<2^r-1.                                     (5)

This is an explicit all-order modular obstruction obtained by a finite
integer loop. No source enumeration needs to be completed to certify it;
C_N is an upper bound on the whole base universe. The missing residue can be
found by enumerating the actual bases and their r-term residue orbits, if a
literal forbidden class rather than its certified existence is desired.

For example H=1,I=0 gives N=1,C_N=2; r=4 already gives2*5<15.
The full source class at this budget has scalar outputs
`{1} union {2^k+1:k>=1}`, so modulo15 it occupies only
`{1,2,3,5,9}`. The counting bound does not claim to be sharp.

## 3. A quantitative local-resource lower bound

Since `C_N<=N*2^(N(N+1)/2)`, universality at M_r requires

    N(N+1)/2+log_2 N+log_2(r+1) >= log_2(2^r-1).        (6)

Consequently any sequence of such budgets universal at M_r as r tends to
infinity has

    N(H_r,I_r) >= (sqrt(2)-o(1))*sqrt(r).                (7)

For I_r>=1, a convenient necessary resource inequality is therefore

    3 log_3 H_r + log_2 I_r
        >= (sqrt(2)-o(1))*sqrt(r).                      (8)

For I_r=0 replace the left side by3 log_3 H_r. More generally if
`log H_r + log(I_r+1)=o(sqrt(r))`, the corresponding actual-source bank cannot
be universal at2^r-1 for all sufficiently large r.

This does not assert that the resources in(7)-(8) suffice, nor that the full
odd spectrum is incomplete. It precisely separates local scalar coverage
from bounded actual-source resources: the width-three controller's local
universality must pay unbounded height or interior cost, even though its
quotient width stays3.

## 4. Verification boundary

`budget.py` independently computes K(H) by both unrestricted 3/4-block search
and its three-candidate form, uses exact integer powers for all logarithmic
floors, certifies explicit obstructing moduli for sample growing budgets, and
checks the stated H1,L0 residue language. Its numerical instances corroborate
(2)-(5); the proof is the inherited actual graph normal form plus the exact
period2^r=1 modulo2^r-1. It does not import or infer any all-order source class
from a finite census.
