# A left-reset paired-source mechanism

Status: proved candidate, independent audit requested. This gives persistent paired rays and a broader binary synchronization criterion. It does not establish eventual coverage or F(8,137).

Use the inherited rooted profile (h,s,e,L,R,d,C,b,T). Delete C,T to obtain the closed seven-state vector (h,s,e,L,R,d,b). Here h,d are whole/deletion Hamilton-path counts; s,e are root-start/root-end counts; L,R are two-bin covers with a distinguished root-end/root-start path; b is the deletion two-bin cover count. Empty ordinary bins are allowed.

Let Q adjoin a vertex z with root -> z and z -> every old nonroot, retaining the old root. Let P adjoin a global source retaining the root. The old O adds a sink; Z adjoins a vertex beating only the root. For a nontrivial old source the exact transfers are

O: (h,s,0,L,2R,d,2b),
Z: (h+R-e,s,d,L+b,2R-d,d,2b),
P: (h,0,e,2L,R,d,2b),
Q: (h+L-s,d,e,2L-d,R+b,d,2b).

O,Z are inherited. Reverse all edges, exchanging s/e and L/R, to obtain P,Q. Thus Q resets the start to d. The old start may be zero or divisible by 3; the reset is feasible for the width8 compiler whenever d=r+3j with r in {1,2}, 0<=j<=8.

## Reduced binary criterion

After one Q, put

X=h+L-d+2j, a=L-d, w=R+b-e, t=d-e.

The five-state vector (X,a,b,w,t) evolves independently of d under Q,Z:

Q: (X+2a,2a,2b,w+2b,t),
Z: (X+w,a+b,2b,2w-2t,0).

These identities follow by substitution in the seven-state transfers. Two physical sources with different residues, equal five-state vectors, and s=d therefore produce equal adjusted bases after every common Q/Z word. Their deletion counts may differ, unlike the earlier common O/Z criterion. Their compiled Hamilton counts are 3X+2 and 3X+4. Existence of such a fully binary pair remains open in this packet.

Pure Q requires only equal (X,a): X_k=2a*2^k+X-2a. This weaker criterion has explicit physical realizations.

## Literal persistent pair

Adjacency rows encode outgoing bits, with root vertex0. Take

A=(108,1,2,6,15,30,62),
B=(34,4,9,19,7,30).

Their full nine-profiles are respectively

(33,1,1,37,55,1,660,64,589),
(31,5,3,41,43,5,328,64,219).

Both have (X,a)=(69,36), with j=0 and j=1. In particular their deletion counts 1 and5 differ. For every k>=0, applying Q exactly k times gives

X_k=72*2^k-3,
H(compiled A_k)=216*2^k-7,
H(compiled B_k)=216*2^k-5.

This is a genuine new source mechanism, not a claim that these numerical values were previously unattained. It is a thin ray, not an interval-cover theorem.

`verify_left_reset.py` checks 444 literal/reduced transfer identities on all labeled tournaments of orders2..4 and independently constructs and counts both compiler outputs for k=0..4. `LEFT_RESET_RESULT.json` contains rows, full profiles and counts; optimized and normal Python receipts agree. `RESET_PAIRED_RAYS.json` records329 matched reset (X,a) states, comprising203 distinct eventual ray signatures after powers of2 are removed from the coefficient. They are not329 disjoint rays.

## A conserved defect for the binary problem

After an additional Z, s=e=d. Put rho=R-d, a=L-d. The reduced four-coordinate action is

Q: (h,a,rho,b) -> (h+a,2a,rho+b,2b),
Z: (h,a,rho,b) -> (h+rho,a+b,2rho,2b).

Consequently

kappa=(h-a-rho+b-3d)/2

is invariant under both moves. Moreover kappa is a nonnegative integer for every physical source with s=e=d.

Proof of nonnegativity: in the squarefree variables of the deletion tournament, let f be the path series including the empty path, E the series of paths that can precede the root (including empty), and S the series of paths that can follow it (including empty). The full-support coefficients of ES, Ef, Sf, f^2 are h,L,R,b. Hence K=h-L-R+b is the coefficient of (f-E)(f-S). It counts ordered two-path covers whose first path ends in a root-outneighbor and whose second begins in a root-inneighbor. Since s=e=d, every Hamilton path of the deletion starts in the outneighbor set and ends in the inneighbor set. Split each such path at its first consecutive outneighbor-to-inneighbor transition. This injects the d paths into the covers counted by K, giving K>=d. Rooted polynomial integrality/reciprocity gives L congruent s and R congruent e modulo2. Since h,d are odd and b even, K-d is even. Therefore kappa=(K-d)/2 is a nonnegative integer. Direct substitution proves invariance.

The adjusted base is

X=2a+rho-b+3d+2kappa+2j.

Thus sources with equal (a,rho,b) synchronize under every Q/Z word exactly when their conserved quantities 3d+2kappa+2j agree. This identifies a precise evolving-state route rather than demanding equality of the old nine-profile. The obstruction is not merely that starts or deletion counts differ: different nonnegative defects could compensate.

For the candidate deletion pairs (d_A,d_B)=(1,5),(19,5),(25,5), equality requires respectively kappa_A-kappa_B=7,-26,-37. The finite structured tests in FINITE_SEARCHES.md do not find these, but do not prove an all-length obstruction.
