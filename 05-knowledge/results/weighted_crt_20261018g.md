# Weighted compatible clocks: the cross-term price and the hidden cycle space

**PROVED:** the weighted conditional-projection criterion, rational cross-term
bound, exact energy identity, cycle-space decoder, and all-length native
hostile. **FINITE-EXACT:** five specified original source-charge tables and
the controls below. **OPEN:** a uniform native cancellation estimate, the
unpaid interaction/arrival debt, and universal ROOT entry. No novelty claim
is made for conditional expectation, Cauchy--Schwarz, or graph cycle spaces.

The useful change is precise: compatible mod15/mod35 registers have **45
marginal-visible coordinates and60 joint interaction coordinates** on the
full105-cell space. The latter are five copies of the cycle space of K3,7.
Small marginal coupling does not remove them. Even an exactly independent
four-cell measure supports a signed native packet that cancels nothing.

## 1. Inheritance and the connection being tested

The nearest proved mechanism is [W1, weighted reflection](105_research_atlas_20261018f.md),
which requires vanishing off-diagonal weighted blocks. Incoming
[P4, nested subgroup layers](digit_125_nested_105_20261018q.md) supplies
unweighted annular orthogonality and preserves a full distance profile;
its named-runner row hostile already prevents arbitrary reweighting.
The native payment obligation is the actual opposite-sign terminal overlap
in [paired energy, equations10--11](collatz_paired_energy_20261015g.md).
[Refinement debt](collatz_refinement_debt_20261016a.md) retains physical rows,
terminal vectors and a common native map; its Gram remainder is not a
replacement for the unpaid coarse term. The [source-clock producer](../../04-computation/experiments/root_clock_channels_20261018o.py)
provides the literal charge used below, not an independence hypothesis.

Corrected near miss: CRT compatibility is not independence of weighted
registers. The canonical hostile below strengthens that warning: exact
independence is still not native terminal cancellation. The least-used
sidecar is the circulation remaining after both marginal projections vanish.

Anchor: source-weighted Collatz payment. Niche: bipartite cycle flows.
Wildcard: what survives the mod15/mod35 compression. The live concepts are
original source, positive measure, conditional projection, interaction flow,
signed endpoint and literal clock. No metric realization follows from the
bipartite support graph.

## 2. A rational certificate for the weighted mixed term

Let S be a finite set of **named original sources**, with weights mu_s>=0.
Zero-weight sources are retained in the receipt but ignored in the weighted
Hilbert space. Set A(s)=s mod15, B(s)=s mod35 and C(s)=s mod5. A compatible
pair (A,B) is equivalent to s mod105. Within one C=c fibre, use i=s mod3
and j=s mod7 and put

    M_ij = sum_{s: (s mod5,s mod3,s mod7)=(c,i,j)} mu_s,
    R_i = sum_j M_ij,   C_j = sum_i M_ij,   T = sum_ij M_ij.

Here C_j denotes a column mass; it is not the register C. Skip a fibre with
T=0, and omit zero rows/columns from divisions. For T>0 define

    D_ij = M_ij - R_i C_j/T,
    Delta = sum_ij D_ij^2/(R_i C_j)
          = sum_ij M_ij^2/(R_i C_j) - 1.                 (1)

The second identity follows by expansion. All quantities are rational when
the original weights are rational; no eigenvalue approximation is needed.

**W2 (PROVED).** For real or Hilbert-valued functions f_i,g_j satisfying
sum_i R_i f_i=0 and sum_j C_j g_j=0, let

    E_f=sum_i R_i ||f_i||^2,   E_g=sum_j C_j ||g_j||^2,
    K=sum_ij M_ij <f_i,g_j>.

Then

    |K|^2 <= Delta E_f E_g.                              (2)

Proof: the product term in M vanishes by centering. Apply Cauchy--Schwarz
to the arrays D_ij/sqrt(R_i C_j) and
sqrt(R_i C_j)<f_i,g_j>. Their second squared norm is at most E_f E_g.
This proves (2), including every zero-energy case.

Consequently replacing the **mixed term only** by its conditional-product
value changes ||f+g||_M^2 by at most2 sqrt(Delta E_f E_g). Flipping one
detail changes this energy by at most4 sqrt(Delta E_f E_g). These are
bounds on the original weighted form; they do not authorize replacing
its source weights. Sum the bounds over the five c-fibres, or use their
maximum Delta and another Cauchy--Schwarz inequality.

The sharp operator constant is the squared spectral norm of the matrix
D_ij/sqrt(R_i C_j); Delta is its squared Frobenius norm, a convenient exact
upper certificate. If a sharper rational certificate is needed, a proposed
lambda>=0 is valid whenever

    D diag(C)^(-1) D^T <= lambda diag(R).                (3)

This is a rational positive-semidefinite matrix test, not a fitted SVD.

**Sharp positive control.** Take f=(1,-1,0), g=(1,-1,0,0,0,0,0), and
M_ij=2+f_i g_j. All entries are1,2 or3. The row masses are14 and column
masses6; Delta=1/21, K=4, E_f=28 and E_g=12. Thus equality holds in (2).

**Exact zero criterion.** Let P15,P35,P5 be actual weighted conditional
expectations on the three registers. Then

    P15 P35 = P5  iff every active M table has rank one
                    iff every Delta is zero.           (4)

Indeed the left side applied to every function of B says that its
conditional law in each A-row is the common C-fibre law. This is exactly
M_ij=R_i C_j/T. The converse follows by direct averaging. This proof also
works when some rows or columns have zero mass. Merely saying the two
projections **commute** is insufficient with missing cells: on the two
diagonal cells of a2x2 table both projections are the identity, while P5
is the mean projection and Delta=1.

## 3. The exact energy decomposition keeps the unpaid part

For any Hilbert-valued observable F on the same physical sources put

    m=P5 F,  x=P15 F-m,  y=P35 F-m,
    z=F-m-x-y.

**W3 (PROVED).**

    ||F||^2 = ||m||^2+||x||^2+||y||^2+||z||^2-2<x,y>.   (5)

Proof: conditional expectation is an orthogonal projection. Thus m is
orthogonal to x,y,z; also <F,x>=||x||^2 and <F,y>=||y||^2. Therefore
<z,x>=<z,y>=-<x,y>, and expansion proves (5). Bound (2) prices only its
last term. If Delta=0, z is orthogonal to both marginal detail spaces.
For nonzero Delta, z need not have zero marginal projections: calling it
an orthogonal interaction would be false.

The squared norm ||z||^2 remains in (5), together with the mean and the
two detail energies. If several physical sources occupy one105-cell,
there is also an orthogonal within-cell zero-mean component. Replacing
F by its cell mean loses that component; the receipt does not do so.

### The hidden component has a cycle-space decoder

For one c-fibre form the bipartite support graph: active i-rows on the
left, active j-columns on the right, edge ij iff M_ij>0. A cell function h
has P15h=P35h=0 exactly when

    J_ij=M_ij h_ij

has zero sum on every row and every column. Orient edges from rows to
columns: this is precisely a real circulation. Multiplication by the
positive edge masses is invertible, so the joint kernel has dimension

    number_of_edges - number_of_vertices + number_of_components. (6)

A spanning forest supplies a complete fundamental-cycle decoder. On full
K3,7 a particularly simple basis consists of the12 rectangles

    J_00=J_ij=1,  J_0j=J_i0=-1,
    i in{1,2}, j in{1,...,6}.

Every basis vector has its own interior entry ij, so they are independent;
equation (6) proves completeness. There are12 hidden and9 marginal-visible
cell coordinates per mod5 fibre. Thus the full105-cell space has60 hidden
and45 visible dimensions. This statement is valid for arbitrary positive
cell masses; the functions are J_ij/M_ij, not unweighted signs unless the
masses agree. Missing edges use (6), not the full-support count.

For p positive-weight physical sources occupying e cells, the additional
within-cell kernel has dimension p-e. In the five tested tables e=105.
The bipartite graph here records compatible observations; it is neither
the Collatz trajectory graph nor a planar unit-distance graph.

## 4. Zero coupling still permits zero native cancellation

The following one-row signed packet uses actual positive odd sources:

| Sign | Source | mod15 | mod35 | Two actual odd-map steps |
|---|---:|---:|---:|---|
| + | 231 | 6 | 21 | 231 ->347 ->521 |
| - | 351 | 6 | 1 | 351 ->527 ->791 |
| - | 511 | 1 | 21 | 511 ->767 ->1151 |
| + | 631 | 1 | 1 | 631 ->947 ->1421 |

Every source is1 modulo5. With source weights1, the active2x2 M table
is all ones, so Delta=0. Both signed marginals vanish. The packet is
exactly the four-cycle interaction; its source squared energy is4.
Every actual valuation letter is1, both before and after the first step.
The four endpoints remain distinct and grow. Thus native l1 norm remains4.
The common shortcut clock is2; the unaccelerated clock is4. There is no
earlier ROOT hit.

**W4 (PROVED all-length hostile).** For every L>=1 this construction can
be made with an actual length-L all-one odd word. Choose one positive
source at each of the same four CRT corners and additionally impose

    n=-1 mod2^(L+1).

CRT applies because105 is odd. Along the first L odd steps,

    U^j(n)+1=3^j(n+1)/2^j,

whose2-adic valuation is at least2 before each of those steps. Therefore
each exact letter is1 and every value increases. The endpoint is
3^L(n+1)/2^L-1, an injective affine function of n. Both signed marginals
remain zero at the source but no terminal coefficients cancel.

This refutes a universal discount for **arbitrary signed native packets**
based only on compatible registers, marginal cancellation and Delta=0.
It does not refute a theorem for the particular annular ROOT coefficients.
Such a theorem must exploit additional structure of those coefficients
and their actual signed arrival fibres.

## 5. Five original ROOT-charge measures, without an arrival oracle

For each requested N, retain every canonical source

    1 <= s <4N, s odd, s mod8 !=5.

The inherited producer defines a three-real-row charge by summing over
m=s,4s+1,4(4s+1)+1,... below4N. Its residue weight vectors are
(0,-3,3) below2N and(-1,-1,2) above. Writing the resulting row as
(2a-b,-a-b,-a+2b) gives A_N(s)=(a,b) and

    mu_s=a^2-ab+b^2 in{0,1,3}.

Only `charge` and `norm` are imported. No route function, discovered
terminal, orbit cutoff, or ROOT oracle is called. A separate implementation
sums the inverted shell charges directly: below2N these are
(1,2),(-2,-1),(1,-1) by m mod3; above they are
(0,1),(-1,-1),(1,0). It agrees at **every retained source**.
The JSON records all source strings `s:a:b:mu`, their canonical hash,
all105 cell masses, each rational Delta, and the original source counts.
Zero-weight source records remain present.

| N | Retained sources | Positive-weight sources | Sum mu | Maximum Delta, exact upper bound |
|---:|---:|---:|---:|---:|
| 4096 | 6144 | 6047 | 9949 | <=1/2368 |
| 5997 | 8996 | 8852 | 14562 | <=1/2539 |
| 11002 | 16503 | 16242 | 26722 | <=1/16925 |
| 16384 | 24576 | 24186 | 39790 | <=1/38703 |
| 65536 | 98304 | 96743 | 159157 | <=1/475082 |

All five tables have full cell support. The bounds are exact rational
comparisons; no asymptotic rate is inferred from their decrease. The code
also verifies (5) on the **diagnostic observable F(s)=s**, including its
physical within-cell variance. This is a test of the source-preserving
projection algebra, explicitly not a ROOT-amplitude estimate.

The actual ROOT energy additionally uses the direction of A_N(s), signed
integer arrival clocks, equality of actual endpoint labels, and the joint
interaction. None is determined by mu_s or its small Delta. In particular,
these finite tables do not establish the open source-weighted overlap in
the paired-energy note. They price one precise marginal cross term and
identify the component which that price leaves untouched.

## 6. Exact controls and reproduction

Run from the repository root:

    python -B 04-computation/experiments/weighted_crt_20261018g.py
    python -B -O 04-computation/experiments/weighted_crt_20261018g.py --verify

The declared controls include all512 support graphs on3+3 vertices, with
independent rational incidence ranks and conditional-projection products;
the complete12-vector K3,7 rectangle basis; the sharp positive table;
zero-mass and disconnected-support boundaries; native hostile lengths
1,2,3,8,16,32,64; literal Boolean/float/negative/ragged rejection; and the
independent charge check at every source in the five scales. Tests use
active exceptions and remain enabled under optimized Python.

Normal and optimized replay agree with the saved receipt: **311,876 checks**,
LF SHA256 `77e06b6dc8c4e7fe6b75db192463e7b8495e4dcab96a4d32db78b9d01d6e4df1`.
An independent peer read the full proofs and source interfaces, then ran
the optimized saved comparison: **PASS**. The source-table consumer explicitly
accepts only the five declared scales; three malformed/out-of-universe
controls keep that finite scope visible. All files are owned by this
packet; no shared theorem or navigation file was changed.
