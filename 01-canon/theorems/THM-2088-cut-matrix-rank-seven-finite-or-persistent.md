---
id: THM-2088
title: "Rank-seven cut matrix: finite terminal or persistent two-parameter template"
status: >
  PROVED. In THM-2087's no-pair branch, choose one height-57 relation on
  each edge of the forced complete bipartite cut. Every cut forest gives
  independent coefficient rows, so the matrix rank is six or seven. Rank
  seven makes the primitive last-guard/terminal tuple a maximal-minor vector
  with every coordinate at most 91421508108581. Rank six leaves exactly a
  persistent two-parameter cut template; every chord relation is identically
  satisfied on it. This is a finite-or-persistent reduction, not LRC(14).
source: codex-2026-07-22-LRC14-cut-matrix-rank
depends_on:
  - THM-2087
related:
  - THM-2052
  - THM-2062
  - THM-2065
  - THM-2069
  - THM-2086
---

# THM-2088 -- the cut matrix is finite or persistent

Assume the no-two-term-relation branch of THM-2087. Thus a nontrivial
partition

```text
Q=A disjoint-union B,             |Q|=7,                (1)
```

has the property that for every `i in A`, `j in B` one can choose a row

```text
rho_ij=(a_ij,b_ij,c_ij)
```

encoding

```text
a_ij h+b_ij q_i+c_ij q_j=0,
b_ij c_ij!=0,
max(|a_ij|,|b_ij|,|c_ij|)<=57.                         (2)
```

Regard `rho_ij` as a row of `Z^(1+7)`, supported on the guard coordinate and
the two endpoint-speed coordinates. Let `M_(A,B)` be the matrix of all
`|A||B|` chosen rows.

## 1. Forest rows are independent

If `F` is any forest in the bipartite graph `K_(A,B)`, then

```text
{rho_e:e in F} is linearly independent over Q.          (3)
```

Indeed, in a putative linear dependence choose a leaf speed vertex `q_i` of
`F`. Its coordinate occurs in only its incident edge row, with nonzero
coefficient by (2), so that row's multiplier is zero. Delete the leaf and
continue. The guard coordinate never enters the induction. QED.

In particular, any spanning tree of `K_(A,B)` has six edges and proves

```text
rank_Q M_(A,B)>=6.                                     (4)
```

On the other hand the positive vector

```text
v=(h,(q_i)_(i in Q)) in Z_(>0)^8                      (5)
```

lies in the kernel, so

```text
rank_Q M_(A,B)<=7.                                     (6)
```

Consequently the cut matrix has exactly one of the ranks

```text
6 or 7.                                                (7)
```

This proof upgrades THM-2087's particular double-star elimination: every cut
spanning tree is a valid six-row basis candidate.

## 2. Rank seven is an explicit finite terminal

Suppose `rank M_(A,B)=7` and select seven independent original cut rows. They
form a `7 x 8` integer matrix `K`. Its kernel is one-dimensional and contains
`v`. Since the terminal core `Q` is primitive,

```text
gcd(h,q_1,...,q_7)=1,                                  (8)
```

so `v` is the primitive kernel generator. Hence its coordinates are the
signed maximal minors of `K`, divided by their common gcd.

Every row of every `7 x 7` maximal-minor matrix has at most three nonzero
entries, each of magnitude at most `57`. Its Euclidean norm is at most
`sqrt(3)*57`. Hadamard therefore gives

```text
max(h,max Q)
 <=floor((sqrt(3)*57)^7)
 =floor(sqrt(3^7*57^14))
 =91421508108581.                                      (9)
```

There are only finitely many coefficient matrices and finitely many primitive
positive terminal tuples in this box. Thus rank seven is a genuine finite
terminal. The bound is far too large for the existing THM-2078 height-24
replay; (9) is a structural finiteness result, not a completed enumeration.

## 3. Rank six is exactly the persistent cut branch

Suppose `rank M_(A,B)=6`. Choose the six rows of any cut spanning tree `T`
and put

```text
V_T=ker_Q{rho_e:e in T} subset Q^8.                    (10)
```

By Section 1, `dim V_T=2`. Every chord row belongs to the six-dimensional row
span, and therefore

```text
rho_e.u=0             for every chord e and every u in V_T. (11)
```

Thus all cut relations persist identically on the same two-parameter template
`V_T`. Conversely, if every chord relation persists on `V_T`, the full matrix
has rank six. Hence

```text
rank six
 iff every fundamental-cycle chord is persistent on V_T. (12)
```

After saturating the lattice `V_T intersect Z^8`, this is precisely a
two-anchor coefficient-row family for the terminal block plus last guard.
When the cut has type `1+6`, its complete bipartite graph is already a tree,
so rank six is automatic. For cut types `2+5` and `3+4`, the respective
numbers of chord tests are

```text
|E|-|T|=4 and 6.                                       (13)
```

One rank-raising chord freezes the terminal tuple by Section 2; otherwise all
of them are persistent template identities.

## 4. Complete rank-seven containment trichotomy

Combining THM-2087 with Sections 2--3, every rank-seven terminal containment
lies in one of three explicit branches:

```text
I.   bounded guard ratio q/h=r/s, coprime r,s<=57;
II.  rank-seven cut matrix, max(h,max Q)<=91421508108581;
III. saturated rank-six persistent cut template.       (14)
```

In the dyadic application `h` is odd, so the denominator `s` in Branch I is
odd. Branch II is finite in principle. Branch III is the only unbounded
no-pair residual produced by the short-relation cut.

This still does not parameterize the earlier dyadic guards or the two original
seam tails. To apply THM-2062/2065 to the full thirteen-speed row, Branch III
must be spliced to those five outer coordinates, preserving owner addresses
and translated-grid incidence. THM-2086's claimed lacunary cone is a
complementary attempt to remove part of that splice.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that many short relations merely enlarge a
hyperplane ledger. On a complete cut they form a representable matroid with a
graphic core: every forest is independent. Six tree edges leave two
parameters; the first algebraically new chord kills one; all other chords are
then redundant. The correct binary observable on a chord is therefore

```text
"raises coefficient-row rank" versus "persistent on the tree template". (15)
```

Possible vertices for Tournament Analysis were runners, cut edges, spanning
trees, fundamental cycles, and proof obligations. Runners lose the row
coefficients; cut edges retain them. Orienting chord obligations by the first
rank-raising event gives a transitive scheduling tournament whose SCCs are
singletons and whose tie Hamiltonian path is arbitrary among persistent
chords. These fingerprints add no proof information. The faithful carrier is
the coefficient-row matroid over `K_(A,B)`, with the terminal positive kernel
and dyadic owner labels as sidecars. QED.
