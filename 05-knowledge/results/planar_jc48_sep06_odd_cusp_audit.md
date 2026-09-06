# Independent audit: odd-cusp divisor passports

**Status: PASS — analytic proof, source, and normal/optimized/frozen replay.
The geometric spectrum is necessary; the displayed abstract passports do
not assert Keller realization. Local topology and low-degree exclusions
are CITED.** September 6, 2026.

Audited primary: [odd-cusp passport](planar_jc48_sep06_odd_cusp.md),
[standalone source](../../04-computation/planar_jc48_sep06_odd_cusp.py), and
[frozen output](planar_jc48_sep06_odd_cusp.out). This is an independent
reader's proof audit and replay, not an independent proof of the cited
low-degree theorems. No producer source or output was changed.

## 1. Typed conclusion and inherited geometric information

The domain is an actual nonautomorphic polynomial Keller map
`C^2 -> C^2` whose **whole** nonproperness set is irreducible, has
normalization `A1`, and has exactly one analytic `(2,m)` cusp with odd
`m>=3`, together with `N>=1` ordinary nodes and no other singularities.
The mapping degree `d`, generic actual fibre count `a`, missing degree
`delta=d-a`, cusp fibre count `n`, and cusp exponent `m` are distinct.
The inherited finite-envelope arguments supply `1<=a<d`, constant `a`
on the smooth support, pure deleted boundary, and exact node overlap
counts. These are actual geometric inputs, not consequences of an
arbitrary permutation representation.

Under these hypotheses the primary conclusion is correct:

```text
q>=3 an odd divisor of m,   n in {0,1},
d=q+n,   a=(q-1)/2+n,   delta=(q+1)/2.
```

For `N>=2`, necessarily `n=1`, so `d=q+1` and `a=delta=d/2`.
The cited mapping-degree-three/four exclusions remove `q=3`. The
two-label formal boundary is excluded separately by `delta=1` purity.
The primary's repaired table explicitly reports only the degrees left
by its cited degree-two-through-four filter; it makes no claim that
every listed degree is unresolved elsewhere in the literature.

The actual-page and purity inputs were independently checked in the
[ordinary-cusp audit](planar_jc48_sep06_cusp_passport_audit.md) and
[infinity audit](planar_jc48_sep06_infinity_audit.md). No withdrawn
stronger connectedness claim is used.

## 2. Marked meridians and the actual intersection

I reread the primary author text of Shimada's
[*Lectures on Zariski van-Kampen theorem*](https://www.math.sci.hiroshima-u.ac.jp/shimada/LectureNotes/LNZV.pdf),
Section 6, pages 24–25, in the downloaded PDF's extracted text. Its
presentation retains the fibre meridians. At `p=m=2r+1`, `q=2`, set
`sigma=ell_1`, `tau=ell_0`, and `w=sigma tau`. The relation
`ell_(2r+1)=ell_0` reads

```text
g sigma g^-1=tau,   g=w^r.
```

These are marked normal meridians generating the local complement
group. They are not arbitrary generators chosen only after passing to
an abstract torus-knot presentation. The second parity relation follows
from this one: conjugation by `g` fixes `w`, hence sends `tau` to
`tau^-1 w`; consequently
`w^(r+1) tau w^-(r+1)=sigma`. Thus the positive cycle controls also
define the full two-meridian local action.

The geometric re-access formula is correctly oriented. If `c` accesses
a smooth cusp transversal and `A_eta` denotes its retained endpoint
labels, then `A=T_c^-1 A_eta`. The declared rightmost-first convention
and the new access path `c'=c g^-1` give

```text
T_c'^-1 A_eta=gA=B,   B subset Fix(tau).
```

This re-accesses the **same smooth point**. It does not identify `g`
with a preferred path along the singular curve or assume a peripheral
centralizer property.

A common label in `A intersect B` is fixed by both local generators,
hence gives a singleton local cover component. Its finite normal
extension is a copy of the target neighborhood. Its generic cusp
point is retained; removing only its origin would contradict pure
deleted boundary. Conversely an actual cusp preimage supplies an etale
local inverse and a common label. Therefore
`n=|A intersect B|` is exact. Replacing actual retained subsets by all
inertia-fixed labels, or imposing transitivity on the local action,
would lose this conclusion.

## 3. Euler budget and complete sparse-cycle argument

The unibranch cusp changes the smooth-stratum Euler count exactly as
in the ordinary-cusp case: `chi(S_smooth)=-2N` and
`chi(C^2\S)=N`. With the inherited exact node formula, Euler integration
gives

```text
n+sum omega_p=1,
omega_p>=max(0,d-2a),
h=|Omega\(A union B)|=d-2a+n<=1.
```

For two or more nodes, positivity of the integer `d-2a` would already
spend at least two units, so `d<=2a` and `h<=n`.

I independently reconstructed the sparse permutation classification.
Write `I=A intersect B`, `U=A\B`, `V=B\A`, and `|U|=|V|=k`.
If `h=0`, both generators preserve `A`, so `w^r A=A=B`; then their
union is `A=Omega`, contradicting `a<d`. This argument needs no Artin
relation.

If `h=1`, let `o` be the outside label. A `sigma`-cycle contained in
`V` is fixed pointwise by `tau` and therefore invariant under `w` and
all its powers. It cannot be contained in `w^r U=V` while disjoint
from `U`. Hence the one cycle through `o` contains all of `V`.
The same invariant-set argument applies to `tau` and `U`. After
labelling their two cycles,

```text
sigma=(o,v_1,...,v_k),  tau=(o,u_1,...,u_k),
w=(o,u_1,...,u_k,v_1,...,v_k).
```

For `k>=1`, put `q=2k+1`. The cyclic interval `U` has a unique entry
boundary, so `w^r U=V` holds exactly when `r=k mod q`. This is
equivalent to `q|m`. Conversely the displayed cycles and this
congruence give the subset relation and the marked Artin relation;
the classification is exact as an abstract statement. It retains all
divisors, rather than only the largest possible cycle.

If `k=0`, `1<=a=|I|<=1` leaves only `d=2,a=n=1` and identity
permutations. Its missing degree is one. The weighted generic
boundary length then has a single summand with ramification and
residue degree one. The retained divisors are unramified as well.
Purity makes the finite normal envelope etale everywhere, and a
connected finite etale cover of `C^2` has degree one. This excludes
the formal boundary geometrically, not by permutation algebra.

For `N>=2`, `h=1<=n<=1` gives `n=1` and all node overlaps zero.
For `N=1,n=0`, the unique node overlap is one and its actual fibre
count is zero. Both statements agree with the exact formula; no
unproved nonemptiness of all special fibres is being assumed.

## 4. Scope, hostiles, and the remaining test

The exponent-nine controls with both `q=3` and `q=9` refute a
maximal-degree-only interpretation of the local lemma. The five
producer hostiles correctly isolate lost information: actual
re-access, the outside-label budget, the common-label budget,
orientation, and the geometric identity boundary. The `q=5,n=1`
positive control is a five-label orbit plus one common singleton;
local transitivity would incorrectly remove it.

For an actual whole support with one `(2,5)` cusp and at least two
nodes, this argument leaves the necessary passport
`d=6,a=delta=3,n=1`, with zero actual fibres at the nodes. It supplies
neither a global complement action nor a Keller realization. The
suggested curve `(t^4+t^2,t^6+t^5+t^2)` is a separate curve-classification
target supplied by the geometry lane, not a classification proved by
this audit or permutation source. The next missing condition concerns
global realization of the degree-six local data.

The low-degree inputs remain CITED with their actual meaning: generic
preimage degree, not polynomial coordinate degree. The primary
Orevkov and Domrina statements were already independently read for the
ordinary-cusp audit; their full original proofs are not re-audited here.
The ambient planar Jacobian conjecture remains OPEN.

## 5. Source and replay

I read the entire standard-library producer. Its first universe is all
`400` independent supported permutation pairs for `k=1,2` and
`r=0,...,9`, without an Artin prefilter. It finds `11` re-access
survivors. Its second universe has `882` canonical rows, with odd
`3<=q<=19`, odd `3<=m<=99`, and `n=0,1`; exactly `116` satisfy the
divisor condition. Tuple composition is compared with cyclic-position
translation. All named controls and the declared equality cases were
checked. These finite controls support, but do not replace, the
all-size proof.

From the worktree root, I independently ran

```sh
python3 -B 04-computation/planar_jc48_sep06_odd_cusp.py
python3 -B -O 04-computation/planar_jc48_sep06_odd_cusp.py
```

Both passed **4,072 explicit gates**. A separate subprocess comparison
confirmed that both complete stdout byte strings equal the frozen
output exactly (`2,657` bytes). Recomputed pins:

```text
source SHA-256:
91f487bac1217da20f235ec0746fd7d587b8d5fa3e4ce8f304d19472df790879
output SHA-256:
18368ecf7ee29d02c11d5daca31d1879582ac5d5bd30c1aecc19e6191aac3cbc
semantic SHA-256:
20873175031fe168cff89a75e43a8fc7f2c1068a64c8a888e3246846702526dd
```

Final acceptance: **PASS**. The scope clarification about the limited
low-degree filter was applied and reread. No mathematical correction
or source/output change was required. This audit is frozen for the
root's integration.
