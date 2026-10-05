# A rooted refuel tree with one finite base proof

2026-10-05. **PROVED** for the explicitly constructed family and its exact
membership decoder. **FINITE-EXACT** for the controls below. This is an
all-height closure theorem from one checked seed, not a proof that every
positive integer belongs to the family.

The base is the literal first-hit route \(3\to5\to1\). A single guarded
constructor takes a previously certified member to infinitely many larger
members. Its inverse recovers the unique smaller parent and the chosen
branch. Consequently its apparent recursive proof loop has a decreasing
integer parameter; it is not an unsupported circular ROOT assumption.

## 1. Inheritance and the question being answered

The closest proved mechanism is the common-future dependency
\[
 H(n)=\frac{729n+669}{1024},\qquad n\equiv155\pmod{2048},
\]
from [collatz_recursive_dependency_kernel_20261004.md](collatz_recursive_dependency_kernel_20261004.md),
sections 2–4. That note already transports a supplied terminal proof,
repeats H while binary fuel lasts, and constructs completed families at
every requested repeat depth. H is not an actual Collatz edge.

The other inherited operation is \(S(u)=4u+1\), which preserves the next
accelerated odd successor. Its complete inverse-ray addresses and retained
ROOT proofs are in
[inverse_ray_ternary_addresses_20261004.md](inverse_ray_ternary_addresses_20261004.md).
The general guard-refuelling principle is already proved in
[collatz_affine_guarded_lifts_20261004.md](collatz_affine_guarded_lifts_20261004.md),
section 6. No novelty is claimed for inverse Collatz construction.

The inherited hostile is a fixed inverse block which exhausts its ternary
fuel. The precise pumping obstruction is
[collatz_guarded_pumping_memory_20261004.md](collatz_guarded_pumping_memory_20261004.md):
every regular language of actual first-hit inverse words has bounded
odd-step depth. The corrected near miss is that an available decreasing
dependency is not a root proof without a grounded terminal. The least-used
sidecar here is the **positive refuel exponent**, recovered from a valuation
rather than discarded as a finite residue. The working board is finite
base / recursive parent / source phase / ternary fuel / first-hit suffix /
source identity.

The present construction uses a new refuel at each generation. It gives
a particularly simple injective rooted family and a total membership
algorithm, with every logical base assumption explicitly discharged.

## 2. One universally guarded constructor

Let \(U(n)=(3n+1)/2^{v_2(3n+1)}\). Work on positive parents
\(u\equiv3\pmod8\); this includes the seed 3 and all constructed children.
In particular \(u\ge3\) and \(v_2(3u+1)=1\).

There is exactly one \(\kappa(u)\in\{1,\ldots,729\}\) such that
\[
 4^{\kappa(u)}(3u+1)\equiv334\pmod{2187}.                 \tag{1}
\]
Indeed \(4\) generates the principal units modulo \(3^7\), with order
\(3^6=729\), and both \(334\) and \(3u+1\) are principal units. For
completeness, cubing shows
\(v_3(4^{3^j}-1)=j+1\); this gives the claimed order and, by counting,
the whole principal-unit subgroup. Only \(u\bmod729\) is needed for (1).
The positive representative 729 is used when the residue is zero.

For any branch \(t\ge0\), put
\[
 k=\kappa(u)+729t,\qquad
 c=S^k(u)=\frac{4^k(3u+1)-1}{3},\qquad
 n=\frac{1024c-669}{729}
   =\frac{1024\,4^k(3u+1)-3031}{2187}.                 \tag{2}
\]
Equation (1) says \(c\equiv111\pmod{729}\), so n is integral.
Since c and 111 are odd, this is also \(c\equiv111\pmod{1458}\).
Writing \(c=111+1458j\) gives \(n=155+2048j\).
As \(c>u\ge3\), the integer j is nonnegative. Moreover
\[
 n-c=\frac{295c-669}{729}>0,\qquad u<c<n.              \tag{3}
\]
Thus n has exactly the native H guard and \(H(n)=c\).
This does not infer a root proof from a residue: the proof for u is
transported next.

Let \(V=(1,2,1,1,1,2)\). On the native guard these are six actual
valuations, with
\[
 U^6(n)=F_V(n)=\frac{729n+925}{256}=4c+1=S^{k+1}(u).
\]
For a supplied first-hit word \((a,\tau)\) of u,
\[
 \boxed{\quad (a,\tau)\quad\longmapsto\quad
              V\,(a+2k+2)\,\tau \quad}               \tag{4}
\]
is a first-hit word of n: the seventh edge has the same endpoint \(U(u)\).
Here a is always 1. Every nonempty prefix of V has multiplier greater
than 1 and positive carry, so all its states exceed n. The next state
\((3u+1)/2<c<n\) is positive and not the root. The remaining supplied
suffix is already first-hit. This proves all intermediate legality and
the absence of root padding, not merely affine equality.

In particular every nonseed in this family already has its **first
descent at the seventh accelerated odd step**. The result below is recursive
rooted closure and lossless storage; no extra generic first-descent
coverage is claimed.

The construction is specific to H. It does not work by replacing H with
an arbitrary paid rule. For example \(S^k(u)\equiv1\pmod4\) when \(k\ge1\),
whereas native L outputs lie \(11\pmod{16}\). This refuel therefore cannot
serve the native inverse L guard.

## 3. A finite proof schema with a well-founded recursive call

Define \(N(())=3\). For a finite tuple of nonnegative branch integers,
define \(N((t_1,\ldots,t_d,t))\) by applying (2) to
\(u=N((t_1,\ldots,t_d))\). Let \(\mathcal F\) be this image.

**PROVED.** Every member of \(\mathcal F\) has an actual first-hit root
certificate. Its only base proof is \(3\to5\to1\), with valuation word
\((1,4)\); its only inductive rule is (4).

This is induction on tuple length. Equivalently, it is strong induction
on the source integer, because (3) makes each recursive parent smaller.
The same proof schema can be depicted by one recursive constructor node
and a base node, but that finite diagram is justified by the retained
decreasing parameter. Deleting that parameter would delete the proof.
There is no finite list of unproved ROOT labels waiting to be trusted.

At depth d, the exact odd first-hit rank is
\[
 r_d=6d+2,
\]
the halving cost is \(5+\sum_{i=1}^d(2k_i+10)\), and the ordinary rank
(counting \(3n+1\) and every halving separately) is
\[
 7+\sum_{i=1}^d(2k_i+16).
\]
For the all-zero branch address, the first exponents \(k_i\) are
675, 60, 174, 605, 604, 723, 316, 498. The first child already has 1353
binary digits. This is an explicit family, not an efficient small-source
enumeration.

The raw D/R inverse language of these certificates is **not regular**:
its odd depth \(6d+2\) is unbounded, so the inherited first-hit pumping
theorem applies. A finite program with an unbounded integer source,
branch address, valuation, and guard calculation is compatible with
that obstruction. A finite proof schema is not a finite-state language.
Likewise the bounded-forward/inverse-depth no-go in
[collatz_partition_cover_20261004.md](collatz_partition_cover_20261004.md)
does not assert impossibility of this explicitly parameterized family.

## 4. Source-only recognition and lossless addresses

The family has a total membership decoder on supplied positive odd integers.
Accept n=3 as the base. Otherwise:

1. Require \(n\equiv155\pmod{2048}\) and compute \(c=H(n)\).
2. Let \(b=v_2(3c+1)\). Require b odd and \(b\ge3\), and set \(k=(b-1)/2\).
3. Recover
   \[
   u=\frac{(3c+1)/4^k-1}{3}.
   \]
   Require a positive integer \(u\equiv3\pmod8\), and
   \(k=\kappa(u)+729t\) for an integer \(t\ge0\).
4. Recurse on u, accepting only if the recursion eventually reaches 3.

Reject at the first failed requirement. Every retained step has
\(u<c<n\), so this algorithm terminates even on nonmembers.
For a constructed source the valuation is exactly \(2k+1\), so every
step recovers its actual parent and branch. Conversely every accepted
step is equation (2), so reaching 3 proves membership. Thus the decoder
is an **iff**, and \(N\) is injective across all addresses and depths.

Native guard membership alone is weaker: n=155 gives c=111 and b=1,
so the required positive refuel is absent. Even adding 2048 to the first
constructed source keeps its native H residue but fails this family
decoder. These rejected inputs may have other completed root proofs;
rejection makes no assertion about their convergence.

The preserved object is the supplied integer together with its exact
root certificate; the address is an alternative lossless representation.
The source-to-address map preserves the complete recursive parent chain.
Keeping only a finite residue loses the branch quotient and the source
identity. The integer t, or its unconsumed high digits, is the required
sidecar.

## 5. All ternary branch addresses, with the source gap retained

For a fixed parent u, write \(n_t\) for (2). For distinct s,t,
\[
 v_3(n_t-n_s)
 =v_3(4^{729(t-s)}-1)-7
 =v_3(t-s).                                         \tag{5}
\]
The other numerator factors in their difference are 3-units.
More precisely, \(4^{729}\equiv1+3^7\pmod{3^8}\), so for all
\(a\ge0\), integers \(t\ge0\), and \(e\in\{0,1,2\}\),
\[
 n_{t+e3^a}\equiv n_t+e3^a\pmod{3^{a+1}}.            \tag{6}
\]
The zero increment is included separately. Equation (6) follows by
binomial expansion; at higher a the first nonzero coefficient stays 1
modulo 3.

Consequently \(t\bmod3^a\mapsto n_t\bmod3^a\) is a bijection. The next
target digit is selected by one carry subtraction, as in the inherited
inverse-ray codec. Every certified parent has certified children in
every ternary residue at every finite precision. Every such child still
lies \(155\bmod2048\). This is coverage of **constructed residue
representatives**, not of each integer in those residue classes.

There is a concrete finite-observer boundary even inside the rooted
family. The parents with addresses (0) and (729) have the same residue
366 modulo729 and the same selected \(\kappa=60\). Their branch-zero
children nevertheless have residues 675 and 253 modulo729. Therefore
the observed pair \((u\bmod729,\kappa)\) is not a closed state system.

On a fixed selected-k branch, changing u changes n by
\(1024\,4^k\,\Delta u/729\); the quotient consumes six ternary digits.
The exact unexpanded reader retains parent u modulo \(729m\) to obtain
the child modulo m. It evaluates \((4^k-1)/3\) with one additional
factor 3 of modular precision before division. Iterating this reader
retains the required precision at every level.

The tuple \((10^{100},7,10^{80})\) is a small stored example with enormous
exponents. Without expanding the source, its residues modulo
(8,19,729,2048) are (3,4,99,155), its odd rank is 20, and its halving
cost has 104 decimal digits. This is a proved symbolic integer and
route, not a claim that its literal orbit was printed.

## 6. Reproduction, controls, and consequence

Run from the repository root:

    python -B 04-computation/experiments/finite_seed_refuel_tree_20261005.py
    python -B -O 04-computation/experiments/finite_seed_refuel_tree_20261005.py

The [script](../../04-computation/experiments/finite_seed_refuel_tree_20261005.py)
and [saved output](finite_seed_refuel_tree_20261005.out) contain:

- all 729 parent phases, checked independently against all 729 candidate
  exponents;
- every branch word over \(\{0,1,3\}\) of depths 0–4: 121 distinct sources,
  2798 literal odd edges, exact source recognition, and 847 modular checks;
- independent inherited ROOT-codec exports for 13 addresses;
- first-descent-at-seven checks for every nonseed in that finite universe;
- the first eight canonical branch-zero generations;
- complete ternary branch tables through precision 5 for three rooted
  parents: 1089 target-address checks, plus 108 derivative checks;
- the rooted observer collision, source-identity and native-guard hostiles,
  exact integer-type failures, and the unexpanded large-address control.

Normal and optimized execution agree; all checks survive removal of Python
assertions. The LF output SHA256 is
\( \texttt{32fa4a0abacc72a444c0ced4fd2be9a9fc5100714e6fe8621acd001d7f634524}\).

The finite proof obligation has succeeded for this family: check the base,
check one guarded rule, and retain the well-founded recursive call.
Universal Collatz would require an additional theorem routing every supplied
source into an already grounded family or another well-founded obligation.
Neither the branch residue bijections nor this infinite generated tree
supply that missing universal entry theorem.
