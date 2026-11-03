# Exact path-cover transfer for a triangle wrapper

Status: exact elementary derivation and finite checks, not a claim of novelty. Inherited work already uses path-cover polynomials for tournament substitution, including `05-knowledge/results/path_cover_transfer_opus_S446.out`, THM-3181, and `05-knowledge/results/session26_random_bridges_20261022/session26-attachments/proof.md`. Those sources use their own normalizations. Here every cover's paths are ordered, so divide by k! to compare with unordered covers.

Let c_k(T) be the number of ordered k-tuples of disjoint nonempty directed paths covering V(T), with each path viewed as its ordered vertex sequence. Thus c_1(T)=H(T), and k! divides c_k(T). For tournament Q on vertices 1,...,r, substitute nonempty tournaments T_i, with all cross arcs oriented according to Q. Let W_Q(k_1,...,k_r) count words with k_i occurrences of i, in which every consecutive pair is a distinct forward arc of Q. Then

H(Q[T_1,...,T_r]) = sum_{k_i>=1} W_Q(k_1,...,k_r) product_i c_{k_i}(T_i).

Proof: uniquely split any Hamilton path at transitions between modules. The resulting word is counted by W_Q and the runs inside each module, in their appearance order, form an ordered path cover. Conversely, concatenate the chosen runs in the word order. No factorial factors arise: equal module letters are indistinguishable in W and the cover tuple orders their actual runs. This is an exact compositional state, not a claim H alone is sufficient. A genuinely nonuniform interface must retain endpoint information, because permitted cross arcs then depend on endpoints.

## Explicit two-vertex wrapper

Let F(T) replace one vertex of a directed triangle by T; write cross orientations T -> x -> y -> T. The only admissible words with one T-run are Txy, xyT, yTx, and with two T-runs is TxyT. There are no more. Therefore

H(F(T)) = 3 c_1(T) + c_2(T).   (1)

In particular, no general affine map H -> aH+b exists for this wrapper, indeed no scalar function of H exists. For vertices 0,...,4, enumerate pairs (01,02,03,04,12,13,14,23,24,34); bit j=0 or 1 means the earlier-to-later or later-to-earlier arc, respectively. Masks 40 and 76 both have H=15, but their ordered profiles (c_1,...,c_5) are respectively (15,78,198,240,120) and (15,90,210,240,120). Equation (1) gives 123 and 135. The wrapper preserves strong connectivity, for every nonempty T: any two module vertices connect via x,y, and all three blocks lie on a quotient directed cycle.

## Linear transfer on the whole profile

For k>=1 consider all words w on k indistinguishable T letters and one x and one y. Let b(w) be the number of adjacent pairs not in {(T,x),(x,y),(y,T)}. This includes T,T, because two adjacent T-runs must be separated by a path-cover cut. For q>=1 set

M(q,k) = sum_w binom(k+1-b(w), q-1-b(w)),

with binomial coefficient zero outside 0 <= lower <= upper. Then

c_q(F(T)) = sum_{k=1}^{|T|} M(q,k)c_k(T).   (2)

Indeed, an ordered q-cover upstairs is a single ordered vertex string with q-1 cuts. Split at cuts and module changes. These produce an ordered T-cover of k runs and a word w; every bad adjacent letter pair must be cut, while each other letter boundary is optionally cut. There are k+1 letter boundaries, yielding the displayed binomial. Conversely, a T-cover, a word and a legal selection of q-1 boundary cuts reconstruct the upstairs ordered cover uniquely. Cuts inside a T-run are forbidden by construction and would instead increase k. This establishes a concrete integer transfer matrix, with finite state dimension |T| at each step; iterating F grows the dimension by two per step. The direct-product construction does not by itself establish an affine recurrence in H, a spectrum filling theorem, or new missing values.

A simple structural constraint for any n-vertex tournament is c_n=n!, c_{n-1}=n!(n-1)/2 for n>=2: the former consists of singletons; the latter has exactly one two-vertex directed path and all remaining paths singleton. Moreover each c_k is a multiple of k!. In particular (1) implies H(F(T)) ≡ 3H(T) (mod 2), consistent with Rédei's oddness theorem, but does not offer a stronger claimed universal obstruction.

## Reproduction

Run `python3 check.py`. The script computes profiles by enumerating permutations and cut positions, verifies (1) independently using subset DP on every labelled tournament of orders 1 through 5 (1,099 cases total), and verifies all coordinates of (2) against direct permutation/cut profile enumeration for all labelled tournaments of orders 1 through 4 (75 cases total). Results are recorded in `results.txt`. The upper-triangular-pair mask convention and direct orientations are implemented in the script. The finite tests verify implementation, not the general proof.
