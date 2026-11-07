# Common right-amalgamation maximality and marked-module certificates

Status: PROVED candidate, awaiting independent audit. Finite search results are labeled separately. The unrestricted Hamiltonian-path spectrum and F(8,137) remain OPEN. No improved coverage interval or eventual threshold is claimed.

## 1. Inheritance and live concepts

The closest constructive mechanism is the Session215 sink synchronizer: matching (X,V,d), followed by O, produces matching (X,e,V,d). The hostile example is Session214's permanent separation after the first common J on each of its four fixed pairs. The least-used sidecars tested here are arbitrary rooted right factors and marked-module run counts. Matching a scalar X alone is not persistent synchronization. The concept board comprises common right products, asymmetric right products, marked-module substitutions, ordered run covers, and source-bank replacement.

The inherited exact arbitrary right-product law is Section4 of frozen `dependencies/ROW_NORMALIZATION.md`. Its quantities are B=h+L-s, V=R+b, and d=H(U-root). A rooted right input W of at least two vertices, with q=d_W, acts as

B'=d*h_W+(V-d)*L_W+(e-d)*s_W+(B-e-V+d)*q,
e'=d*e_W,
V'=b_W*V+d*(R_W-b_W),
d'=d*q, s'=s*q.

The two rooted inputs identify their roots and all left nonroot vertices beat all right nonroot vertices. The singleton right input is the identity. All formulas concern physical tournaments, not freely assigned integers.

## 2. Exact maximality theorem for a common right factor

Let A and B be any two persistent two-residue carriers: their (X,e,V,d) agree, their starts are s_A=1+3j_A and s_B=2+3j_B, and X=B_A+2j_A=B_B+2j_B. No assumption about which source bank supplied them is made.

Apply the SAME arbitrary physical rooted right factor W to both, and reassign the residue labels and j values from the new starts whenever these residues are nonzero. Then the resulting adjusted bases are equal if and only if q=H(W-root)=1. This remains necessary even if the j bound is dropped. For q divisible by3 the starts both have residue0 and hence cannot supply the required residue pair. For other q, let r'_A=(q mod3) and r'_B=(2q mod3), as elements of {1,2}. Then exactly

X'_A-X'_B = 2(-q-r'_A+r'_B)/3.

Thus the difference is 2(1-q)/3 for q=1 mod3, and -2(q+1)/3 for q=2 mod3. A subsequent common O/Z word cannot repair this difference: the new e,V,d agree, so every such continuation adds equal increments to the two X coordinates.

Proof. Write delta for A minus B. Initially delta B=-2 delta j and delta s=3 delta j-1. The displayed law gives delta B'=q delta B because e,V,d agree. Reassigning j gives 3 delta j'=q delta s-r'_A+r'_B. Substitution proves the difference formula. Positivity of q then proves the iff assertion. The propagated coordinates e,V,d agree, proving the continuation statement. Equivalently, the two physical four-core readouts, initially separated by2, become separated by exactly2q under this common right multiplication. This is a real consequence of the inherited law; it is not an additional independent counting formula.

Furthermore q=1 right factors are EXACTLY common O/Z words, up to relabeling of newly adjoined vertices. A tournament has one Hamiltonian path exactly when it is transitive: this follows, for example, by strong-component factorization and the fact that a nontrivial strong tournament has at least three Hamilton paths (or the inherited minimal strong count theorem). Thus W-root is transitive. List its vertices in their transitive order and adjoin them successively. Every old nonroot vertex beats the next one. If the root also beats it use O; otherwise use Z. This produces precisely the ordered root product with W. Conversely every O/Z word leaves the new nonroot vertices transitive, so its right factor has q=1.

For completeness, the cited strong-component argument is: the condensation of a tournament is a transitive tournament, and Hamilton paths concatenate independent Hamilton paths of its strong components. Camion's theorem gives a Hamilton cycle in any nontrivial strong component; cutting each of its at least three arcs gives at least three paths. Hence H=1 forces every strong component to be a singleton. Camion's theorem is already a cited input in repository `01-canon/theorems/THM-343-H7-impossible.md` (A10); only that standard cited fact, not its H=7 argument, is used here.

Scope: this classifies common right factors on an already synchronized pair. It does not classify asymmetric right factors, root changes, different initial bases, or replacement of either source.

## 3. Exact marked-module interface

Replace a NONROOT vertex v in a rooted tournament by a nonempty module M, preserving all outside incidence. Let c_k(M) be the number of ordered k-tuples of nonempty directed paths partitioning M. Every designated path/cover decomposes uniquely into its maximal nonempty visits to M. Collapsing each visit to a placeholder gives a run skeleton. Consecutive placeholders in the same path are forbidden; external vertices occur once in total; at most n-1+q placeholders occur in a q-path cover. Thus every required profile coordinate is an exact finite nonnegative integer linear combination of the c_k(M). Empty designated rooted paths are forbidden. Empty ordinary bins are allowed.

`probe.py` computes these coefficients by subset DP. For one path its state is (external subset, last external vertex or placeholder, number of placeholders). For two ordered paths, convolve their external subsets and placeholder counts. Filtering starts or ends at the root supplies s,e,L,R. Omitting the root supplies d,b. A fixed ordering of bins and traversals assigns the k module paths to the skeleton's k visits, so no division by k! occurs.

For a transitive m-vertex module, c_k=k! S(m,k), where S is a Stirling number of the second kind. These functions are linearly independent as functions of positive integer m, by their triangular evaluations at m=1,2,... . In particular c1=1 and c2=2^m-2. These are NOT ordinary polynomials in m. Replacing the marked vertex by a transitive chain preserves a fixed start count for every m precisely when all its coefficients for k>=2 vanish. Therefore coefficient matching is a proof for all module sizes, not interpolation from a numerical cutoff.

## 4. Complete finite marked-position result for six fixed pairs

FINITE-EXACT source scope: the four Session203 pairs, plus the two Session215 pairs at4225 and4357 after their stated JO or ZJO prefixes. Every nonroot position in each member is tested, allowing different positions on the two sides. The retained conditions are: start stays constant and in its original allowed residue/budget, and (X,e,V,d) agrees for every transitive module size. Precisely four ordered position pairs pass. In the root-first labeling retained in `matches.json`, they are:

-759, positions1 and1: (X,e,V,d)=(759+449c2,29,927+449c2,29).
-821, positions1 and1: (821+523c2,33,1079+523c2,33).
-4225, positions9 and9: (4225,0,8910+4455c2,195).
-4357, positions9 and9: (4357,0,9836+4918c2,251).

These are NOT new infinite scalar families. In the first two cases the marked vertex beats only the root, so replacing it by a transitive m-module is exactly appending m-1 copies of Z. In the last two cases it is a global sink, so the replacement is appending m-1 copies of O. Thus this complete finite marked-position probe gives no escape from the old binary language. This finite conclusion does not constrain other source pairs, simultaneous substitution at multiple vertices, nontransitive modules, or changing the root.

`verify.py` directly recounts all eight selected substitutions for modules of sizes1,2,3 and a directed triangle. The latter checks the arbitrary-module interface beyond the transitive family. Arithmetic controls separately verify the exact right-factor gap equation; neither finite control replaces the proofs above.

## 5. Asymmetric right-factor exploration

FINITE-EXACT search only. `BANK.txt` is the unchanged Session215 supplied bank; no class-completeness is assumed. There are12,893 bounded-start rooted sources and420 rooted right inputs through order6. The search retains q>1, s*q<=26, and s*q nonzero modulo3, giving266,768 right products. Hashing their exact (X,e,V,d) finds no two-residue matches, even when the original unmultiplied source bank is available as an alternative counterpart. Hashing after a sink reset likewise finds none. Counts are in the named text receipts. Extending each left source by every O/Z/J word of length at most2 gives167,609 presentations. The same nontrivial-right filter yields3,467,984 candidates and1,588,044 sink-reset keys after also inserting the unmultiplied prefix bank. There are still zero matches involving a nontrivial right product. This does NOT prove nonexistence beyond this finite bank and operation menu.

Reproduction: run `python3 probe.py`; run `python3 verify.py` and `python3 -O verify.py` and compare receipts; compile the C++ searches with `g++ -O2 -std=c++17` and pass BANK.txt. No published or inherited file is modified.
