# Independent audit: PASS

Frozen subject: session203-positive-profile-defects/PROOF.md, verify.py, counting.py, compiler.py, and the two inherited dependencies listed in SOURCE_SHA256SUMS. Candidate bytes were preserved.

## Theoretical review

1. All elementary cover inequalities have the stated domains. The f0 overlap is exactly the root-ending spanning paths. The q0 endpoint partition retains the empty-first-bin contribution in L. The a0 injection inserts the root immediately before the first outneighbor in a deterministically chosen labeled bin, with a legal predecessor if present, and has a deletion left inverse. The w0 families are disjoint exactly on nontrivial inputs.
2. The universal forest polynomial expands correctly. Its only negative forest type consists of two paths both of type (1,1). Splitting at the first outneighbor-to-inneighbor transition creates uniquely typed (1,0), (0,1), (1,1) components. Their coefficient is +2. Concatenating the uniquely typed first two components recovers the split path, so cancellation is injective. Singleton extraction gives tau=1. The theorem is all-input.
3. Independently expanded every displayed g transfer and both gate scalar-defect identities. All are exact. Gate g positivity uses only inequalities valid on every rooted input. J/Z positivity uses a root outneighbor; every actual grammar input has one. O preserves the inductive nonnegative g cone.
4. The J lower bound is valid: in a nontrivial tournament no spanning path both starts and ends at the root, so h>=s+e and R-e>=h-e>=s. Thus g>=0 implies B_J-2B>=2s. The other nonneutral defects are positive on grammar inputs. Their parity, the q-step lower bound, and the Z-ray formula and invariants are correct.
5. The improved K1 estimate is L_out=(C-2L-b)+2T>=2T>=2b+2s>=2^n. Its vertex cost is two, giving the same floor(log2 B)+2 bound as the other three nonneutral operations. Removing terminal O operations is valid on the physical word and preserves the unique gate for start two. The proof does not delete a non-grammar source, assert inverse existence at an arbitrary target, or generalize g>=0 to every positive-start tournament.
6. Singleton Z and K1 start exceptions are excluded from grammar inputs and are not silently used. All-input gate g statements remain valid at the singleton, independently checked below.

## Code and finite controls

- The compiler differs from the frozen predecessor only in its docstring and order cap, from max(2,bound.bit_length()+2) to max(2,bound.bit_length()+1). No unadvertised algorithmic change.
- Reviewed verify.py: direct rooted counts, forest cancellation including collision and left-inverse tests, literal children, identities, ray checks, and cap-support comparisons have appropriate domains. Its requires do not disappear under Python optimization.
- Reran verify.py normally and with python3 -O. NORMAL.json and OPTIMIZED.json are byte-identical to each other and to the author RESULT.json.
- Totals: 1,099 rooted inputs through order five; 5,490 literal children; 18,803 forests; 192 negative forests and 192 distinct cancellation images. Complete old/new cap comparisons agree at bounds 139, 997, and 9999, with respectively 29, 302, and 4009 attained start-two bases. Positive, absent-start-two, and resource-limit controls pass.
- Independent implementation independent.py enumerates all vertex permutations rather than using subset endpoint DP, then computes labeled-bin counts by convolution. It checks all labeled inputs through order four and their five literal children, 450 profiles total. It agrees with the candidate counts and g transfer identities, and explicitly verifies singleton Z start 0 versus formal 1 and singleton K1 start 1 versus formal 2.

No mathematical or implementation gap found. PASS is confined to the stated grammar, universal tau inequality, halving, and search cap. It supplies no staggered occupancy or eventual-spectrum theorem.
