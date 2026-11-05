# Independent review: complete uniquely Hamiltonian tournament class

**Result: PASS.** No correction was required to the frozen proof or verifier. The elementary all-order argument establishes that every strong tournament with exactly one directed Hamiltonian cycle has a strong bridge. It therefore removes the structural qualification in the earlier rooted/Fibonacci packets. This review makes no literature-priority claim.

## Bound source snapshot

The reviewed files are in `session197-lowcycle/`:

- `C1_COMPLETE.md`: `d38a53b1c475df3578d58e83d9e9f97785d9455084c0a094aaed752577d8b289`
- `verify_c1_complete.py`: `645ddfc551fd5ce2760f2d1f71318c08f57e1115e54edf9d28d7db3e6ad91d3b`
- `c1_complete_certificate.json`: `6b039e86e066ecdfbdacba4ab92100e991adfbaa5d006c95fd711c079ed1cd30`

Its consequences use the separately frozen `ROOTED.md` (`ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e`) and `FIBONACCI_CLASS.md` (`d95a133e360b96589009ebd8a1c8449525c4a2847a88b33a4b3f2b839739e63e`). Those files were not changed by this completion.

## Symbolic proof checks

1. **Deletion/insertion.** If T-v is strong, Camion's standard theorem gives a spanning directed cycle D there. Strongness of T supplies both neighbor types of v. Following D from an in-neighbor to the first out-neighbor gives a genuine adjacent insertion gap. Inserting v gives the unique cycle of T, and deleting v recovers D. This proves both that C-v is a cycle and that it is unique. There is no assumption that every deletion is strong.

2. **2-strongness forces second-step arcs under uniqueness.** Applying that lemma separately at each vertex forces i to point to i+2 in the unique cycle's indexing. A 2-strong tournament has both minimum semidegrees at least two, so its order is at least five. This handles orders three and four without trying to use inconsistent second-step orientations on those small cycles.

3. **Odd parity.** For odd order, step two is relatively prime to the order, hence the prescribed arcs form a spanning directed cycle different from the original one.

4. **Even parity, chord 4 to 1.** The cycle 0,2,4,1,3,5,6,...,n-1,0 uses each vertex once and only existing step-two arcs, the tested chord, and original cycle arcs. In particular at n=6 its last edge is the original 5 to 0, so the boundary case works.

5. **Even parity, chord 1 to 4.** Cutting 2 to 4 from the even step-two cycle and 1 to 3 from the odd one, then adding 2 to 3 and 1 to 4, reconnects the two directed paths into one cycle rather than two. The explicit order is 4,6,...,n-2,0,2,3,5,...,n-1,1,4. At n=6 this becomes 4,0,2,3,5,1,4. It is distinct from the original in every allowed even order.

6. **Separating-vertex product.** Removing a vertex from a directed Hamiltonian cycle leaves one Hamiltonian path. In a tournament that path must run through the SCCs of T-r consecutively in their intrinsic order. Fixing r as the cyclic cut removes rotational multiplicity. The displayed product of initial rooted-start choices, middle path counts, and terminal rooted-end choices is therefore a genuine bijection, including singleton components and k=2.

7. **The actual bridge.** Strongness guarantees nonempty initial/terminal neighbor sets. Each relevant rooted path count is positive by Camion, or by the singleton convention. A product of positive integers equal to one forces exactly one initial out-neighbor a of r. Across the cut D1 versus the remaining components together with r, every arc points out of D1 except r to a. Removing that arc leaves no path back into D1. This proves strong-bridge status directly, not merely that the arc lies on the unique Hamiltonian cycle.

8. **Consequences.** The established rooted equivalences now give an outdegree-one vertex, a root-start-one representation, and the complete arbitrary-gadget grammar for every strong C=1 tournament. The previously proved unlabelled enumeration therefore becomes a_3=1 and a_n=F_(2n-6) for all n>=4. The rigidity statement gives n! F_(2n-6) labelled graphs at n>=4, while the triangle contributes only two labelled orientations. The prior closed bounded-image calculation can consequently be interpreted as the complete C=1 spectrum through 5,001, without introducing an order bound or claiming a conductor. No conclusion about the full C=2 class or the unrestricted {7,21} spectrum follows.

The only imported existence theorem in the new structural argument is the explicitly identified standard Hamiltonian-cycle theorem for nontrivial strong tournaments.

## Reproduction and independent counterexample search

The source verifier was copied before execution because it writes its certificate beside itself. Normal and optimized runs both match the frozen certificate byte-for-byte, and no source file was modified.

The independent checker `independent.py` does not import the author's verifier. It uses bitset reachability and a capped directed-cycle DFS, fixing the least vertex as the cycle start. It exhausts all **33,867 labelled tournaments through order six**.

- It finds 6,146 uniquely Hamiltonian tournaments: 2 at order three, 24 at order four, 360 at order five, and 5,760 at order six.
- All 2,664 2-strong tournaments found have at least two cycles.
- All 22,728 strong-deletion inheritance checks pass, including explicit insertion gaps and uniqueness of the smaller cycle.
- All 13,734 separating-vertex checks expose exactly the claimed initial reverse arc; deleting it destroys strong connectivity. The analogous terminal arc is also independently checked as a bridge.
- Every uniquely Hamiltonian case also has an outdegree-one vertex, as the now-complete rooted classification predicts.
- 1,920 additional cycle-surgery tests through order 100 randomize every unused arc, checking both orientations of the crucial chord at even orders. Both even cases receive 480 tests; odd step-two receives 960.

The unique-cycle/bridge-data digest is `bf9db6824205bd21e08770b8b7cec13b457e3c2d4ded4fb17ebe5728c9515807`.

The finite searches are controls and counterexample checks. The symbolic argument above, not the census, establishes the all-order theorem. There are no unresolved scope or implementation concerns in this snapshot.
