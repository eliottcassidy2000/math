# Independent review: sharp self-converse positive-unit ternary clock

**PASS. No correction or blocker.** The frozen source packet proves exactly its stated marked self-converse transitive-cloning theorem, including eventual identities, global least positive-unit amplitude 72, explicit prime members, and the unique 7p member 217.

## Frozen source

Reviewed without alteration in `session197-orbit-wildcard/positive-clock/`:

- `PROOF.md`: b725136569ee84ddd46c5ef4f0b45b40e216eec9c314ffb604609e16a1e96f6f
- `verify.cpp`: ec9956ba32fd101105c8e27208e4f056d44da61a755d20902dc536a3f2ce7b6b
- `receipt.txt`: 229e78dbfe4523b3dbf4dc2c19ab8bb0080671296add086185749d1854a4832d

The copied verifier was compiled both with `-O3` and with `-O0 -g -D_GLIBCXX_ASSERTIONS`. Both outputs reproduce the frozen receipt byte-for-byte. Copies and outputs are in this audit directory.

## Mathematical review

The nonnegative root-run expansion makes the eventual-to-exact step valid: the largest nonzero run index supplies the leading exponential base, and distinct bases are independent. A ternary identity therefore has rank two and the exact three coefficient equations. The even-clone self-converse congruence forces A divisible by four for the +1 sign.

The root-SCC reduction is sound. A singleton root SCC would give a constant sequence. Every other component factor divides the odd gcd of 3A+1 and 9A+1, which is one. Nontrivial strong components cannot have count one, so all exterior components are singletons. The root component is invariant under the chosen anti-involution, which has no second fixed vertex, hence has odd order.

For A<72, A<=68 and H<=205. The displayed exact inequality 5^10>205^3 gives strong order at most ten, and oddness gives at most nine. The normalized universe proof uses a fixed Hamilton word of reversal composed with the specified anti-involution; odd Hamilton-path parity supplies such a word. It preserves the marked root as the middle vertex, so it does not omit a different root orbit.

The source run profile (217,216,72), matching rank two, physical self-converse clone construction, and scalar recurrence 3H-2 are correct. The proof does not confuse the +1 exponential constant with a +2 affine step. Strongness of the quotient lifts to every clone vertex, using a directed path out of and back into the root module when two vertices lie in that module.

The all-order 7p obstruction is exact. The powers of 3 modulo 7 have period six, giving m=1 mod6. Writing m=1+6k produces x=6*9^k and H_m=(x+1)(x^2-x+1). At k=0 the factors are 7 and31. At every k>=1 both exceed7, which is impossible for a number 7p with p prime. This does not require either factor to be prime, nor a numerical scan. The three asserted prime values pass independent sieve-based trial division.

## Independent implementation

`check_positive_clock.cpp` does not import the author's code. It uses a different complete source parametrization:

- label the out-neighbors of the fixed root;
- choose every tournament on that side;
- choose every symmetric binary crossing matrix;
- force the opposite side by the anti-involution.

This directly represents every marked anti-involutive source and needs no Hamilton-word normalization. Through order nine it has 66,067 physical presentations, with 4,146 rank-two inputs. The rank histograms are:

- order1: 1;
- order3: 1,1;
- order5: 2,4,10;
- order7: 8,24,168,312;
- order9: 64,256,3968,13056,48192.

For every rank-two input, a direct repeated-root word DP counts a_1,a_2,a_3. It never extracts them by clone differences and never uses the square theorem. Exactly amplitudes72 and200 satisfy the positive-unit coefficient equations, with24 physical presentations of each, all at order9. These multiplicities differ from the author's11 and19 because the parametrizations differ; the amplitude support is identical.

The checker constructs clones last rather than first and uses a reverse-endpoint Hamilton DP. It independently verifies all eight displayed first values, strongness, the physical anti-involution, and the prime cases m=5,6,8. Optimized and debug outputs are identical.

As an additional check, all24 amplitude72 inputs are rooted-isomorphic to the displayed core, with exactly one rooted isomorphism each. That core is rigid as an unmarked tournament, and its only ternary +1 amplitude72 root is vertex4. A separate literal enumeration of the384 possible anti-fixed Hamilton words finds11 valid words. These equality findings are recorded in the separate `EQUALITY.md` extension and are not needed for this PASS decision.

## Scope

This audit does not infer that arbitrary tournament spectra are covered, that the family contains infinitely many primes, or that any other self-converse mechanism has a similar 7p restriction. The all-order minimality applies to a source equipped with an anti-involution fixing its transitive-clone root. The source files remain unchanged.
