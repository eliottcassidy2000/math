# Exact target-budget numeration for physical one-apex joins

Status: PROVED arithmetic reorganization conditional only on the inherited physical join law; novelty of exact decoding is NOT claimed. Chan–Pak OPEN. No all-target coverage, asymptotic density, or new numerical cutoff is claimed.

## Inheritance and non-obvious imports

The physical input is THM-4885-one-apex-join-representation-mass-pi2-over-12.md, section A/C, and THM-4875: a rooted tournament with triple (a,b,c)=(H,E,R) joined to the one-apex tournament T(v) gives

    Phi(v)=a g(v)-b e(v)+c(v+1),
    g(0)=1, g(2w+1)=g(w),
    g(2w)=g(w)+w+(w mod 2),
    e(w)=1 if w is even, 0 otherwise.

Here a,c positive, 0<=b<=a, and c>=a. The graph is constructed by the inherited one-arc join. We use that claim conditionally, not re-audit it.

Session246's descent/PROOF.md controls D/U by binary words but loses exact integer scale. This packet instead transports the requested target at every bit. The older niche source is 05-knowledge/results/mahler_carry_interface_20261014i.md, sections 1–2: branch words encode *two* distinct source coordinates, and the source must travel with its code. The companion recursion_carry_bridges_20261014i.md explains why a lost reflection sheet must be retained. We transfer this source-retaining design, not any theorem about Collatz or Mahler to tournaments. The older Berggren comparison in 07-reflections/berggren-parabolic-golden-boolean-owner-synthesis-codex-20260812.md distinguishes positive words from group completion and retains the primitive normalizer. Here positivity is automatic in the dual coefficients and the two-state endpoint bit is the necessary normalizer.

These sources motivate a concrete map, not a numerical resemblance. Board: physical join / exact residual / two-state carry / noncommuting coefficients / adaptive seed / finite terminal certificate.

### Priority boundary

Session197 ONE_CORE_BINARY_INVERSE and Session198 DIGIT_BLOCK already give exact one-core inverse arithmetic; Session228 recursive/PROOF.md gives a complete target compiler for a different recursive physical grammar. This packet does not supersede them and does not claim that an exact membership algorithm is new. Its contribution is a source-retaining nonnegative-budget presentation, the explicit child-set gluing identity, and the elementary source bound 2R+H<=N for nonterminal representations. These are derived here, but not certified as absent from the broader corpus. No operative all-target source selector was found.

## 1. Dual digit maps

A state (a,b,c;N) means solve Phi(v)=N. Peeling the low bit gives exactly

    v=2w+1: (a,0,2c;N),
    v=2w:   (a,a,a+2c;N+b+c).

Proof: the odd case is immediate. For the even case, parity(w)=1-e(w), including w=0, and

    Phi(2w)=a g(w)-a e(w)+(a+2c)(w+1)-b-c.

The terminal v=0 is N=a-b+c. Leading zero digits may repeat this terminal description but are omitted from the canonical word.

The affine coefficient maps c->2c and c->2c+a do not commute: their two orders differ by a. Unlike inverse commutators of THM-5096, both are positive, source-valid dual transitions. They are not new graph extension letters; they describe exact inverse arithmetic for an existing physical family.

## 2. The target becomes a nonnegative budget

Define

    T=N+b-a-c,    q=c+b.

The terminal condition is exactly T=0. The same transitions become

    bit 0: (b,c,T) -> (a,2c+a,T),
    bit 1: (b,c,T) -> (0,2c,T-c-b).

Equivalently the two-state q-carry rules are

    bit 0: (b,q,T) -> (a,2(q-b+a),T),
    bit 1: (b,q,T) -> (0,2(q-b),T-q).

After the first step b is either 0 or a. Thus a source-specific exact numeration has a free zero digit which enlarges future denomination, and a one digit which spends precisely the current denomination. All arithmetic is integral and no target is rounded. The word is read least-significant bit first.

**Budget theorem.** Given positive a,c and 0<=b<=a, Phi(v)=N has a solution v>=0 iff the initial T is nonnegative and a finite legal word of these budget transitions reaches T=0. Every one step must have q<=T; zero steps are unrestricted before pruning. The word produces v by its literal binary digits and hence a physical join whenever the supplied source triple is physical.

Proof: the dual identities are equivalences for each specified parity. At T=0 the quotient v=0 solves the current equation. Reversing the finite parity sequence reconstructs v and all earlier equalities. Conversely any solution supplies its finite binary expansion, ending at v=0. For every w, g(w)>=1 and e(w)<=1, so Phi(w)>=a-b+c. Therefore every intermediate target budget is nonnegative. This proves necessity of each budget test. QED.

**Finite termination and pruning.** If T<0 reject. If T=0 accept. If T>0 and c+b>T, reject: a one digit overspends, while any string of zero digits increases c and cannot change T. At depth k, c>=2^k c_initial while T<=T_initial. Consequently every successful canonical word has at most floor(log2(T_initial/c_initial))+1 digits when T_initial>=c_initial; if T_initial<c_initial and T_initial>0 it fails immediately. This yields a fully terminating exact target solver with at most O(1+T_initial/c_initial) explored tree nodes, without enumerating graphs.

This is not a finite automaton uniform in N: the budget and coefficient magnitudes are unbounded. The endpoint state b is finite, but throwing away T or c would destroy the theorem.

## 3. Adaptive physical-source criterion

Let B(N) be any supplied, finite, target-dependent collection of literal rooted tournaments, not a fixed bank. Run the budget decoder from each of their triples. A successful terminal path returns an explicit tournament with H=N. Conversely failure of every decoder means exactly that the specified one-apex joins over B(N) miss N, not that N is absent from the full spectrum.

For each source the initial budget is N+b-a-c. Positive-budget success requires

    c+b <= N+b-a-c, hence 2c+a <= N.

The other possibility is immediate terminal N=a-b+c. Thus *every nontrivial one-apex join representation automatically has a source satisfying 2R+H<=N*. This is a useful target-dependent source pruning rule. It is arithmetic, not a bound on tournament order; transitive padding can leave counts unchanged.

An exact sufficient completion obligation is now: for every required odd N, find an actual rooted source with either N=H-E+R, or 2R+H<=N whose budget tree reaches zero. This avoids projective or CRT surjectivity claims. It does not itself produce that source. It makes adaptive physical-source selection the only missing ingredient, and gives a finite certificate for each selection.

## 4. Two-sided gluing and obstruction location

At any prefix, the remaining target belongs to the recursively defined set

    A(b,c) = {0} union A(a,2c+a) union (c+b+A(0,2c)).

The union includes the zero branch even though that branch can also lead to terminal zero. This set identity is exact, with a fixed a throughout. It is an integer-scale self-similarity, not a real limiting interval. A target budget T is represented iff T lies in A(b,c). A full integer interval can be certified by covering it with these two child sets and translating the second by c+b. Thus future interval or sieve work can operate on exact budget sets while retaining their source and endpoint state.

For a finite cutoff L, use A_L with values in [0,L]. Every positive contribution from c>L vanishes; the recursion then terminates after O(log L) levels. This gives a bottom-up exact interval certificate. No claim is made that any one such set contains a cofinal interval: the inherited finite-bank obstruction is compatible with this recursion.

The exact fresh S21 input at origin/main, THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md, supplies actual whole rays 1+c2^s on specified odd-index lines. A line seed is a valid success certificate, but its mere existence does not fill this decoder's other budgets. The two encodings should be joined by their literal source triples, not by identifying their two uses of c. Here c=R; in S21 c is an odd line index. No conflation or all-line coverage is asserted.

## 5. Useful reorganization and its boundary

The operative reformulation turns nonlinear g-evaluation into positive, noncommuting, exact target-budget dynamics with a two-state carry. It supplies a canonical finite witness, a finite obstruction certificate for a chosen family, source pruning, and exact gluing of integer sets. No projective precision or stochastic multiplicity must be upgraded to equality. It does not prove that any prescribed target has a physically available source passing the test. The next useful theorem is a physical source family whose budget sets overlap without gaps; a larger fixed-bank census is not the intended next step.

The small independent arithmetic check tests only these newly derived identities and witness equivalence, not inherited tournament classifications.
