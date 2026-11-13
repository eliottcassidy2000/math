# Incoming mathematics and the changed construction obligation

The session started at06847e4c6e6709585bf1b4cfa6f84713081034ea and read
incoming work through ebd71e68d8. The claims below retain the scope of
their source proofs. We did not rerun their broad computations; their
role here is to identify which physical interfaces new work must supply.

## 1. Bounded signed descent is available; its graph execution is missing

[Session305 CAP4](../session305_collatz_core_changes_20261113/CAP4.md)
supplies a mixed-sign odd descent outside{1,5,7,21}: at most26 signed
steps, with exact halving exponents1 through4. Its finite dyadic
certificate has1068 disjoint cylinders, precision at most43, an affine
threshold at65, and separate finite boundary controls. This is a theorem
about a chosen-sign arithmetic system, not ordinary positive Collatz.

The associated reverse ladders avoid7 and21. This removes one arithmetic
obstruction to a tournament compiler, but a ladder's intermediate count
does not automatically come with the required root or module state.
Induction on its smaller endpoint alone does not construct those states.
Our finite target selectors return actual modules and a smaller source
where they succeed; they do not yet execute every ladder in that cover.

**Concrete new question:** can the ordered determinant condition
e=j(6p+h)−6hq=±1 be refreshed while a core-change macro follows one of
these bounded ladders? Our fixed-partner and padding bounds rule out a
finite stock of nontrivial pairs with transitive padding as the whole
answer. They leave changing strong modules and new observation interfaces
available.

## 2. There is another physical minus-Collatz edge, with a different cost

[Session306 MATCHED_COLLATZ](../session306_relative_isolators_20261113/MATCHED_COLLATZ.md)
has a crossing matched core X(a,b,c), with

    H=1+2^(a+c+1)3^b+2(2^a+2^c)2^b.

For u,b≥1, a single adjacent vertex slide / one arc reversal gives

    X(u+1,b,u+1) -> X(u,b+1,u+1), H -> (3H−1)/2.

The resulting outer widths are unequal, so it is not immediately another
entry of the same rule. The paper's full fixed-outer-width classification
also excludes a fixed signed macro of length at least2 in that chart
family. These are chart results, not exclusions of changing cores.

Our Π chart instead spends one dyadic tail vertex per transfer and can
take several consecutive steps. Its exact zero-site deficit18H(M), and
the native refill, identify the physical cost at its own boundary. Both
examples point to a measurable resource, but their resources differ:
outer-width equality versus a nonempty module site. A general compiler
must explicitly transport or restore that resource.

## 3. Higher visits explain the new response coordinate

[Session305 APEX_FUSION](../session305_collatz_core_changes_20261113/APEX_FUSION.md)
uses equal normalized two-cover responses to move a real peelable factor3
at fixed total and order. This transfers an existing resource; it does
not create it. Its third-visit observation separates some modules that
the two-visit observation identifies.

[Session306 RELATIVE](../session306_relative_isolators_20261113/RELATIVE.md)
shows that a uniform module context with m exterior vertices observes
the first floor(m/2)+1 normalized ordered-cover coordinates. Matching
only the first two is therefore not universal interchangeability.
Our explicit handoff root has

    F=186P_K(3)+48P_K(2)−53H(K),

an actual continuation statistic that needs the third coordinate. This
is a concrete instance of that incoming context principle, and mirrors
the older square-sum requirement to retain named boundary pairings.

## 4. Why this core change is outside the fixed-frame obstruction

[Prime ternary CLONE_PARITY](../prime_ternary_lift_20261113/CLONE_PARITY.md)
shows that a fixed multitransitive-block frame's count modulo4 depends
only on which blocks are singletons. Repeated pure-growth+4 transitions
consume singleton data. Its one-module finite-frame results also exclude
a nonnegative-length+4 cycle in that language.

The new handoff deletes two fixed vertices, inserts three, and changes
the core and its module incidences. It deliberately changes the root
star. The clock refill creates a fresh site at an explicitly changed
count. These moves are outside the pure-growth hypothesis; they do not
refute the obstruction. The unresolved task is controlling their cost
for a prescribed target and proving entry into the next required chart.

## 5. Finite arithmetic registers are still available

[Session306 COLLATZ_LANGUAGE](../session306_relative_isolators_20261113/COLLATZ_LANGUAGE.md)
studies the exact prefix-slope survivor language for ordinary shortcut
Collatz. Its nonregularity and finite-horizon state lower bounds concern
a finite-state language model. They do not exclude finitely many
unbounded integer registers with guarded arithmetic updates, nor a
tournament construction whose changing module is such a register.

The present debt theorem is a sufficient positive-word descent criterion,
not a recognition theorem for all surviving ordinary Collatz words.
Our reusable clock circuit changes the arithmetic map at its boundary.
Keeping these quantifiers prevents either result from becoming an
unjustified convergence or finite-automaton claim.
