# Contextual recoverability: hidden state, readable state, and missing winding

**Status: proved synthesis of the sibling packets.** The general quotient and information facts below are elementary and are not claimed as new information theory. The new concrete achievement is an interface between full-profile-blind tournament storage, controlled Hamilton readouts, holonomy codes, and physically measured return lengths. Each connection specifies the state, allowed experiment, and information that remains inaccessible.

## 1. The organizing theorem

Let X be a nonempty finite set of source states. Each allowed context c supplies a deterministic response R_c:X→Y_c; different contexts can have different output alphabets. A context includes its prescribed attachment or continuation, retained labels, and final observer. Define the **context response** R(x)=(R_c(x))_c and identify x∼x' exactly when all these responses agree. An arbitrary context family is allowed; only finitely many different equivalence classes can occur on X.

**Theorem (recoverability through allowed contexts).**

1. A task V:X→Z is exactly recoverable from all allowed responses if and only if V is constant on every response class. The quotient X/∼ is the coarsest state record preserving every allowed response: any record through which every R_c factors determines this quotient.
2. Adding contexts can only split response classes. An old record P is sufficient for every new context precisely when each new response is constant on every P-fiber. If this fails, processing P alone cannot reproduce the new response; the original state or a richer retained record is needed.
3. For uniform X, the remaining Shannon uncertainty after observing R is

       H(X|R)=Σ_F (|F|/|X|) log₂|F|,

   where F runs over response classes. For any full-support prior, zero conditional uncertainty is equivalent to singleton classes. For a full-support prior, H(V|R)=0 is likewise equivalent to exact recoverability of the specified task on all X; without full support it requires recovery only on the prior support.
4. Suppose X is a K-dimensional binary vector space and the responses are linear measurements collected as E. Let N=ker E and a=rank E. Then the response classes are the affine cosets of N,

       H(X|E(X))=K−a,       I(X;E(X))=a

   under the uniform law. Adding a linear bank D replaces the blind space by N∩ker D. Exactly the linear questions in im E* are answerable. For nonlinear or richer responses, the equivalence-class statement remains valid but rank need not describe it.
5. If every query has at most b outputs, deterministic exact identification of all |X| states requires worst-case query depth at least ceil(log_b|X|), for b≥2. More specifically, fewer than K adaptive binary linear queries cannot certify the zero vector in the worst case: their all-zero branch has a nonzero common kernel. Neither bound treats an arbitrary integer-valued count as one bit.

**Proof.** A function of the responses must agree wherever the responses agree. Conversely, a fiber-constant task defines a well-defined function on response classes. Any record supporting all responses refines those classes, proving (1). Appending response coordinates intersects equality relations, proving (2). A uniform source is uniform conditionally on each response class, giving (3); with full support a class has zero conditional uncertainty exactly when it contains one possible state, and the same argument applies to V. Rank–nullity and affine fibers prove (4); a linear functional factors through E exactly when it vanishes on ker E, equivalently when it lies in im E*. A b-ary depth-q decision tree has at most b^q leaves. On an all-zero binary-linear branch, fewer than K homogeneous equations leave a nonzero solution, proving (5). ∎

This theorem does not claim that forgotten data spontaneously reappears. A new physical experiment may access structure still present in the original object but absent from an old statistic. If only the old statistic was retained and the source destroyed, the missing response cannot be manufactured.

For a linear dynamical system with state space V, observation C, and allowed continuation operators A_i, the sibling framework strengthens this into a compositional theorem: the blind subspace is

    N=intersection_w ker(C A_w),

not merely ker C. It is the largest invariant subspace inside ker C. In dimension D the closure stabilizes by length D−rank C. These statements require all finite words in the declared operators, including the empty word. A finite-state restriction must first be incorporated into a finite-dimensional controller-augmented system, with its resulting dimension used in the bound; arbitrary language restrictions are not covered. The finite-source theorem above does not assume a tournament invariant is itself linear or that all contexts preserve X.

## 2. The exact tournament implementation

The [storage and read-head packet](../information/NOTE.md) constructs, from the S170 holonomy code

    C_r=RM(d−r,d),       M=2^d,       K=Σ_(j=0)^(d−r) binom(d,j),

with d≥2 and 2≤r≤d, actual tournaments T_h on 5M vertices. Each is an ordinal sum of M five-vertex strong seeds. The seeds are nonisomorphic and differ by one specified arc, but have identical labelled-slot polynomial

    F(u)=u⁵+6u³+2u.

The resulting 2^K tournaments are pairwise nonisomorphic, have common complete profile F(u)^M and common H=9^M, and satisfy the exact labelled arc-distance identity

    distance(T_h,T_g)=weight(h+g).

Thus the bulk profile observer has one response class and loses all K uniform bits. This is not because the graphs are identical: their ordered strong-component types recover h.

An allowed context adds one vertex in addressed block y, with fixed local attachment: it beats only local vertex 0, is beaten by earlier blocks, and beats later blocks. The same rule is used without first reading the seed bit. Its Hamilton response is

    H(Q_y(T_h))=(17+2h(y))9^(M−1).

Consequently the addressed context is exactly a binary evaluation h↦h(y), after known rescaling. The theorem now gives the entire recovery contract: a bank E leaves precisely ker E hidden; the K addresses of weight at most d−r recover all h; fewer binary linear reads cannot universally certify identity.

For d=4,r=2, the source has 80 vertices, 2,048 nonisomorphic states, and 11 hidden bits. Reading the five addresses of weight at most one leaves six bits, namely the unconstrained quadratic coefficients. Reading the six additional weight-two addresses removes that blind space. This is a concrete recovery of information absent even from the *complete* original path-cover profile.

The [Fisher packet](../fisher/NOTE.md) makes the task sensitivity exact. Under the uniform code law, a nonzero character χ_ℓ(h)=(-1)^ℓ(h) survives the bank precisely when ℓ∈im E*. Otherwise it is independent of the transcript and its exponential-tilt Fisher information is completely lost. For exponential tilts about the uniform law, at parameter θ=0, general scores lose exactly their squared Fourier mass outside that answerable character space. Statistical evidence for identity can increase before a deterministic identity certificate becomes possible.

### Rich outputs change the accounting

Adding one read vertex in each of q distinct selected blocks yields

    H(Q_S(T_h))=9^(M−q)17^(q−w)19^w,
    w=Σ_(y∈S)h(y).

For known q, this response recovers the ordinary subset weight w. It has at most q+1 outputs and carries at most log₂(q+1) bits per such query, not necessarily one bit and not necessarily the full upper bound. It generally merges different words of the same subset weight. The theorem applies through its actual fibers, but the binary rank argument does not apply to that response as if it were a single linear bit.

## 3. Winding memory has the same quotient logic, but different state

The [winding-code packet](../winding_codes/NOTE.md) studies S170's different six-vertex seed torus. Every N-block state has H=29^N and the same complete profile. Here a logical loop returns to a pinned-sheet error e∈F₂^N. A syndrome A e hides exactly a chosen linear code C=ker A. The first nonzero hidden return takes exactly

    3 d(C)

literal arc flips, where d(C) is minimum nonzero codeword weight; for C={0}, set d(C)=∞ and there is no hidden return. This is a physically attained threshold, not just an entropy analogy. Syndrome rank prices retained sheet information; code distance prices the shortest return mistake it can conceal. The states here are pinned return sheets, not the nonisomorphic storage tournaments of the previous section.

The [knot packet](../knots/NOTE.md) proves a different and stronger failure for finite winding observations. On a marked standard torus in S³, for any fixed integer m≥1, the primitive slopes (1+2mk,2), k≥1, all have the same mod-m deck signature and isotropic-cover lift count, but their torus-knot genera are mk. Thus a finite deck response leaves unbounded marked knot complexity hidden. No amount of repeated querying of that same residue response separates this family. One needs additional slope information or richer specified experiments, and still needs the ambient marking to interpret a slope as a knot.

These knots, the tournament flip torus, and the knot complement are distinct spaces. In particular, the knot packet proves that every representation of a torus-knot group into a finite nilpotent group has cyclic image. An isolated signed-unit defect may occur as cyclic holonomy; the whole noncommuting signed-unit transport system cannot. A common word such as “holonomy” does not make their allowed continuation algebras interchangeable.

## 4. New concepts and their concrete payoff

The useful new synthesis is not a new claim of general observability theory. It is a small collection of experimentally specified interfaces:

- **Contextual blind space:** the defects invisible after every permitted continuation, rather than only now. It says exactly when a compressed record remains valid under composition.
- **Hamilton read head:** a fixed one-vertex attachment with the exact 17/19 response. It turns a previously inaccessible code coordinate into an actual tournament count without assuming the source bit is already known.
- **Two different code resources:** a code family stores messages with arc-distance protection; a kernel code hides deck returns with a sharp 3d(C) failure threshold. Message space and hidden-error space must not be conflated.
- **Answerable character space:** im E* identifies which binary questions survive, not merely how many bits remain.
- **Marked winding blindness:** a fixed finite transport response can have unbounded knot genus in one fiber. It prevents finite-sheet results from being promoted into knot classification claims.

The [response framework](../framework/NOTE.md) also supplies an exact square-zero boundary response for Hamilton paths, with order-independent elimination. This is a constructive interface retaining visited-vertex data, not merely a scalar count. Its state size may be exponential; no polynomial compression follows from calling it a boundary response.

## 5. Actionable targets for returning to the Hamilton-count spectrum

These are research targets, not consequences already proved.

1. **Measure the missing context space.** Start with the verified equal-profile five-vertex seeds. Enumerate fixed one-vertex attachments and their H-response vectors, then extend only where a rank or fiber test shows genuinely new separation. Compare those observations with the framework's boundary response. The deliverable is a minimal declared context bank, not an unspecified “complete invariant.”
2. **Find a genuinely coupled strong construction.** The present storage family uses ordered strong components and its local read formulas multiply by 9^(M−1). Search for a uniformly strongly connected family with an authenticated bit readout and audited profile behavior. A fixed strong wrapper, new feedback arcs, or altered local seed must have its own exact count formula; ordinal multiplication alone does not provide this coupling.
3. **Separate arithmetic realization from return certification.** The winding carrier amplifies an existing count q only to q·29^N; it does not solve missing prime lanes. A proposed spectrum advance must specify which new target values are obtained, preserve a realizable construction language, and show that its context or composition law escapes the already known multiplicative family. Improved source-return verification is useful infrastructure, not a new spectrum value by itself.
4. **Attach a verified return certificate to a useful arithmetic construction.** Identify a natural operation family in which phase or sheet errors can invalidate a claimed endpoint count. Compute its actual response quotient and test whether the code-distance threshold applies to that family, instead of selecting an arbitrary code and calling it intrinsic. A positive result would protect a real spectrum-building process rather than merely exhibit hidden memory.

The session's strongest established conclusion is therefore operational: a scalar or even complete bulk profile may discard all stored state, while a small authenticated family of allowed contexts recovers it exactly. Whether that capability creates new Hamilton counts remains a separate constructive problem.

## Source and audit boundary

This synthesis uses the exact sibling files `framework/NOTE.md`, `information/NOTE.md`, `fisher/NOTE.md`, `knots/NOTE.md`, and `winding_codes/NOTE.md` read in this session. Its new theorem is proved above; numerical examples are the proved and checked instances in those packets. It does not import an unexamined historical theorem or claim an independent all-repository audit.
