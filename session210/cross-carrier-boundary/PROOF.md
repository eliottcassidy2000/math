# The cover polynomial does not determine the projective carrier

Status: EXACT EXPLICIT COUNTEREXAMPLE, with a direct structural proof and independent review. No global-spectrum or literature-priority claim. This is a boundary between the two audited Session210 correspondences.

Let A and B be the order-six tournaments with outgoing-neighbor bitmask rows

    A=(58,4,1,6,14,30),
    B=(38,8,2,5,15,30).

Bit j of row i means i→j. Their forward-adjacency histograms, indexed by f=0,...,5, are both

    (17,105,238,238,105,17).

By the exact binomial transform in [the gate packet, Section6.1](../graded-spectrum-gates/PROOF.md), both therefore have

    P(q)=sum_f m_f binom(q+f,6)=q^6+8q^4+8q².

Consequently H=17, all nested-gate Hamilton counts agree, and their permutohedron face–sink incidence polynomials agree. This is equality of the scalar geometric incidence polynomials, not an isomorphism of the complete oriented permutohedra.

Their directed odd-cycle packing polynomials differ:

    I_A(z)=1+8z,
    I_B(z)=1+6z+z².

Literal cycle enumeration gives A four directed triangles and four directed five-cycles, with no disjoint pair. B has four triangles and two five-cycles, with precisely one disjoint pair. At order six no three nontrivial odd cycles can be disjoint. Thus these enumerations give the complete displayed polynomials, not truncated approximations.

The [projective-cell correspondence](../projective-cycle-carrier/PROOF.md) makes Q_A eight isolated points. It makes Q_B a circle, represented by two vertices joined by two parallel edges, together with four isolated points. Their cell vectors are (8) and (6,2); their rational and mod-two Betti vectors are (8,0) and (5,1). Both have eight cells, as required by (H−1)/2=8. Their maximum packing sizes are one and two, respectively.

Therefore even the full cover polynomial, equivalently every nested-gate response and the full face–sink incidence polynomial, does not determine the projective carrier's topology or packing number. The failure first becomes possible at order six: below six there are no disjoint nontrivial odd cycles, so Q_T is a discrete set of (H(T)−1)/2 points, already determined by P_T(1).

The standalone verifier enumerates all 720 permutation words of each explicit tournament; expands the binomial transform exactly over rational coefficients; independently counts ordered-bin covers for q=0,...,4; and enumerates every odd directed cycle and compatible pair. It reports the cycles themselves for inspection. The checks use explicit exceptions and remain active under Python optimization. Run:

    python3 verify.py > replay.json
    cmp replay.json RESULT.json
    python3 -O verify.py > replay.optimized.json
    cmp replay.optimized.json RESULT.json

The infinite equality of nested-gate responses follows from the proved cover-operator formula, not from testing finitely many closure depths. The minimal-order statement follows from the six-vertex packing bound, not a claim to have censused all higher-order structures.
