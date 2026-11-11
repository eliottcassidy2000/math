# Fifth wave: cycle resonance, controlled chart changes, and exact return debts

**PROVED structural constraints and a controlled constructor; OPEN complete binary controller.** This wave uses the sharp spectral normal form to classify recurrent finite-control possibilities, rather than enlarging the word census.

## 1. Historical draw and mechanism

Seed261110 draws10 files from the original11,638-file pool. The recorded sources read are:

* `05-knowledge/results/session254_four_frontiers_20261111/jacobian/PROOF.md`: a scalar repair can pay one response while leaving a precise quotient of unpaid memory. Transfer: separate a scalar digit identity from all recurrent spectral debts.
* `05-knowledge/results/session74_eleven_mechanisms_20261027/session74-inheritance/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`: periodic word action becomes multiplication in a retained orbit carrier. Transfer: follow an actual controller cycle before comparing its exponential scale.
* `05-knowledge/hypotheses/HYP-9276-k6-minus-p6-terminal-is-the-adjoint-of-a-weight-3-level-32-newform.md`: tensor-product traces multiply eigenvalues. Its modular identity remains an inherited OPEN claim and is not used as a theorem. Transfer: two-ended spine contexts carry products of left and right eigenvalues.

The new operation is a controlled converse rebase; the new exact condition is a resonance test plus a finite return-debt kernel. Both retain literal physical presentations.

## 2. Append-only finite controllers cannot recur through longer digit blocks

Consider any finite controller for a fixed k-apex spine grammar. Each binary digit appends a specified nonempty word, possibly depending on the controller state, and changes that state. Suppose all digit sequences from a physical seed with H>=3 obey H->2H+1 or H->2H-1 exactly.

Take any reachable directed controller cycle with d>=1 digits. Let its concatenated physical word have length L. Each digit inserts at least one vertex, so L>=d. The composed scalar operation is

    h -> 2^d h + c,       |c|<=2^d-1.

Repeating the cycle m times gives

    h_m = 2^(dm)*(h+c/(2^d-1)) - c/(2^d-1).

The coefficient of2^(dm) is strictly positive because h>=3. By the fixed-k periodic spectral theorem, its visible exponential bases must be j^L, j in{1,...,k+1}. Thus

    2^d=j^L.

Unique factorization gives j=2^a and d=aL for a>=1. Combined with L>=d this forces a=1 and L=d.

**Theorem.** Every edge on a reachable recurrent controller cycle appends exactly one vertex per digit. Longer append-only blocks can occur only on transient edges, never as the recurring mechanism of an exact binary constructor.

This holds for every fixed k, not merely k=2, and permits arbitrary state-dependent choices. It does not assume stationary digits or a common affine line. An operation that physically changes the core, reselects a deletion set, or changes a marked interface is outside this theorem unless it is itself just a fixed appended word in the same grammar.

This explains why increasing append-block length in a stationary scan is not a promising recurrent escape: the spectral obstruction rules it out in advance, even if a few first steps work.

## 3. Exact cycle-debt filters at dimension9 for two apices

For an append-only control cycle, let M be its word matrix, alpha=2^d, and c its signed binary carry. The desired cycle identity is

    L M v = alpha L v + c.

Augment by a constant coordinate and put B=diag(M,1), r=(L M-alpha L,-c). All repetitions are correct iff

    r B^t (v,1)=0       for every t>=0.

For k=2, B is annihilated by the same degree9 polynomial

    t^3(t-1)^3(t-2^L)^2(t-3^L)

as M. Therefore checking t=0,...,8 suffices exactly. This is a compact all-size cycle filter, not a nine-step empirical heuristic.

When resonance gives L=d, the generic periodic shape is

    C*3^(Lm) + (A*m+B)*2^(Lm) + D*m^2+E*m+F.

An exact binary return requires C=A=D=E=0, together with the prescribed constant/leading relation. These are the unpaid return debts hidden by checking only the first scalar digit. They can be imposed as linear constraints on the physical state and intersected with the integral reachable-lattice gate from `../core/INTEGRAL_RETURN_GATES.md`.

Checking a selected finite collection of cycles is only a necessary pruning test for a proposed controller. The full simultaneous row closure in `FOURTH_WAVE.md` remains the sufficient-and-necessary all-word certificate.

## 4. Two-ended cycles have a different resonance law

A composition of append, prepend, and converse operations always has one of the forms

    w -> u w v,        or        w -> u w* v,

for fixed words u,v. Let r=|u|, s=|v|, L=r+s.

For an orientation-preserving cycle w->u w v, its context operator is

    X -> M_v X M_u.

Its eigenvalues are products j^s l^r, with j,l in{0,1,2,3} for k=2, interpreting a zero-length factor as the identity. If the cycle performs d binary digits on an orbit of count>=3, the positive mode2^d must be present. Consequently

    d is one of r, s, r+s,

omitting zero possibilities when d>0. To see this, no positive factor involving3 can multiply to a power of2, so each nonempty-side factor is1 or2.

Thus a longer two-ended cycle can evade append-only resonance only if one side is spectrally neutral (base1) on its growing count mode. When both sides contribute their base2 mode, every inserted vertex must again pay for one digit.

For an orientation-reversing cycle w->u w* v, squaring gives

    w -> u v* w u* v.

Both new side lengths are L. Applying the preserving-cycle rule to its2d digits gives

    d=L/2 or d=L.

This is a precise necessary condition for a finite controller using converse rebases. It is not sufficient: eigenvalue resonance does not eliminate polynomial pollution or prove a physical seed exists.

The marked reselection reported by the constructors lane can escape these fixed-word restrictions more fundamentally: if it changes the sufficient interface by a map not preserving the common spine flag, its cycle operator is not a product of the old letter matrices. This is exactly where a genuinely different physical return should enter, rather than as an unlabelled equality of scalar q.

## 5. A concrete controlled chart-return constructor

Let w_m=(1,2,0,2^m), m>=1. The earlier certified ray gives

    H(w_m)=15*2^m-1.

Use two presentation states:

    state0: literal word w_m;
    state1: literal word w_m*.

In state0, append2 and then apply converse, entering state1 with word w_(m+1)*. In state1, apply converse and then append2, entering state0 with word w_(m+1). Both are specified physical operations and both implement

    H -> 2H+1.

This is an all-size finite-state return with alternating chart choice, not a formal scalar recurrence. Its two-digit cycle inserts exactly two vertices and returns to append22, in agreement with the resonance theorem. It proves the controlled-chart machinery can carry a real constructor, while making its limitation explicit: it supplies only the plus digit. There is no new count coverage beyond the existing ray and no minus digit hidden in the converse.

An inert padding example further shows why resonance alone is weak: prepending0 is a universal-source insertion, so w->0w2 on this ray still gives2H+1 while adding two vertices. The extra left side is base1, exactly the exceptional two-ended resonance case. It changes a presentation but supplies no missing count operation.

## 6. What the next successful certificate must contain

A viable recurrent binary controller must do at least one of the following:

* use a state-dependent one-vertex append menu and satisfy every cycle's mode-debt kernel;
* use two-sided/converse cycles with the stated neutral-channel resonance, while paying all their nonconstant spectral debts;
* use a genuine core/mark/interface rebase outside the common fixed-spine flag, with its full sufficient tensor retained.

For each candidate cycle, test resonance before any expensive seed search, then the nine-row mode kernel where applicable, then the integral orbit-lattice condition, and finally literal physical reachability. These are logically separate gates. The parent/core lane's half-integral obstruction shows why rational kernel feasibility is insufficient.

## 7. Transfer ledger

Source: periodic equivariant orbit carriers. Target: a recurrent binary tournament controller. Map: follow one actual control cycle, compare its exact scalar multiplier with the block spectrum. Preserved predicate: a nonzero growing2^d mode. Destroyed information: transient states and nonvisible modes. Sidecar: physical word length and the positive leading coefficient h+c/(2^d-1). Cheapest test: a purported two-letter append digit repeated on a cycle must fail j^2=2. Outcome: the all-k recurrent length theorem.

Source: a complete response quotient with unpaid memory. Target: the cycle's exact affine return. Map: retain all residual rows rB^t, rather than the first output only. Preserved predicate: the cycle is correct at every repeat. Destroyed information: vector coordinates outside the observable cycle Krylov space. Sidecar: exact annihilator and integral/physical reachability. Cheapest test: nine exact residuals for k=2. Outcome: a finite all-size return-debt gate.

Source: tensor-product spectral multiplication. Target: prepend/append/converse cycles. Map: X->M_v X M_u and the corresponding Kronecker product. Preserved predicate: the allowed exact exponential modes. Destroyed information: coefficients and Jordan debts. Sidecar: both boundary matrices, controller orientation, and literal seed. Cheapest test: the controlled plus-ray and the universal-source padding example. Outcome: a constructive chart return and precise two-sided resonance condition, not complete binary coverage.
