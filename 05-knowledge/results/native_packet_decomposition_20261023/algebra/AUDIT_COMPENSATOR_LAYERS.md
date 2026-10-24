# Independent audit: one-active-helper layers

**PASS.** Independently reviewed `../interface/COMPENSATOR_LAYERS.md` after its
physical-row ledger included the diagonal target-degree contribution.

The state bijection is exact under the stated independent binary supply and
authenticated complete source-current criterion. With no helper at V, the
old V-set has cardinality k and charge c. With helper h at V, its cardinality
is k-1 and its charge is c-lambda_h. These are separate constraints: equality
of source charge alone does not preserve the target margin.

The changed-slot metric is the symmetric-difference size of old V-sets plus
the Hamming distance of the helper occupancy vectors. Those vectors are zero
or standard unit vectors, so their distances are exactly0,1,2. A swap between
active helpers is an atomic two-helper change, not double occupancy at arrival.
Actual path-clock distances may be used as nonnegative per-slot weights.

For a retained character X with equal target values t, an old V-set change
from I0 to I changes the U-owned current by

    2[(sum_I X - sum_I0 X) - t(|I|-|I0|)].

This follows directly from moving paired conductance2 away from U:
removing the edge to source s contributes +2(X_s-t) to the U row. The helper
U-to-V move contributes +2(X_h-t). The total vanishes when the translated
charge equation and the k-1 cardinality both hold. V has the opposite debt,
and each source-row debt is zero because X_U=X_V. The formula is valid over
integers and in the corresponding finite modulus, with the paired factor
retained. Omitting the t terms from a PARTIAL row ledger would be incorrect.

The theorem only represents the one-active-helper restriction. Independence
does not create new source paths or capacities, and the supplied unused
helpers remain inventory records throughout. The conclusion that the actual
two named helper layers connect an eighteen-source fibre is delegated to the
native finite certificate; it is not inferred from the abstract layer theorem.
