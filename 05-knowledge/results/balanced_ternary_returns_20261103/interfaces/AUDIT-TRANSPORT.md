# Independent audit of the named 59 to 181 source exchange

**VERIFIED exact finite certificate and scoped minimum; no global shortest
native-route claim.** The audited object is the eight-name source
`202032`, its H-preserving endpoint exchange, and the subsequent actual
two-birth reset to a balanced H=181 child.

`audit_transport.py` imports no producer or native-compiler code. It
replays each of the six native moves by complementing exactly the edges
outside the supplied Hamilton path. It verifies those paths, temporary
birth restrictions, deletion of only the temporary name, each intended
old-edge reversal, and all fifteen saved stage counts by independent
subset endpoint DP. The producer instead enumerates simple paths.

The audit verifies the physical map from the preceding session's C_2:

    (c_2,x,y,c_1,birth_a,birth_b,reset_root,reset_child)
       -> (0,1,4,5,6,7,2,3).

The named endpoint profiles are `(59,17,17,59)` before the exchange,
`(59,5,5,61)` after it at physical root 2, and `(181,59,59,181)` at the
new reset root 8. All 68 saved general two-port birth queries are rebuilt
as literal graphs and checked independently. Their operation is not the
two-apex word suffix operation.

The physical names persist, and the edge 2->3 is restored at the exchange
endpoint. Some intermediate native moves reverse that edge. The certificate
explicitly checks this hostile, so the result cannot be read as preserving
the port edge or H at every microstep. The exchanged endpoint and final
balanced child are weak: physical name 6 remains a global sink. Strong
intermediates are not asserted.

## The exact scope of the three-edge minimum

The first template audit found one missing admissible image: the protected
child could be mapped to the padded sink as well as to the four ready-core
out-neighbors. The producer now includes that case. Independent exhaustion
of all 3600 maps gives minimum distances 3,4,6,8,9 for those five images,
so the minimum remains three.

A stronger endpoint minimum follows from a separate exhaustive check.
The audit recomputes all 28 one-edge and 378 two-edge reversals from the
actual source. None of the one-edge changes retains H=59. Exactly five
two-edge changes do. The only weak one reverses the protected edge 2->3
and is excluded by the endpoint guard. Of the four strong guarded ones,
two are explicitly isomorphic to the original source and one to its
converse. Those cannot have an induced H=181 extension by the separately
audited original-source exclusion. The fourth has its complete 131072
two-birth cube checked independently by the direct Hamilton-path deletion
polynomial, with no H=181 result.

To pass from the finite birth bound to arbitrary extension order, use the
packet's inherited strong-component argument: a strong H=59 source in an
H=181 extension lies inside its H=181 component, since 181 is prime.
The cited strong-order minimum at order 11 is 225, so this component has
order at most ten, allowing at most two added vertices. This same argument
applies to each strong guarded mutant. It is the additional mathematical
cutoff; the audit does not infer all-order exclusion from a bounded search.

Consequently three is the minimum old-edge distance to an **eight-name
H=59 endpoint retaining 2->3 and admitting an induced H=181 extension**.
This minimum permits a weak endpoint but fixes its order and H. It is not
a minimum number of native moves, a minimum over all temporary sizes, or
a statement about endpoints of different orders or counts.

Run `python3 -B audit_transport.py` and `python3 -B -O audit_transport.py`.
The two receipts must agree byte for byte. Their input hashes pin the
specific producer and interior receipts audited.
