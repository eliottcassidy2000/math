# Proof and implementation audit

Status: **same-agent review with independent computational paths**. No
independent-agent review is claimed. The audit programs import no producer.

## Mathematical gates

1. **P1 is a synthesis, not rediscovery.** The 2026-10-19 simultaneous-pair
   theorem already proves exact nonedges and collisions. The new use adds
   one common fictitious anchor for all singleton queries, an empty-set
   candidate, and the retained actual tournament. The empty label is
   LH+2, not the old construction's label2: placing it after every source
   label makes the increasing-order tournament convention faithful.
2. **Source types.** Arbitrary old T is a sidecar unless its supplied labels
   already obey the square orientation rule. A birth changes only new
   incidences. No information theorem replaces old directed arcs by an
   undirected square graph. In the literal H31 witness the entire old and
   new graphs do obey the square rule in their retained name order.
3. **Exact response includes two-component covers.** Deleting the newborn
   need not leave an old Hamilton path. The exposed-step Q matrix accounts
   for these missing contributions. The quadratic identity, its cut-field
   decomposition and its Walsh interpretation are inherited from
   THM-4104 and THM-4115 with the opposite cut convention stated explicitly.
4. **The free probe uses a real constraint.** R(V)=R(empty)=H gives the one
   coefficient relation. It saves one pair probe, not an arbitrary source
   bit or a port. The linear-dimension result is for the unrestricted cut
   field; actual tournament full rank is certified only at n=4,5 here.
5. **Common source and ports.** All three known branch endpoints are
   reconstructed from the C3 command program, recounted, checked strong,
   and checked against every pair of integer labels. Their seven original
   numeric labels are equal across the three final witnesses. The all-word
   next-level failure and marked reentry classification remain inherited.
6. **MRB scope.** The raw nth-root alternating series is not used: its terms
   do not tend to zero. Every series here uses f(n)=n^(1/n)-1. The period-four
   filters have mean zero. The difference in (9) is absolutely convergent;
   it is not an illicit rearrangement. The derivative of log(x)/x makes
   every omitted term negative. Rational root enclosures, not floating
   point values of MRB, certify the sign. Equation (10) shows exactly how
   much information this analytic selector retains and no more.
7. **Clocks and closure.** MRB's n is the declared query clock. Square
   dilations act on different integer labels. No substitution identifies
   them. Opposite C3 phases and the AP triple anchor remain hostile;
   repeated rescaling supplies no countable embedding or universal spectrum.

## Independent computational paths

`verify.py` enumerates old permutations and their unique exposed bad step
to obtain the cut coefficients. A separate new-permutation algorithm
enumerates every legal incoming mask. They agree on all33,866 cuts from
all1,099 labelled tournaments of orders1..5. The omitted pair is recovered
in every source. All637 nonempty subsets of1..10 of size at most5 have the
full arithmetic bank checked, including candidate-candidate nonedges.

`independent_audit.py` uses direct full Hamilton-permutation counts for the
supplied rank witnesses and branch graphs; the producer uses subset DP on
the branch graphs. Rational Gaussian elimination checks the ranks1,2,4,10,15
of the producer's modularly selected minors. The prime rank modulus is
trial-divided in the producer. Native command replay checks the inherited
physical ports. Integer powers independently check all six root intervals.

Both programs use explicit exceptions, not assertions removed by `-O`.
Normal and optimized executions reproduce their respective frozen JSON
outputs byte for byte. `git diff --check` and local-link checks are part of
close-out. The script hash is retained in the producer receipt.

No demonstrated false canon claim was found; no correction entry is
necessary. The three candidate rediscoveries (scaled pair universality,
quadratic ear response, and one-birth triangle blindness) were identified
before promotion and are credited as inherited throughout.
