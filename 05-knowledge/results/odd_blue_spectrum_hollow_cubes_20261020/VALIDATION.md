# Validation receipt

Status: PASS, independently audited; 2026-10-20.

## Uniform mathematical result

For every odd n>=3, the self-converse classes have equal blue-degree
spectra on the two cyclic-triangle-parity shores. The proof supplies an
actual class involution, and more strongly partitions the classes into
free hollow-frame cube packets of constant blue degree. At n=1 mod4,
n>=5, each blue-degree multiplicity on each shore is even.

Both proof audits approve the anti-involution marking loss, exact
centralizer quotient, odd hollow-frame groups, equivariant orbit
selection, class-level involutivity and free full-cube strengthening.
The exact integer matching-cube c3 law was also independently checked.

## Producer's checks

Command:

    python3 verify_odd_pairing.py --inherited <full-census-directory>

All normalized reflection-frame completions at n=3,5,7,9, respectively
2,16,512,65536, were classified. The SC class counts are 2,8,88,2752.
A separately coded fixed-Hamiltonian-path blue-cube enumeration checks
every class's blue count. The inherited full census is cross-checked on
all SC classes through n=7 when supplied.

The complete auxiliary cube counts are:

- n3: one Q1.
- n5: two Q2.
- n7: four Q1 and ten Q3.
- n9: sixteen Q2 and 168 Q4.

The main verifier checks every selected partner and every full cube,
including stabilizer preservation, symmetric-path preservation,
B preservation, parity reversal, free action, exact triangle affine
law and positive odd hollow blue factor.

## Independent audits

`audit.md` is a separate proof audit with the explicit n7 nonadjacency
and purity-change counterexample.

`audit_independent/AUDIT_RECEIPT.md` is the second proof and independent
computation audit. Its separate bit convention and three programs import
no producer/census classifier. It reproduces all SC classes and cube
counts through n9, checks every one of 98 SC representatives through n7
by all full vertex permutations, and independently checks an n9 B=39
witness by all 9! permutations. That witness also refutes an odd-only
repair of the old 3-power blue cap.

These independent programs run from their own directory:

    python3 independent_pairing_audit.py --max-m 4
    python3 independent_cube_audit.py
    python3 independent_n9_cap_witness.py

## Preserved boundaries

- B is h_rho divided by its rho-centralizer, not h_rho itself.
- Auxiliary cube edges need not be antipodal blue edges.
- The construction need not preserve purity, mass, full Aut, H or scores.
- n1 is excluded and even-order generalization is false.
- The finite n9 data are an SC blue-spectrum census, not a full
  tournament/antipodal-incidence census.
- No pure-blue classification, 3-power cap, LRC statement, or full
  metagraph reconstruction is assumed or deduced.
