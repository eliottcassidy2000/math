# Post-freeze inverse-mask advances

Status: FINITE-EXACT restricted-class constructions and hostile probes.
Chan–Pak completeness and eventual coverage remain OPEN.

The separately frozen `session216-inverse-release` is unchanged by this work.
Generic run transfers, homogeneous-block modes, and saturation are inherited
as credited in PROOF.md. These are not new generic two-apex counting formulas.

## Sharp finite-order coverage boundary

`census_depths.cpp` exhausts all four-letter mask words of length at most14,
with apex arc0->1 and the increasing chain. It records shortest mask witnesses
for every output up to4,000,000. All4^14=268,435,456 depth14 words are actually
visited in this run (no cap pruning occurred above the leaves).

Every odd integer through2,327,673 except7 and21 has a witness on at most16
vertices in this class. The least additional omission at this order is
2,327,675. First additional omissions by maximum chain length0,...,14 are:

3,5,9,17,35,61,297,727,2593,7455,24107,77533,230521,776161,2327675.

This improves only the restricted-class certificate, not Session203's much
larger unrestricted finite range. The omission is finite-order, not all-order.
The literal word120021002033010 has length15 and count2,327,675, independently
verified by subset Hamiltonian-path DP. Consequently its minimum chain length
in this exact two-apex class is15. This is not a claim about its minimum order
among arbitrary tournaments.

The raw record is `census14full.out`. Its missing_odds list gives only the
first100 omissions; missing_count=139132 counts all odd omissions up to4m.
It records1,860,868 attained odds in that full range. `coverage2327673.bin.gz` (reconstructed losslessly from the parts in
`certificate/`, see README.md) contains the complete initial-range certificate in little-endian32-bit format.
There is one record for each h=0,...,2327673. The shipped certificate is
explicitly little-endian. The census C++ binary writer uses host byte order;
the supplied binary-comparison replay therefore requires a little-endian host. Value0xffffffff means absent.
Otherwise the high four bits store length n; successive low two-bit pairs
store the mask letters in increasing chain order. Hence every claimed odd
has a literal graph, not merely a numerical coverage bit.

Correctness/overflow: for each of the two apex permutations, each coordinate
counts assignments of the n chain vertices to3 gaps, so it is at most3^n.
At n<=14 this is4,782,969, safely below uint32_t capacity. Final totals use
uint64_t. At most28 bits encode the word, and lengths<=14 cannot collide with
the absent sentinel. The positive transfer is proved in PROOF.md; the depth
pruning H>cap is its proved induced-monotonicity pruning. Recording the least
visited length gives the exact per-target minimum in the enumerated universe.

`audit_large.py` independently checks1,200 deterministically selected saved
witnesses from the earlier one-million cap, including first/last100 and1,000
seed216 random targets, with the separate subset-DP audit.cpp. All passed.
This is a sampled implementation control, not a claim to have individually
subset-DP checked all1.16m graphs. The first release independently checked
all4,998 witnesses through10,000. The sharp boundary witness is also fully
independently checked. Mathematical exhaustiveness comes from the complete
word traversal, not sample extrapolation.

Reproduce the census from this directory on a little-endian host
(full certificate checks and temporary-file cleanup: run `sh replay.sh`):

    g++ -O3 -std=c++17 census_depths.cpp -o /tmp/session216-census-depths
    /tmp/session216-census-depths 14 4000000 /tmp/session216-census14full.bin > /tmp/session216-census14full.out
    cmp census14full.out /tmp/session216-census14full.out
    python3 reconstruct_certificate.py /tmp/session216-coverage2327673.bin.gz
    gzip -cd /tmp/session216-coverage2327673.bin.gz > /tmp/session216-coverage.bin
    g++ -O3 -std=c++17 audit.cpp -o /tmp/session216-inverse-audit
    echo '2327675 120021002033010' | /tmp/session216-inverse-audit

The search finding the boundary repair is reproducible, but unnecessary to
verify it: compile find_target.cpp and run with arguments2327675 15. Its
result is a literal independently checkable witness. A failed search would
report INCOMPLETE, never absence.

## Global-map digit probes

A context means P w S with P,S fixed. Tested actions preserve these fixed
outer contexts and replace central w by u phi(w) v. The global letter map
phi need not be invertible. Tests preserve a chosen central alphabet so the
result remains in the same language. This does not exhaust finite-state
systems with changing outer contexts.

- Full four-letter alphabet, all256 maps, total wrapper length1..2, all
  contexts of combined length4..6 and baseline H>=23:13,262 binary first-value
  candidates,2 compatible with every singleton central word, no candidate
  survives all central words of length<=3. Code:global_map_wrappers.py.
- Closed two-letter alphabets, all maps into that alphabet, wrapper length
  1..3, same context range:37,490 first-value candidates,4,684 singleton-
  compatible candidates, none survives length<=4. Code:two_letter_wrappers.py.
- Centered ternary targets3H-2,3H,3H+2, full alphabet, wrapper length1..2,
  same context range:18,441 first-value candidates,4 singleton-compatible
  candidates, none survives length<=3. Code:ternary_global_wrappers.py.

The first-value and singleton-compatible counts above count context/wrapper
choices before enumerating their compatible letter maps, not individual map
actions. All compatible maps are then tested.

These are bounded hostile tests. They neither rule out larger contexts nor
larger wrappers, context changes, richer control, higher apex count, or a
non-regular language. No exact infinite-word equivalence was promoted, since
all candidates already failed finite tests.

## Pure ternary modes are physically available

The inherited Session215 homogeneous-block theorem says that for m>=1,
H(P c^m S) is a linear combination of1,2^m,3^m. Scanning contexts of combined
length<=6 finds578 literal context/letter triples with nonzero3-mode and
zero2-mode (70,538 have a nonzero3-mode at all). This includes both small
constant signs:

H(1122 0^m 1)=60*3^m-1,
H(11022 0^m 1)=120*3^m+1,
H(220 1^m 322)=72*3^m+1.

The formulas hold for every m>=0. For m>=1 the inherited three-mode theorem
reduces proof to the three exact evaluations at m=1,2,3; the m=0 check is
separate. The first values including m=0 are respectively:
59,179,539,1619; 121,361,1081,3241; 73,217,649,1945.
They were independently checked by subset DP. The probe computes all578
triples in ternary_modes.out. These are sparse positive families illustrating
actual mode cancellation, not a ternary digit compiler. Switching between
contexts while preserving an arbitrary requested value remains unresolved.
