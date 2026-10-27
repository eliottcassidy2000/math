# Symmetry-first censuses, forced structure rather than forced symmetry, and labelled versus class distance

*klein-2026-10-27, the night after [order-ten-types-fall-out](order-ten-types-fall-out-and-what-eleven-still-withholds-klein-20261027.md).
A REFLECTION: it records why things worked; every claim points to
[THM-4818](../01-canon/theorems/THM-4818-anti-automorphism-cycles-are-2-mod-4-order-11-prime-types-by-symmetry-order-12-ledger-dipole-census-condensation-criterion-and-p7-class-distance-six.md),
the [packet](../05-knowledge/results/order11_prime_types_n12_ledger_dipole_census_20261027/README.md) or the
[external note](../05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md) and carries
that file's status.*

## 1. Enumerate the symmetry, not the population

The order-10 types took an eight-part census of 9.7 million classes; the order-11 types were expected to take hours
of a 904-million-class run (HYP-9422 (4)). They took seven seconds, because the question only concerns the
symmetric classes, and a symmetric class is an invariant tournament of a permutation of a known cycle type: one bit
per orbit on pairs, then one canonical form per candidate. Two lemmas cut the list of permutations (a lone prime
cycle is a module; an anti-automorphism's cycles have length `2 (mod 4)`), and the second lemma was found by the
generator itself: its consistency assertion fired on the types `8 2 1` and `6 4 1`. The rule is older than this
repository (Burnside counts of self-converse tournaments use the same cycle condition), but the repository had only
"even", and the sharper law removed four of seven anti-types. **When a census is wanted for the classes with a
property that forces a symmetry, enumerate the symmetry's fixed points and deduplicate; the population is never
needed.** The validation then comes from identities that see the whole population without listing it: the
self-converse decomposition fixes the self-converse count from one OEIS term, and the labelled Burnside defect
`P(n) − P_lab(n)/n! = Σ(1 − 1/|Aut|)` is a rational number that any miscounted group order breaks.

## 2. "Forced symmetry" is the wrong name for most non-exactness

The reflection of the night before says that where a distance bound fails to be exact, a forced symmetry explains
it. On every class of orders 4–8 the picture is three different pictures. For `λ` and the non-strong locus there is
no failure. For `δ⁰` and the singleton locus the failures are exactly the classes in which the condensation forbids
the needed endpoint (small score must start in the initial strong component, large score must end in the terminal
one); one failing class at order 8 is rigid and not self-converse. That is a forced **structure** (the strong
components) and the criterion is a one-line consequence of the dipole law plus Rédei's observation that a
Hamiltonian path crosses the condensation once. For `L` and the transitive class the failures are generic: 748 of
the 904 local minima at order 8 are rigid, non-self-converse and non-strong. No symmetry explains a generic
failure; what explains it is that the transitive class is a point, not a locus, and the best linear bound (the
fractional dipole LP, whose dual certificates are all the linear potentials at once) already has gaps at order 6
on rigid classes. The one constraint the LP cannot see, that the kept path be directed in the current tournament, is
the whole difficulty. **Ask first whether the target is a locus with a descent certificate or a point; symmetry
explains exceptions only in the first case, and even there "structure" may be the better word.**

## 3. Say *labelled* or *class*

`P_7` is seven moves from the transitive tournament in THM-4810 and six in THM-4818; both are right. The seven is
the distance from one labelling of `P_7` to one labelling of `TT_7`; half of the 5,040 transitive labellings are six
away. A potential that is an isomorphism invariant can only be judged against the class distance, and the table
that called `L` inexact at `P_7` (MISTAKE-682) had judged it against a labelled one. The same confusion is available
for every "distance to X" statement in the metagraph work, where the labelled graphs have eccentricities 13 and 17
and the class graphs 6 and 7.

## 4. External preprints annotate, they do not close

Three preprints of the openai/math collection claim `N(6) = 3`, Sogge's conjecture with Kakeya in every dimension,
and the refutation of ETH. The repository's practice from the 2026-10-06 reading holds: type every contact, change no
status. The useful contacts were small and exact: the qutrit overlap `(1+2ω)/3` is the Eisenstein unit of the
directed 3-cycle; the dipole `e_a − e_b` is a root of `A_{n−1}` and `L` is the root-word metric, which is what the
QCD/root analogy had said in other words; nothing in the three papers touches a finite object of the repository.
Two readers stalled on 67-page PDFs until told to save notes after every ten pages; the cheap fix is to make the
reader write as it reads.
