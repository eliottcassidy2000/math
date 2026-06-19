# LRC14 Correction Rank Scaling

The correction-rank prompt sharpened the additive-energy route, but it also
knocked out a tempting simplification.

The useful split is not "visible relations are dangerous, hidden relations are
safe."  The exact scan shows that relation rank behaves like capacity, while
the correction coefficient is a signed/coset phase.  AP rows are still the top
envelope, but the k=8 row

```text
E = (0,1,3,5,7,9,11,13)
```

has visible rank `0`, hidden rank `5`, and `Corr_y=0.215709575`.  It is an
odd-coset AP shell: the pair sums are hidden, but the hidden shell phase is
still positively aligned.  This is the same "two bounded signs can oppose"
phenomenon seen in the reciprocal-lift work, now on the sector-correction
side.

The KPS third-pocket rows supply the other guardrail.  Their pair-sum ranks are
not large enough to explain their relation coverage, but their bounded weighted
relation rank is full.  Ordinary pair sums are therefore only the mass-2 slice;
the proof needs weighted summand fibres before the signed phase is assigned.

The current theorem-shaped target is:

```text
Corr_y(E) <= sum_P phase_coeff(P) * rank_capacity(P) + peel_error(E),
```

where `P` ranges over typed relation packets: total Freiman rank, visible fold
rank, hidden shell rank, bounded weighted rank, parity/coset sign, and
multiplicand visibility.  Rank supplies the number of possible correction
channels.  Phase decides whether those channels reinforce, cancel, or feed a
finite AP/GAP pocket.

This also explains why the Fourier offset-lattice rank from the earlier
analytic route cannot be the scalar answer: for one-parameter integer rows its
rank is essentially fixed.  Its covolume and short vectors matter, but the
variable structural object is the Freiman/weighted summand rank packet.

The post-rebase KPS S12 note adds an important reconciliation: the mod-27
antipodal summand graph is inert for the small binding clusters because their
elements sit below `27/2`.  The "hidden shell" in this scout is not that
antipodal residue graph.  It is the integer pair-sum or weighted-fibre shell
inside `E+E`.  Mod 27 diagnoses why n=14 is the hard row; the small-cluster
correction is carried by genuine integer relations and their signed/coset
phase.

The next proof obligation is a signed/coset packet lemma, not another scalar
energy theorem:

1. low weighted rank gives independent-limit or peel;
2. high rank plus AP/odd-coset phase routes to the finite near-AP table;
3. high rank plus cancelling phase gives the dimension-penalty tail;
4. low pair rank but high weighted rank goes through HYP-2637's weighted-fibre
   GAP split.
