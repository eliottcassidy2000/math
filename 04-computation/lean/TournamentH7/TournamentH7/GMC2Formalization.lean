import TournamentH7.GMC2Reduction
import TournamentH7.GMC2ChargeGeometry
import TournamentH7.GMC2WickChannels
import TournamentH7.GMC2MomentRelations
import TournamentH7.GMC2MomentTransport
import TournamentH7.GMC2ConstantTermRelations
import TournamentH7.GMC2IntegralRelations
import TournamentH7.GMC2NormalizedMoment
import TournamentH7.GMC2AlgebraicDescent
import TournamentH7.GMC2TorusDescent
import TournamentH7.GMC2NullconeDescent
import TournamentH7.GMC2SupportDescent
import TournamentH7.GMC2IntegralSpecialization
import TournamentH7.GMC2IntegralTorusSpecialization
import TournamentH7.GMC2LowestFaceExistence
import TournamentH7.GMC2LowestFaceBridge
import TournamentH7.GMC2LowestFacePackage
import TournamentH7.GMC2FrobeniusFace
import TournamentH7.GMC2FaceDictionary
import TournamentH7.GMC2FaceHeightFloor
import TournamentH7.GMC2DvdKInterface
import TournamentH7.GMC2FaceSeed
import TournamentH7.GMC2FaceSeedChannel
import TournamentH7.GMC2FaceReferenceChannel
import TournamentH7.GMC2FaceSeedDescent
import TournamentH7.GMC2IntegralFaceSeedDescent
import TournamentH7.GMC2ChannelDilation
import TournamentH7.GMC2GoodReduction
import TournamentH7.GMC2FrobeniusResidue
import TournamentH7.GMC2ResidueAssembly

/-!
# Formalization spine for NC2 and GMC(2)

This aggregation module exposes the checked components of `THM-2022` through
one stable import.  It deliberately does **not** assert an end-to-end `NC2`
theorem until the last concrete channel assembly has been proved.

The present spine has seven layers.

1. `GMC2Reduction` and `GMC2ChargeGeometry` define `NC2`, prove both strict
   charge branches, and prove `NC2 -> GMC(2)`.
2. `GMC2WickChannels`, the rational/integral relation modules, and
   `GMC2NormalizedMoment` give the exact Wick expansion and the integral
   minimum-factorial cancellation used before reduction.
3. The lowest-face modules construct a rational exposing line and exact
   equality face, identify polynomial coordinates with charge/radial data,
   prove the balanced height floor and strict integral off-face gap, and
   extract a concrete reference channel from a nonzero seed.
4. `GMC2DvdKInterface` records the one-variable Duistermaat--van der Kallen
   input as an explicit premise; `GMC2FaceSeed` derives the nonzero lowest-face
   constant term without installing a custom axiom.
5. The algebraic modules provide two seed-preserving descents: a number-field
   route and a shorter direct finite-field route through a finite-type
   `Z`-algebra.  The direct route universally preserves every integral zero
   relation, so its residue characteristic may be learned before selecting
   the characteristic-dependent normalized moment.
6. `GMC2ChannelDilation` formalizes extension by zero, componentwise
   dilation, all mass/charge/height/coefficient scaling laws, the exact image
   in `piAntidiag`, and finite-sum reindexing.
7. The Frobenius modules prove non-dilated multinomial isolation, multinomial
   Lucas, whole-face `Q^p`, the strict factorial-gap divisibility, and an
   abstract three-case finite-sum assembly.  `GMC2GoodReduction` retains the
   independent number-field reduction route, although it is no longer on the
   shortest critical path.

The remaining composition task is to instantiate the abstract three-case
assembly with the normalized Wick channels and the exact dilation image, then
feed its contradiction through `GMC2IntegralFaceSeedDescent` to obtain the
conditional theorem `DvdK1 -> NC2`.  Formalizing DvdK itself is a separate
published-theorem project; here it remains visible in theorem hypotheses.

Tournament-assumption audit: orienting channels pairwise preserves the
lower-face/off-face comparison but loses coefficient magnitudes and tied-face
cancellation.  The proof must retain the entire equality face as a finite sum;
it cannot replace that hyperedge by a tournament winner.
-/
