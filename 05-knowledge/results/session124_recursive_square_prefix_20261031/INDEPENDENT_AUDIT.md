# Dependency-minimal core: independent mathematical audit

**PASS.** The unscaled-prefix cycle theorem needs only:

1. The11 existing compressed B-base chunks and their manifest, covering every index53 through2620;
2. The new98 recursive paired B certificates;
3. The new98 direct terminal-cycle certificates;
4. Their integer replay/decoding code and the elementary induction below.

No old B transition, A transition, final-cycle transition, or original Gerbicz program is required for this core. The original independent audit checked additional older machinery because the initial target included simultaneous-copy corollaries; those extra replays are not logical dependencies of the minimal core.

## Closed recursive constructor

Define B(n), for n>=53, as follows:

- If53<=n<=2620, decode the certified base pair.
- Otherwise put `m=floor((n-24)/49)` and `d=n-49m`, recursively construct B(m), and apply the new paired certificate at `(m mod2,d)`.

For n>=2621, one has53<=m<n and24<=d<=72. The recursion therefore terminates at a certified base. Every base is a weak B pair: square-sum Hamiltonian members on consecutive intervals, common positional parity, prefix1,3, and terminal8 for even member length or10 for odd member length.

Crucially, the **new** paired table alone maps every weak B pair to another weak B pair. Its longer output prefix begins1,3; all other weak B invariants are exactly those checked by the independent universal replay. Thus strong induction on n proves B(n) exists for every n>=53. There is no assumption of an old all-size supply and no circular dependency.

## Guaranteed depth and prefix

Let `T_l=(107*49^l-1)/2`, with T0=53 and `T_(l+1)=49*T_l+24`. The same fixed constructor B(n), not a different existential family for each l, has common prefix P_l whenever n>=T_l.

For l=0 this is the weak B invariant. For l>=1, n>=T_l implies n>=2621, so the recursive branch is used and its input index satisfies m>=T_(l-1). By induction in l the recursively constructed input has prefix P_(l-1). The exact five-atom opening of the new paired table yields `(1,3,6,19),D4(P_(l-1))=P_l`. This argument proves the minimum number of recursive layers from the size bound; no uniform exact depth is assumed.

For an arbitrary N>=T_l, l>=1, decompose N=49m+d in the same way. Then m>=T_(l-1), so the fixed constructor supplies B(m) with prefix P_(l-1). The direct terminal table yields a Hamiltonian cycle with the original unscaled P_l. This proves the core theorem and the eventual extension of every finite prefix of the explicit ray.

## Audit evidence inherited by this reduced proof

- All2568 base pairs were independently decoded and fully checked in the main audit, with compressed hashes matched to the certified manifest and Git baseline.
- All98 new paired certificates were independently replayed with universal ownership, endpoint, symbolic-length and parity semantics; complete two-layer concrete expansions were checked separately.
- All98 direct-cycle certificates were independently replayed with terminal semantics, plus196 full concrete cycles and hostile mutations. Normal/optimized receipts match.
- These checks do not rely on correctness of any old transition to establish the new certificates' implications.

The existing main audit scripts also perform extra checks for the wider theorem package; that is an executable convenience, not a mathematical need for those old inputs. A dependency-minimal publication constructor may omit every old transition object without weakening this core argument.

## Explicit exclusions

The2303 simultaneous-copy theorem and the48*49^(r-1) depth/width amplification corollary still depend on the identified old B/A/final transitions. They must be presented outside the minimal core, with their separate dependencies. The new direct-cycle pipeline does not by itself assert the482-based seed or those copies on the same cycle. Any additional common-intact-offset bonus is likewise not needed for this core audit.
