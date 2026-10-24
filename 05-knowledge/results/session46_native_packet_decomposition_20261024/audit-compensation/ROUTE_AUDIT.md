# Independent audit: restored native charged compensators

**FINITE-EXACT / PASS within the stated partial-packet scope.** The audit imports no producer code. The complete native domains are independently reconstructed in `audit.py`; `audit_routes.py` then checks the producer's actual physical matrices, owner bijections, all four K/converse presentations, entire accepted state fibres and every direct image equation against the original graph.

## Results

| Base circuit | Added restored buffers | Complete binary states | Legal states | Exact slot bottleneck | Exact ordinary-clock bottleneck |
|---|---|---:|---:|---:|---:|
| 16 | none | 65,536 | 2 | 16 | 32 |
| 16 | 3015, cost 3 | 131,072 | 6 | 6 | 13 |
| 16 | 1394, cost 2 | 131,072 | 12 | 8 | 16 |
| 16 | 1394 and 2905, both cost 2 | 262,144 | 25 | 6 | 12 |
| 18 | none | 262,144 | 2 | 18 | 43 |
| 18, original cost 43 | 1394, cost 2 | 524,288 | 19 | 6 | 17 |
| 18, original cost 37 | none | 262,144 | 2 | 18 | 37 |
| 18, original cost 37 | 2071, cost 2 | 524,288 | 38 | 6 | 12 |

Here a legal state has the same target margins and every original nodewise dyadic-current signature as the enlarged starting packet. It need not be harmonic on its own. Exact bottlenecks come from exhaustive reachability at every smaller threshold in the complete authenticated menu, not merely the largest edge on a chosen route. Both slot and actual owner-clock edge weights are checked independently. The advertised routes attain these minima except that the last producer route uses four moves and has clock maximum 13; the independent audit finds a five-move route attaining the exact clock minimum 12.

For the original 16-source base, bit positions use the producer's source order, independently authenticated against the original matrices. Bit zero means the source's unique pinned U=1965 path; bit one means its unique pinned V=2234 path. Buffers are appended.

- Buffer 3015: 65280 → 85840 → 92372 → 100534 → 255. Slot costs 6,6,6,6; actual owner-clock costs 13,12,12,13.
- Buffer 1394: 65280 → 73058 → 66471 → 255. Slot costs 8,8,6; actual owner-clock costs 16,16,12.
- Buffers 1394,2905: 65280 → 260104 → 210993 → 86071 → 255. Slot costs 6,6,6,6; actual owner-clock costs 12,12,12,12.

For the 18-source base, the appended buffer is 1394:

- 261632 → 437825 → 274217 → 272671 → 511. Slot costs 6,6,6,6; actual owner-clock costs 13,17,14,15.

The distinct cost-37 eighteen-source base, with appended buffer 2071, has the advertised route 261632 → 351784 → 479592 → 332150 → 511, with clock costs 12,12,12,13 and total cost 49. Exhaustive weighted reachability finds the sharper clock-bottleneck route:

- 261632 → 453144 → 428849 → 360765 → 360700 → 511.
- Slot costs 6,6,6,4,6; clock costs 12,12,12,9,12; total cost 57.

The cost-three source 2835 changes in the four-slot step. This is a genuine trade-off between fewer moves/less total churn and lower maximum simultaneous clock support. All states of this second path are covered by the complete literal-current audit. The original twenty direct image vectors certify the initial advertised routes. Nine additional vectors now certify the new minimum-work and minimum-peak variants; those exact equations and the additional weighted path are independently verified in VARIANT_AUDIT.md.

Every buffer starts at its U option and returns to exactly that same actual resource, physical owner lift and literal closure. Returning to a resource later is not duplicate simultaneous use.

## Physical rows, not only compressed labels

The producer's physical owner permutation is checked to be an actual relabeling of the independently reconstructed canonical source matrix. All common and private fixed clocks agree with the exhaustive-domain result. The audit reconstructs each K partner and each converse matrix, verifies their Hamiltonian paths, and computes their four words directly. It never substitutes a charge label for an actual native path.

Every accepted state has one distinct normalized resource per slot and exactly four different closed ports per slot. The entire selected 28-source two-target palette contains 56 different resources, or 224 ports. The route union sizes are 136,136,144,152 and 152 ports for the five advertised compensated routes, respectively. No hidden automorphism lift is counted as capacity; the independent foundation checks all source automorphism groups are trivial.

For each of the five authentic order-256 generators, the audit constructs a full 3,528-entry current vector from each state's literal black incidences, using the genuine multiplicity two. It compares every physical metagraph row with the baseline. Each step is then checked again by independently assembling its changed conductances and verifying that this produces exactly the difference of the endpoint current vectors.

The complete state enumeration is driven by actual target-row current contributions, not by the producer's subset-sum code. Each margin-correct rejected state has a nonzero residue at the actual target row V for an authentic original harmonic generator. This establishes the complete accepted fibre, not only validity of the advertised route.

The chosen ordinary per-slot X phase modulo eight is retained because X[U]=X[V]=146. No statement is made about every ordinary individual-port harmonic character.

## Exact images and the whole dyadic tower

All **20** supplied full image vectors are multiplied anew through the **438,216 original edge rows**. At each of 3,528 physical rows the audit verifies the exact integer equation

    L z = charge + 256 correction.

These include the common target charge 2(eU−eV) and the source charge of every compensated step. Grounding, integer ranges and zero sums are checked. The step charges are also summed explicitly to recover the original unbuffered endpoint charge for the 16-source circuit and both 18-source circuits.

The separate independent Session 45 filtration audit was copied without mathematical changes as `independent_dyadic.py` and rerun against these same original bytes. Both ordinary and optimized Python pass **530,454 checks** with byte-identical receipts. It independently proves the complete original 2-primary group and its doubled character coordinates, including saturation beyond modulus 256. Thus the current checks and row-image certificates preserve the full original **nodewise dyadic bank at every dyadic modulus**, not just five unexplained tests. They do not claim odd-primary preservation or equality of the old and new kernels. Evaluating one unchanged integer X at every larger modulus would be a different, generally false demand.

## Which minimum extends to the larger supply?

The selected-menu minimum eight for the one-buffer-1394 16-source packet is not a universal minimum: adding buffer 2905 lowers it to six.

In contrast, six is a true lower bound throughout the complete 28-source common-clock supply with the original target margins. The foundation independently enumerates all 27 four-slot trades, and every such trade preserves

    I = b1638 + b1974.

Distinct source codes exclude two-slot trades. Target margins force even changed support. The original 16-source endpoint values of I are 0 and 2; the 18-source endpoint values differ by one (0 and 1 within that base). Any other supplied sources restored at the endpoints add equal amounts at both ends. The sole extra raw pinned path, 1923→1559, cannot occur because the required target-1559 margin is zero.

Consequently neither endpoint pair is connected by moves of at most four slots anywhere in this authenticated 28-source universe, with any restored buffer initialization or any number of steps. The exhibited six-slot routes make the slot bound sharp there. This is not a lower bound after changing pins, target margins, source supply, or capacities; the displayed ordinary-clock minima are asserted only in their respective selected menus.

## Hostile controls and receipts

Normal and optimized route audits each pass **459,554 explicit checks**, and their JSON receipts are byte-identical. Four additional hostile controls, in ordinary and optimized Python, are rejected at the expected guard:

1. Corrupt a physical owner lift.
2. Rename one actual resource as another supposed choice.
3. Insert an illegal intermediate state.
4. Corrupt one coordinate of a full Laplacian image vector.

Only temporary copies are changed by these controls. The foundation independently passes 841 checks in each mode. A separate exhaustive single-buffer audit passes 3,436 checks in both modes: all twelve outside sources and both starting bits are tested. Each of the four cost-two buffers has slot minimum eight and clock minimum sixteen; the minimum clock bottleneck over all twenty-four one-buffer cases is thirteen, attained by 3015 and 3360. The final manifest pins the code, reports and all receipts; the receipts also pin all original and producer inputs.

The separate `audit_contract.py` passes 3,629 checks in each mode and pins the root guard/reuse interpretation note. Direct multiplication verifies the sparse countercharacter chi=4e1965 modulo eight: U has 406 nonloop incidences, all even, with total weighted degree 926. It changes every separately marked ordinary U/V phase between four and zero. This precisely certifies the stronger-per-port obstruction rather than merely disclaiming untested characters.

The separate event/work optimization and nine extra image-vector audit are recorded in VARIANT_AUDIT.md.

## Reproduction

Let REPO be the original repository and BUNDLE its result directory `05-knowledge/results/session46_native_packet_decomposition_20261024`. After reconstructing `audit-compensation/normal.json` with the foundation command in AUDIT.md, run from BUNDLE:

    python audit-compensation/audit_routes.py --repo REPO --packet compensated-decomposition/packet.json --foundation audit-compensation/normal.json --images compensated-decomposition/image_witnesses.json --extra-packet compensated-decomposition/packet_eighteen.json --extra-images compensated-decomposition/image_witnesses_eighteen.json --extra37-packet compensated-decomposition/packet_eighteen37.json --extra37-images compensated-decomposition/image_witnesses_eighteen37.json --out audit-compensation/routes.normal.json

Repeat with `python -O` and the optimized output name, then compare the JSON bytes. To reproduce the independent filtration, run `independent_dyadic.py --repo REPO --packet REPO/05-knowledge/results/session45_injective_multisource_moves_20261023/dyadic-constraints --out OUTPUT`. The hostile checker accepts `--repo`, `--packet`, `--images`, `--foundation` and `--out` with the corresponding paths above.

## Strict temporal nonreturn is a different requirement

The advertised routes have independently reconstructed resource-reentry event counts 7,5,6,5,5, in the order above. A reentry means reacquiring a normalized resource that was occupied earlier and then released. Consequently **none of the compensated routes satisfies strict temporal nonreturn**, even though every individual state is capacity-one injective and all buffers return exactly.

If “without reusing ports” means never reacquiring a released resource, the primitive obstruction survives: each restored auxiliary binary slot must remain constant, and each original binary slot can change at most once. Every intermediate original state would then be a conformal subset of the endpoint changes. Because each of the three original complete pinned domains has only its two guarded endpoints, the first nontrivial move must change all sixteen or eighteen original slots together. This proof applies throughout the complete authenticated supply with restored buffers and fixed target margins; changing pins or supplying new path alternatives changes its hypotheses.

Thus the positive result is reversible unit-capacity borrowing, not strict no-reacquisition. The negative result under strict temporal nonreturn is equally substantive and must accompany any user-facing conclusion.

## Scope and unresolved completion

These are enlarged selected partial observer packets. An unchanged global remainder must avoid **every normalized resource visited by the route**, retain the declared compatible owner coordinates, and supply the same complementary current signature for the enlarged baseline. This audit constructs no such remainder in the fully occupied canonical atlas. It also does not assert that all induced arcs of different source tournaments agree in one common host.

The mathematical distinction is substantive: neutral copies cannot pay an electrical debt, while actual compensators with new charge signatures can split a primitive circuit and return. Existence of such injective native supply at arbitrary order remains an additional theorem to prove, not a consequence of these finite examples.
