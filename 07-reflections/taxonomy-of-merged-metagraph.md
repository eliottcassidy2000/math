# Taxonomy of the Merged Meta-Graph G_n/Z_2

*opus-2026-03-23-S236*

## The Five Edge Colorings

Every edge in G_n/Z_2 carries five independent binary labels:

| Labeling | Categories | What it measures |
|----------|-----------|------------------|
| **SC-type** | BLUE / BLACK / GREEN | SC↔SC / SC↔NS / NS↔NS (structural position) |
| **Score** | AMBER / VIOLET | Same score sequence / different (invariant preservation) |
| **Fiber** | TEAL / MAGENTA | Same parent at n-1 / different (Mode A recursion) |
| **Step size** | STEP / LEAP / LEVEL | \|dH\|=2 / \|dH\|>2 / dH=0 (gradient magnitude) |
| **c3 parity** | CROSSING / PRESERVING | c3 parity changes / stays same (3-cycle layer) |

## Structural Laws Revealed

### Law 1: TEAL = 100% at n=5,6
Every edge in the merged meta-graph connects nodes sharing at least one parent class at n-1. No MAGENTA edges exist. This means the meta-graph is entirely contained within the fiber overlap structure of Mode A recursion.

### Law 2: BLUE edges are c3-PRESERVING at n=6
At n=6, 100% of BLUE (SC↔SC) edges preserve c3 parity. At n=5, 75% cross (the blue bipartition). This parity law reverses between even and odd n.

### Law 3: AMBER is rare on BLACK edges at n=5 (0%) but present at n=6 (20%)
At n=5, no SC↔NS edge preserves the score sequence. At n=6, 1 in 5 do. The score sequence is more "stable" under the black ribs as n grows.

### Law 4: SC nodes have NO GREEN edges
SC↔SC is BLUE, SC↔NS is BLACK. GREEN only connects NS↔NS. This is definitional but creates a clean three-zone partition: nodes are either "all BLUE+BLACK" (SC) or "GREEN+BLACK" (NS).

### Law 5: Dominant edge profile evolves
- n=5: Most common = BLACK + VIOLET + TEAL + LEAP + THICK (3/21)
- n=6: Most common = GREEN + VIOLET + TEAL + LEAP + c3_CROSS + THICK (31/143 = 22%)

The typical edge at large n is: GREEN (NS↔NS), VIOLET (score changes), TEAL (same parent), LEAP (|dH| > 2), c3-CROSSING (parity flips), THICK (high multiplicity). This is the "generic" edge — connecting two non-symmetric tournaments that differ in score, cross c3 parity, but share an ancestor.

### Law 6: SC nodes are all c3_EVEN at n=6
At n=6, every SC node has even c3 count. At n=5, SC splits 50/50. This relates to the bilateral structure at odd n.

## Node Categories (27 types identified)

### Primary axes (3):
- SC / NS (self-complementary vs not)
- RIGID / SYMMETRIC (trivial vs nontrivial automorphism)
- HAS_SELF_LOOP / NO_SELF_LOOP (class is self-adjacent)

### Connectivity axes (5):
- HAS_BLUE_EDGE / NO_BLUE_EDGE (connected to SC backbone)
- HAS_BLACK_EDGE / NO_BLACK_EDGE (bridges SC-NS boundary)
- HAS_GREEN_EDGE / NO_GREEN_EDGE (part of NS sea)
- HAS_AMBER_EDGE / NO_AMBER_EDGE (has score-preserving neighbor)
- HAS_TEAL_EDGE / NO_TEAL_EDGE (has same-parent neighbor — always YES at n≤6)

### H-structure axes (3):
- FLOOR / BULK / CEILING (H-value extremal or not)
- SINGLETON_LEVEL / MULTI_LEVEL (unique H or shared)
- c3_EVEN / c3_ODD (3-cycle count parity)

### Recursion axes (2):
- SINGLE_PARENT / MULTI_PARENT (unique or multiple parents at n-1)
- UNIQUE_SCORE / SHARED_SCORE (score sequence distinguishes this node)

## Key Compound Categories

| Category | n=5 | n=6 | Meaning |
|----------|-----|-----|---------|
| SC + NO_BLACK (isolated SC) | 2/10 | 0/34 | SC nodes with no NS neighbors — vanish at even n |
| SC + RIGID | 5/10 | 9/34 | Most SC nodes are rigid |
| NS + NO_BLUE | 2/10 | 22/34 | NS nodes only in GREEN+BLACK zone — grows to dominate |
| SC + HAS_AMBER | 5/10 | 9/34 | SC nodes with score-preserving neighbors |
| MAX_DEGREE | 1/10 | 3/34 | Hub nodes |

## Boolean Node Profiles

At n=5: **10 distinct profiles for 10 nodes** — every node is distinguishable by its boolean flags alone.

At n=6: **17 distinct profiles for 34 nodes** — profiles start collapsing. The most common profile (5 nodes) is: SC + RIGID + LOOP + BLUE + BLACK + AMBER + TEAL + c3_EVEN.

## Edge Profile Distribution

At n=6, 30 distinct edge profiles exist among 143 edges. The top 3 profiles account for 47% of all edges:
1. GREEN + VIOLET + TEAL + LEAP + c3_CROSS + THICK (31 edges, 22%)
2. GREEN + VIOLET + TEAL + LEAP + THICK (19 edges, 13%)
3. BLACK + VIOLET + TEAL + LEAP + c3_CROSS + THICK (17 edges, 12%)

Profile #1 is the "generic" edge of the meta-graph — the bulk of the NS sea.

## What Changes from n=5 to n=6

| Feature | n=5 | n=6 | Trend |
|---------|-----|-----|-------|
| SC fraction | 80% | 35% | Decreasing (SC drowns in NS) |
| GREEN edge fraction | 5% | 59% | Exploding (NS sea fills) |
| AMBER edge fraction | 14% | 18% | Slowly increasing |
| LEVEL edges | 1 (5%) | 5 (4%) | Stable low fraction |
| c3_CROSSING fraction | 67% | 52% | Decreasing toward 50% |
| Distinct profiles/V | 100% | 50% | Nodes becoming less distinguishable |
| SELF_LOOP fraction | 50% | 56% | Slowly increasing |

The merged meta-graph transitions from an "SC-dominated small world" at n=5 to an "NS-dominated sea with SC islands" at n=6.
