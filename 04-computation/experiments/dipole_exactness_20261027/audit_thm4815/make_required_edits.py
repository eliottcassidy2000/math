"""make_required_edits.py -- builds required_edits.json for the THM-4815 audit and checks that every "old" string
occurs exactly once in its file, also after the earlier edits to the same file are applied (dry run on copies in
memory; no repository file is modified).  Run: python3 -I make_required_edits.py"""
import json
import os

ROOT = '/Users/e/Documents/GitHub/math-wt-pyth-cuboids/'
THM = ROOT + '01-canon/theorems/THM-4815-dipole-principle-native-move-distance-bounds-endpoint-relation-cut-or-colour-exactness-and-the-p7-parity-correction.md'
NOTE = ROOT + '05-knowledge/results/dipole_principle_and_prime_cycle_index_tower_20261027.md'
REFL = ROOT + '07-reflections/why-the-advances-work-ii-potentials-symmetry-budgets-positive-coordinates-macmini-20261027.md'
HYP = ROOT + '05-knowledge/hypotheses/HYP-9420-strong-untouchables-1-7-21-complement-graph-components-paley-depth-quadratic-bound.md'

THOMASSEN = ('Thomassen, Proc. London Math. Soc. (3) 45 (1982) 151-168 (Theorem 5.5 of Bang-Jensen and Gutin, '
             '"Paths, trees and cycles in tournaments")')

E = []
def edit(f, old, new, why):
    E.append({'file': f, 'old': old, 'new': new, 'why': why})

# ---------------- THM-4815 ----------------
edit(THM, "strong failures at orders <= 9 all have a cut vertex",
     "strong L-local minima through order 10 all have a cut vertex, though 2-connected failures above them occur from order 7",
     "False as stated: strong 2-connected L-failures exist (1, 5, 143 classes at orders 7, 8, 9); only the local minima have cut vertices.")
edit(THM, "every tournament has a Hamiltonian path from a maximum-score to a minimum-score vertex, so at odd orders",
     "every tournament has a Hamiltonian path from a maximum-score to a minimum-score vertex (classical: a corollary of Thomassen 1982; the case of a minimum vertex beating a maximum one is Lu 2018), so at odd orders",
     "B7 is prior art (Thomassen's characterisation; Lu's first theorem is case (i)).")
edit(THM, "  the reversal lemma for max-to-min paths (given Lu 2018); (C) the forward-arc parity law and its consequences.",
     "  (B7) the max-to-min Hamiltonian path theorem (self-contained re-proof of a classical result: Thomassen 1982, case (i)\n"
     "  = Lu 2018) and (B8) the regular-class distance at odd orders; (C) the forward-arc parity law and its consequences.",
     "No 'reversal lemma' exists in the body; B7/B8 were not listed; 'given Lu' contradicted 'proved here in full'.")
edit(THM, "orders 4-8 (order 9 in the note), the endpoint-relation census",
     "orders 4-8 (order 9 in exact_n9_summary.out, completed after filing), the endpoint-relation census",
     "The note has no order-9 class-graph data (it says 'in progress at filing').")
edit(THM, "Lu, Graphs Combin. 34 (2018) 613-618 (min-max pairs; prior art for B7, abstract only)",
     "Lu, Graphs Combin. 34 (2018) 613-618 (min-max pairs: case (i) of B7 is its first theorem; abstract only); " + THOMASSEN + " (implies B7)",
     "Credit Thomassen 1982 and state precisely what Lu proves.")
edit(THM, "Session74 (source/sink\n  trap)",
     "Session74 (source/sink\n  pruning lemma)",
     "Session74 calls it a 'safe pruning lemma'; its 'trap' is a different (strong, order-11) state.")
edit(THM, "| regular class (this theorem, B6) |",
     "| regular class (this theorem, B8) |",
     "The regular-class distance is B8, not B6.")
edit(THM, "Session74 source/sink trap |",
     "Session74 source/sink pruning lemma |",
     "Session74's name for the statement.")
edit(THM, "the endpoint principle and cut confinement (B1-B3)",
     "the endpoint criterion and the exact form of cut confinement (B1-B3; THM-4818 (C2) already uses the endpoint form of the dipole law and that every word starts in the initial and ends in the terminal component)",
     "Fair credit: THM-4818 (C2) contains B1 and the easy half of B3.")
edit(THM, "the max-to-min theorem and the regular-class distance (B7, B8), the alternating half of the parity\nlaw (C)",
     "the max-to-min theorem (classical, see B7) and the regular-class distance (B7, B8), the explicit alternating half of\nthe parity law (C; implicit in THM-4810 Theorem P's odd toggle count and used in THM-4732 Lemma 1)",
     "B7 is classical; the alternating parity is a one-line consequence already used in the repository.")
edit(THM, "The \"only if\" half is Session74's source/sink trap.",
     "The \"only if\" half is Session74's source/sink pruning lemma (supplement, section 4).",
     "Session74's name for the statement.")
edit(THM, "Lu (Graphs Combin. 34 (2018)\n613-618) studies Hamiltonian paths between min-max pairs and proves them with a few exceptions (published abstract;\n"
          "the full text was not consulted). Case (i) and the existential statement are proved here in full, and no priority over\nLu's text is claimed.",
     "Prior art: B7 is classical. Thomassen (Proc. London Math. Soc. (3) 45 (1982) 151-168; Theorem 5.5 of Bang-Jensen and\n"
     "Gutin's survey *Paths, trees and cycles in tournaments*) characterised when a specified vertex `x` of maximum and `y` of\n"
     "minimum out-degree are joined by a Hamiltonian path from `x` to `y`: always, except (a) when the terminal component of\n"
     "`T - x` is one vertex `z` with `z -> x`, preceded by the component `{y}`, and `x -> y`; (b) the mirror configuration in\n"
     "`T - y`; (c) for `R(1,2)` and two small exceptional tournaments. In (a) the vertex `z` also has minimum score and beats\n"
     "`x`, in (b) dually, so the existential statement follows ((c) by direct check). Lu (Graphs Combin. 34 (2018) 613-618; published\n"
     "abstract, full text not consulted) proves case (i) for every pair (a minimum vertex beating a maximum vertex always\n"
     "gives a Hamiltonian path from the latter to the former) and the case `x -> y` up to a few exceptions. The proof above is\n"
     "a self-contained re-proof; the endpoint relation agrees with the characterisation for every (max, min) pair at orders\n"
     "3-9 (dipole_exactness_20261027/audit_thm4815/audit_thomassen_n3_9.out).",
     "Lu was misdescribed (exceptions only in the x -> y case) and Thomassen 1982, which implies B7, was not cited.")
edit(THM, "(two at\norder 6 and one at order 7, all for the regular target)",
     "(two at\norder 6 for the almost-regular target defined below, one at order 7 for the regular target)",
     "Order 6 is even: the target there is the almost-regular class with potential ceil(R'/2).")
edit(THM, "orders 4-8 (order 9 in the note). The endpoint relation",
     "orders 4-8 (order 9 in exact_n9_summary.out, completed after filing). The endpoint relation",
     "The note has no order-9 class-graph data.")
edit(THM, "(THM-4730: palindromic stacks of `K1`, `C3`, `T5`)",
     "(THM-4730 with THM-4729 L4: palindromic stacks of `{K1, C3}` halves with centre in `{none, K1, C3, T5}`)",
     "Precise form of the pure-blue classification (equivalent at orders <= 9).")
edit(THM, "The smallest rigid failure is the strong 6-tournament with scores `1,1,2,3,4,4` (`111011111111101` in gentourng\n"
          "order): `|Aut| = 1`, not self-converse, `L = 1`, distance 3 to the transitive class. Vertex 3 has the single\n"
          "out-neighbour 4, and 5 the single out-neighbour 3, so every Hamiltonian path through them contains `5 -> 3 -> 4`. As a\n"
          "result no path joins a minimum-score vertex (3 or 5) to a maximum-score vertex (0 or 1). These are the only pairs in\n`D_L`.",
     "The smallest rigid failures are non-strong: `K1 => S4` and `S4 => K1` at order 5 (`1111100111`, `1101111111`;\n"
     "`|Aut| = 1`, not self-converse, `L = 1`, distance 2; a source or a sink, B5). The smallest strong rigid ones occur at\n"
     "order 6 (10 of the 20 rigid `L`-failures there). All three strong `L`-local minima of order 6 have scores `1,1,2,3,4,4`\n"
     "and `|Aut| = 1`: `111011111111101` (in gentourng order; not self-converse, `L = 1`, distance 3 to the transitive\n"
     "class), its converse `101011111111110`, and the self-converse `111011111111110`. In the first, vertex 3 has the single\n"
     "out-neighbour 4 and 5 the single out-neighbour 3, so a Hamiltonian path starting at 3 must end at 5, and the only one\n"
     "starting at 5 is `5, 3, 4, 0, 1, 2`. As a result no path joins a minimum-score vertex (3 or 5) to a maximum-score vertex\n"
     "(0 or 1). These are the only pairs in `D_L`.",
     "Rigid failures exist at order 5; and only 5 of the 15 Hamiltonian paths contain 5 -> 3 -> 4 (the others end at 3 or 5).")
edit(THM, "near-regular and planted-tie samples, 960,000 tournaments in all (endpoint_random.c).",
     "near-regular and planted-tie samples, 1,120,000 tournaments in all (endpoint_random_n11_16.out).",
     "The output file lists 1,120,000 distinct samples (3 x 100k at n = 11-13, 3 x 50k at 14, 3 x 20k at 15, 10k at 16).")
edit(THM, "| of which symmetric (n <= 8) |",
     "| of which symmetric, `|Aut| > 1` or self-converse (n <= 8) |",
     "Define 'symmetric' (the failures' symmetric counts 2, 4, 13, 52, 173 are in the note).")
edit(THM, "read the labelled 7 as a class distance.",
     "read the labelled 7 as a class distance (the reflection and the body of HYP-9420 (4) are now corrected; Session74's\nsupplement is left as filed).",
     "The quoted sources have since been corrected.")
edit(THM, "For `n ≡ 0, 3 (mod 4)` it alternates, so the labelled move graph is bipartite",
     "For `n ≡ 0, 3 (mod 4)` it alternates (THM-4810's toggle count is odd there; THM-4732 Lemma 1 uses this for detours), so\nthe labelled move graph is bipartite",
     "Credit: the alternation is a one-line consequence of THM-4810's toggle count.")
edit(THM, "and rigid examples exist from order 6.",
     "and rigid examples exist from order 5 (strong ones from order 6); 2-connected failures that are not local minima occur\n   from order 7.",
     "K1 => S4 (order 5) is a rigid non-self-converse failure; 2-connected strong L-failures exist at orders 7-9.")

# ---------------- note ----------------
edit(NOTE, "regular class is exactly `R = (1/2) sum |s_v - (n-1)/2|` | PROVED |",
     "regular class is exactly `R = (1/2) sum |s_v - (n-1)/2|` | PROVED (B7 classical: Thomassen 1982; case y -> x is Lu 2018) |",
     "B7 is prior art.")
edit(NOTE, "Lu (Graphs Combin. 34 (2018) 613-618) studies\nHamiltonian paths between min-max pairs (abstract only consulted).",
     "Lu (Graphs Combin. 34 (2018) 613-618) proves\ncase (i) of B7 (abstract only consulted); B7 itself follows from " + THOMASSEN + ".",
     "Credit Thomassen 1982 and state what Lu proves.")
edit(NOTE, "The source/sink trap and the `P_11`\ncertificate are Session74.",
     "The source/sink pruning lemma and the `P_11`\ncertificate are Session74.",
     "Session74's name for the statement.")
edit(NOTE, "The identity is THM-4811 (A) (this machine's previous session) and Session75 section 6 (`K -> -K + 2L`).",
     "The identity is Session75 section 6 (`K -> -K + 2L`), the dyadic-module-spectrum session and THM-4811 (A) (this machine's previous session), as THM-4811 records.",
     "THM-4811 credits the matrix form to Session75 (6.1) and the dyadic-module-spectrum session.")
edit(NOTE, "The smallest rigid example is the strong\n6-tournament with scores `1,1,2,3,4,4` (gentourng `111011111111101`, `|Aut| = 1`, not self-converse). Its arcs\n"
           "`5 -> 3 -> 4` are forced, so no Hamiltonian path joins a minimum-score vertex to a maximum-score vertex.",
     "The smallest rigid examples are non-strong\n(`K1 => S4` and `S4 => K1` at order 5); strong ones start at order 6, e.g. the 6-tournament with scores `1,1,2,3,4,4`\n"
     "(gentourng `111011111111101`, `|Aut| = 1`, not self-converse). There vertex 3 can only be followed by 4 and 5 only by 3,\n"
     "so a path from 3 ends at 5 and the only path from 5 ends at 2: no Hamiltonian path joins a minimum-score vertex to a\nmaximum-score vertex.",
     "Order-5 rigid failures exist; the arcs 5 -> 3 -> 4 lie on only 5 of the 15 Hamiltonian paths.")
edit(NOTE, "**Theorem (THM-4815 B7).**",
     "**Theorem (THM-4815 B7; classical, a corollary of Thomassen 1982; the case `y -> x` is Lu 2018).**",
     "B7 is prior art.")
edit(NOTE, "This matches Lu's min-max pair theorem with exceptions.",
     "This matches Lu's min-max pair theorem (no exceptions when `y -> x`, a few when `x -> y`) and Thomassen's characterisation\n(its exceptional configurations need `x -> y`).",
     "State Lu's result precisely.")
edit(NOTE, "non-adjacency, Session74's source/sink trap);",
     "non-adjacency, Session74's source/sink pruning lemma);",
     "Session74's name for the statement.")
edit(NOTE, "two at order 6, one at order 7 for the regular\ntarget)",
     "two at order 6 for the almost-regular target, one at order 7\nfor the regular target)",
     "Order 6 is even.")

# ---------------- reflection ----------------
edit(REFL, "The first is the common one, and it has rigid, non-self-converse examples from order 6.",
     "The first is the common one, and it has rigid, non-self-converse examples from order 5 (strong ones from order 6).",
     "K1 => S4 and S4 => K1 (order 5) are rigid non-self-converse cut failures.")
edit(REFL, "module\nshattering and the source/sink trap.",
     "module\nshattering and Session74's source/sink pruning lemma.",
     "Session74's name for the statement.")
edit(REFL, "**The regular target never fails.** Every tournament has a Hamiltonian path from a maximum-score to a minimum-score\nvertex.",
     "**The regular target never fails.** Every tournament has a Hamiltonian path from a maximum-score to a minimum-score\nvertex (a classical fact: Thomassen 1982; Lu 2018 for a minimum vertex beating a maximum one).",
     "B7 is prior art.")
edit(REFL, "fails at sources, sinks and cut vertices (section 5)",
     "its local minima sit at sources, sinks and cut vertices (section 5)",
     "Failures themselves can be 2-connected (orders 7-9); the statement holds for local minima.")

# ---------------- HYP-9420 (title still carries the labelled 7) ----------------
edit(HYP, "exact data are d(P_7) = 7 and d(P_11) = 15",
     "exact class distances are d(P_7) = 6 and d(P_11) = 15, both equal to L (THM-4810 (D)'s 7 is the labelled distance to i -> j iff i < j, MISTAKE-682)",
     "MISTAKE-682's addendum says HYP-9420 is corrected, but the title still gives d(P_7) = 7.")

# dry run: apply sequentially per file on in-memory copies; each old must occur exactly once at its turn
texts = {}
ok = True
for e in E:
    f = e['file']
    if f not in texts:
        texts[f] = open(f, encoding='utf-8').read()
    c = texts[f].count(e['old'])
    if c != 1:
        ok = False
        print('NOT UNIQUE (count %d): %s :: %r' % (c, os.path.basename(f), e['old'][:80]))
    texts[f] = texts[f].replace(e['old'], e['new'], 1)
# also check against the pristine files (each old must occur once in the file as it stands now)
for e in E:
    c = open(e['file'], encoding='utf-8').read().count(e['old'])
    if c != 1:
        ok = False
        print('NOT UNIQUE in pristine file (count %d): %r' % (c, e['old'][:80]))
print('edits:', len(E), 'all unique:', ok)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'required_edits.json')
with open(out, 'w', encoding='utf-8') as fh:
    json.dump(E, fh, ensure_ascii=False, indent=1)
print('wrote', out)
