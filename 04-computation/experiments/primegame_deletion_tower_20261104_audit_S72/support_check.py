import sys
D = "data"
for n in range(3, 10):
    sups = set(); decks = set(); wsup = set(); N = 0
    for line in open(f"{D}/audit_nodes_{n}.txt"):
        if line.startswith("#"): continue
        a, W, md = line.split("|")
        md = tuple(int(x) for x in md.split())
        Wd = tuple(sorted(int(t.split(":")[0]) for t in W.split()))
        decks.add(md); sups.add(tuple(sorted(set(md)))); wsup.add(Wd); N += 1
    print(f"n={n}: nodes {N}, distinct merged decks {len(decks)}, distinct deck support sets {len(sups)}, distinct W supports {len(wsup)}")
