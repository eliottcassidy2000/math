import importlib.util
spec = importlib.util.spec_from_file_location('m','/Users/e/Documents/GitHub/math/04-computation/hampath_cycle_correspondence_mac-mini-2026-06-15-S6.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

# For each n, find the HC=1 graphs where spine+wrap NOT findable, and HC=0 where it IS,
# and the HC vs exists-strong mismatches. Probe a cleaner cycle correspondence.
for n in range(3,8):
    reps = m.all_noniso_graphs(n)
    sw_false_neg=[]   # HC but no spine+wrap order
    sw_false_pos=[]   # spine+wrap order but no HC
    es_false_neg=[]   # HC but not exists-strong
    es_false_pos=[]   # exists-strong but not HC
    for es in reps:
        adj=m.edges_to_adj(n,es)
        hc=m.has_hamiltonian_cycle(n,adj)
        sw=m.exists_order_spine_and_wrap(n,adj)
        exS=m.exists_order_tournament_strong(n,adj)
        if hc and not sw: sw_false_neg.append(sorted(es))
        if sw and not hc: sw_false_pos.append(sorted(es))
        if hc and not exS: es_false_neg.append(sorted(es))
        if exS and not hc: es_false_pos.append(sorted(es))
    print(f"n={n}: spine+wrap  FN(HC,no-sw)={len(sw_false_neg)} FP(sw,no-HC)={len(sw_false_pos)}")
    print(f"      exists-strong FN(HC,no-exS)={len(es_false_neg)} FP(exS,no-HC)={len(es_false_pos)}")
    # show the FN graphs (these break HC=>spine+wrap)
    if sw_false_neg:
        print(f"      HC-but-no-spine+wrap examples (<=3): {sw_false_neg[:3]}")
    if es_false_neg:
        print(f"      HC-but-not-exists-strong examples (<=3): {es_false_neg[:3]}")
