import importlib.util
spec = importlib.util.spec_from_file_location('m','/Users/e/Documents/GitHub/math/04-computation/hampath_cycle_correspondence_mac-mini-2026-06-15-S6.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

# Characterize exactly when exists-order-T_G-strong is FALSE.
for n in range(2,8):
    reps=m.all_noniso_graphs(n)
    full=n*(n-1)//2
    false_list=[]
    for es in reps:
        adj=m.edges_to_adj(n,es)
        if not m.exists_order_tournament_strong(n,adj):
            ne=len(es)
            label = 'EMPTY' if ne==0 else ('COMPLETE' if ne==full else f'OTHER({ne}e)')
            false_list.append(label)
    print(f"n={n}: exStrong=FALSE graphs -> {false_list}")

# Now the cleanest cycle correspondence. T_G edges are ALWAYS acyclic (oriented by order).
# Strongness must use backward arcs from NON-edges. 
# Test: G has Ham cycle  <=>  exists order with T_G containing a directed Ham cycle?
# T_G always contains a Ham cycle (Camion: strong tournament => Ham cycle; but T_G strong
#   requires non-complete). Actually EVERY tournament has a Ham PATH (Redei). 
# The right object: does T_G contain a directed Ham cycle using the SPINE pattern?
# Let's instead directly test the candidate the task names: 
#   "all consecutive arcs forward AND wrap arc" under exists-order  == spine+wrap (done, fails).
# And report the spine+wrap false-neg structure: they are circulant/self-complementary-ish.
print()
print("spine+wrap false-negatives (HC but no order with spine-forward+wrap):")
for n in range(3,8):
    reps=m.all_noniso_graphs(n); full=n*(n-1)//2
    fn=[]
    for es in reps:
        adj=m.edges_to_adj(n,es)
        if m.has_hamiltonian_cycle(n,adj) and not m.exists_order_spine_and_wrap(n,adj):
            ne=len(es)
            fn.append('COMPLETE' if ne==full else f'{ne}e')
    print(f"  n={n}: {fn}")
