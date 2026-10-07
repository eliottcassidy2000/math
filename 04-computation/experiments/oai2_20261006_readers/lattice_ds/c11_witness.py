# Dump exact witnesses: N-point subsets of Z[w] (pairs (a,b) = a + b*w, w = e^{i pi/3}) with u(N) pairs at norm D.
import random, math, json, importlib.util
spec = importlib.util.spec_from_file_location("s", "/private/tmp/claude-501/-Users-e-Documents-GitHub-math/f1b41b3c-5f00-4184-8c88-059bb9352d19/scratchpad/oai2/lds/calc/c9_sections_search.py")
src = open(spec.origin).read().split("Ds = [1, 3, 7, 13, 49, 91]")[0]
ns = {}; exec(src, ns)
anneal, shell, count_edges = ns['anneal'], ns['shell'], ns['count_edges']
uN = ns['uN']
W = {}
sh7 = shell(7)
for N in (9, 10, 11, 12, 21):
    best = None
    for seed in range(20):
        e, S = anneal(N, sh7, iters=4000, seed=seed)
        if best is None or e > best[0]: best = (e, sorted(S))
        if e == uN[N]: break
    e, S = best
    assert count_edges(S, sh7) == e
    W[N] = {"D": 7, "edges": e, "u(N)": uN[N], "points": S}
    print(f"N={N}: D=7 witness with {e} edges (u(N)={uN[N]}): {S}")
json.dump(W, open("/private/tmp/claude-501/-Users-e-Documents-GitHub-math/f1b41b3c-5f00-4184-8c88-059bb9352d19/scratchpad/oai2/lds/calc/witnesses_D7.json", "w"))
