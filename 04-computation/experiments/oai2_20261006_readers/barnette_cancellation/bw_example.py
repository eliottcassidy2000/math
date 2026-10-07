import random
import board_paths as BP
exec(open('board_weave_paths.py').read().split("rng = random.Random")[0])
rng = random.Random(7)
found = []
for trial in range(400):
    P = rand_path(9, rng, within=inner)
    if P is None: continue
    EP = edges_of(P)
    if not BP.closure_ok(nV, E, adj, eid, EP): continue
    r = BP.find_hc(E, nV, EP)
    found.append(([divmod(x, N) for x in P], r))
    if len(found) >= 4: break
for P, r in found:
    print("locally consistent 9-move path inside rings 1-2:", P, "in a closed tour:", r)
# sanity: an 8-move path inside rings 1-2 passing closure IS extendable?
cnt = [0, 0]
for trial in range(300):
    P = rand_path(8, rng, within=inner)
    if P is None: continue
    EP = edges_of(P)
    if not BP.closure_ok(nV, E, adj, eid, EP): continue
    r = BP.find_hc(E, nV, EP); cnt[0 if r else 1] += 1
print("8-move inner paths passing closure: extendable / not:", cnt)
