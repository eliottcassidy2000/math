# Dump the 37 genuine dense seven-vertex survivors of the Paley TT7 census, with their signed multiplicities, from the
# audit's symbolic elimination engine (audit/audit_census_symbolic.py), in its canonical labelling with the order
# orientation (arc i -> j for i < j).  Output: genuine37.json.  mac-mini-2026-10-25-pyth.
import json, sys, time
sys.path.insert(0, '../../audit')
import audit_census_symbolic as acs
t0 = time.time()
res = acs.census(7)
out = []
for key, v in res.items():
    if key == 'one' or key[1] != 7: continue
    mask = key[2]; deg = acs.degrees(7, mask)
    E = [list(acs.EDGES[7][n]) for n in acs.bits(mask)]
    out.append({'mask': mask, 'edges': E, 'm': v[0] if len(v) == 1 else v, 'degrees': list(deg),
                'all_even': all(d % 2 == 0 for d in deg)})
json.dump(out, open('genuine37.json', 'w'), indent=0)
print(len(out), 'seven-vertex survivors;', sum(1 for o in out if not o['all_even']), 'genuine;', f'{time.time() - t0:.0f}s')
