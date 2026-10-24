#!/usr/bin/env python3
"""Match independently obtained score and mass results to the frozen producer."""
import collections
import json
import re
from pathlib import Path
import independent as a
HERE = Path(__file__).resolve().parent

def main():
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    a.verify_manifest(HERE/'producer',pins['producer_manifest_sha256'])
    a.verify_manifest(HERE/'nonrigid',pins['nonrigid_manifest_sha256'])
    local = json.loads((HERE/'normal.json').read_text())
    producer = json.loads((HERE/'producer'/'normal.json').read_text())
    for lhs,rhs in [('producer/sources/nonrigid-receipt.json','nonrigid/normal.json'),
                    ('producer/sources/nonrigid-theorem.md','nonrigid/THEOREM.md'),
                    ('producer/sources/nonregular-audit.md','nonrigid/inherited/nonregular-audit.md')]:
        a.require((HERE/lhs).read_bytes() == (HERE/rhs).read_bytes(),'source byte identity')
    pairs = [('labelled_regular11','labelled_regular11'),
             ('nonrigid_labelled_mass','nonrigid_regular11_labelled_mass'),
             ('rigid_regular11_classes','rigid_regular11_classes'),
             ('regular11_classes','regular11_classes'),
             ('almost_regular10_classes','almost_regular10_classes'),
             ('regular11_error','regular11_error')]
    for ours,theirs in pairs:
        a.require(local['closure'][ours] == producer[theirs],'independent closure comparison: '+ours)
    for item in local['regular_labelled_counts']:
        a.require(item['labelled_count'] == producer['labelled_regular_counts'][str(item['order'])],'regular labelled crosscheck')
    a.require(local['fixed_assigned_parent_count'] == producer['fixed_profile_almost_regular10_labelled_count'],'fixed assigned parent')
    a.require(local['score_assignment_multiplier'] == producer['all_profiles_factor'],'score profile multiplier')
    hist = collections.Counter()
    for entry in local['nonrigid']['group_and_orbit_spectrum']:
        hist[11*(entry['vertex_orbits']-1)] += entry['classes']
    hist[110] += local['closure']['rigid_regular11_classes']
    a.require([list(x) for x in sorted(hist.items())] == producer['regular11_error_histogram'],'complete regular defect histogram')
    # The proof is a pinned inherited input; parse its two explicit formulas.
    proof = (HERE/'nonrigid'/'inherited'/'nonregular-audit.md').read_text()
    strong_offset = int(re.search(r'E_s,11 = E_regular11 \+ (\d+)',proof)[1])
    total_offset = int(re.search(r'E_11   = E_regular11 \+ (\d+)',proof)[1])
    strong = local['closure']['regular11_error'] + strong_offset
    total = local['closure']['regular11_error'] + total_offset
    a.require(strong == producer['strong11_error'] and total == producer['total11_error'],'inherited nonregular transfer')
    a.require(producer['unrestricted_order11_count_formula'] == f'(U11-{total})/11','unrestricted residual formula only')
    print(json.dumps(dict(status='PASS',producer_manifest_sha256=pins['producer_manifest_sha256'],
                         independently_matching_closure_fields=len(pairs),labelled_orders_checked=6,
                         regular_error_histogram=[list(x) for x in sorted(hist.items())],
                         inherited_strong_error=strong,inherited_total_error=total,
                         unrestricted_order11_count_formula=f'(U11-{total})/11',
                         no_full_parent_or_target_atlas_used=True),sort_keys=True,indent=2))
if __name__ == '__main__':
    main()
