"""Independent literal map and compiler-state check of the one-vertex gap bridge."""
from pathlib import Path
import json
import module_lift as m
import intertwiner_probe as i
from audit_defect_two import construct

HERE=Path(__file__).resolve().parent


def main():
    cert=json.loads((HERE.parent/'gap_bridge.json').read_text())
    k=17;D=1<<k;bits=(1,)*k+(0,0,0,1)
    x=sum(b<<j for j,b in enumerate(bits))
    old,p,r=construct((2,4,1),(0,),x)
    # Independently apply the full six-coordinate word matrix to the C3 core.
    core=(3,1,3,12,3,7)
    parent=tuple(sum(a*b for a,b in zip(row,core)) for row in i.matrix(bits))
    H=parent[0]+parent[4]-1;q=parent[0]
    z=(H,q-1,q,H+1)
    assert z==(135*D+23,27*D+26,27*D+27,135*D+24)
    for _ in range(2):
        h,s,q,F=z
        z=(h+F-s,h-s,h,h+F)
        old=m.birth(old,r);r=len(old)-1
    assert tuple(old)==tuple(cert['source_rows'])
    assert z==tuple(cert['source_state'])==(405*D+71,135*D+24,243*D+21,513*D+68)
    new,np,nr=construct((0,),(0,),45157534)
    previous=json.loads((HERE/'defect_two_audit_receipt.json').read_text())['adjacent_ready_control']
    assert new==tuple(cert['target']['rows'])==tuple(previous['rows'])
    assert tuple(cert['target']['state'])==tuple(previous['state'])
    assert cert['target']['state'][0]==3*z[0]+2
    mapping=dict(cert['retained_map'])
    assert set(mapping)==set(range(28)) and len(set(mapping.values()))==28
    assert set(range(29))-set(mapping.values())=={26}
    assert (r,cert['source_predecessor'],mapping[r],mapping[26])==(27,26,nr,np)
    flips=[]
    for a in range(28):
        for b in range(a+1,28):
            if ((old[a]>>b)&1)!=((new[mapping[a]]>>mapping[b])&1):flips.append((a,b))
    assert flips==list(map(tuple,cert['changed_old_pairs'])) and len(flips)==75
    assert all(r not in pair for pair in flips)
    module=tuple(range(3,20))
    assert tuple(mapping[v] for v in module)==tuple(range(1,18))
    assert m.induced(old,module)==m.induced(new,tuple(mapping[v] for v in module))==m.tt(17)
    assert [v for v in range(28) if old[v]>>r&1]==[26]
    assert [v for v in range(29) if new[v]>>nr&1]==[np]
    changed_module_ports=[v for v in module if ((old[26]>>v)&1)!=((new[np]>>mapping[v])&1)]
    assert changed_module_ports==[4,5,6,7,10,13,14,19]
    result=dict(status='PASS',scope='Independent reconstruction of both actual graphs; six-state source compiler and two N births; all retained old pairs and exact target state; no large inverse or large subset census',
                source_order=28,target_order=29,source_state=z,target_state=previous['state'],
                retained_vertices=28,inserted_target_vertex=26,changed_old_pairs=75,
                retained_root=27,target_root=nr,retained_predecessor=26,target_predecessor=np,
                retained_transitive_module_order=17,changed_predecessor_module_arcs=changed_module_ports,
                all_old_root_star_pairs_preserved=True,all_module_internal_pairs_preserved=True)
    (HERE/'gap_bridge_audit_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
