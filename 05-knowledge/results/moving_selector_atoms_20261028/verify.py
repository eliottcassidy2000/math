"""Replay all finite controls and verify pinned inputs and artifacts."""
from pathlib import Path
from hashlib import sha256
import json,contextlib,io
import probes,frames,folded_minimum,tree_actions
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
need=probes.need


def digest(path):
    data=path.read_bytes()
    if path.suffix!='.gz':data=data.replace(b'\r\n',b'\n')
    return sha256(data).hexdigest()


def seals():
    for name,base in [('sources.json',ROOT),('manifest.json',HERE)]:
        for path,wanted in json.loads((HERE/name).read_text(encoding='utf-8')).items():
            need(digest(base/path)==wanted,('changed pinned input',path))


def main():
    seals();outputs={}
    with contextlib.redirect_stdout(io.StringIO()):
        for name,job in [('probe-receipt.json',probes.main),('frame-receipt.json',frames.main),
                         ('folded-receipt.json',folded_minimum.main),('tree-receipt.json',tree_actions.main)]:
            saved=json.loads((HERE/name).read_text(encoding='utf-8'));fresh=json.loads(json.dumps(job()))
            need(saved==fresh,('fresh replay changed',name));outputs[name]=fresh
    pr=outputs['probe-receipt.json'];fr=outputs['frame-receipt.json'];dp=outputs['folded-receipt.json']
    need(pr['padded'][1]['first_separating_position']==2,'nine-vertex cursor separator')
    need(pr['run_hostile']['all_strong_children']['H']==[2899675,2906091],'actual strong composition targets')
    need(pr['wreath_scan']['presentations']==pr['wreath_scan']['cursor_profile_cells']==39,'full one-source presentation census')
    need([dp[k]['cells'] for k in ('minimal','selected_vertex','all_last_vertices')]==[16,35,64],'exact guarded quotient costs')
    need(dp['minimal']['count']==33 and dp['baseline']['orbits']==24,'folded mass and inherited orbit minimum')
    need([v['group'] for v in fr['frames']['atoms'][-2:]]==[81,27],'hostile full groups')
    need(fr['paley27']['alpha_orbit']==[3,5,4] and fr['paley27']['fixed_ports']==[0,1,2],'literal field phase')
    seals()
    report=dict(status='PASS',optimized=not __debug__,cursor_separator_at9=2,
                strong_partner_H=[2899675,2906091],forest_controls=sum(v['forests'] for v in fr['forests']),
                folded_states=[64,24,16],guarded_folded_states=[35,64],small_lifts=dp['small_lifts'],
                tree_controls=sum(v['trees'] for v in outputs['tree-receipt.json']['orders']),
                hostile_groups=[81,27],Paley27_Frobenius_cycle=fr['paley27']['alpha_orbit'])
    print(json.dumps(report,indent=2));return report


if __name__=='__main__':main()
