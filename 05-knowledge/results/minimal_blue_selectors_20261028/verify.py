"""Replay the complete bounded universe and reject changed pinned inputs."""
from hashlib import sha256
from pathlib import Path
import contextlib,io,json
import blue_minimum as m
import selector_probes,bridges,incoming_alignment
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(path):
    data=path.read_bytes()
    if path.suffix!='.gz':data=data.replace(b'\r\n',b'\n')
    return sha256(data).hexdigest()


def seals():
    for name,base in [('sources.json',ROOT),('manifest.json',HERE)]:
        for path,wanted in json.loads((HERE/name).read_text(encoding='utf-8')).items():
            m.need(digest(base/path)==wanted,('changed pinned input',path))


def main():
    seals()
    jobs=[('blue-receipt.json',m.main),('selector-receipt.json',selector_probes.main),('bridge-receipt.json',bridges.main),
          ('incoming-receipt.json',incoming_alignment.main)]
    with contextlib.redirect_stdout(io.StringIO()):
        results={}
        for name,job in jobs:
            saved=json.loads((HERE/name).read_text(encoding='utf-8'))
            fresh=json.loads(json.dumps(job()))
            m.need(fresh==saved,('fresh replay differs',name));results[name]=fresh
        bank=m.r.rooted_bank(7)
        m.need(bank==json.loads((HERE.parent/'singleton_pair_reconstruction_20261028/rooted-seven.json').read_text(encoding='utf-8')),'independent full rooted-bank rebuild')
    selector=results['selector-receipt.json']
    m.need(selector['carrier']['minimal_cells']==32272 and selector['carrier']['successor_port_cells']==32768,'minimal and selected carrier census')
    reset=selector['selectors']['successor_start_reset'][1]
    m.need(reset['legal_path_counts']==[23,21] and
           [row.get('29',0) for row in reset['partner_hamilton_histograms']]==[1,0],'normalized H29 event')
    seals()
    result=dict(status='PASS',optimized=not __debug__,all_presentations=sum(v['presentations'] for v in results['blue-receipt.json']['universes']),
                minimal_blue_counts=[v['blue_states'] for v in results['blue-receipt.json']['universes']],
                exact_minimal_cells_at7=32272,successor_port_cells_at7=32768,
                successor_counterexample_order=9,hostile_groups=[81,27],borwein_denominators=[15,113,80499],
                incoming_wreath_checks=results['incoming-receipt.json']['inherited']['active_checks'])
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':main()
