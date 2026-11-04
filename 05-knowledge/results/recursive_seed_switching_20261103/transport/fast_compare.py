"""Compare the compiled counter with the full Python prefix-port receipt."""
from pathlib import Path
import hashlib, json, sys
import fast_gap

HERE=Path(__file__).resolve().parent; RESULTS=HERE.parent.parent
OLD=RESULTS/'zero_ready_component_returns_20261103/interface'
sys.path.insert(0,str(OLD))
import gap
import phase


def main():
    inp=HERE.parent/'selector/named_repair.json'; verified=HERE/'verify.json'
    target=json.loads(inp.read_text())['target']; proof=json.loads(verified.read_text())
    rows=tuple(target['rows']);names=target['names'];fb=tuple(names.index(v) for v in (4,7,14,39))
    gap.count=fast_gap.count; phase.count=fast_gap.count
    r,s=names.index(50),names.index(51)
    fields=gap.pair_fields(rows,fb,r,s);observer=phase.observer17(rows,fb,(r,s))
    expected=next(x for x in proof['joint_pairs'] if x['ports']==[50,51])
    if fields != [expected['H'],expected['P2'],*expected['port_data'][0],*expected['port_data'][1],*expected['ordered_W2']]:
        raise RuntimeError('Compiled full named pair differs from Python receipt')
    if observer!=expected['observer17']:raise RuntimeError('Compiled17 differs from Python receipt')
    pins=[inp,verified,HERE/'fast_gap.cpp',HERE/'fast_gap.py',Path(__file__),OLD/'gap.py',OLD/'phase.py']
    print(json.dumps(dict(status='PASS',scope='Independent implementation of the same gap algorithm, compared with the completed original Python computation.',
        H=target['H'],ports=[50,51],fields=fields,observer17=observer,
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}),indent=2,sort_keys=True))


if __name__=='__main__':main()
