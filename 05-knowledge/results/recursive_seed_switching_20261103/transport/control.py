"""Fresh literal subset-DP replay of the inherited named 81/27 hostile."""
from pathlib import Path
import hashlib, json, sys

HERE=Path(__file__).resolve().parent; RESULTS=HERE.parent.parent
OLD=RESULTS/'zero_ready_component_returns_20261103'
sys.path.insert(0,str(OLD/'inverse'))
import reader
sys.path.insert(0,str(OLD/'interface'))
from check import phase


def main():
    source=RESULTS/'return6325_completion_20261103/halving_response.json'
    hostile=json.loads(source.read_text())['hostile'];records=[]
    reader.need(reader.tables((2,4,1))[0][-1]==3,'Triangle positive control')
    for rows,response in zip(hostile['rows'],hostile['responses']):
        rows=tuple(rows);fields=reader.all_fields(rows)
        reader.need(fields[8][0]==hostile['baseline'],'Same actual source Hamilton count')
        charts={}
        for port in (8,6):
            h,e,S2,F,p,q,P2,S3=fields[port]
            reader.need((S2-e)%2==0,'Integral golden coordinates')
            golden=[h,(S2-e)//2,(S2+e)//2]
            charts[str(port)]=dict(eight=list(fields[port]),golden=golden,charts=phase(golden))
        new=reader.birth(rows,((1<<len(rows))-1)^hostile['cut'])
        counted=reader.tables(new)[0][-1]
        reader.need(counted==response,'Actual same-cut extension differs')
        records.append(dict(rows=list(rows),names=list(range(len(rows))),cut=hostile['cut'],
                            H=fields[8][0],new_H=counted,named_port_charts=charts))
    reader.need(records[0]['new_H']!=records[1]['new_H'],'Hostile really separates')
    pins=[source,Path(__file__),OLD/'inverse/reader.py',OLD/'interface/check.py']
    print(json.dumps(dict(status='PASS',scope='Actual 81/27 pair and distinguished cut; automorphism labels inherited, not recomputed.',
        records=records,source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}),indent=2,sort_keys=True))


if __name__=='__main__':main()
