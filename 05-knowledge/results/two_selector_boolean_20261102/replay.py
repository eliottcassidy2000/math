"""Pinned-source replay, optimization check, and literal transport audit."""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import subprocess
import sys
from constructor import construct
from experiment import need, rows

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def named_rows(tile,order):
    current=rows(7,tile,True);out=[0]*7
    for i in range(7):
        for j in range(7):
            if current[i]>>j&1:out[order[i]]|=1<<order[j]
    return out


def main():
    manifest=json.loads((HERE/'sources.json').read_text(encoding='utf8'))
    for entry in manifest['files']:
        need(sha256((ROOT/entry['path']).read_bytes()).hexdigest()==entry['sha256'],
             'source integrity: '+entry['path'])
    runs=[]
    for script,receipt in (('experiment.py','experiment.json'),('independent.py','independent.json'),
                           ('inspiration_audit.py','inspiration_audit.json')):
        hashes=[]
        for optimization in (False,True):
            command=[sys.executable,'-B']+(['-O'] if optimization else [])+['-X','utf8',str(HERE/script)]
            result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
            need(result.returncode==0,script+': '+result.stderr)
            hashes.append(sha256((HERE/receipt).read_bytes()).hexdigest())
        need(hashes[0]==hashes[1],'normal and optimized receipt agreement')
        runs.append({'script':script,'receipt':receipt,'sha256':hashes[0]})
    negative=construct(868,17,1);positive=construct(868,17,2)
    need(not negative['reachable'] and positive['minimum_commands']==4,'second-selector inverse control')
    for step in positive['trace']:
        before=named_rows(step['source'],step['source_vertex_order'])
        after=named_rows(step['target'],step['target_vertex_order'])
        if step['command']!='C':
            need(before==after,'re-marking retains every named source arc')
        else:
            order=step['source_vertex_order']
            kept={tuple(sorted(pair)) for pair in zip(order,order[1:])}
            for i,j in combinations(range(7),2):
                changed=(before[i]>>j&1)!=(after[i]>>j&1)
                need(changed==((i,j) not in kept),'actual complement keeps exactly its named path')
    report={'status':'PASS','pinned_files':len(manifest['files']),'runs':runs,
            'constructor_negative':negative,'constructor_positive':positive,
            'transport':'every named original-label arc checked at every step'}
    (HERE/'replay.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':'PASS','pinned_files':len(manifest['files']),
                      'byte_agreeing_receipts':len(runs),'inverse_word':positive['word']}))


if __name__=='__main__':main()
