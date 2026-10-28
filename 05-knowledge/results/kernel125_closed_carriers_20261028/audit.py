"""Recheck receipts, pinned dependencies, and small alternate executions."""
import contextlib
import gzip
import hashlib
import io
import json
import re
import subprocess
import sys
from pathlib import Path
import root_minimum as root

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CHECKS=0
def need(condition,message):
    global CHECKS
    CHECKS+=1
    if not condition: raise RuntimeError(message)
def read(name): return json.loads((HERE/name).read_text(encoding='utf-8'))
def semantic(row): return {k:v for k,v in row.items() if k!='seconds'}

def main():
    pins=read('source_pins.json')
    for path,sha in pins['sha256'].items():
        need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,'source pin '+path)
    for script,receipt in (('kernel125.py','kernel125.json'),
                           ('strong_ports.py','strong_ports.json'),
                           ('legacy_audit.py','legacy_repairs.json')):
        before=(HERE/receipt).read_bytes()
        result=subprocess.run([sys.executable,'-O','-B',str(HERE/script)],
                              cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
        need(result.returncode==0,script+' optimized: '+result.stderr)
        need((HERE/receipt).read_bytes()==before,script+' optimized receipt identity')
    kernel=read('kernel125.json'); ports=read('strong_ports.json')
    need(kernel['unmarked_hub_quotient']==25 and kernel['states']==125,'quotient')
    need(kernel['diagonal_guard_table']==[[i==j for j in range(5)] for i in range(5)],'native guards')
    need(ports['full_chart_group']==21 and ports['four_atom_relative_gain_states']==9261,'actual ports')
    degree=read('root_degree.json'); plain=read('root_minimum.json')
    reference=json.loads((root.PRIOR/'root_switch.json').read_text(encoding='utf-8'))
    for a,b,c in zip(degree['results'],plain['results'],reference['results']):
        for key in ('classes','pointed_states','automorphism_order_histogram'):
            need(a[key]==b[key]==c[key],'root receipt '+key)
        need(a['rejected_candidates']+a['classes']==a['pointed_states'],'candidate partition')
    need(degree['results'][-1]['reroots_examined']==278461,'degree certificate work')
    need(degree['exact_class_bijection_to_prior_bank'],'full exact reference coverage')
    raw=gzip.decompress((root.PRIOR/'tournament_bank.json.gz').read_bytes())
    need(hashlib.sha256(raw).hexdigest()==plain['bank_sha256'],'plain exact bank identity')
    for report in (degree,plain):
        h=report['hostile_81_27']
        need(h['automorphism_orders']==[81,27] and h['labelled_cyclic_triples']==67,'hostile groups')
        need(h['Hamilton_paths']==[462915,462915] and h['identity_order_legal']==[True,False],'hostile guards')
    # The saved order-nine run is not repeated. Reproduce its short prefix
    # with optimization enabled by running this audit itself under -O.
    parents=[()]
    with contextlib.redirect_stdout(io.StringIO()):
        for n in range(1,9):
            parents,row=root.generate(parents,n,65536,True)
            # JSON roundtrip normalizes integer histogram keys.
            row=json.loads(json.dumps(row))
            need(semantic(row)==semantic(degree['results'][n-1]),'small degree-first replay')
    channels=read('sequence_channels.json')
    orders=[row['m'] for row in channels['rows']]
    need(orders==list(range(31))+[50,100,200],'exact executed channel universe')
    for row in channels['rows']:
        w=row['work']
        need(w['separate_odd_factor_evaluations']==3*w['shared_odd_factor_evaluations'],'shared factors')
    # Only this packet's documentation is in the link universe.
    links=0
    for file in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',file.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'): continue
            target=target.split('#',1)[0]
            need((file.parent/target).exists(),'local link '+str(file.name)+':'+target)
            links+=1
    report={'status':'PASS','checks':CHECKS,'source_pins':len(pins['sha256']),
            'optimized_byte_identical_receipts':['kernel125.json','strong_ports.json','legacy_repairs.json'],
            'degree_generation_replayed_through':8,'order9_full_receipt_checked':True,
            'order9_reference_audit':'full exact bijection was executed in root_minimum.py --degree-first',
            'local_document_links_checked':links}
    (HERE/'audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))

if __name__=='__main__':main()
