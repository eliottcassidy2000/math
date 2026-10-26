#!/usr/bin/env python3
"""Compare independent output, instrumented production output, and full n8 ledgers."""
import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

ap=argparse.ArgumentParser()
ap.add_argument('production_output',type=Path)
ap.add_argument('--audit-dir',type=Path,default=Path(__file__).resolve().parent)
a=ap.parse_args();p=a.production_output;d=a.audit_dir
ind=json.loads((d/'independent_results.json').read_text())['rows']
prod=json.loads((p/'census.json').read_text())
instrumented=json.loads((d/'production_direct_audit.json').read_text())
require(prod==instrumented,'instrumented production summary mismatch')
require([r['n'] for r in ind]==list(range(1,8)),'independent order coverage')
require([r['n'] for r in prod]==list(range(1,9)),'production order coverage')
for row,target in zip(ind,prod):
    for key in ('n','tilings','ordinary_classes','merged_classes'):
        require(row[key]==target[key],('summary',key))
    for key,pk in [('self_converse_classes','self_converse'),('right_legal_incidences_including_collapses','legal_right_relocations'),('ordinary_fibre_sizes','tc_histogram')]:
        require(row[key]==target[pk],('summary',key))
    if row['n']>=3:
        for key,pk in [('blue_lines','blue_lines'),('merged_loops_including_multiplicity','loop_lines')]:
            require(row[key]==target[pk],('summary',key))
        require(row['complement_lines']==target['blue_lines']+target['black_lines'],'line count')
with (p/'nodes_8.tsv').open() as f:
    nodes=list(csv.DictReader(f,delimiter='\t'))
merged_mass=Counter()
for r in nodes:
    merged_mass[int(r['merged_representative'])]+=int(r['tilings_tc'])
degree=Counter();colour=Counter();loops=0;rows=0
with gzip.open(p/'lines_8.tsv.gz','rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        u,v,b,k=[int(r[j]) for j in ('u','v','blue','multiplicity')]
        require(k>0 and b in (0,1) and u<=v,'invalid line record')
        require(u in merged_mass and v in merged_mass,'unknown endpoint')
        degree[u]+=k;degree[v]+=k;colour[b]+=k;loops+=(u==v)*k;rows+=1
require(degree==merged_mass,'degree/mass mismatch')
require(sum(merged_mass.values())==1<<21,'total mass')
require(len(nodes)==6880 and len(merged_mass)==3528 and sum(int(r['self_converse']) for r in nodes)==176,'n8 node totals')
require((colour[1],colour[0],loops,rows)==(2048,1046528,622,438374),'n8 line totals')
receipt={'status':'PASS','direct_arc_LUT_and_legality_audit_through_n':8,'direct_audit_summary_matches_production':True,'independent_complete_node_and_line_table_match_through_n':7,'n8_node_line_accounting':{'ordinary_nodes':len(nodes),'merged_nodes':len(merged_mass),'half_edges':sum(merged_mass.values()),'blue_lines':colour[1],'black_lines':colour[0],'loops':loops,'coloured_endpoint_rows':rows,'each_node_degree_equals_cube_fibre_mass':True}}
(d/'production_comparison.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
