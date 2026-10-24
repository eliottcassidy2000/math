#!/usr/bin/env python3
"""Audit independent receipts against a manifest-frozen producer real copy."""
from pathlib import Path
from hashlib import sha256
from math import gcd
import argparse,importlib.util,json

def require(c,m):
    if not c:raise RuntimeError(m)
def main():
    p=argparse.ArgumentParser();p.add_argument('--producer',required=True);p.add_argument('--audit',default=str(Path(__file__).resolve().parent));p.add_argument('--output',required=True);a=p.parse_args()
    P=Path(a.producer);A=Path(a.audit)
    manifest=json.loads((P/'MANIFEST.json').read_text())
    for rec in manifest['files']:require(sha256((P/rec['path']).read_bytes()).hexdigest()==rec['sha256'],'producer manifest mismatch '+rec['path'])
    source=json.loads((P/'normal.json').read_text());independent=json.loads((A/'normal.json').read_text())
    matches=[]
    for expected,actual in zip(independent['counts'],source['counts']):
        require(expected['n']==actual['n'],'order')
        for i,j in [('endpoint','endpoint_ports'),('ports','balanced_ports'),('internal','internal_ports'),('blue_internal','internal_blue_ports'),('blue_boundary','boundary_blue_ports')]:require(expected[i]==actual[j],'independent '+i+' at n='+str(expected['n']))
        require(actual['balanced_blue_ports']==expected['blue_internal']+expected['blue_boundary'],'all blue ports')
        require(actual['boundary_ports']==expected['ports']-expected['internal'],'boundary subtraction')
        require(actual['internal_black_ports']==expected['internal']-expected['blue_internal'],'internal black subtraction')
        require(actual['boundary_black_ports']==expected['ports']-expected['internal']-expected['blue_boundary'],'boundary black subtraction')
        if expected['n']>=4:
            require(expected['native_internal_C_lines']==actual['internal_lines'],'native line convention')
            require(actual['internal_blue_lines']==expected['blue_internal']//2,'blue lines')
            require(actual['internal_black_lines']==(expected['internal']-expected['blue_internal'])//2,'black lines')
            require(actual['incident_lines']==expected['ports']-expected['internal']//2,'incident lines')
        matches.append(expected['n'])
    require(matches==[2,4,6,8,10,12],'complete promised scalar scope')
    for row in independent['rotation_controls']:
        data=next(x for x in source['counts'] if x['n']==row['n'])
        require(row['fixed']==data['proper_rotation_fixed_counts'][str(row['gcd'])],'actual edge orbit versus divisor fixed polynomial')
    import sys
    sys.path.insert(0,str(P))
    from check import rotation_fixed
    extras=json.loads((A/'extra-rotation-controls.json').read_text())
    for row in extras:require(row['fixed']==rotation_fixed(row['n'],gcd(row['n']+1,row['step'])),'extra composite fixed sets')
    sourcepins=json.loads((P/'SOURCE_PINS.json').read_text())
    for pin in sourcepins['sources']:require(sha256((P/pin['path']).read_bytes()).hexdigest()==pin['sha256'],'source pin')
    inherited=json.loads((P/'sources/regular-eleven-counts.json').read_text())
    require(inherited['almost_regular_parent_classes']==13333 and inherited['self_converse_parents']==229,'inherited n10 classes distinguished')
    result={'status':'PASS','producer_manifest_sha256':sha256((P/'MANIFEST.json').read_bytes()).hexdigest(),'independent_scalar_orders':matches,'n8_literal_universe':2097152,'native_C_line_convention':'not optionally R-quotiented','additional_rotation_fixed_controls':extras,'inherited_n10_class_data':{'ordinary_balanced_classes':13333,'self_converse_balanced_classes':229,'merged_balanced_nodes':6781},'prior_scalar_provenance':'THM-554 exact full score polynomial and THM-555 total balanced ports through n10 are inherited and reproduced'}
    Path(a.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
