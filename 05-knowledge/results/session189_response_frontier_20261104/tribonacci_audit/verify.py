#!/usr/bin/env python3
"""Independent finite root-response coverage with one strong witness per odd value."""
from pathlib import Path
import json,hashlib,subprocess,sys
P=Path(__file__).resolve().parent

def need(x,s):
 if not x:raise RuntimeError(s)
def main():
 incoming=json.load(open(P/'INCOMING.json'))
 for f,info in incoming.items():
  b=(P/f).read_bytes();need(hashlib.sha256(b).hexdigest()==info['sha256'],'incoming SHA');need(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==info['git_blob'],'incoming exact Git blob')
 codes=json.load(open(P/'CODES.json'));data=b''.join((P/f).read_bytes()for f in codes['files']);(P/'all_codes.txt').write_bytes(data)
 for source,exe in [('coverage.cpp','coverage'),('witness_check.cpp','witness_check')]:subprocess.run(['g++','-O3','-std=c++17',str(P/source),'-o',str(P/exe)],check=True)
 out=subprocess.check_output([str(P/'coverage'),str(P/'all_codes.txt'),str(P/'witnesses.hex')]);(P/'coverage.json').write_bytes(out);r=json.loads(out);need((r['rooted_total'],r['distinct_total'])==(1782468,1728779),'full rooted totals');need(r['height1_only_additions']==[1,11,13,47,67],'exact omitted height-one values');F={x['family']:x for x in r['families']};expected={'two':(15449,8637),'two_strong':(15449,8594),'nearjoin':(382687,674533),'tower3to10':(77,801304),'tower1to10':(432485,801327),'source_union':(702189,832898),'source_union_strong':(702189,832195)}
 for name,(endpoint,count)in expected.items():need((F[name]['maximum_complete_odd_interval'],F[name]['distinct_values'])==(endpoint,count),'exact finite coverage values')
 out=subprocess.check_output([str(P/'witness_check'),str(P/'witnesses.hex')]);(P/'witness_receipt.json').write_bytes(out);w=json.loads(out);need(w['verified_strong_constructions']==351061 and w['maximum_output_order']==19,'all strong physical construction witnesses')
 lines=(P/'witnesses.hex').read_text().splitlines();parts=[]
 for j in range(0,len(lines),45000):
  f=f'witnesses.part{j//45000}.hex';(P/f).write_text('\n'.join(lines[j:j+45000])+'\n');parts.append(f)
 (P/'WITNESSES.json').write_text(json.dumps({'parts':parts,'first_value':69,'value_step':2,'witness_count':len(lines),'encoding':'Unsigned hexadecimal. Bits0..35 column-major code;36..39 physical root;40..43 base order;44 converse;45..46 operation(0 C2,1 nearjoin,2 tower);47..50 parameter k;51 base strong flag. Line i certifies H=69+2i.','all_strong':True,'maximum_output_order':19},sort_keys=True,indent=2)+'\n')
 opt=['-O']if sys.flags.optimize else[]
 for script,name in [('controls.py','controls'),('samples.py','samples')]:
  out=subprocess.check_output([sys.executable,*opt,str(P/script)]);(P/(name+'.replay.stdout')).write_bytes(out)
 sources=['verify.py','coverage.cpp','witness_check.cpp','controls.py','ham.cpp','samples.py'];artifacts=['coverage.json','witness_receipt.json','WITNESSES.json','controls_receipt.json','sample_receipt.json','actual_graphs.in','actual_values.txt','membership_controls.json','tower_controls.json','sample_graphs.in','sample_counts.txt','sample_controls.json']+codes['files']+parts
 receipt={'status':'PASS','rooted_bases':1782468,'distinct_root_response_rows':1728779,'strong_witnesses':351061,'maximum_output_order':19,'scope':'Finite construction coverage only. No all-order selector, joint-port grammar closure, or ready-generator completeness claimed.','source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in sources},'artifact_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in artifacts},'incoming_pins_sha256':hashlib.sha256((P/'INCOMING.json').read_bytes()).hexdigest()};(P/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
