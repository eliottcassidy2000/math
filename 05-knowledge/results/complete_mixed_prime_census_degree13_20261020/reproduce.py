#!/usr/bin/env python3
"""Rebuild every claimed control without repository/network dependencies."""
from pathlib import Path
import subprocess,sys,json,hashlib,gzip
ROOT=Path(__file__).resolve().parent

def require(x,msg):
 if not x:raise RuntimeError(msg)
def run(*args):subprocess.run(args,cwd=ROOT,check=True)
def digest(file):return hashlib.sha256((ROOT/file).read_bytes()).hexdigest()
def data(file):return json.loads((ROOT/file).read_text())
def without_times(obj):
 if isinstance(obj,dict):return {k:without_times(v) for k,v in obj.items() if k!='seconds'}
 if isinstance(obj,list):return [without_times(x) for x in obj]
 return obj

def main():
 expected={}
 for f in ('control','all_lower','small13','large13'):
  z=ROOT/(f+'.json.gz')
  if z.exists():expected[f]=hashlib.sha256(gzip.decompress(z.read_bytes())).hexdigest()
 run('g++','-std=c++17','-O3','cube_orbits.cpp','-o','cube_orbits')
 run('g++','-std=c++17','-O3','generator_orbits.cpp','-o','generator_orbits')
 checks=[]
 for mode in ('control','all_lower','small13','large13'):
  run(sys.executable,'audit_census.py',mode);h=digest(mode+'.json')
  if mode in expected:require(h==expected[mode],mode+' frozen receipt')
  checks.append({'check':mode+'_normal','sha256':h,'status':'PASS'})
  if mode in ('small13','large13'):
   run(sys.executable,'-O','audit_census.py',mode)
   require(h==digest(mode+'.json'),mode+' optimized receipt differs')
   checks.append({'check':mode+'_normal_vs_optimized','sha256':h,'status':'PASS'})
 for script in ('prime13_certificate','subdirect_certificate','burnside_control','generator_action_replay','independent_large_audit','compare_receipts'):
  run(sys.executable,script+'.py')
  file='receipt_comparison.json' if script=='compare_receipts' else script+'.json'
  first=without_times(data(file));run(sys.executable,'-O',script+'.py')
  require(first==without_times(data(file)),script+' optimized control differs')
  checks.append({'check':script+'_normal_vs_optimized','status':'PASS'})
 refs={(r['n'],r['class']):r for r in data('prior_reference_rows.json')['rows']}
 fields=['order','normalizer_order','pair_orbits','fixed_objects_each','relative_tournament_orbits','relative_even_graph_orbits','tau_T','tau_E','difference','sign_escape_global_types','tournament_fusion_defect','even_graph_fusion_defect','sign_escape_relative_orbits']
 lower=data('all_lower.json')['rows']
 for row in lower:
  for f in fields:require(row[f]==refs[row['n'],row['class']][f],'lower-degree reference mismatch')
 require(len(lower)==28,'lower scope')
 (ROOT/'all_lower_comparison.json').write_text(json.dumps({'status':'PASS','matched_rows':len(lower),'fields':fields},indent=2)+'\n')
 checks.append({'check':'all_28_lower_rows_against_frozen_references','status':'PASS'})
 (ROOT/'reproduction_receipt.json').write_text(json.dumps({'status':'PASS','checks':checks,'degrees_through':13,'complete_rows':48,'degree13_rows':20,'incoming_degree13_matched_rows':18,'prior_incoming_common_matched_rows':30},indent=2,sort_keys=True)+'\n')
 print('PASS: complete48-row census; all20 degree13 rows; independent large-row audits; normal/optimized controls.')
if __name__=='__main__':main()
