"""Separate1537271 observer audit. Does not modify all-clock artifacts."""
import json,hashlib
from pathlib import Path
from candidate_phases import chart,need
HERE=Path(__file__).parent

def main():
 joint_path=HERE.parent/'ports/continuation_joint.json'
 candidate_path=HERE.parent/'returns/continuation/candidate.json'
 j=json.loads(joint_path.read_text());c=json.loads(candidate_path.read_text())
 need(j['H']==c['target_H']==1537271)
 need(j['fields_order']==['E1','S1','deletion_H','E2','S2'])
 # Revalidate precisely the four pinned historical chart outputs.
 old=json.loads((HERE.parent/'ports/input_audit.json').read_text())
 for s in old['states']:need(json.loads(json.dumps(chart(s['golden'])))==s['charts'])
 states=[]
 for name,fields in j['fields'].items():
  e1,s1,deletion_h,e2,s2=fields;need((s2-e1)%2==0)
  g=[j['H'],(s2-e1)//2,(s2+e1)//2]
  states.append({'physical_name':int(name),'endpoint_fields':fields,'golden':g,'charts':chart(g)})
 out={'status':'EXACT_CURRENT_OBSERVER_CODECS_FROM_INDEPENDENT_ENDPOINT_DP','H':j['H'],'source_H':c['source_H'],'source_ref':'4ad839b20dd280cdcfbadf86cb62b2be3c9615ad','codec_source':'05-knowledge/results/return6325_completion_20261103/transport.py','field_modulus':[1,0,1,1],'historical_codec_controls':len(old['states']),'input_sha256':{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (joint_path,candidate_path)},'states':states,'assumptions':['Endpoint fields are the independently recomputed physical-name fields supplied by ports/continuation_joint.json; the chart audit does not independently enumerate tournament paths.','The candidate1537271 is a separate supplied relation. No equality of phases is inferred from the preceding512423 relation.','Current F343 vector codec is not the ancestral spectral clock; these observations do not establish tournament automorphisms or authorize re-marking.','Native guards and deletions are audited separately; this file establishes observer values at the stated endpoint only.']}
 (HERE/'continuation_phases.json').write_text(json.dumps(out,indent=2)+'\n')
 for s in states:print(s['physical_name'],s['golden'],s['charts'])
 print('PASS four historical controls and six continuation physical-port charts')
if __name__=='__main__':main()
