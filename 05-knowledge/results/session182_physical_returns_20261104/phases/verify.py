"""Actual named endpoint and observer audit, inherited algorithm with new inputs."""
import json,subprocess,hashlib
from pathlib import Path
from codec import chart,need
from automorphism import automorphisms
P=Path(__file__).resolve().parent

def run(rows):
 n=len(rows);indices=list(range(n));inp=' '.join(map(str,[n,n,*rows,*indices]))
 text=subprocess.run([str(P/'joint')],input=inp,text=True,capture_output=True,check=True).stdout
 h=None;fields={};joints=[]
 for line in text.splitlines():
  a=line.split()
  if a[0]=='H':h=int(a[1])
  elif a[0]=='PORT':fields[int(a[1])]=list(map(int,a[2:]))
  elif a[0]=='J':joints.append(dict(zip(['end_index','start_index','value'],map(int,a[1:]))))
 return h,fields,joints

def main():
 old_path=P/'input_audit.json';old=json.loads(old_path.read_text())
 for s in old['states']:need(json.loads(json.dumps(chart(s['golden'])))==s['charts'])
 out={'scope':'Recomputed named endpoint fields and observer codecs; not native legality or field actions as automorphisms. Independent of discovery counts, using inherited S180 endpoint DP and codec.','codec_source_ref':'4ad839b20dd280cdcfbadf86cb62b2be3c9615ad','historical_codec_controls':len(old['states']),'field_modulus':[1,0,1,1],'states':[],'hostile':[],'input_sha256':{}}
 for filename in ['source.json','candidate2.json','candidate4.json']:
  path=P/filename;c=json.loads(path.read_text());names=c['names'];h,f,j=run(c['rows']);need(h==c['target_H']);records=[]
  for i in range(len(names)):
   e1,s1,deletion,e2,s2=f[i];need((s2-e1)%2==0);g=[h,(s2-e1)//2,(s2+e1)//2]
   records.append({'physical_name':names[i],'fields':f[i],'golden':g,'charts':chart(g)})
  for rec in j:rec['end']=names[rec.pop('end_index')];rec['start']=names[rec.pop('start_index')]
  root=names.index(39);need(f[root][3]-h==c['root_defect']);inports=[names[i] for i,r in enumerate(c['rows']) if r>>root&1];need(sorted(inports)==sorted(c['in_ports']))
  out['states'].append({'file':filename,'H':h,'fields_order':['E1','S1','deletion_H','E2','S2'],'root_defect':f[root][3]-h,'in_ports':inports,'ports':records,'joints':j})
  out['input_sha256'][filename]=hashlib.sha256(path.read_bytes()).hexdigest()
  print(filename,h,'root defect',f[root][3]-h)
  for name in [39,50]:
   r=records[names.index(name)];print(name,r['golden'],r['charts'])
 for rows,order,response in zip(old['hostile']['rows'],[81,27],[1061775,1395279]):
  h,_,_=run(rows);need(h==462915);n=len(rows);cut=256
  # cut bits mean old vertices beating the newborn, matching inherited control.
  child=[r|(((cut>>i)&1)<<n) for i,r in enumerate(rows)]+[((1<<n)-1)^cut]
  hc,_,_=run(child);need(hc==response)
  autos,nodes=automorphisms(rows);need(len(autos)==order);need(len(set(autos))==order)
  out['hostile'].append({'H':h,'cut':cut,'response':hc,'full_automorphism_order':order,'search_nodes':nodes})
 out['input_sha256']['input_audit.json']=hashlib.sha256(old_path.read_bytes()).hexdigest()
 (P/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS all named fields, directed joints, historical codecs and full 81/27 hostile controls')
if __name__=='__main__':main()
