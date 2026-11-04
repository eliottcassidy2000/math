import json
from model import D,P,inflate,profile,need
from local_checks import decode

def main():
 records=json.loads((P/'five_child_records.json').read_text());q=[2,4,1];op=[4,1,2];cases=[];lines=[]
 for r in records:
  if r['code'] not in [2,5,8,11]:continue
  target=next(t for t in r['targets'] if t['profile']!=r['profile'])
  for a in range(8):
   for c in range(8):
    row=dict(code=r['code'],family=target['family'],target_code=target['code'],a=a,c=c,indices=[])
    for qq,local in [(q,r['rows']),(op,r['rows']),(op,decode(target['code'],5))]:
     graph,_=inflate(qq,{0:local},a,c);idx=len(lines);row['indices'].append(idx);lines.append(' '.join(map(str,[idx,0,len(graph),1]+graph)))
    cases.append(row)
 (D/'paired_physical.in').write_text(str(len(lines))+'\n'+'\n'.join(lines)+'\n');(D/'paired_records.json').write_text(json.dumps(cases,separators=(',',':'))+'\n')
if __name__=='__main__':main()
