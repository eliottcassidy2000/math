#!/usr/bin/env python3
from pathlib import Path
from itertools import product
import hashlib,json,subprocess
P=Path(__file__).resolve().parent

def need(b,s):
 if not b:raise RuntimeError(s)
def main():
 banks=json.load(open(P/'physical_banks.json'));controls=json.load(open(P/'physical_contrasts.json'));requests=[];maps=[];records=[]
 for j,ctrl in enumerate(controls):
  m=ctrl['m'];source,target=banks[2*j:2*j+2];n=len(source['rows']);sides=[[i for i in range(m,n)if bool(source['rows'][i]&1)==bool(side)]for side in [0,1]]
  for side,group in enumerate(sides):
   if not group:continue
   for r in sorted({x['r']for x in ctrl['checks']}):
    start=len(requests)
    for sign,B in [(1,source),(-1,target)]:
     for bits in range(16):
      a=(1 if bits&1 else 0)|(2 if bits&2 else 0);c=((1<<r)if bits&4 else 0)|(sum(1<<i for i in group)if bits&8 else 0);value=sum(z for ai,A in enumerate(B['masks'])if A&a==A for ci,C in enumerate(B['masks'])if C&c==C for z in [B['coefficients'][ai*len(B['masks'])+ci]])
      requests.append({'rows':B['rows'],'a':a,'c':c,'full_polynomial_value':value});maps.append(sign*(-1)**bits.bit_count())
    expected=sum(x['difference']for x in ctrl['checks']if x['r']==r and x['i']in group);records.append({'case':j,'m':m,'r':r,'physical_group':group,'grouped_quartic':expected,'query_start':start,'query_count':32})
 text=str(len(requests))+'\n'+''.join(str(len(x['rows']))+'\n'+' '.join(map(str,x['rows']))+'\n'+str(x['a'])+' '+str(x['c'])+'\n'for x in requests);(P/'grouped_queries.txt').write_text(text);subprocess.run(['g++','-O3','-std=c++17',str(P/'queries.cpp'),'-o',str(P/'queries')],check=True);out=subprocess.check_output([str(P/'queries')],input=text,text=True);(P/'grouped_query_values.txt').write_text(out);values=list(map(int,out.split()));need(len(values)==len(requests),'all actual queries')
 for request,value in zip(requests,values):need(value==request['full_polynomial_value'],'physical Boolean query agrees with complete literal polynomial')
 for rec in records:
  start=rec['query_start'];value=sum(values[k]*maps[k]for k in range(start,start+32));need(value==rec['grouped_quartic'],'literal grouped finite difference')
 # Seven cases all have at least one usable grouped certificate.
 need({r['case']for r in records if r['grouped_quartic']}==set(range(len(controls))),'each grouped physical battery detects')
 (P/'grouped_certificates.json').write_text(json.dumps(records,sort_keys=True,indent=2)+'\n');receipt={'status':'PASS','literal_physical_queries':len(values),'grouped_finite_differences':len(records),'physical_cases':len(controls),'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['query_audit.py','queries.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['physical_banks.json','physical_contrasts.json','grouped_queries.txt','grouped_query_values.txt','grouped_certificates.json']}};(P/'query_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
