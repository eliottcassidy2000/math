"""Independent literal controls and large-target decoder tests."""
from decode import *
from math import comb
import tempfile,subprocess,hashlib

def run():
 d=load();need(len(d['even_records'])==24 and len(d['records'])==46,'all70 charts')
 receipts=[];cases=[];families=[]
 for kind,rs in [('even',d['even_records']),('odd',d['records'])]:
  for i,r in enumerate(rs):
   E=profile(r,kind);gc=[sum(c*(b-1)*comb(b,j) for b,c in E if b>=j) for j in range(1,max(b for b,c in E)+1)]
   need(gc==r['growth_Stirling'] and all(x>=0 for x in gc) and gc[0]>0,'recheck exact positive growth coefficients')
   for m in (1,2):
    c=construct(d,kind,i,m,1)
    cases.extend([(c['source_rows'],c['selected'],1,c['U'],None),(c['rows'],0,None,None,c['H'])])
    hits=recognize(c['H'],d);need({'kind':kind,'family':i,'m':m,'t':1} in hits,'target compiler recovers literal construction')
    trace=descent({'kind':kind,'family':i,'m':m,'t':1},d);need(trace[-1]['kind']=='even' and trace[-1]['m']==trace[-1]['t']==1,'finite base descent')
   for m,t in ((19,7),(101,31),(1000,173)):
    N=1+(1<<t)*val(E,m);hits=recognize(N,d)
    need({'kind':kind,'family':i,'m':m,'t':t} in hits,'large exact target recognition')
    trace=descent({'kind':kind,'family':i,'m':m,'t':t},d)
    families.append({'kind':kind,'family':i,'m':m,'t':t,'target_bits':N.bit_length(),'steps':len(trace)-1})
 for N in (1,2,3,5,7,21):need(not recognize(N,d),'small outside-language control')
 with tempfile.TemporaryDirectory() as td:
  exe=str(Path(td)/'literal');subprocess.run(['g++','-O2','-std=c++17',str(HERE/'literal.cpp'),'-o',exe],check=True)
  inp=''.join(str(len(rr))+' '+str(mask)+' '+' '.join(map(str,rr))+'\n' for rr,mask,_,_,_ in cases)
  proc=subprocess.run([exe],input=inp,text=True,capture_output=True,check=True)
  vals=[list(map(int,line.split())) for line in proc.stdout.splitlines()];need(len(vals)==len(cases),'literal response count')
  for case,v in zip(cases,vals):
   rr,mask,delta,u,h=case
   if delta is not None:need(v[1:]==[delta,u],'literal source endpoint packet')
   if h is not None:need(v[0]==h,'literal tail HP count')
   receipts.append({'order':len(rr),'selected':mask,'H':v[0],'delta':v[1],'U':v[2]})
 out={'scope':'70 fixed-free-port charts; exact all-length profiles independently tested on280 literal source/tail graphs and210 large target states. No universal spectrum claim.','literal':receipts,'large':families,'profile_sha256':hashlib.sha256((HERE/'profiles.json.gz').read_bytes()).hexdigest()}
 name='verification'+('.optimized' if not __debug__ else '')+'.json.gz';(HERE/name).write_bytes(gzip.compress((json.dumps(out,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0));print(json.dumps({'literal_graphs':len(cases),'large_targets':len(families)}))
if __name__=='__main__':run()
