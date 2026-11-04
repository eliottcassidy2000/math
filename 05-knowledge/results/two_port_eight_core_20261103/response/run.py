"""Exact12-core converse hostile: all endpoint tensors and full two-port battery."""
import hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROWS=(3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def need(ok,message):
    if not ok:raise RuntimeError(message)
def main():
    pins={name:sha((HERE/name).read_bytes()) for name in ('run.py','check.cpp')}
    controls=[ROWS,tuple(4095^(1<<i)^r for i,r in enumerate(ROWS))]
    with tempfile.TemporaryDirectory(prefix='twelve-full-response-') as tmp:
        tmp=Path(tmp);inp=tmp/'input.txt';exe=tmp/'check.exe'
        inp.write_text('2\n'+'\n'.join('12 '+' '.join(map(str,rows)) for rows in controls)+'\n',encoding='ascii')
        compiler=shutil.which('g++');need(compiler,'existing g++ required')
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'check.cpp'),'-o',str(exe)],check=True,capture_output=True)
        result=json.loads(subprocess.run([str(exe),str(inp)],check=True,stdout=subprocess.PIPE,text=True,encoding='utf8').stdout)
    need(pins=={name:sha((HERE/name).read_bytes()) for name in pins},'source changed')
    a,b=result['controls'];n=12
    need(a['S']==b['E'] and a['E']==b['S'],'converse endpoint involution')
    need(all((x-y)%2==0 for x,y in zip(a['S'],a['E'])),'endpoint-parity hostile')
    skew=[]
    for i in range(n):
        for j in range(n):
            need(a['M'][i*n+j]==b['M'][j*n+i] and a['C'][i*n+j]==b['C'][j*n+i],'ordered pair involution')
            if i<j and (a['C'][i*n+j]-a['C'][j*n+i])%4:skew.append([i,j,a['C'][i*n+j],a['C'][j*n+i]])
            for k in range(n):
                need(a['U'][(i*n+j)*n+k]==b['V'][(k*n+j)*n+i],'two-path endpoint involution')
                need(a['V'][(i*n+j)*n+k]==b['U'][(k*n+j)*n+i],'other two-path endpoint involution')
                for l in range(n):need(a['D'][((i*n+j)*n+k)*n+l]==b['D'][((l*n+k)*n+j)*n+i],'full three-path reversal involution')
    masks=[m for m in range(1<<n) if m.bit_count()<=2]
    differences=[i for i,(x,y) in enumerate(zip(a['battery'],b['battery'])) if x!=y]
    need(differences and len(skew)==28,'full-response hostile and modular boundary')
    at=differences[0];witness={'old_to_x_mask':masks[at//len(masks)],'old_to_y_mask':masks[at%len(masks)],'port_arc':'x->y','H_T_plus_ports':a['battery'][at],'H_converse_plus_same_ports':b['battery'][at]}
    result.update(source_sha256=pins,endpoint_parity_balanced=True,nonzero_mod4_skew_pairs=skew,distinguishing_battery_entries=len(differences),first_actual_query_witness=witness,scope='One specified12-core and its converse; complete tensors and6241-query interpolation batteries, plus9 high-weight direct checks each. Not a complete order12 census.')
    print(json.dumps(result,sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
