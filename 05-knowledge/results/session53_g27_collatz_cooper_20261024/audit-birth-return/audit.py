import json,hashlib,shutil,pathlib,subprocess,sys,importlib.util
BASE=pathlib.Path(__file__).resolve().parent
ROOT=BASE.parent

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(cmd,cwd):
    p=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True)
    need(p.returncode==0,repr(cmd)+p.stderr)
    return p.stdout

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

summary={}
for folder,manifest in [('session53-g27-growth','MANIFEST.json'),('session53-collatz27-return','manifest.json')]:
    src=BASE/'inputs'/folder; obj=json.loads((src/manifest).read_text()); files={x['path']:x['sha256'] for x in obj['files']} if 'files' in obj else obj
    dst=BASE/'relocated'/folder
    if dst.exists():shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for name,h in files.items():
        need(sha(src/name)==h,'hash mismatch '+name)
        target=dst/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src/name,target)
    shutil.copy2(src/manifest,dst/manifest)
    outputs=[]
    for flags in ([],['-O']):outputs.append(run([sys.executable,*flags,'check.py'],dst))
    need(outputs[0]==outputs[1],'normal optimized mismatch')
    expected='results.json' if 'g27' in folder else 'normal.json'
    need(outputs[0]==(dst/expected).read_text(),'saved output mismatch')
    item={'manifest_sha256':sha(src/manifest),'file_count':len(files),'normal_optimized_saved_equal':True}
    if 'g27' in folder:
        tests=[run([sys.executable,*flags,'test_rejections.py'],dst) for flags in ([],['-O'])]
        need(tests[0]==tests[1]==(dst/'test.out').read_text(),'test output mismatch');item['general_tests']=json.loads(tests[0])
    summary[folder]=item
# A separately written direct ordinary and shortcut orbit certificate for 27.
n=27;j=s=0;b=0;cert=[];x=n
while True:
    odd=x%2;s+=odd;b=3*b+2**j if odd else b
    x=(3*x+1)//2 if odd else x//2;j+=1
    cert.append({'shortcut_time':j,'ordinary_time':j+s,'odd_count':s,'slope_numerator':3**s,'denominator':2**j,'carry':b,'baseline_value':x,'slope_difference':3**s-2**j})
    if x<n:break
need(j==59 and j+s==96 and x==23,'independent first descent')
need(all(c['slope_difference']>0 for c in cert[:-1]) and cert[-1]['slope_difference']<0,'all height slope certificate')
for t in [0,1,2,7,100,2**200+1]:
    n=27+2**59*t;x=n;i=0
    while True:
        x=3*x+1 if x%2 else x//2;i+=1
        if x<n:break
        need(i<=96,'ordinary descent late')
    need(i==96 and x==23+3**37*t,'ordinary lifted descent')
(BASE/'prefix_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
# Mutation tests must fail under normal and optimized Python in relocated copies.
mutations=[
 ('graph_same_component_join','session53-g27-growth','check.py','sum(len(E[i])*len(E[j]) for i,j in combinations(range(len(E)),2))','sum(len(E[i])*len(E[j]) for i in range(len(E)) for j in range(i,len(E)))','test_rejections.py'),
 ('graph_singleton_ports','session53-g27-growth','check.py','for a in sorted(set(p)&N):','for a in sorted(a for a in p if a in N):','test_rejections.py'),
 ('graph_missing_join','session53-g27-growth','check.py','for i,j in combinations(range(len(state)),2):','for i,j in []:','test_rejections.py'),
 ('graph_false_middlecase','session53-g27-growth','test_rejections.py','m if emin or jnext else m+1','m if emin else m+1','test_rejections.py'),
 ('collatz_wrong_odd_slope','session53-collatz27-return','check.py','if e:a*=3;','if e:a*=5;','check.py'),
 ('collatz_missing_final_bit','session53-collatz27-return','check.py','n=27+2**71*h','n=27+2**70*h','check.py'),
 ('collatz_universal_orientation','session53-collatz27-return','check.py',"chi(A)==(-1)**sum(w)","chi(A)==1",'check.py'),
 ('collatz_wrong_source_clock','session53-collatz27-return','check.py','[96,59,37]','[70,59,37]','check.py'),
 ('collatz_equal_total_stopping','session53-collatz27-return','check.py',"prefix(27+2**70,70)[0][-1]!=1","prefix(27+2**70,70)[0][-1]==1",'check.py')]
mut=[]
for label,folder,file,old,new,entry in mutations:
    dst=BASE/'mutants'/label;shutil.copytree(BASE/'relocated'/folder,dst,dirs_exist_ok=True)
    p=dst/file;text=p.read_text();need(text.count(old)==1,'mutation target missing '+label);p.write_text(text.replace(old,new))
    for flags in ([],['-O']):
        r=subprocess.run([sys.executable,*flags,entry],cwd=dst,capture_output=True,text=True)
        need(r.returncode!=0,'mutant survived '+label);mut.append({'name':label,'mode':'optimized' if flags else 'normal','rejection':r.stderr.splitlines()[-1]})
summary['mutations']=mut;summary['independent_prefix_certificate_steps']=len(cert)
(BASE/'audit_results.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
