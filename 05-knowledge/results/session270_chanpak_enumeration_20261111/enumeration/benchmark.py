"""Independent controls and interleaved same-machine benchmark, pinned baselines."""
import hashlib,json,platform,re,statistics,subprocess,sys,tempfile,time,types
from pathlib import Path
from math import factorial,gcd
from counter import count
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
BASE='f2e03e29358d7a6e566f27f75794ec601d0449f5'
DP='05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py'
C='04-computation/a000568_gmp_enum_v4.c'
def source(path): return subprocess.check_output(['git','show',BASE+':'+path])
def direct(n):
    f=factorial(n); total=0
    def rec(rem,lo,parts):
        nonlocal total
        if not rem:
            e=sum(x//2 for x in parts)+sum(gcd(x,y) for i,x in enumerate(parts) for y in parts[:i])
            z=1
            for x in set(parts): z*=x**parts.count(x)*factorial(parts.count(x))
            total+=(f//z)<<e
        for x in range(lo,rem+1,2): rec(rem-x,x,parts+[x])
    rec(n,1,[])
    a,r=divmod(total,f)
    if r: raise ArithmeticError('direct division')
    return a

def main():
    raw=source(DP); module=types.ModuleType('inherited');exec(compile(raw,DP,'exec'),module.__dict__)
    for n in range(31):
        a=direct(n)
        for sw in (3,9,15,999):
            if count(n,switch=sw)[0]!=a: raise AssertionError((n,sw))
        if count(n,gauge=False)[0]!=a: raise AssertionError(n)
    report={'base_commit':BASE,'python':sys.version,'platform':platform.platform(),
      'sources':{DP:hashlib.sha256(raw).hexdigest(),C:hashlib.sha256(source(C)).hexdigest()},
      'new_sources':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest() for p in ['counter.py','counter.cpp']},
      'controls':'all targets 0..30 equal independent direct Burnside; switches3,9,15,999 and raw',
      'scope':'single coefficient engines vs inherited full prefix; compiled C is single coefficient, one thread', 'benchmarks':[]}
    with tempfile.TemporaryDirectory() as tmp:
        cpath=Path(tmp)/'baseline.c';cpath.write_bytes(source(C)); exe=Path(tmp)/'baseline'
        subprocess.run(['gcc','-O3',str(cpath),'-o',str(exe),'-lgmp','-lpthread'],check=True)
        fresh=Path(tmp)/'new'
        subprocess.run(['g++','-O3','-std=c++17',str(Path(__file__).with_name('counter.cpp')),'-o',str(fresh),'-lgmpxx','-lgmp'],check=True)
        for n in range(31):
            for sw in ('9','-1'):
                result=subprocess.run([str(fresh),str(n),sw],capture_output=True,text=True,check=True)
                if int(result.stdout)!=direct(n):raise AssertionError(('cpp',n,sw))
        for n in (80,120,160,200,250):
            times={name:[] for name in ['projective','raw_single','inherited_prefix','cpp_projective','cpp_raw']}
            if n<=120: times['compiled_single']=[]
            stats={}; ref=None
            for repeat in range(1 if n==250 else 3):
                # Rotate order to reduce monotone shared-host bias.
                methods=list(times);methods=methods[repeat%len(methods):]+methods[:repeat%len(methods)]
                for name in methods:
                    t=time.perf_counter()
                    if name.startswith('cpp_'):
                        result=subprocess.run([str(fresh),str(n),'9' if name=='cpp_projective' else '-1'],capture_output=True,text=True,check=True)
                        a=int(result.stdout);r=json.loads(result.stderr)
                    elif name=='compiled_single':
                        out=subprocess.check_output([str(exe),str(n),'1'],text=True)
                        a=int(re.search(r'a\(\d+\) = (\d+)',out)[1]);r={}
                    elif name=='inherited_prefix':
                        arr,r=module.local_dp(n);a=int(arr[n])
                    else: a,r=count(n,gauge=name=='projective',switch=9)
                    times[name].append(time.perf_counter()-t)
                    if ref is None: ref=a
                    if a!=ref: raise AssertionError((n,name))
                    stats[name]={k:v for k,v in r.items() if k!='rows'}
            row={'n':n,'seconds':times,'medians':{k:statistics.median(v) for k,v in times.items()},
                 'statistics':stats,'answer_sha256':hashlib.sha256(str(ref).encode()).hexdigest(),'answer_bits':ref.bit_length()}
            report['benchmarks'].append(row);print(json.dumps(row),flush=True)
            Path(__file__).with_name('receipt.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
