"""Check final source provenance, cache integrity, and small independent controls."""
import ast,hashlib,json,subprocess,tempfile
from pathlib import Path
from benchmark import direct
from counter import count
from cache import load
P=Path(__file__).resolve().parent
receipt=json.loads((P/'receipt.json').read_text())
final=(P/'counter.py').read_text()
timed=final.replace('# Gauge coordinates b_j=(j*q_1-q_j)/2; b_1 vanishes identically.',
                    '# Gauge coordinates b_j=j*q_1-q_j; b_1 vanishes identically.')
if hashlib.sha256(timed.encode()).hexdigest()!=receipt['new_sources']['counter.py']:
    raise AssertionError('exact timed source reconstruction failed')
if ast.dump(ast.parse(timed))!=ast.dump(ast.parse(final)):
    raise AssertionError('executable Python statements changed after timing')
values=load()
for n in range(31):
    ref=direct(n)
    if values[n]!=ref:raise AssertionError(('cache',n))
    for opts in ({'gauge':False},{'switch':3},{'switch':9},{'switch':15},{'switch':999},{'terminal3':False}):
        if count(n,**opts)[0]!=ref:raise AssertionError((n,opts))
with tempfile.TemporaryDirectory() as tmp:
    exe=Path(tmp)/'counter'
    subprocess.run(['g++','-O3','-std=c++17',str(P/'counter.cpp'),'-lgmpxx','-lgmp','-o',str(exe)],check=True)
    for n in list(range(31))+[120]:
        for switch in ('9','-1'):
            r=subprocess.run([str(exe),str(n),switch],capture_output=True,text=True,check=True)
            if int(r.stdout)!=values[n]:raise AssertionError(('cpp',n,switch))
report={'status':'PASS','python_final_source_sha256':hashlib.sha256(final.encode()).hexdigest(),
 'timed_python_source_sha256':hashlib.sha256(timed.encode()).hexdigest(),
 'exact_timed_source_reconstructed':True,'python_executable_AST_unchanged':True,
 'post_timing_change':'One comment corrected b_j to (j*q_1-q_j)/2. No executable statements changed.',
 'final_python_controls':'Independent direct Burnside all n=0..30, six recurrence modes, cache equality.',
 'final_cpp_controls':'Compiled current source with -O3; raw/projective both match checked cache at all n=0..30 and120.',
 'cache':'SHA256 verified; all251 terms previously independently checked by two full-prefix recurrences.',
 'final_script_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(P.glob('*.py'))},
 'final_cpp_sha256':hashlib.sha256((P/'counter.cpp').read_bytes()).hexdigest()}
(P/'FINAL_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
