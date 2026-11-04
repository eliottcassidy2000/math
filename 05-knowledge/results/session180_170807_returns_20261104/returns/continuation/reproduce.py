from pathlib import Path
import tempfile,subprocess,json,sys
D=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='continuation-') as t:
 t=Path(t);subprocess.run(['g++','-O3',str(D/'search.cpp'),'-o',str(t/'search')],check=True)
 r=subprocess.run([str(t/'search'),'2'],cwd=t,input=(D/'input.txt').read_text(),text=True,capture_output=True,check=True)
 if json.loads(r.stdout)!=json.loads((D/'search_summary.json').read_text()):raise RuntimeError('summary mismatch')
 if (t/'hits.jsonl').read_text()!=(D/'hits.jsonl').read_text():raise RuntimeError('menu mismatch')
 subprocess.run([sys.executable,*(['-O'] if sys.flags.optimize else []),str(D/'build.py'),'--check'],check=True)
 print('PASS deletion-allowed continuation menu and prefix classification')
