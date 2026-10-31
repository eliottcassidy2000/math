"""Run after compiling adjacency_dp.cpp; optionally pass its executable path."""
from pathlib import Path
import subprocess,sys,tempfile
p=Path(__file__).resolve().parent;lines=[]
if len(sys.argv)!=2:raise SystemExit('Usage: mutation.py ADJACENCY_DP_EXECUTABLE')
workspace=tempfile.TemporaryDirectory(prefix='spectrum-mutation-')
exe=sys.argv[1]
for f in sorted(p.glob('adjacency_*.txt')):
 data=f.read_text().splitlines();a=[list(x) for x in data[1:]]
 for i,j in ((0,1),(1,0)):a[i][j]='1' if a[i][j]=='0' else '0'
 m=Path(workspace.name)/'mutated.txt';m.write_text(data[0]+'\n'+'\n'.join(''.join(x) for x in a)+'\n')
 r=subprocess.run([exe,str(m)],capture_output=True,text=True)
 if r.returncode!=1:raise RuntimeError(('edge mutation not rejected',f,r.stdout))
 lines.append(f.name+' edge(0,1) flipped: correctly rejected; '+r.stdout.strip())
workspace.cleanup();print('\n'.join(lines))
