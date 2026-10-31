"""Require malformed inputs and out-of-scope constructor targets to be rejected."""
from pathlib import Path
import subprocess
import sys
import tempfile
if len(sys.argv)!=2:raise SystemExit('Usage: hostile_cli.py ADJACENCY_DP_EXECUTABLE')
exe=Path(sys.argv[1]).resolve();root=Path(__file__).resolve().parent.parent
cases=['','0 1\n','23 1\n','2 1\n0\n10\n','2 1\n0x\n10\n','2 1\n10\n10\n','2 1\n01\n10\n','2 1\n01\n00\nextra\n']
with tempfile.TemporaryDirectory(prefix='spectrum-hostile-') as directory:
 p=Path(directory)/'input.txt'
 if subprocess.run([str(exe)],capture_output=True).returncode!=2:raise RuntimeError('Missing arguments accepted')
 for case in cases:
  p.write_text(case)
  if subprocess.run([str(exe),str(p)],capture_output=True).returncode!=2:raise RuntimeError(('Malformed input accepted',case))
 for target in (0,-1,2,7,21,1000001):
  cmd=[sys.executable,*(['-O'] if sys.flags.optimize else []),str(root/'construct_witness.py'),str(target),'--adjacency',str(p)]
  if subprocess.run(cmd,capture_output=True).returncode!=2:raise RuntimeError(('Invalid constructor target accepted',target))
print('PASS: 9 malformed/missing adjacency inputs; 6 out-of-scope constructor targets rejected')
