from pathlib import Path
import shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
for case in ('corrupt','missing','extra'):
    with tempfile.TemporaryDirectory() as t:
        root=Path(t)/'packet';shutil.copytree(ROOT,root)
        if case=='corrupt':(root/'sources'/'clock-fractal.md').write_text('corrupt')
        elif case=='missing':(root/'check.py').unlink()
        else:(root/'sources'/'unexpected.txt').write_text('extra')
        p=subprocess.run([sys.executable,str(root/'replay.py'),'--verify-only'],capture_output=True)
        if p.returncode==0:raise RuntimeError('accepted '+case)
print('PASS corruption, missing member, unexpected nested member rejected')
