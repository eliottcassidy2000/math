#!/usr/bin/env python3
"""Run the unchanged prime-witness portable suite with temporary build outputs."""
from pathlib import Path
import subprocess,sys,tempfile
with tempfile.TemporaryDirectory(prefix='session142-prime-build-') as work:
 subprocess.run([sys.executable,'-B',str(Path(__file__).resolve().parent/'prime80557/run_checks.py'),'--work-dir',work],check=True)
