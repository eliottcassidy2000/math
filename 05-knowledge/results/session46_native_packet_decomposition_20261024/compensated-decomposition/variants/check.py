#!/usr/bin/env python3
"""Run the parent independent verifier including every certified route variant."""
from pathlib import Path
import runpy,sys
sys.argv.append('--include-variants')
runpy.run_path(str(Path(__file__).resolve().parents[1]/'check.py'),run_name='__main__')
