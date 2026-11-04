"""Inherited five actual Q4 presentations; no replacement scalar chart."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
import return6325
from golden_transport import golden,phase
def five_charts(g):
    return {name:dict(golden=golden.step(g,golden.BANK[name]),phase=phase(golden.step(g,golden.BANK[name]))) for name in ('q0','q1','q2','q4','q5')}
