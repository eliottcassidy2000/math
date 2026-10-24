"""Replay frozen native supply, controller, compatibility and fibre certificates."""
from pathlib import Path
import subprocess,sys
HERE=Path(__file__).resolve().parent
flags=['-B']+(['-O'] if sys.flags.optimize else [])
for relative in ('supply/check.py','compatibility/resilience.py','supply/shared_owner.py','controller/check.py','controller/padding_probe.py','controller/partner_gluing.py','compatibility/check.py','fibres/reproduce.py'):
    subprocess.run([sys.executable,*flags,str(HERE/relative)],check=True)
