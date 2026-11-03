"""Replay the inherited actual 81/27 control; no new result claimed."""
from pathlib import Path
import importlib.util
import json
import sys

base=Path(__file__).resolve().parent.parent/'primegame_port_spectrum_20261103'
sys.path.insert(0,str(base))
spec=importlib.util.spec_from_file_location('prior_hostile_verifier',base/'verify.py')
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
print(json.dumps(mod.hostile(),indent=2,sort_keys=True))
