#!/usr/bin/env python3
import importlib.util, sys
spec = importlib.util.spec_from_file_location(
    "tm", r"C:\Users\Eliott\Documents\GitHub\math\04-computation\lrc14_tourmap_switch-parity_kps-S2-wf.py")
# prevent the module's __main__ from running
import types
src = open(spec.origin, encoding="utf-8").read()
# strip the __main__ blocks
cut = src.split('if __name__ == "__main__":')[0]
mod = types.ModuleType("tm")
exec(compile(cut, spec.origin, "exec"), mod.__dict__)
sys.modules["tm"] = mod

probe_map = mod.probe_map if hasattr(mod,'probe_map') else None
# probe_map/report_probe live AFTER first __main__; re-grab from full source minus mains
# Actually they are defined between the two __main__ blocks; our cut dropped them.
