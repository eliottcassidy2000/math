#!/usr/bin/env python3
"""Explicit inherited theorem replay, separate from the fresh Session49 audit."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent
namespace={'__name__':'inherited_independent_filtration'}
source=P/'archives/inherited_independent_filtration.py.txt'
exec(compile(source.read_text(),str(source),'exec'),namespace)
class OriginalInputAdapter:
    def __truediv__(self,name):
        if name!='05-knowledge/results':raise RuntimeError('unexpected original path')
        return P/'inputs'
filtration,labels=namespace['independent_filtration'](OriginalInputAdapter(),P/'inputs/prior-filtration')
expected=json.loads((P/'inputs/prior-filtration/guard-completeness.normal.json').read_text())['filtration']
if json.dumps(filtration,sort_keys=True)!=json.dumps(expected,sort_keys=True):raise RuntimeError('inherited complete filtration changed')
print(json.dumps({'status':'PASS','checks':namespace['COUNT'],'paired_charge_order':filtration['paired_charge_order'],'critical_two_primary':filtration['critical_two_primary'],'scope':'Earlier independent filtration re-executed, not a newly derived Session49 theorem.'},sort_keys=True))
