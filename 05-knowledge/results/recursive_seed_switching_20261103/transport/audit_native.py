"""Independent native interpreter: no compiler or producer helper imported."""
from pathlib import Path
import hashlib, importlib.util, json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent
ELEMENTARY = RESULTS/'zero_ready_component_returns_20261103/component/audit_native.py'
spec = importlib.util.spec_from_file_location('elementary_audit', ELEMENTARY)
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)


def main():
    inp = HERE.parent/'selector/named_repair.json'
    program_path = HERE/'native.json'
    data = json.loads(inp.read_text()); target = data['target']; source = data['source']
    program = json.loads(program_path.read_text())['program']
    names = source['names']; rows = source['rows']; A.valid(rows)
    A.need(names == target['names'] and len(names)==23, 'All physical names identical')
    A.need(program['permanent_births']==[], 'No permanent births')
    A.need(A.differences(rows,target['rows'],names)==target['old_arc_differences'], 'Complete change list')
    for z in program['changes']:
        A.need(z['helper'] not in names, 'Fresh helper')
        extended = names+[z['helper']]
        state = A.add_vertex(rows,z['helper_out_mask'])
        for word in z['paths']:
            state = A.native(state,extended,word)
        after = [row & ((1 << len(names))-1) for row in state[:-1]]
        A.need(A.differences(rows,after,names)==[z['changed_pair']], 'One arc per gadget')
        A.need(hashlib.sha256(json.dumps(after,separators=(',',':')).encode()).hexdigest()==z['after_delete_sha256'], 'Intermediate hash')
        for a,b in A.differences(rows,after,names):
            A.need(a not in target['prefix'] and b not in target['prefix'], 'Prefix incidences retained at gadget endpoints')
        rows = after
    A.need(rows==target['rows'], 'Exact target rows')
    for hostile in ([2,1,0],[0,1,1]):
        rejected=False
        try:A.native([6,4,0],[0,1,2],hostile)
        except RuntimeError:rejected=True
        A.need(rejected,'Illegal guard hostile rejected')
    A.need(A.native([6,4,0],[0,1,2],[0,1,2])==[2,4,1], 'Positive complement control')
    pins=[inp,program_path,ELEMENTARY,Path(__file__)]
    print(json.dumps(dict(status='PASS',parent_H=source['H'],target_H=target['H'],
        order=23,peak_order=24,native_complements=2*len(program['changes']),
        helper_births=len(program['changes']),deleted_old_names=[],prefix=target['prefix'],
        old_mark=31,new_mark=target['mark'],source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}),indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
