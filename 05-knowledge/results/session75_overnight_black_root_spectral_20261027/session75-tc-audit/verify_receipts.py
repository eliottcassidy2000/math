import hashlib,json,pathlib
R=pathlib.Path(__file__).resolve().parent

def need(c,m):
 if not c:raise RuntimeError(m)
expected='''Z3xZ1 sets=2 min_tc=1 small 1:2
Z5xZ1 sets=4 min_tc=3 small 3:4
Z7xZ1 sets=8 min_tc=9 small 9:2
Z9xZ1 sets=16 min_tc=39 small
Z3xZ3 sets=16 min_tc=39 small
Z11xZ1 sets=32 min_tc=1729 small
Z13xZ1 sets=64 min_tc=94949 small
Z15xZ1 sets=128 min_tc=162165 small
Z17xZ1 sets=256 min_tc=792878221 small
'''
need((R/'cayley_independent.out').read_text()==expected,'Cayley counts')
need((R/'vt15_independent.out').read_text()=='15 H=197030475 Aut=1215 tc=162165\n'*8+'15 H=193215375 Aut=375 tc=515241\n'*8,'VT15 counts')
need((R/'vt15.normal.out').read_bytes()==(R/'04-computation/experiments/tc_fibres_n12_20261027/vt15_noncayley.out').read_bytes(),'VT15 source replay')
for stem in ['analytic','group_checks']:
 need((R/(stem+'.normal.json')).read_bytes()==(R/(stem+'.optimized.json')).read_bytes(),stem+' optimization')
a=json.loads((R/'analytic.normal.json').read_text());g=json.loads((R/'group_checks.normal.json').read_text())
need(a['greedy_checks']==192 and a['tail_exact_base'] and len(a['configs'])==24,'analytic')
need(g['subgroups']==50 and g['transitive_subgroups']==13 and g['remaining_regular_C3xC3']==4,'degree9')
W=pathlib.Path('/tmp/session75-tc-window-audit')
s=json.loads((W/'SUMMARY.json').read_text())
for name,expected in [('normal',dict(jobs=83,strong=185642,autcalls=48908,checks=93646,failures=0,hits=0)),('n12',dict(jobs=45,strong=46968,autcalls=10104,checks=45384,failures=0,hits=0))]:
 for key,value in expected.items():need(s[name][key]==value,name+' '+key)
print('PASS: analytic/group normal/-O;526 Cayley sets;16 VT15 substitutions;83-window plus45-n12 receipts.')
