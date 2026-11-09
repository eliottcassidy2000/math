import ast,json
from pathlib import Path
root=Path(__file__).resolve().parent
for name in ('check.py','check_kernel.py'):
 tree=ast.parse((root/name).read_text())
 if any(isinstance(n,ast.Assert) for n in ast.walk(tree)):
  raise RuntimeError('assert found in '+name)
 node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='require')
 ns={};exec(compile(ast.Module(body=[node],type_ignores=[]),name,'exec'),ns)
 try:ns['require'](False,'hostile sentinel')
 except RuntimeError as e:
  if str(e)!='hostile sentinel':raise
 else:raise RuntimeError('require failed to reject false condition')
print(json.dumps({'correctness_assert_count':0,'require_rejects_false_in_both_checkers':True,'status':'PASS'}))
