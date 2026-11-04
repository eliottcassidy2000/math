"""Normal/-O deterministic replay and exact receipt comparison."""
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
def run():
    results=[]
    for script,receipt in [('check.py','certificate.json'),('clocks.py','clocks.json'),('native_reset.py','native_reset.json'),('audit.py','audit.json')]:
        subprocess.run([sys.executable,'-X','utf8','-B',str(HERE/script)],check=True)
        a=(HERE/receipt).read_bytes()
        subprocess.run([sys.executable,'-O','-X','utf8','-B',str(HERE/script)],check=True)
        b=(HERE/receipt).read_bytes()
        if a!=b:raise RuntimeError('optimized mismatch: '+script)
        results.append({'script':script,'receipt':receipt,'sha256':hashlib.sha256(a).hexdigest(),'normal_optimized_identical':True})
    (HERE/'reproduction.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':run()
