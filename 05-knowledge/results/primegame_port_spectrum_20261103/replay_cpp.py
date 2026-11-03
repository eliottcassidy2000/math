"""Portable full C++ replay, with binaries kept in a temporary directory."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json,shutil,subprocess
HERE=Path(__file__).resolve().parent
def main():
    compiler=shutil.which('g++')
    if compiler is None:raise SystemExit('A C++17 g++ compiler is required.')
    with TemporaryDirectory(prefix='math-primegame-port-') as scratch:
        binaries={}
        for name in ('tree_scan','verify_tree','first_jet'):
            target=Path(scratch)/(name+'.exe')
            subprocess.run([compiler,'-O3','-std=c++17',str(HERE/(name+'.cpp')),'-o',str(target)],check=True,capture_output=True)
            binaries[name]=str(target)
        for n,filename in ((10,'tree_scan.json'),(11,'tree_scan11.json')):
            result=subprocess.run([binaries['tree_scan'],str(n)],check=True,capture_output=True,text=True).stdout
            json.loads(result);(HERE/filename).write_text(result,encoding='utf-8',newline='\n')
        data=json.loads((HERE/'tree_scan11.json').read_text())
        inputs=''.join(str(h)+' '+str(len(w)+1)+' '+' '.join(map(str,w))+'\n' for h,w in data['witnesses'].items())
        result=subprocess.run([binaries['verify_tree']],input=inputs,check=True,capture_output=True,text=True).stdout
        (HERE/'tree_independent.json').write_text(result,encoding='utf-8',newline='\n')
        result=subprocess.run([binaries['first_jet'],'7'],check=True,capture_output=True,text=True).stdout
        (HERE/'first_jet7.json').write_text(result,encoding='utf-8',newline='\n')
    print('Complete tree censuses, every saved witness, and all first jets through seven replayed.')
if __name__=='__main__':main()
