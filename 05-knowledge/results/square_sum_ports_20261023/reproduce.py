"""Rebuild the exact universe, then replay all independent frozen consumers."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
BYTE_PROTOCOL='UTF-8 text bytes with CRLF replaced by LF before comparison or SHA-256; no other normalization.'

def portable_bytes(path):
    return path.read_bytes().replace(b'\r\n',b'\n')

def run(command):
    subprocess.run(command,cwd=HERE,check=True)

def main():
    with tempfile.TemporaryDirectory(prefix='square_sum_ports_') as tmp:
        exe=Path(tmp)/('paths.exe' if sys.platform=='win32' else 'paths')
        candidate=Path(tmp)/'paths.json'
        run(['g++','-O3','-std=c++17',str(HERE/'paths.cpp'),'-o',str(exe)])
        run([str(exe),'32',str(candidate)])
        if portable_bytes(candidate)!=portable_bytes(HERE/'paths.json'):raise RuntimeError('exact path universe changed')
        for stem in ('tight_cut_search','marked_cut_search'):
            compiler_output=Path(tmp)/(stem+('.exe' if sys.platform=='win32' else ''))
            result=Path(tmp)/(stem+'.json')
            run(['g++','-O3','-std=c++17',str(HERE/(stem+'.cpp')),'-o',str(compiler_output)])
            run([str(compiler_output),str(result)])
            if portable_bytes(result)!=portable_bytes(HERE/(stem+'.json')):raise RuntimeError(('cut universe changed',stem))
    for optimize in ([],['-O']):
        for script in ('verify.py','recurrence.py','cuts.py','guarded_forbidden.py'):
            run([sys.executable,*optimize,str(HERE/script)])
    manifest=json.loads((HERE/'manifest.json').read_text())
    if manifest['byte_protocol']!=BYTE_PROTOCOL:raise RuntimeError('manifest byte protocol changed')
    for name,digest in manifest['sha256'].items():
        data=portable_bytes(HERE/name)
        if hashlib.sha256(data).hexdigest()!=digest:raise RuntimeError(('manifest mismatch',name))
        # Hostile checkout control: materialize CRLF and canonicalize again.
        if data.replace(b'\n',b'\r\n').replace(b'\r\n',b'\n')!=data:raise RuntimeError(('line-ending round trip',name))
    print('PASS: exact C++ universe, independent edge-subset census, endpoint convolution, cut witnesses, normal/-O receipts, manifest')

if __name__=='__main__':main()
