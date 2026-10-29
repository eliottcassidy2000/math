"""Read-only evidence replay; all checks stay active with python -O."""
from pathlib import Path
from hashlib import sha256
import contextlib,io,json,re
import check,replay_inherited
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def need(value,message):
    if not value:raise RuntimeError(message)

def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def seals():
    for filename,base in [('sources.json',ROOT),('manifest.json',HERE)]:
        for name,wanted in json.loads((HERE/filename).read_text(encoding='utf-8')).items():
            need(digest(base/name)==wanted,('changed pinned input',name))

def main():
    seals();fresh={}
    for name,job in [('receipt.json',check.main),('inherited.json',replay_inherited.main)]:
        saved=json.loads((HERE/name).read_text(encoding='utf-8'))
        with contextlib.redirect_stdout(io.StringIO()):value=job(write=False)
        value=json.loads(json.dumps(value))
        need(saved==value,('complete replay disagrees',name));fresh[name]=value
    native=fresh['receipt.json'];inherited=fresh['inherited.json']['controls']
    need([x['minimal_states'] for x in native['prime_clocks']['minimal_joint_banks']]==[60,24,24,2520],'joint clock minima')
    need([x['transporter_0_to_4'] for x in inherited['hostile']]==[27,0],'hostile transporters')
    need([x['binary_attachment_orbits'] for x in inherited['hostile']]==[176,512],'hostile joint attachments')
    need(native['holonomy']['maximal_path_independent_vertex_orbits']==[[0],[1],[2],[3,4,5],[6,7,8]],'five-cell vertex quotient')
    links=0
    for doc in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf-8')):
            if target.startswith(('https://','http://','#')):continue
            need((doc.parent/target.split('#')[0]).resolve().exists(),('broken local link',doc.name,target));links+=1
    seals()
    report={'status':'PASS','optimized':not __debug__,'native_checks':native['active_checks'],
            'named_clock_minima':[60,24,24,2520],'hostile_groups':[81,27],
            'holonomy_orbit_sizes':[1,1,1,3,3],'local_links':links,
            'scope':'exact finite replay plus explicitly numerical Decimal controls; same-agent review'}
    filename='verification.optimized.json' if not __debug__ else 'verification.json'
    (HERE/filename).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2));return report

if __name__=='__main__':main()
