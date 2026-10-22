from pathlib import Path
import sys,json
BASE=Path(__file__).resolve().parent
SOURCE=BASE.parent/'session26-runners'
sys.path.insert(0,str(SOURCE))
from probe import curve
rows=json.loads((BASE/'independent-census.json').read_text())
for row in rows:
 _,a=curve(row['vs'])
 if [[str(h),str(c)] for h,c in a.items()]!=row['atoms']:raise ValueError(row['vs'])
for name in ['probe','family']:
 if (SOURCE/(name+'-output.json')).read_bytes()!=(BASE/(name+'-rerun.json')).read_bytes():raise ValueError(name)
 if (SOURCE/(name+'-optimized.json')).read_bytes()!=(BASE/(name+'-rerun-optimized.json')).read_bytes():raise ValueError(name+' optimized')
print('PASS: independent full atom maps match all 164 original census/control rows; normal and optimized artifacts byte-match.')
