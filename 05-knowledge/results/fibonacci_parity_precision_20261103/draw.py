"""Draw verified n=4 incidence with Pillow, emitting the same SVG geometry."""
from html import escape
from pathlib import Path
import json
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
data = json.loads((HERE/'results.json').read_text())['diagrams']
labels = {'TT4':'Transitive','Q4':'Strong Q4','1_to_C3':'1 → C3','C3_to_1':'C3 → 1'}
blue,black,bg = '#2166ba','#303b49','#f6f8fc'
image = Image.new('RGB',(1300,730),bg)
draw = ImageDraw.Draw(image)
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1300" height="730" viewBox="0 0 1300 730">',
       f'<rect width="1300" height="730" fill="{bg}"/>']


def font(size):
    for candidate in ['C:/Windows/Fonts/segoeui.ttf','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf']:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate,size)
    return ImageFont.load_default(size=size)


def text(x,y,value,size=20):
    draw.text((x,y),value,fill=black,font=font(size),anchor='mm')
    svg.append(f'<text x="{x}" y="{y}" text-anchor="middle" dominant-baseline="central" '
               f'font-family="Segoe UI,DejaVu Sans,sans-serif" font-size="{size}" fill="{black}">{escape(str(value))}</text>')


def line(a,b,color):
    draw.line([a,b],fill=color,width=3)
    svg.append(f'<line x1="{a[0]}" y1="{a[1]}" x2="{b[0]}" y2="{b[1]}" stroke="{color}" stroke-width="3"/>')


def ellipse(box,color,fill=None):
    draw.ellipse(box,outline=color,width=3,fill=fill)
    a,b,c,d=box
    svg.append(f'<ellipse cx="{(a+c)/2}" cy="{(b+d)/2}" rx="{(c-a)/2}" ry="{(d-b)/2}" '
               f'stroke="{color}" stroke-width="3" fill="{fill or "none"}"/>')


text(650,38,'Complement pairs presentations; the quotient retains incidence',29)
text(300,110,'8 marked presentations',25)
text(975,110,'4 ordinary isomorphism classes',25)
for j,b in enumerate(range(4)):
    y=185+j*105
    color=blue if b in data['n4_blue_charts'] else black
    line((160,y),(440,y),color)
    for x,key in [(160,b),(440,b^7)]:
        ellipse((x-23,y-23,x+23,y+23),color,'white')
        text(x,y,str(key),19)
        text(x,y+42,labels[data['n4_chart_classes'][key]],18)
pos={'TT4':(975,185),'1_to_C3':(735,375),'C3_to_1':(1215,375),'Q4':(975,375)}
for edge in data['n4_lines']:
    a,b=edge['ends']
    color=blue if a==b or 'TT4' in (a,b) else black
    if a==b:
        ellipse((925,385,1025,510),color)
    else:
        line(pos[a],pos[b],color)
for name,(x,y) in pos.items():
    ellipse((x-49,y-30,x+49,y+30),black,'white')
    text(x,y,'Q4' if name=='Q4' else 'TT4' if name=='TT4' else labels[name],19)
    text(x+(76 if name=='Q4' else 0),y-54,'degree '+str(data['n4_fibres'][name]),19)
text(300,577,'Four disjoint complement pairs.',20)
text(975,577,'Loops count twice. Every degree is odd.',20)
text(650,625,'Blue and black are actual path-line colours. The quotient star has no perfect matching.',20)
text(650,674,'Adding within-class re-markings gives class pairs: TT4 ↔ Q4 and (1 → C3) ↔ (C3 → 1).',20)
text(650,706,'Actual chart trails: 0—7 and 3—4—5—2—1—6; the latter retains two re-markings inside Q4.',17)
svg.append('</svg>')
image.save(HERE/'complement_n4.png')
(HERE/'complement_n4.svg').write_text('\n'.join(svg)+'\n',encoding='utf-8',newline='\n')
print('Wrote complement_n4.png and complement_n4.svg from the checked receipt.')
