"""Generate the original vector diagram; optional Pillow PNG uses same primitives."""
from pathlib import Path
from itertools import combinations
from math import sqrt
from html import escape

ROOT=Path(__file__).resolve().parent
W,H=1040,535
BG='#f7f8fb'; INK='#203040'; BLUE='#2466a0'; RED='#c14343'; MUTED='#526579'
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
     f'<rect width="{W}" height="{H}" fill="{BG}"/>']
try:
    from PIL import Image,ImageDraw,ImageFont
    img=Image.new('RGB',(W,H),BG)
    draw=ImageDraw.Draw(img)
    def font(size):
        return ImageFont.truetype('C:/Windows/Fonts/arial.ttf',size)
except ImportError:
    draw=None

def line(a,b,color=INK,width=2):
    svg.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" x2="{b[0]:.3f}" y2="{b[1]:.3f}" stroke="{color}" stroke-width="{width}"/>')
    if draw: draw.line([a,b],fill=color,width=width)

def circle(p,color,r=8):
    svg.append(f'<circle cx="{p[0]:.3f}" cy="{p[1]:.3f}" r="{r}" fill="{color}"/>')
    if draw: draw.ellipse((p[0]-r,p[1]-r,p[0]+r,p[1]+r),fill=color)

def text(p,value,size=19,color=INK):
    svg.append(f'<text x="{p[0]}" y="{p[1]}" text-anchor="middle" font-family="Arial,sans-serif" font-size="{size}" fill="{color}">{escape(value)}</text>')
    if draw: draw.text(p,value,font=font(size),fill=color,anchor='ms')

text((W/2,42),'One collision forces the A2 flower',28)
left={m:(100+160*((m>>0)&1)+65*((m>>2)&1),
         350-145*((m>>1)&1)-65*((m>>2)&1)) for m in range(8)}
edges=[(a,b) for a,b in combinations(range(8),2) if (a^b).bit_count()==1]
for a,b in edges: line(left[a],left[b],MUTED,2)
for m,p in left.items():
    circle(p,RED if m in (0,7) else BLUE)
    offset={0:(-25,26),1:(23,26),2:(-28,-15),3:(26,0),4:(-25,-15),
            5:(28,7),6:(-28,-15),7:(26,-15)}[m]
    label=''.join(str((m>>i)&1) for i in range(3))
    text((p[0]+offset[0],p[1]+offset[1]),label,18,RED if m in (0,7) else INK)
text((215,425),'8 labels · 12 edges · χ=2',22)
text((215,453),'Abstract cube shown schematically',16,MUTED)

line((405,245),(555,245),INK,3)
line((555,245),(540,235),INK,3)
line((555,245),(540,255),INK,3)
text((480,205),'one fibre',18)
text((480,225),'000 ≡ 111',20,RED)

unit=[(1,0),(-.5,sqrt(3)/2),(-.5,-sqrt(3)/2)]
coords={m:(sum(((m>>i)&1)*unit[i][0] for i in range(3)),
           sum(((m>>i)&1)*unit[i][1] for i in range(3))) for m in range(8)}
right={m:(795+125*coords[m][0],245-125*coords[m][1]) for m in range(8)}
for a,b in edges: line(right[a],right[b],MUTED,2)
for m in range(7):
    circle(right[m],RED if m==0 else BLUE)
    if m:
        p=right[m]; d=coords[m]
        label=''.join(str((m>>i)&1) for i in range(3))
        text((p[0]+28*d[0],p[1]-27*d[1]+6),label,18)
text((795,276),'000 = 111',18,RED)
text((795,425),'7 points · 12 unit edges · χ=3',22)
text((795,453),'Unit spokes and unit hexagon sides',16,MUTED)
text((W/2,506),'Every planar unit-edge drawing of Q3 with seven distinct images has this form.',18)
svg.append('</svg>')
(ROOT/'cube_to_flower.svg').write_text('\n'.join(svg)+'\n',encoding='utf-8')
if draw: img.save(ROOT/'cube_to_flower.png')
print('Wrote cube_to_flower.svg' + (' and cube_to_flower.png' if draw else ''))
