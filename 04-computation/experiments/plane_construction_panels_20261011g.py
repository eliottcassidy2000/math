"""Illustrative vector/raster diagrams; proofs use the exact companion scripts.

Requires Pillow. Drawing coordinates are floating point and are NOT certificates.
"""
from pathlib import Path
from math import sqrt
from html import escape
from PIL import Image, ImageDraw, ImageFont


def main():
    width, height = 1500, 690
    image = Image.new('RGB', (width, height), '#fafbfc')
    draw = ImageDraw.Draw(image)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="#fafbfc"/>']
    fonts = {}

    def font(size):
        if size not in fonts:
            fonts[size] = ImageFont.truetype('C:/Windows/Fonts/segoeui.ttf', size)
        return fonts[size]

    def text(x, y, value, size=20, fill='#182b3a'):
        draw.text((x, y), value, font=font(size), fill=fill)
        svg.append(f'<text x="{x}" y="{y+size}" font-family="Segoe UI,Arial,sans-serif" font-size="{size}" fill="{fill}">{escape(value)}</text>')

    def line(points, fill='#182b3a', weight=2):
        draw.line(points, fill=fill, width=weight)
        svg.append(f'<polyline points="{" ".join(f"{x},{y}" for x,y in points)}" fill="none" stroke="{fill}" stroke-width="{weight}"/>')

    def polygon(points, fill, stroke='#d2d9df', weight=1):
        draw.polygon(points, fill=fill)
        line(points+[points[0]], stroke, weight)
        svg.append(f'<polygon points="{" ".join(f"{x},{y}" for x,y in points)}" fill="{fill}" stroke="{stroke}" stroke-width="{weight}"/>')

    def dot(x, y, radius=6, fill='#182b3a'):
        draw.ellipse((x-radius,y-radius,x+radius,y+radius),fill=fill)
        svg.append(f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{fill}"/>')

    text(30, 20, 'What the constructions preserve', 30)
    text(30, 65, 'Exact constraints, full-plane distance rules, and a restricted periodic model', 21, '#526776')
    for x in (510, 1000):
        line([(x,115),(x,650)], '#d9e0e5', 1)

    text(30, 120, 'Moser spindle', 25)
    text(30, 155, 'Two equality diamonds + one unit edge', 18)
    a = complex(.5, sqrt(3)/2)
    b = complex(5/6, sqrt(11)/6)
    zs = [0j,1+0j,a,1+a,b,a*b,(1+a)*b]
    points = [(72+215*z.real,540-215*z.imag) for z in zs]
    for i in range(7):
        for j in range(i):
            if abs(abs(zs[i]-zs[j])-1) < 1e-9:
                line([points[i],points[j]], '#c94b35' if {i,j}=={3,6} else '#526776', 4)
    for x,y in points:
        dot(x,y)
    for i,label,dx,dy in ((0,'o',-22,0),(3,'u',8,-5),(6,'v',8,-23)):
        x,y=points[i]
        text(x+dx,y+dy,label,22)
    text(30, 570, 'In 3 colours: o = u = v, but uv is an edge.', 18)
    text(30, 600, 'The same equality fails in larger palettes.', 18)

    text(535, 120, 'Stechkin triangle-region residual', 25)
    text(535, 155, 'A unit triangle lies inside one colour class', 18)
    def uv(u,v):
        return (545+145*(u+v/2), 525-145*sqrt(3)*v/2)
    # Draw the same translated triangular gaps used in the exact packet.
    for i in range(3):
        for j in range(3):
            up=[(i+.5,j),(i,j+.5),(i+.5,j+.5)]
            down=[(i+.5,j+1),(i+1,j+.5),(i+.5,j+.5)]
            for triangle,colour in ((up,'#c1d9f0'),(down,'#efcfbd')):
                points0=[uv(*p) for p in triangle]
                if all(530 <= x <= 975 and 195 <= y <= 535 for x,y in points0):
                    polygon(points0,colour)
    tri=[uv(1/3,1/3),uv(4/3,1/3),uv(1/3,4/3)]
    line(tri+[tri[0]],'#c94b35',4)
    for x,y in tri:
        dot(x,y,6,'#a73525')
    text(535, 570, 'Keeping the four hexagon colours cannot', 18)
    text(535, 600, 'leave a two-colour repair for these regions.', 18)

    text(1025, 120, 'Two-size square pattern', 25)
    text(1025, 155, 'A valid 8-colour control at a = 1/sqrt(8)', 18)
    colours=['#99c6ed','#d3e8f8','#edae91','#f8dcca',
             '#9dceb1','#d5eddf','#bea8d8','#e4dcef']
    def square(x,y,w):
        p,q=x//2,y//2
        c=4*((p//2+q//2)%2)+2*(p%2)+(q%2)
        left,top=1060+42*x,535-42*(y+w)
        polygon([(left,top),(left+42*w,top),(left+42*w,top+42*w),(left,top+42*w)],colours[c],'#82919a')
        text(left+21*w-6,top+21*w-13,str(c+1),18)
    for i in range(2):
        for j in range(2):
            w=2 if (i+j)%2==0 else 1
            for x in range(0,4,w):
                for y in range(0,4,w):
                    square(4*i+x,4*j+y,w)
    for x in (1060,1228,1396):
        line([(x,199),(x,535)],'#526776',2)
    for y in (199,367,535):
        line([(1060,y),(1396,y)],'#526776',2)
    text(1025, 570, 'Whole-cell, 8a-period colouring needs 8.', 18)
    text(1025, 600, 'This does not rule out finer 6-colour designs.', 18)
    text(30, 655, 'Illustrations only. Exact coordinates, boundary conventions and proofs are in the companion research packets.', 17, '#526776')
    svg.append('</svg>')
    out=Path(__file__).resolve().parents[2]/'05-knowledge/results'
    image.save(out/'plane_construction_panels_20261011g.png')
    (out/'plane_construction_panels_20261011g.svg').write_text('\n'.join(svg)+'\n',encoding='utf-8')
    print('Rendered PNG and SVG; illustrative floating-point geometry, not a proof certificate.')


if __name__=='__main__':
    main()
