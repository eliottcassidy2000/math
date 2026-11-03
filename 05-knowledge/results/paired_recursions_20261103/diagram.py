"""Render the exact n4 commuting-pairing example; Pillow is only for this figure."""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

HERE=Path(__file__).resolve().parent
im=Image.new('RGB',(1000,580),'#f7f8fb');d=ImageDraw.Draw(im)
fontpath=Path('C:/Windows/Fonts/arial.ttf')
def font(n):return ImageFont.truetype(str(fontpath),n) if fontpath.exists() else ImageFont.load_default(size=n)
def text(x,y,s,size=23,color='#172b46',anchor='mm'):
    d.text((x,y),s,font=font(size),fill=color,anchor=anchor)
text(500,45,'Four classes; only one free pairing axis commuting with converse',28)
text(500,88,'a(4) = 4        self-converse classes = 2        maximum rank = 1',21)
for y in (230,390):
    d.line([(290,y),(710,y)],fill='#cb6819',width=5)
    text(500,y-25,'chosen pairing J',20,color='#a94d0b')
for x in (250,750):
    d.arc((x-55,150,x+55,220),190,530,fill='#3269a8',width=4)
    text(x,128,'fixed by converse R',18,color='#3269a8')
for y,labels in ((230,('Transitive T4  |  H = 1','Strong Q4  |  H = 5')),
                 (390,('Source + C3  |  H = 3','C3 + sink  |  H = 3'))):
    for x,label in zip((250,750),labels):
        d.ellipse((x-28,y-28,x+28,y+28),fill='#ffffff',outline='#172b46',width=3)
        text(x,y+53 if y==230 else y-58,label,22)
d.arc((250,340,750,475),0,180,fill='#3269a8',width=4)
text(500,495,'converse R exchanges the bottom pair',20,color='#3269a8')
text(500,532,'J changes H on the top row. No free H-preserving pairing can cover all four.',21)
im.save(HERE/'pairing_n4.png')
print(HERE/'pairing_n4.png')
