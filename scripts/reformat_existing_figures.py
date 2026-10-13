"""Re-layout existing PDF vector paths; no numerical values are fitted, sampled or aggregated."""
from pathlib import Path
import fitz, json, subprocess, shutil, re
O=Path(__file__).resolve().parents[1]; M=O/'manuscript'
R=O/'audit/figure_previews'; R.mkdir(parents=True,exist_ok=True)
SRC=O/'figures_source'
OUT=M/'figures';OUT.mkdir(exist_ok=True)
W=119/25.4*72
FONT='helv'; BLACK=(0,0,0)
records=[]
def text(p,xy,s,size=8.5,align='left',rot=0,font=FONT):
 s=s.replace('\u2212','-')
 width=fitz.get_text_length(s,fontname=font,fontsize=size)
 x,y=xy
 if align=='center':x-=width/2
 elif align=='right':x-=width
 p.insert_text((x,y),s,fontname=font,fontsize=size,rotate=rot,color=BLACK)
def transform(src,dst):
 sx=dst.width/src.width;sy=dst.height/src.height
 def f(a):return fitz.Point(dst.x0+(a.x-src.x0)*sx,dst.y0+(a.y-src.y0)*sy)
 return f

def pathcopy(p,q,tr,dash=None):
 shape=p.new_shape()
 for it in q['items']:
  if it[0]=='l':shape.draw_line(tr(it[1]),tr(it[2]))
  elif it[0]=='c':shape.draw_bezier(*(tr(v) for v in it[1:]))
  elif it[0]=='re':shape.draw_rect(fitz.Rect(tr(it[1].tl),tr(it[1].br)))
  elif it[0]=='qu':shape.draw_quad(fitz.Quad(*(tr(v) for v in it[1])))
  else:raise ValueError(it[0])
 shape.finish(color=q.get('color'),fill=q.get('fill'),width=max(.45,q.get('width') or .6),dashes=dash or q.get('dashes'),closePath=q.get('closePath',False),even_odd=q.get('even_odd',False),fill_opacity=q.get('fill_opacity') or 1,stroke_opacity=q.get('stroke_opacity') or 1)
 shape.commit()

def embed(raw,dest):
 subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite','-dCompatibilityLevel=1.4','-dEmbedAllFonts=true','-dSubsetFonts=true','-dAutoRotatePages=/None','-sOutputFile='+str(dest),'-c','<< /NeverEmbed [] >> setdistillerparams','-f',str(raw)],check=True)
 # EPS preferred by journal. Keep PDF for pdfLaTeX.
 subprocess.run(['pdftops','-eps','-level3',str(dest),str(dest.with_suffix('.eps'))],check=True)
 raw.unlink()

for fn,outfn,label in [('kkt_by_gradients.pdf','Fig1.pdf','Cumulative gradient evaluations'),('kkt_by_updates.pdf','FigS1.pdf','Accepted updates')]:
 srcdoc=fitz.open(SRC/fn);sp=srcdoc[0];draw=sp.get_drawings();doc=fitz.open();p=doc.new_page(width=W,height=456)
 curves=[]
 for k in range(3):
  src=fitz.Rect(49.82+312*k,35.68,310.5+312*k,267.92)
  dst=fitz.Rect(47,20+151*k,W-5,126+151*k);tr=transform(src,dst)
  p.draw_rect(dst,color=BLACK,width=.6)
  panelcurves=[]
  for q in draw:
   if q['type']=='s' and len(q['items'])>100 and src.x0-1<=q['rect'].x0<src.x1:
    ua=q['color'][0]>.5
    pathcopy(p,q,tr,dash='[3 2] 0' if ua else '[] 0');panelcurves.append(q)
    curves.append({'panel':k,'method':'UA' if ua else 'FA','segments':len(q['items'])})
  # Original tick locations and original labels, relocated at final-size 8.5 pt.
  spans=[s for b in sp.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
  for s in spans:
   x0,y0,x1,y1=s['bbox'];cx=(x0+x1)/2;cy=(y0+y1)/2
   if 275<y0<290 and src.x0-3<cx<src.x1+2:
    pt=tr(fitz.Point(cx,src.y1));p.draw_line(pt,pt+(0,3),color=BLACK,width=.6);text(p,(pt.x,pt.y+13),s['text'],align='center')
   if s['text'].startswith('\u2212') and src.x0-25<cx<src.x0:
    pt=tr(fitz.Point(src.x0,cy));p.draw_line(pt,pt+(-3,0),color=BLACK,width=.6);text(p,(pt.x-5,pt.y+3),s['text'],align='right')
  thresh=tr(fitz.Point(src.x0,259.16)).y;p.draw_line((dst.x0,thresh),(dst.x1,thresh),color=BLACK,width=.6,dashes='[3 3] 0')
  text(p,(dst.x0,12+151*k),'('+chr(97+k)+')',size=10,font='hebo')
  text(p,((dst.x0+dst.x1)/2,151*k+151),label,align='center')
  # One rotated y label at every panel, without subscript-size reduction.
  text(p,(12,108+151*k),'log10(KKT + budget)',rot=90)
  colors={('UA' if q['color'][0]>.5 else 'FA'):q['color'] for q in panelcurves}
  for a,method in enumerate(['FA','UA']):
   yy=33+151*k+a*12; p.draw_line((W-67,yy),(W-47,yy),color=colors[method],width=.9,dashes='[3 2] 0' if method=='UA' else '[] 0');text(p,(W-42,yy+3),method)
 raw=OUT/('raw-'+outfn);doc.save(raw,deflate=True);embed(raw,OUT/outfn)
 records.append({'source':fn,'output':outfn,'width_mm':119,'height_mm':456/72*25.4,'smallest_font_pt':8.5,'curve_paths':curves,'data_change':'none; existing vector segments affinely relocated without dropping points'})

# CV panels: preserve the original vector data paths, all 48 positions, missing path and selection marks.
sp=fitz.open(SRC/'cv_paths_local.pdf')[0];draw=sp.get_drawings();doc=fitz.open();p=doc.new_page(width=W,height=444)
spans=[s for b in sp.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
for k in range(6):
 col=k%2;row=k//2;src=fitz.Rect(54.91+257.86*col,79.75+173.75*row,248.49+257.86*col,171.82+173.75*row)
 dst=fitz.Rect(38+W/2*col,41+row*137,W/2-7+W/2*col,131+row*137);tr=transform(src,dst)
 for q in draw:
  qr=q['rect'];cx=(qr.x0+qr.x1)/2;cy=(qr.y0+qr.y1)/2
  if src.x0-3<cx<src.x1+3 and src.y0-3<cy<src.y1+3:
   if q['fill']==(1.,1.,1.) and qr.width>src.width*.8:continue
   # Copy all original line and point artwork, including diamond selections.
   pathcopy(p,q,tr)
 for s in spans:
  x0,y0,x1,y1=s['bbox'];cx=(x0+x1)/2;cy=(y0+y1)/2
  if s['text'] in ['.01','.025','.06','.15'] and src.x0<cx<src.x1 and src.y1<y0<src.y1+20:
   pt=tr(fitz.Point(cx,src.y1));text(p,(pt.x,dst.y1+12),s['text'],size=8,align='center')
  elif re.fullmatch(r'[\u2212-]?\d\.\d{1,2}',s['text']) and src.x0-40<cx<src.x0 and src.y0<cy<src.y1:
   pt=tr(fitz.Point(src.x0,cy));text(p,(dst.x0-5,pt.y+2.5),s['text'],size=8,align='right')
  elif s['text']=='*' and src.x0<cx<src.x1+10 and src.y0<cy<src.y1:
   pt=tr(fitz.Point(cx,cy));text(p,(min(pt.x,dst.x1-2),pt.y), '*',size=10)
 text(p,(dst.x0,33+137*row),'('+chr(97+k)+')',size=10,font='hebo')
 text(p,((dst.x0+dst.x1)/2,dst.y1+27),'Penalty',size=8.5,align='center')
 # Use one common vertical axis label in left margin, ticks are panel-specific.
 if k==3:text(p,(dst.x0+2,dst.y0+9),'UA .01 unavailable',size=8)
text(p,(10,300),'Three-fold validation score',size=8.5,rot=90)
# Legend taken from original data-path colours, with shapes and dashes retained.
fa=next(q['color'] for q in draw if q['type']=='s' and q['color'] and .4<q['color'][1]<.5 and q['color'][0]<.1)
ua=next(q['color'] for q in draw if q['type']=='s' and q['color'] and q['color'][0]>.8 and q['color'][1]<.5)
for xx,lab,c,dash in [(42,'FA',fa,'[] 0'),(102,'UA',ua,'[4 2] 0')]:
 p.draw_line((xx,9),(xx+18,9),color=c,width=1,dashes=dash);text(p,(xx+23,12),lab,size=8.5)
# A diamond key; retain explicit border contrast.
p.draw_polyline([(167,5),(171,9),(167,13),(163,9),(167,5)],color=BLACK,width=.65)
text(p,(178,12),'Selected',size=8.5);text(p,(235,12),'* Upper boundary',size=8.5)
raw=OUT/'raw-Fig2.pdf';doc.save(raw,deflate=True);embed(raw,OUT/'Fig2.pdf')
records.append({'source':'cv_paths_local.pdf','output':'Fig2.pdf','width_mm':119,'height_mm':444/72*25.4,'smallest_font_pt':8,'data_change':'none; original paths and markers affinely relocated; original missing candidate kept'})
(O/'audit/figure_formatting.json').write_text(json.dumps(records,indent=2))
for fn in ['Fig1','Fig2','FigS1']:
 d=fitz.open(OUT/(fn+'.pdf'));d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(R/(fn+'.png'))
print('formatted',len(records),'existing vector figures')
