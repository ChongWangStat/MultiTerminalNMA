#!/usr/bin/env python3
"""Figure 2 (flow_polygon_overview.pdf): original vector figure drawn with matplotlib; no external images used.
From the v15 source (figures_draw.py); output path changed to this folder and the font pinned."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family']='DejaVu Sans'  # pinned: matplotlib's bundled font, so output does not depend on system fonts
matplotlib.rcParams['mathtext.fontset']='dejavusans'
matplotlib.rcParams['pdf.fonttype']=3
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Circle
from matplotlib.lines import Line2D
from pathlib import Path
root=Path(__file__).resolve().parent
fig,axes=plt.subplots(1,2,figsize=(10,3.2), gridspec_kw={'width_ratios':[1.05,1]})
for ax in axes:
 ax.set_axis_off(); ax.set_xlim(0,10); ax.set_ylim(0,7.5)
ink='#1b3445'; teal='#187985'; gray='#69757c'
def node(ax,x,y,s,kind='t'):
 if kind=='t':
  ax.add_patch(Circle((x,y),0.43,fill=False,linewidth=1.4,edgecolor=ink));ax.text(x,y,s,ha='center',va='center',fontsize=10,weight='bold',color=ink)
 else:
  ax.add_patch(plt.Rectangle((x-.43,y-.27),.86,.54,fill=False,linewidth=1.2,edgecolor=teal));ax.text(x,y,s,ha='center',va='center',fontsize=9,color=teal)
def arr(ax,x1,y1,x2,y2,lab=None,off=(0,0),wid=1.9):
 a=FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=11,linewidth=wid,color=teal,shrinkA=1.5,shrinkB=1.5);ax.add_patch(a)
 if lab:ax.text((x1+x2)/2+off[0],(y1+y2)/2+off[1],lab,fontsize=8,color=ink,ha='center',va='center')
a=axes[0]
a.text(.2,7,'A   Electrical currents fixed',fontsize=11,weight='bold',color=ink)
node(a,1.15,3.7,'A');node(a,5,3.7,'C');node(a,9,3.7,'B')
for x,y,s in [(3,5.25,'k1'),(3,2.15,'k2'),(7,5.25,'k3'),(7,2.15,'k4')]:node(a,x,y,s,kind='s')
for x1,y1,x2,y2,off in [(1.7,4.15,2.65,5.0,(0,.25)),(1.7,3.25,2.65,2.35,(0,-.25)),(3.45,5,4.55,4.1,(0,.3)),(3.45,2.35,4.55,3.35,(0,-.3)),(5.45,4.15,6.55,5,(0,.3)),(5.45,3.25,6.55,2.35,(0,-.3)),(7.45,5.0,8.55,4.15,(0,.3)),(7.45,2.35,8.55,3.25,(0,-.3))]:arr(a,x1,y1,x2,y2,'0.5',off)
a.text(5,.52,'Four routes;  t ranges from 0 to 0.5',ha='center',fontsize=9,color=ink)
a.text(5,.08,'Proportional routing selects t = 0.25',ha='center',fontsize=9,color=gray)
b=axes[1]
b.text(.2,7,'B   One covariance-coupled study element',fontsize=11,weight='bold',color=ink)
# Neutral wires: connections to terminals do NOT prescribe input/output pairings.
center=(5.15,3.72)
vertices=[(5.15,5.75),(7.02,3.72),(5.15,1.69),(3.28,3.72)]
# Draw exactly as a hidden star construction: rays end at the four future
# polygon vertices (terminals). Join those vertices into the polygon and
# cover the entire interior. No junction, star spokes, or internal flow is shown.
# Mapping: A -> N, B -> E, C -> W, D -> S.
import math
port_wires=[((1.25,5.65),vertices[0]), ((9.05,5.65),vertices[1]),
            ((1.25,1.68),vertices[3]), ((9.05,1.68),vertices[2])]
for (tx,ty),(px,py) in port_wires:
 dx,dy=px-tx,py-ty
 norm=math.hypot(dx,dy)
 start=(tx+.43*dx/norm,ty+.43*dy/norm)
 b.plot([start[0],px],[start[1],py],color=ink,linewidth=1.7,zorder=1)
b.add_patch(Polygon(vertices,closed=True,edgecolor=ink,facecolor='#1B3449',linewidth=1.6,zorder=2))
# Open engineering terminals are placed precisely at the polygon vertices.
for px,py in vertices:
 b.add_patch(Circle((px,py),radius=.105,facecolor='white',edgecolor=ink,
                    linewidth=1.35,zorder=6))
for x,y,lab in [(1.25,5.65,'A'),(9.05,5.65,'B'),(1.25,1.68,'C'),(9.05,1.68,'D')]:node(b,x,y,lab)
b.text(5.15,4.86,'Study $k$',ha='center',va='center',fontsize=11,weight='bold',color='white',zorder=4)
# Conventional resistor zigzag shown as a GLYPH, not a route through the block.
zx=[4.15,4.32,4.43,4.56,4.69,4.82,4.95,5.08,5.21,5.34,5.47,5.60,5.73,5.89,6.11]
zy=[4.16,4.16,4.40,3.93,4.40,3.93,4.40,3.93,4.40,3.93,4.40,3.93,4.16,4.16,4.16]
b.plot(zx,zy,color='white',linewidth=2.8,solid_joinstyle='round',zorder=5)
b.text(5.15,3.14,r'$\Sigma_k^{-1}$',ha='center',va='center',fontsize=17,color='white',zorder=5)
# the polygon has no internal routing; the matrix label is its only model data
b.text(5.15,.68,'Terminal law:  $i_k=K_kv_k$',ha='center',fontsize=9,color=ink)
b.text(5.15,.19,r'$K_k=C_k^{\mathsf{T}}\Sigma_k^{-1}C_k$',ha='center',fontsize=9,color=gray)
fig.tight_layout(pad=0.7,w_pad=1.3)
out=root/'flow_polygon_overview.pdf';fig.savefig(out,bbox_inches='tight');fig.savefig(root/'flow_polygon_overview.png',dpi=220,bbox_inches='tight');print(out)
