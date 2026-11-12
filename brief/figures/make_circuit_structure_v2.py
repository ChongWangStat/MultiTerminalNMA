import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family']='DejaVu Sans'  # pinned: matplotlib's bundled font, so output does not depend on system fonts
matplotlib.rcParams['mathtext.fontset']='dejavusans'
matplotlib.rcParams['pdf.fonttype']=3
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from pathlib import Path as P
import math

root=P(__file__).resolve().parent
ink='#1b3445'; fill='#1B3449'; gray='#5c6970'

pos={
    'A':(7.5,7.1),
    'B':(0.9,2.4),
    'C':(4.8,-2.2),
    'D':(10.2,-2.2),
    'E':(14.1,2.4),  # symmetric with B about A's x=7.5 axis
}
recovery={'A':(4.8,4.2),'B':(2.95,2.4),'C':(4.8,0.6),'D':(6.65,2.4)}
solidarity={'A':(10.2,4.2),'C':(8.35,2.4),'D':(10.2,0.6),'E':(12.05,2.4)}

fig,ax=plt.subplots(figsize=(10.5,6.7))
ax.set_axis_off(); ax.set_xlim(-1,17.2); ax.set_ylim(-4.1,10.1)

def draw_node(k,label):
    x,y=pos[k]
    ax.add_patch(Circle((x,y),0.56,fill=False,lw=1.6,edgecolor=ink,zorder=6))
    ax.text(x,y+0.02,label,ha='center',va='center',fontsize=12,weight='bold',color=ink,zorder=7)

def shorten_line(p1,p2,dr1=0.56,dr2=0.12):
    x1,y1=p1; x2,y2=p2
    dx,dy=x2-x1,y2-y1; L=math.hypot(dx,dy)
    ux,uy=dx/L,dy/L
    return (x1+ux*dr1,y1+uy*dr1),(x2-ux*dr2,y2-uy*dr2)

def draw_polygon(terminals, title, title_shift=0.42, resistor_shift=0.0, title_size=10.7):
    keys=list(terminals.keys())
    pts=[terminals[k] for k in keys]
    ax.add_patch(Polygon(pts, closed=True, facecolor=fill, edgecolor=ink, lw=1.6, zorder=4))
    for p in pts:
        ax.add_patch(Circle(p,0.12,facecolor='white',edgecolor=ink,lw=1.35,zorder=5))
    cx=sum(p[0] for p in pts)/len(pts); cy=sum(p[1] for p in pts)/len(pts)
    ax.text(cx, cy+title_shift, title, color='white', fontsize=title_size, weight='bold', ha='center', zorder=6)
    zx=[cx-0.85,cx-0.68,cx-0.55,cx-0.40,cx-0.25,cx-0.10,cx+0.05,cx+0.20,cx+0.35,cx+0.50,cx+0.67,cx+0.85]
    zy=[cy+resistor_shift,cy+resistor_shift,cy+0.24+resistor_shift,cy-0.18+resistor_shift,cy+0.24+resistor_shift,cy-0.18+resistor_shift,
        cy+0.24+resistor_shift,cy-0.18+resistor_shift,cy+0.24+resistor_shift,cy-0.18+resistor_shift,cy+resistor_shift,cy+resistor_shift]
    ax.plot(zx,zy,color='white',lw=2.3,zorder=6)

def draw_wire(treat, term):
    p1,p2=shorten_line(pos[treat],term)
    ax.plot([p1[0],p2[0]],[p1[1],p2[1]],color=ink,lw=1.5,zorder=2)


def bezier(P0,P1,P2,P3,n=220):
    pts=[]
    for i in range(n+1):
        t=i/n
        x=(1-t)**3*P0[0]+3*(1-t)**2*t*P1[0]+3*(1-t)*t**2*P2[0]+t**3*P3[0]
        y=(1-t)**3*P0[1]+3*(1-t)**2*t*P1[1]+3*(1-t)*t**2*P2[1]+t**3*P3[1]
        pts.append((x,y))
    return pts

def draw_rotated_box_with_terminals(center, ux, uy, length=1.28, height=0.52, terminal_r=0.11):
    vx,vy=-uy,ux
    hl=length/2; hh=height/2
    corners=[]
    for sx,sy in [(-hl,-hh),(hl,-hh),(hl,hh),(-hl,hh)]:
        corners.append((center[0]+ux*sx+vx*sy, center[1]+uy*sx+vy*sy))
    ax.add_patch(Polygon(corners, closed=True, facecolor=fill, edgecolor=ink, lw=1.4, zorder=4))
    left=(center[0]-ux*hl, center[1]-uy*hl)
    right=(center[0]+ux*hl, center[1]+uy*hl)
    for p in [left,right]:
        ax.add_patch(Circle(p,terminal_r,facecolor='white',edgecolor=ink,lw=1.2,zorder=5))
    # internal zigzag in white
    offs=[-0.38,-0.28,-0.18,-0.08,0.02,0.12,0.22,0.32,0.42]
    zig=[]
    for j,o in enumerate(offs):
        amp=0 if j in [0,len(offs)-1] else (0.10 if j%2 else -0.10)
        zig.append((center[0]+ux*o+vx*amp, center[1]+uy*o+vy*amp))
    ax.plot([p[0] for p in zig],[p[1] for p in zig],color='white',lw=2.0,zorder=6)
    return left,right

def draw_two_terminal_study(start_key,end_key,ctrl1,ctrl2,label,lab_xy,box_len=1.28):
    p0=pos[start_key]; p3=pos[end_key]
    pts=bezier(p0,ctrl1,ctrl2,p3)
    mid_idx=len(pts)//2
    mid=pts[mid_idx]
    before=pts[mid_idx-3]; after=pts[mid_idx+3]
    dx,dy=after[0]-before[0],after[1]-before[1]
    L=math.hypot(dx,dy)
    ux,uy=dx/L,dy/L
    left,right=draw_rotated_box_with_terminals(mid,ux,uy,length=box_len,height=0.50,terminal_r=0.11)

    # Crop the Bézier at the circular treatment-node boundaries.  Rendering
    # a transparent treatment circle over a wire would falsely suggest that
    # the conductor passes THROUGH the treatment node.
    def clip_to_node_boundary(points, node_center, radius=.56, from_start=True):
        q=points if from_start else list(reversed(points))
        def dist(t): return math.hypot(t[0]-node_center[0],t[1]-node_center[1])
        j=next(i for i,x in enumerate(q) if dist(x)>=radius)
        if j==0: cut=q[0]
        else:
            inside,outside=q[j-1],q[j]
            lo,hi=0.0,1.0
            for _ in range(30):
                m=(lo+hi)/2
                pt=(inside[0]+m*(outside[0]-inside[0]),inside[1]+m*(outside[1]-inside[1]))
                if dist(pt)<radius: lo=m
                else: hi=m
            cut=(inside[0]+hi*(outside[0]-inside[0]),inside[1]+hi*(outside[1]-inside[1]))
        return [cut]+q[j:]

    pts=clip_to_node_boundary(pts,p0,from_start=True)
    pts=list(reversed(clip_to_node_boundary(pts,p3,from_start=False)))
    # Remove the center of the curve around the black resistor box.
    def nearest_index(target):
        return min(range(len(pts)), key=lambda i:(pts[i][0]-target[0])**2+(pts[i][1]-target[1])**2)
    i1=nearest_index(left)
    i2=nearest_index(right)
    if i1>i2: i1,i2=i2,i1
    seg1=pts[:i1+1]+[left]
    seg2=[right]+pts[i2:]
    ax.plot([z[0] for z in seg1],[z[1] for z in seg1],color=ink,lw=1.35,zorder=2)
    ax.plot([z[0] for z in seg2],[z[1] for z in seg2],color=ink,lw=1.35,zorder=2)
    # Label sits outside the component on a white background, so that
    # neither the resistor body nor adjacent wires obscure study names.
    ax.text(lab_xy[0],lab_xy[1],label,ha='center',va='center',fontsize=9,color=gray,
            bbox=dict(facecolor='white',edgecolor='none',boxstyle='round,pad=0.22',alpha=0.98),zorder=8)

# draw nodes
for k,lab in [('A','A'),('B','B'),('C','C'),('D','D'),('E','E')]:
    draw_node(k,lab)
# wires to polygons
for t,p in recovery.items(): draw_wire(t,p)
for t,p in solidarity.items(): draw_wire(t,p)
# polygons
rec_order={'A':recovery['A'],'D':recovery['D'],'C':recovery['C'],'B':recovery['B']}
draw_polygon(rec_order,'RECOVERY',title_shift=0.54,resistor_shift=-0.20,title_size=10.1)
sol_order={'A':solidarity['A'],'E':solidarity['E'],'D':solidarity['D'],'C':solidarity['C']}
draw_polygon(sol_order,'SOLIDARITY',title_shift=0.54,resistor_shift=-0.20,title_size=10.1)
# two-arm study elements with similar styling
# move labels away enough and use element boxes instead of bare zigzag lines
# CoDEX and REMAP-CAP between A and B; ACTT-1 between A and E
draw_two_terminal_study('A','B',(4.0,6.05),(2.25,4.20),'CoDEX',(3.55,4.16),box_len=1.16)
draw_two_terminal_study('A','B',(3.25,7.05),(1.00,5.25),'REMAP-CAP',(2.20,6.68),box_len=1.40)
draw_two_terminal_study('A','E',(11.0,7.4),(13.9,6.0),'ACTT-1',(12.4,5.32),box_len=1.18)

ax.text(8.0,9.6,'Illustrative five-study circuit before selecting a comparison',ha='center',va='center',fontsize=13,weight='bold',color=ink)
ax.text(8.0,8.95,'Treatment nodes are connected by two-terminal resistive elements and multi-terminal study elements.',ha='center',fontsize=9.5,color=gray)
ax.text(8.0,-3.4,'Any target comparison is induced later by injecting +1 at one treatment node and extracting -1 at another.',ha='center',fontsize=9,color=gray)
ax.text(12.85,-2.86,'Black polygons/boxes: study resistive elements\nWhite rings: explicit study terminals\nWhite zigzag: symbolic resistor glyph',fontsize=8.6,color=gray,ha='left')
out_pdf=root/'COVID19_Circuit_Structure_NoFlow.pdf'
out_png=root/'COVID19_Circuit_Structure_NoFlow.png'
fig.savefig(out_pdf,bbox_inches='tight')
fig.savefig(out_png,dpi=220,bbox_inches='tight')
print(out_pdf)
