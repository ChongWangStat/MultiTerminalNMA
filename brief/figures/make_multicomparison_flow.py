import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family']='DejaVu Sans'  # pinned: matplotlib's bundled font, so output does not depend on system fonts
matplotlib.rcParams['mathtext.fontset']='dejavusans'
matplotlib.rcParams['pdf.fonttype']=3
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, FancyArrowPatch
from pathlib import Path as P
import pandas as pd
import math

root=P(__file__).resolve().parent
ink='#22384A'; fill='#172C40'; gray='#69757c'
study_colors={'RECOVERY':'#B45A14','SOLIDARITY':'#3767AF','ACTT-1':'#23733F','CoDEX':'#198A9A','REMAP-CAP':'#6A48A8'}

# ---------------------------------------------------------------------
# Freeze the comparison-specific flow figure onto the exact same circuit
# layout used in Figure 1.
# ---------------------------------------------------------------------
pos={
    'A':(7.5,7.1),
    'B':(0.9,2.4),
    'C':(4.8,-2.2),
    'D':(10.2,-2.2),
    'E':(14.1,2.4),
}
recovery={'A':(4.8,4.2),'B':(2.95,2.4),'C':(4.8,0.6),'D':(6.65,2.4)}
solidarity={'A':(10.2,4.2),'C':(8.35,2.4),'D':(10.2,0.6),'E':(12.05,2.4)}
study_id_to_name={'S1':'RECOVERY','S2':'SOLIDARITY','S3':'ACTT-1','S4':'CoDEX','S5':'REMAP-CAP'}
display_names={'RECOVERY':'S1','SOLIDARITY':'S2','ACTT-1':'S3','CoDEX':'S4','REMAP-CAP':'S5'}

direct_specs={
    'CoDEX': {
        'nodes': ('A','B'),
        'ctrl1': (4.0,6.05), 'ctrl2': (2.25,4.20),
        'label': (3.10,3.75), 'box_len': 1.16,
    },
    'REMAP-CAP': {
        'nodes': ('A','B'),
        'ctrl1': (3.25,7.05), 'ctrl2': (1.00,5.25),
        'label': (1.55,6.95), 'box_len': 1.40,
    },
    'ACTT-1': {
        'nodes': ('A','E'),
        'ctrl1': (11.0,7.4), 'ctrl2': (13.9,6.0),
        'label': (11.90,4.95), 'box_len': 1.18,
    },
}

edges=pd.read_csv(root/'data'/'covid_CSP_edges_common.csv')
ae=pd.read_csv(root/'data'/'covid_AE_study_polygon_currents.csv')

def currents_from_edges(target):
    sub=edges[edges.target==target]
    out={}
    for (study,id_),g in sub.groupby(['study','id']):
        name=study_id_to_name[study]
        flows={}
        for _,r in g.iterrows():
            w=float(r['weight'])
            flows[r['from']]=flows.get(r['from'],0)+w
            flows[r['to']]=flows.get(r['to'],0)-w
        out[name]=flows
    return out

def currents_AE():
    out={}
    for study,g in ae.groupby('study'):
        flows={}
        for _,r in g.iterrows():
            flows[r['t']]=float(r['signed_current'])
        out[study]=flows
    return out

currents={'A:E':currents_AE(),'A:B':currents_from_edges('A:B'),'B:E':currents_from_edges('B:E')}
for target in currents:
    for s in study_colors:
        currents[target].setdefault(s,{})

# --------------------------- drawing helpers --------------------------
def add_node(ax,k,r=0.56):
    x,y=pos[k]
    ax.add_patch(Circle((x,y),r,facecolor='white',edgecolor=ink,lw=1.5,zorder=12))
    ax.text(x,y+0.02,k,ha='center',va='center',fontsize=12,weight='bold',color=ink,zorder=13)


def draw_wire(ax,p1,p2,lw=1.35,z=1):
    ax.plot([p1[0],p2[0]],[p1[1],p2[1]],color=ink,lw=lw,zorder=z)


def shorten_line(p1,p2,dr1=0.56,dr2=0.12):
    x1,y1=p1; x2,y2=p2
    dx,dy=x2-x1,y2-y1
    L=math.hypot(dx,dy)
    ux,uy=dx/L,dy/L
    return (x1+ux*dr1,y1+uy*dr1),(x2-ux*dr2,y2-uy*dr2)


def bezier(P0,P1,P2,P3,n=220):
    pts=[]
    for i in range(n+1):
        t=i/n
        x=(1-t)**3*P0[0]+3*(1-t)**2*t*P1[0]+3*(1-t)*t**2*P2[0]+t**3*P3[0]
        y=(1-t)**3*P0[1]+3*(1-t)**2*t*P1[1]+3*(1-t)*t**2*P2[1]+t**3*P3[1]
        pts.append((x,y))
    return pts


def clip_to_node_boundary(points,node_center,radius=.56,from_start=True):
    q=points if from_start else list(reversed(points))
    def dist(t): return math.hypot(t[0]-node_center[0],t[1]-node_center[1])
    j=next(i for i,x in enumerate(q) if dist(x)>=radius)
    if j==0:
        cut=q[0]
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


def draw_polygon(ax,terminals,title,sigma,label_y=None,title_shift=0.54,resistor_shift=-0.20,title_size=10.1,sigma_size=15):
    pts=list(terminals.values())
    ax.add_patch(Polygon(pts,closed=True,facecolor=fill,edgecolor=ink,lw=1.6,zorder=4))
    for p in pts:
        ax.add_patch(Circle(p,0.12,facecolor='white',edgecolor=ink,lw=1.35,zorder=5))
    cx=sum(p[0] for p in pts)/len(pts); cy=sum(p[1] for p in pts)/len(pts)
    if title:
        ax.text(cx, cy+title_shift, title, color='white', fontsize=title_size, weight='bold', ha='center', va='center', zorder=6)
    zx=[cx-0.85,cx-0.68,cx-0.55,cx-0.40,cx-0.25,cx-0.10,cx+0.05,cx+0.20,cx+0.35,cx+0.50,cx+0.67,cx+0.85]
    zy=[cy+resistor_shift,cy+resistor_shift,cy+0.24+resistor_shift,cy-0.18+resistor_shift,cy+0.24+resistor_shift,cy-0.18+resistor_shift,
        cy+0.24+resistor_shift,cy-0.18+resistor_shift,cy+0.24+resistor_shift,cy-0.18+resistor_shift,cy+resistor_shift,cy+resistor_shift]
    ax.plot(zx,zy,color='white',lw=2.3,zorder=6)
    sig_y=label_y if label_y is not None else cy-0.56
    ax.text(cx,sig_y,sigma,ha='center',va='center',fontsize=sigma_size,color='white',zorder=6)


def draw_rotated_box_with_terminals(ax, center, ux, uy, length=1.28, height=0.50, terminal_r=0.11):
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
    offs=[-0.38,-0.28,-0.18,-0.08,0.02,0.12,0.22,0.32,0.42]
    zig=[]
    for j,o in enumerate(offs):
        amp=0 if j in [0,len(offs)-1] else (0.10 if j%2 else -0.10)
        zig.append((center[0]+ux*o+vx*amp, center[1]+uy*o+vy*amp))
    ax.plot([p[0] for p in zig],[p[1] for p in zig],color='white',lw=2.0,zorder=6)
    return left,right


def build_direct_geometry(spec):
    start_key,end_key=spec['nodes']
    p0=pos[start_key]; p3=pos[end_key]
    pts=bezier(p0,spec['ctrl1'],spec['ctrl2'],p3)
    pts=clip_to_node_boundary(pts,p0,from_start=True)
    pts=list(reversed(clip_to_node_boundary(pts,p3,from_start=False)))
    mid_idx=len(pts)//2
    mid=pts[mid_idx]
    before=pts[mid_idx-3]; after=pts[mid_idx+3]
    dx,dy=after[0]-before[0],after[1]-before[1]
    L=math.hypot(dx,dy)
    ux,uy=dx/L,dy/L
    left,right=draw_rotated_box_with_terminals(ax_dummy, mid, ux, uy, length=spec['box_len'])
    return None

# We compute geometry without drawing by duplicating the box math.
def direct_geometry(spec):
    start_key,end_key=spec['nodes']
    p0=pos[start_key]; p3=pos[end_key]
    pts=bezier(p0,spec['ctrl1'],spec['ctrl2'],p3)
    pts=clip_to_node_boundary(pts,p0,from_start=True)
    pts=list(reversed(clip_to_node_boundary(pts,p3,from_start=False)))
    mid_idx=len(pts)//2
    mid=pts[mid_idx]
    before=pts[mid_idx-3]; after=pts[mid_idx+3]
    dx,dy=after[0]-before[0],after[1]-before[1]
    L=math.hypot(dx,dy)
    ux,uy=dx/L,dy/L
    vx,vy=-uy,ux
    hl=spec['box_len']/2
    left=(mid[0]-ux*hl, mid[1]-uy*hl)
    right=(mid[0]+ux*hl, mid[1]+uy*hl)
    def nearest_index(target):
        return min(range(len(pts)), key=lambda i:(pts[i][0]-target[0])**2+(pts[i][1]-target[1])**2)
    i1=nearest_index(left); i2=nearest_index(right)
    if i1>i2: i1,i2=i2,i1
    seg1=pts[:i1+1]+[left]
    seg2=[right]+pts[i2:]
    return {'pts':pts,'mid':mid,'ux':ux,'uy':uy,'left':left,'right':right,'seg1':seg1,'seg2':seg2}


def draw_direct_structure(ax,study):
    spec=direct_specs[study]
    geom=direct_geometry(spec)
    ux,uy=geom['ux'],geom['uy']
    draw_rotated_box_with_terminals(ax, geom['mid'], ux, uy, length=spec['box_len'], height=0.50, terminal_r=0.11)
    draw_wire(ax, geom['seg1'][0], geom['seg1'][-1], lw=0)  # no-op to satisfy linter-like style
    ax.plot([z[0] for z in geom['seg1']],[z[1] for z in geom['seg1']],color=ink,lw=1.35,zorder=2)
    ax.plot([z[0] for z in geom['seg2']],[z[1] for z in geom['seg2']],color=ink,lw=1.35,zorder=2)
    ax.text(spec['label'][0],spec['label'][1],display_names[study],ha='center',va='center',fontsize=9,color=gray,
            bbox=dict(facecolor='white',edgecolor='none',boxstyle='round,pad=0.22',alpha=0.98),zorder=8)
    return geom


def arrow_on_segment(ax,p1,p2,color,width,z=7):
    arr=FancyArrowPatch(p1,p2,arrowstyle='-|>',mutation_scale=10+2.2*width,
                        lw=width,color=color,shrinkA=0,shrinkB=0,zorder=z)
    ax.add_patch(arr)
    return arr


def point_along_polyline(points, frac):
    d=[0.0]
    for i in range(1,len(points)):
        d.append(d[-1]+math.hypot(points[i][0]-points[i-1][0], points[i][1]-points[i-1][1]))
    total=d[-1]
    if total==0: return points[0]
    target=frac*total
    for i in range(1,len(points)):
        if d[i] >= target:
            rem=target-d[i-1]
            seg=d[i]-d[i-1]
            t=0 if seg==0 else rem/seg
            return (points[i-1][0]+t*(points[i][0]-points[i-1][0]), points[i-1][1]+t*(points[i][1]-points[i-1][1]))
    return points[-1]


def subsegment(points, f1, f2):
    if f2 < f1:
        rev=subsegment(list(reversed(points)), 1-f1, 1-f2)
        return rev
    p_start=point_along_polyline(points,f1)
    p_end=point_along_polyline(points,f2)
    d=[0.0]
    for i in range(1,len(points)):
        d.append(d[-1]+math.hypot(points[i][0]-points[i-1][0], points[i][1]-points[i-1][1]))
    total=d[-1]
    out=[p_start]
    for i,p in enumerate(points[1:-1], start=1):
        frac=d[i]/total if total>0 else 0
        if f1 < frac < f2:
            out.append(p)
    out.append(p_end)
    return out


def add_line_label(ax,p1,p2,text,color='black',offset=(0,0),fontsize=8.0):
    x=(p1[0]+p2[0])/2+offset[0]; y=(p1[1]+p2[1])/2+offset[1]
    ax.text(x,y,text,ha='center',va='center',fontsize=fontsize,color=color,
            bbox=dict(facecolor='white',edgecolor='none',pad=0.18,alpha=0.95),zorder=14)


def draw_multiterminal_arrows(ax,target,study,terminals):
    flows=currents[target][study]
    color=study_colors[study]
    for t,term in terminals.items():
        val=flows.get(t,0.0)
        if abs(val)<1e-4:
            continue
        w=1.2+10*abs(val)
        node_clip, term_clip = shorten_line(pos[t],term,0.56,0.12)
        if val>0:  # treatment -> study
            p1,p2=node_clip,term_clip
        else:
            p1,p2=term_clip,node_clip
        arrow_on_segment(ax,p1,p2,color,w,z=7)
        dx,dy=term[0]-pos[t][0], term[1]-pos[t][1]
        L=math.hypot(dx,dy)
        nx,ny=-dy/L,dx/L
        if abs(dy)<1e-9:  # horizontal wires (B and E terminals): label below the arrow
            nx,ny=0.0,-1.0
        mid=((p1[0]+p2[0])/2,(p1[1]+p2[1])/2)
        off=0.30+w/40.0
        if abs(dy)<1e-9:  # short horizontal wires: clear the arrowhead (half-width in data units)
            off=0.28+(10+2.2*w)/43.0
        lx,ly=mid[0]+off*nx, mid[1]+off*ny
        ax.text(lx,ly,f'{abs(val):.4f}',fontsize=8.1,color=color,ha='center',va='center',
                bbox=dict(facecolor='white',edgecolor='none',pad=0.16,alpha=0.95),zorder=14)


def draw_direct_arrows(ax,target,study,geom):
    spec=direct_specs[study]
    a,b=spec['nodes']
    va=currents[target][study].get(a,0.0)
    vb=currents[target][study].get(b,0.0)
    w=max(abs(va),abs(vb))
    if w<1e-4:
        return
    color=study_colors[study]
    if va>0 and vb<0:
        left_to_right=True
        label_side='left'
    elif vb>0 and va<0:
        left_to_right=False
        label_side='right'
    else:
        return
    width=1.2+10*w
    seg1=geom['seg1']
    seg2=geom['seg2']
    if left_to_right:
        sub1=subsegment(seg1,0.18,0.88)
        sub2=subsegment(seg2,0.12,0.82)
    else:
        sub1=subsegment(list(reversed(seg1)),0.18,0.88)
        sub2=subsegment(list(reversed(seg2)),0.12,0.82)
    arrow_on_segment(ax,sub1[0],sub1[-1],color,width,z=7)
    ax.plot([p[0] for p in sub1],[p[1] for p in sub1],alpha=0)  # keep polyline object count stable
    arrow_on_segment(ax,sub2[0],sub2[-1],color,width,z=7)
    # value label near the study label, fixed structure but changing value only
    lx,ly=spec['label']
    ax.text(lx,ly-0.44,f'{w:.4f}',ha='center',va='center',fontsize=8.2,color=color,
            bbox=dict(facecolor='white',edgecolor='none',pad=0.18,alpha=0.95),zorder=14)


def draw_source_sink(ax,source,sink):
    # External unit-current source/sink markers.  Treatment B always uses
    # a horizontal arrow placed to its LEFT, never a vertical arrow.
    # A:B  -> sink B: arrow points LEFT, out of B.
    # A:E  -> B is neither source nor sink: no external marker at B.
    # B:E  -> source B: arrow points RIGHT, into B.
    if source == 'A':
        sx,sy=pos['A']
        arrow_on_segment(ax,(sx,sy+1.45),(sx,sy+0.62),'black',2.6,z=6)
        ax.text(sx+0.30,sy+1.02,'+1',fontsize=10,weight='bold',color='black',zorder=14,
                bbox=dict(facecolor='white',edgecolor='none',pad=0.12))
    if source == 'B':
        bx,by=pos['B']
        arrow_on_segment(ax,(bx-1.45,by),(bx-0.62,by),'black',2.6,z=6)
        ax.text(bx-1.03,by+0.34,'+1',fontsize=10,weight='bold',color='black',zorder=14,
                ha='center',bbox=dict(facecolor='white',edgecolor='none',pad=0.12))
    if sink == 'B':
        bx,by=pos['B']
        arrow_on_segment(ax,(bx-0.62,by),(bx-1.45,by),'black',2.6,z=6)
        ax.text(bx-1.03,by+0.34,'-1',fontsize=10,weight='bold',color='black',zorder=14,
                ha='center',bbox=dict(facecolor='white',edgecolor='none',pad=0.12))
    if sink == 'E':
        ex,ey=pos['E']
        arrow_on_segment(ax,(ex+0.62,ey),(ex+1.45,ey),'black',2.6,z=6)
        ax.text(ex+1.03,ey+0.32,'-1',fontsize=10,weight='bold',color='black',zorder=14,
                ha='center',bbox=dict(facecolor='white',edgecolor='none',pad=0.12))

XLIM=(-0.9,15.9); YLIM=(-3.05,9.75)
def draw_panel(ax,target,source,sink):
    ax.set_axis_off(); ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.text(7.5,9.30,f'Target comparison {target}',ha='center',va='center',fontsize=11.5,weight='bold',color=ink)

    # fixed wires from treatments to multi-terminal elements
    for t,p in recovery.items():
        draw_wire(ax,*shorten_line(pos[t],p),lw=1.35,z=2)
    for t,p in solidarity.items():
        draw_wire(ax,*shorten_line(pos[t],p),lw=1.35,z=2)

    # fixed study elements
    rec_order={'A':recovery['A'],'D':recovery['D'],'C':recovery['C'],'B':recovery['B']}
    sol_order={'A':solidarity['A'],'E':solidarity['E'],'D':solidarity['D'],'C':solidarity['C']}
    draw_polygon(ax,rec_order,'S1',r'$\Sigma_1^{-1},label_y=1.55,title_shift=0.52,resistor_shift=-0.04,title_size=7.5,sigma_size=12)
    draw_polygon(ax,sol_order,'S2',r'$\Sigma_2^{-1},label_y=1.55,title_shift=0.52,resistor_shift=-0.04,title_size=7.5,sigma_size=12)
    geom_direct={study:draw_direct_structure(ax,study) for study in ['REMAP-CAP','CoDEX','ACTT-1']}

    # changing arrows and values only
    draw_multiterminal_arrows(ax,target,'RECOVERY',recovery)
    draw_multiterminal_arrows(ax,target,'SOLIDARITY',solidarity)
    for study in ['REMAP-CAP','CoDEX','ACTT-1']:
        draw_direct_arrows(ax,target,study,geom_direct[study])

    # treatment nodes on top
    for k in ['A','B','C','D','E']:
        add_node(ax,k)

    draw_source_sink(ax,source,sink)


# ------------------------------- key panel ----------------------------
def draw_key(ax):
    ax.set_axis_off(); ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.text(7.5,9.30,'How to read the panels',ha='center',va='center',fontsize=11.5,weight='bold',color=ink)
    y=7.95
    for study in ['RECOVERY','SOLIDARITY','ACTT-1','CoDEX','REMAP-CAP']:
        arrow_on_segment(ax,(0.2,y),(2.6,y),study_colors[study],3.2,z=7)
        ax.text(3.0,y,display_names[study],ha='left',va='center',fontsize=10.5,color=study_colors[study],weight='bold')
        y-=0.95
    ax.text(8.6,7.95,'Arrow direction: sign of the\nstudy-terminal current',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,6.15,'Arrow width and number:\nits magnitude',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,4.35,'+1 / \u22121: unit current injected\nat the source, extracted at the sink',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # opaque multi-terminal element
    demo={'A':(2.0,2.1),'D':(3.4,0.8),'C':(2.0,-0.5),'B':(0.6,0.8)}
    draw_polygon(ax,demo,'',r'$\Sigma_k^{-1}$',label_y=0.18,title_shift=0.5,resistor_shift=0.42,title_size=1,sigma_size=12)
    ax.text(4.3,0.8,'Multi-arm study: opaque element;\ncorners are its treatment terminals;\nno internal pairing is drawn',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # two-arm resistor
    l,r=draw_rotated_box_with_terminals(ax,(2.0,-1.9),1.0,0.0,length=1.30,height=0.50,terminal_r=0.11)
    ax.plot([0.2,l[0]-0.11],[-1.9,-1.9],color=ink,lw=1.35,zorder=2); ax.plot([r[0]+0.11,3.8],[-1.9,-1.9],color=ink,lw=1.35,zorder=2)
    ax.text(4.3,-1.9,'Two-arm study: resistor',ha='left',va='center',fontsize=9.8,color=ink)


fig=plt.figure(figsize=(10.2,8.2))
gs=fig.add_gridspec(2,2,wspace=0.02,hspace=0.04)
cells=[gs[0,0],gs[0,1],gs[1,0],gs[1,1]]
for cell,(target,src,snk) in zip(cells[:3],[('A:B','A','B'),('A:E','A','E'),('B:E','B','E')]):
    ax=fig.add_subplot(cell); draw_panel(ax,target,src,snk); ax.set_aspect('equal')
axk=fig.add_subplot(cells[3]); draw_key(axk); axk.set_aspect('equal')
fig.suptitle('Comparison-specific currents on the same fixed study-element circuit',fontsize=14,weight='bold',color=ink,y=0.985)
fig.subplots_adjust(left=0.01,right=0.99,bottom=0.01,top=0.95)
out_pdf=root/'COVID19_MultiTerminal_MultiComparison.pdf'
out_png=root/'COVID19_MultiTerminal_MultiComparison.png'
fig.savefig(out_pdf,bbox_inches='tight')
fig.savefig(out_png,dpi=220,bbox_inches='tight')
print(out_pdf)
,label_y=1.55,title_shift=0.52,resistor_shift=-0.04,title_size=7.5,sigma_size=12)
    draw_polygon(ax,sol_order,'SOLIDARITY',r'$\Sigma_S^{-1}$',label_y=1.55,title_shift=0.52,resistor_shift=-0.04,title_size=7.5,sigma_size=12)
    geom_direct={study:draw_direct_structure(ax,study) for study in ['REMAP-CAP','CoDEX','ACTT-1']}

    # changing arrows and values only
    draw_multiterminal_arrows(ax,target,'RECOVERY',recovery)
    draw_multiterminal_arrows(ax,target,'SOLIDARITY',solidarity)
    for study in ['REMAP-CAP','CoDEX','ACTT-1']:
        draw_direct_arrows(ax,target,study,geom_direct[study])

    # treatment nodes on top
    for k in ['A','B','C','D','E']:
        add_node(ax,k)

    draw_source_sink(ax,source,sink)


# ------------------------------- key panel ----------------------------
def draw_key(ax):
    ax.set_axis_off(); ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.text(7.5,9.30,'How to read the panels',ha='center',va='center',fontsize=11.5,weight='bold',color=ink)
    y=7.95
    for study in ['RECOVERY','SOLIDARITY','ACTT-1','CoDEX','REMAP-CAP']:
        arrow_on_segment(ax,(0.2,y),(2.6,y),study_colors[study],3.2,z=7)
        ax.text(3.0,y,study,ha='left',va='center',fontsize=10.5,color=study_colors[study],weight='bold')
        y-=0.95
    ax.text(8.6,7.95,'Arrow direction: sign of the\nstudy-terminal current',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,6.15,'Arrow width and number:\nits magnitude',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,4.35,'+1 / \u22121: unit current injected\nat the source, extracted at the sink',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # opaque multi-terminal element
    demo={'A':(2.0,2.1),'D':(3.4,0.8),'C':(2.0,-0.5),'B':(0.6,0.8)}
    draw_polygon(ax,demo,'',r'$\Sigma_k^{-1}$',label_y=0.18,title_shift=0.5,resistor_shift=0.42,title_size=1,sigma_size=12)
    ax.text(4.3,0.8,'Multi-arm study: opaque element;\ncorners are its treatment terminals;\nno internal pairing is drawn',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # two-arm resistor
    l,r=draw_rotated_box_with_terminals(ax,(2.0,-1.9),1.0,0.0,length=1.30,height=0.50,terminal_r=0.11)
    ax.plot([0.2,l[0]-0.11],[-1.9,-1.9],color=ink,lw=1.35,zorder=2); ax.plot([r[0]+0.11,3.8],[-1.9,-1.9],color=ink,lw=1.35,zorder=2)
    ax.text(4.3,-1.9,'Two-arm study: resistor',ha='left',va='center',fontsize=9.8,color=ink)


fig=plt.figure(figsize=(10.2,8.2))
gs=fig.add_gridspec(2,2,wspace=0.02,hspace=0.04)
cells=[gs[0,0],gs[0,1],gs[1,0],gs[1,1]]
for cell,(target,src,snk) in zip(cells[:3],[('A:B','A','B'),('A:E','A','E'),('B:E','B','E')]):
    ax=fig.add_subplot(cell); draw_panel(ax,target,src,snk); ax.set_aspect('equal')
axk=fig.add_subplot(cells[3]); draw_key(axk); axk.set_aspect('equal')
fig.suptitle('Comparison-specific currents on the same fixed study-element circuit',fontsize=14,weight='bold',color=ink,y=0.985)
fig.subplots_adjust(left=0.01,right=0.99,bottom=0.01,top=0.95)
out_pdf=root/'COVID19_MultiTerminal_MultiComparison.pdf'
out_png=root/'COVID19_MultiTerminal_MultiComparison.png'
fig.savefig(out_pdf,bbox_inches='tight')
fig.savefig(out_png,dpi=220,bbox_inches='tight')
print(out_pdf)
,label_y=1.55,title_shift=0.52,resistor_shift=-0.04,title_size=7.5,sigma_size=12)
    geom_direct={study:draw_direct_structure(ax,study) for study in ['REMAP-CAP','CoDEX','ACTT-1']}

    # changing arrows and values only
    draw_multiterminal_arrows(ax,target,'RECOVERY',recovery)
    draw_multiterminal_arrows(ax,target,'SOLIDARITY',solidarity)
    for study in ['REMAP-CAP','CoDEX','ACTT-1']:
        draw_direct_arrows(ax,target,study,geom_direct[study])

    # treatment nodes on top
    for k in ['A','B','C','D','E']:
        add_node(ax,k)

    draw_source_sink(ax,source,sink)


# ------------------------------- key panel ----------------------------
def draw_key(ax):
    ax.set_axis_off(); ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.text(7.5,9.30,'How to read the panels',ha='center',va='center',fontsize=11.5,weight='bold',color=ink)
    y=7.95
    for study in ['RECOVERY','SOLIDARITY','ACTT-1','CoDEX','REMAP-CAP']:
        arrow_on_segment(ax,(0.2,y),(2.6,y),study_colors[study],3.2,z=7)
        ax.text(3.0,y,study,ha='left',va='center',fontsize=10.5,color=study_colors[study],weight='bold')
        y-=0.95
    ax.text(8.6,7.95,'Arrow direction: sign of the\nstudy-terminal current',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,6.15,'Arrow width and number:\nits magnitude',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,4.35,'+1 / \u22121: unit current injected\nat the source, extracted at the sink',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # opaque multi-terminal element
    demo={'A':(2.0,2.1),'D':(3.4,0.8),'C':(2.0,-0.5),'B':(0.6,0.8)}
    draw_polygon(ax,demo,'',r'$\Sigma_k^{-1}$',label_y=0.18,title_shift=0.5,resistor_shift=0.42,title_size=1,sigma_size=12)
    ax.text(4.3,0.8,'Multi-arm study: opaque element;\ncorners are its treatment terminals;\nno internal pairing is drawn',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # two-arm resistor
    l,r=draw_rotated_box_with_terminals(ax,(2.0,-1.9),1.0,0.0,length=1.30,height=0.50,terminal_r=0.11)
    ax.plot([0.2,l[0]-0.11],[-1.9,-1.9],color=ink,lw=1.35,zorder=2); ax.plot([r[0]+0.11,3.8],[-1.9,-1.9],color=ink,lw=1.35,zorder=2)
    ax.text(4.3,-1.9,'Two-arm study: resistor',ha='left',va='center',fontsize=9.8,color=ink)


fig=plt.figure(figsize=(10.2,8.2))
gs=fig.add_gridspec(2,2,wspace=0.02,hspace=0.04)
cells=[gs[0,0],gs[0,1],gs[1,0],gs[1,1]]
for cell,(target,src,snk) in zip(cells[:3],[('A:B','A','B'),('A:E','A','E'),('B:E','B','E')]):
    ax=fig.add_subplot(cell); draw_panel(ax,target,src,snk); ax.set_aspect('equal')
axk=fig.add_subplot(cells[3]); draw_key(axk); axk.set_aspect('equal')
fig.suptitle('Comparison-specific currents on the same fixed study-element circuit',fontsize=14,weight='bold',color=ink,y=0.985)
fig.subplots_adjust(left=0.01,right=0.99,bottom=0.01,top=0.95)
out_pdf=root/'COVID19_MultiTerminal_MultiComparison.pdf'
out_png=root/'COVID19_MultiTerminal_MultiComparison.png'
fig.savefig(out_pdf,bbox_inches='tight')
fig.savefig(out_png,dpi=220,bbox_inches='tight')
print(out_pdf)
,label_y=1.55,title_shift=0.52,resistor_shift=-0.04,title_size=7.5,sigma_size=12)
    draw_polygon(ax,sol_order,'SOLIDARITY',r'$\Sigma_S^{-1}$',label_y=1.55,title_shift=0.52,resistor_shift=-0.04,title_size=7.5,sigma_size=12)
    geom_direct={study:draw_direct_structure(ax,study) for study in ['REMAP-CAP','CoDEX','ACTT-1']}

    # changing arrows and values only
    draw_multiterminal_arrows(ax,target,'RECOVERY',recovery)
    draw_multiterminal_arrows(ax,target,'SOLIDARITY',solidarity)
    for study in ['REMAP-CAP','CoDEX','ACTT-1']:
        draw_direct_arrows(ax,target,study,geom_direct[study])

    # treatment nodes on top
    for k in ['A','B','C','D','E']:
        add_node(ax,k)

    draw_source_sink(ax,source,sink)


# ------------------------------- key panel ----------------------------
def draw_key(ax):
    ax.set_axis_off(); ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.text(7.5,9.30,'How to read the panels',ha='center',va='center',fontsize=11.5,weight='bold',color=ink)
    y=7.95
    for study in ['RECOVERY','SOLIDARITY','ACTT-1','CoDEX','REMAP-CAP']:
        arrow_on_segment(ax,(0.2,y),(2.6,y),study_colors[study],3.2,z=7)
        ax.text(3.0,y,study,ha='left',va='center',fontsize=10.5,color=study_colors[study],weight='bold')
        y-=0.95
    ax.text(8.6,7.95,'Arrow direction: sign of the\nstudy-terminal current',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,6.15,'Arrow width and number:\nits magnitude',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    ax.text(8.6,4.35,'+1 / \u22121: unit current injected\nat the source, extracted at the sink',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # opaque multi-terminal element
    demo={'A':(2.0,2.1),'D':(3.4,0.8),'C':(2.0,-0.5),'B':(0.6,0.8)}
    draw_polygon(ax,demo,'',r'$\Sigma_k^{-1}$',label_y=0.18,title_shift=0.5,resistor_shift=0.42,title_size=1,sigma_size=12)
    ax.text(4.3,0.8,'Multi-arm study: opaque element;\ncorners are its treatment terminals;\nno internal pairing is drawn',ha='left',va='center',fontsize=9.8,color=ink,linespacing=1.25)
    # two-arm resistor
    l,r=draw_rotated_box_with_terminals(ax,(2.0,-1.9),1.0,0.0,length=1.30,height=0.50,terminal_r=0.11)
    ax.plot([0.2,l[0]-0.11],[-1.9,-1.9],color=ink,lw=1.35,zorder=2); ax.plot([r[0]+0.11,3.8],[-1.9,-1.9],color=ink,lw=1.35,zorder=2)
    ax.text(4.3,-1.9,'Two-arm study: resistor',ha='left',va='center',fontsize=9.8,color=ink)


fig=plt.figure(figsize=(10.2,8.2))
gs=fig.add_gridspec(2,2,wspace=0.02,hspace=0.04)
cells=[gs[0,0],gs[0,1],gs[1,0],gs[1,1]]
for cell,(target,src,snk) in zip(cells[:3],[('A:B','A','B'),('A:E','A','E'),('B:E','B','E')]):
    ax=fig.add_subplot(cell); draw_panel(ax,target,src,snk); ax.set_aspect('equal')
axk=fig.add_subplot(cells[3]); draw_key(axk); axk.set_aspect('equal')
fig.suptitle('Comparison-specific currents on the same fixed study-element circuit',fontsize=14,weight='bold',color=ink,y=0.985)
fig.subplots_adjust(left=0.01,right=0.99,bottom=0.01,top=0.95)
out_pdf=root/'COVID19_MultiTerminal_MultiComparison.pdf'
out_png=root/'COVID19_MultiTerminal_MultiComparison.png'
fig.savefig(out_pdf,bbox_inches='tight')
fig.savefig(out_png,dpi=220,bbox_inches='tight')
print(out_pdf)
