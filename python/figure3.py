#!/usr/bin/env python3
"""Figure 3 a-e: signed anatomy panels and Yeo region-occupancy radar.

The radar uses recorded vROI values for 27 ordered regions, fixed component
colours and a percentage ROI-occupancy scale.
"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

COLORS={'Comp 1':'#0f4c5c','Comp 2':'#e36414','Comp 3':'#6a994e','Comp 4':'#7b2cbf'}
ORDER=list(COLORS)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel-dir',type=Path,required=True)
    p.add_argument('--yeo-csv',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--width-mm',type=float,default=180)
    p.add_argument('--dpi',type=int,default=600)
    a=p.parse_args()
    d=pd.read_csv(a.yeo_csv)
    required={'network_order','yeo_network','manuscript_label','saved_identifier','network_occupancy_percent'}
    if not required.issubset(d.columns): raise ValueError('Input must be yeo_radar_display.csv exported from vROI records')
    mapping=d[['manuscript_label','saved_identifier']].drop_duplicates().set_index('manuscript_label').saved_identifier.to_dict()
    assert mapping==dict(zip(ORDER,['COMP01','COMP02','COMP04','COMP06']))
    names=d[['network_order','yeo_network']].drop_duplicates().sort_values('network_order').yeo_network.tolist()
    m=d.pivot(index='yeo_network',columns='manuscript_label',values='network_occupancy_percent').loc[names,ORDER]
    if not np.isfinite(m.to_numpy()).all(): raise ValueError('Missing recorded radar values: refusing zero imputation')
    plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':7})
    fig=plt.figure(figsize=(a.width_mm/25.4,250/25.4),facecolor='white')
    for i in range(4):
        path=a.panel_dir/f'comp_{i+1}_source.png'
        if not path.is_file(): raise FileNotFoundError(path)
        ax=fig.add_axes([.045+i*.215,.77,.199,.19])
        ax.imshow(mpimg.imread(path),interpolation='none'); ax.axis('off')
        ax.set_title(f'{chr(97+i)}   Comp {i+1}',fontsize=9,fontweight='bold',loc='left',pad=6)
    cmap=LinearSegmentedColormap.from_list('signed',[(0,'#0000ff'),(.4,'#000000'),(.6,'#000000'),(1,'#ff0000')])
    cax=fig.add_axes([.923,.77,.011,.19])
    cb=fig.colorbar(ScalarMappable(norm=Normalize(-15,15),cmap=cmap),cax=cax,ticks=[-15,-10,-5,0,5,10,15])
    cb.set_label('Voxel weight',fontsize=7); cb.ax.tick_params(labelsize=6,length=2); cb.outline.set_visible(False)
    ax=fig.add_axes([.23,.205,.54,.3888],projection='polar')
    angles=np.linspace(0,2*np.pi,len(names),endpoint=False); closed=np.r_[angles,angles[0]]
    for comp in ORDER:
        row=m[comp].to_numpy(); row=np.r_[row,row[0]]
        ax.plot(closed,row,color=COLORS[comp],lw=1,label=comp)
        ax.fill(closed,row,color=COLORS[comp],alpha=.08)
    ax.set_theta_offset(np.pi/2); ax.set_theta_direction(-1)
    ax.set_xticks(angles); ax.set_xticklabels([])
    for angle,name in zip(angles,names):
        deg=90-np.degrees(angle); align='left'
        if deg < -90: deg+=180; align='right'
        ax.text(angle,47,name,fontsize=6.1,rotation=deg,rotation_mode='anchor',ha=align,va='center',clip_on=False)
    ax.set_ylim(0,45); ax.set_yticks(np.arange(5,46,5)); ax.tick_params(axis='y',labelsize=5.5)
    ax.set_rlabel_position(90); ax.spines['polar'].set_color('#bfbfbf')
    ax.xaxis.grid(True,color='#cccccc',lw=.4,alpha=.7); ax.yaxis.grid(True,color='#cccccc',lw=.4,alpha=.7)
    fig.text(.035,.715,'e',fontsize=10,fontweight='bold')
    handles,labels=ax.get_legend_handles_labels(); fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.055),ncol=4,frameon=False,fontsize=7)
    fig.text(.5,.033,'Yeo network occupancy (% of network ROI voxels)',ha='center',fontsize=7)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf','svg']: fig.savefig(a.output_dir/f'Figure_3_components_anatomy.{ext}',dpi=a.dpi,facecolor='white')
    fig.savefig(a.output_dir/'Figure_3_components_anatomy_preview.png',dpi=150)
    plt.close(fig)
    print(f'PASS: Figure 3 a-e; {len(names)} networks; recorded values; 180 x 250 mm; PDF native bitmap embedding')

if __name__=='__main__': main()
