#!/usr/bin/env python3
"""Render Part D's saved final configurations using OVITO (no simulation rerun).
Run from any directory: QT_QPA_PLATFORM=offscreen python scripts/render_d_ovito.py
Requires ovito 3.14.1 (the rendering version used for this report). C/O labels encode A/B bead types, not chemical species.
"""
from pathlib import Path
import json
import numpy as np
import ovito
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.modifiers import PythonScriptModifier
from ovito.vis import Viewport, OSPRayRenderer
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'
CASES=[('d1','D1_snapshot.png',True),('d_n3_a28','D2_snapshot_N3.png',False),('d_n5_a28','D2_snapshot_N5.png',False),('d_n10_a28','D2_snapshot_N10.png',False)]

def add_caption(path,manifest,version):
    image=Image.open(path).convert('RGB')
    canvas=Image.new('RGB',(image.width,image.height+95),'white')
    canvas.paste(image,(0,0));draw=ImageDraw.Draw(canvas)
    try:font=ImageFont.truetype('DejaVuSans.ttf',19)
    except OSError:font=ImageFont.load_default(size=19)
    text=f'N = {manifest["chain_length"]} | aAB = {manifest["aab"]:g} | t = {manifest["steps"]*manifest["dt"]:g} | OVITO {version} / OSPRay'
    draw.text((25,image.height+5),text,fill='black',font=font)
    draw.text((25,image.height+42),'A beads',fill=(37,119,187),font=font)
    draw.text((180,image.height+42),'B beads',fill=(223,108,49),font=font)
    canvas.save(path)

def render(tag,filename,show_bonds):
    manifest=json.loads((DATA/f'{tag}_manifest.json').read_text())
    L=np.array(manifest['L'],float)
    pipeline=import_file(str(DATA/f'{tag}.pdb'),multiple_frames=True)
    frame=pipeline.source.num_frames-1
    def appearance(frame,data):
        assert data.particles.count==manifest['num_part']
        type_property=data.particles_.particle_types_
        for original_type in type_property.types:
            particle_type=type_property.make_mutable(original_type)
            particle_type.color=(0.145,0.467,0.733) if particle_type.name=='C' else (0.875,0.424,0.192)
            particle_type.radius=0.16 if show_bonds else 0.21
        if show_bonds:
            n=manifest['chain_length']
            starts=np.array([c*n+j for c in range(data.particles.count//n) for j in range(n-1)])
            topology=np.column_stack((starts,starts+1))
            positions=np.asarray(data.particles.positions)
            shift=-np.rint((positions[starts+1]-positions[starts])/L).astype(int)
            bonds=data.particles_.create_bonds(count=len(topology))
            bonds.create_property('Topology',data=topology)
            bonds.create_property('Periodic Image',data=shift)
            bonds.vis.width=0.065
            bonds.vis.color=(0.18,0.18,0.18)
    pipeline.modifiers.append(PythonScriptModifier(function=appearance))
    data=pipeline.compute(frame)
    pipeline.add_to_scene()
    vp=Viewport(type=Viewport.Type.Ortho)
    center=L/2
    direction=np.array((-0.55,1.0,-0.6))
    vp.camera_pos=center-direction*35
    vp.camera_dir=direction
    vp.zoom_all(size=(1000,720))
    vp.render_image(filename=str(DATA/filename),size=(1000,720),frame=frame,
                    background=(1,1,1),renderer=OSPRayRenderer(samples_per_pixel=4,refinement_iterations=1))
    add_caption(DATA/filename, manifest, ovito.version_string)
    pipeline.remove_from_scene()
    result=dict(tag=tag,file=filename,renderer='OVITO OSPRay',ovito_version=ovito.version_string,
                source=f'{tag}.pdb',frame=frame,frames=pipeline.source.num_frames,
                beads=data.particles.count,chain_length=manifest['chain_length'],aAB=manifest['aab'],
                time=manifest['steps']*manifest['dt'],box=list(L),bonds_shown=show_bonds,
                A_color='blue',B_color='orange')
    print(json.dumps(result))
    return result

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--case',type=int)
    args=parser.parse_args()
    if args.case is not None:
        result=render(*CASES[args.case])
        (DATA/f'D_ovito_case_{args.case}.json').write_text(json.dumps(result,indent=2)+'\n')
    else:
        import subprocess,sys
        for i in range(len(CASES)):
            subprocess.run([sys.executable,__file__,'--case',str(i)],check=True)
        results=[json.loads((DATA/f'D_ovito_case_{i}.json').read_text()) for i in range(len(CASES))]
        (DATA/'D_ovito_render_manifest.json').write_text(json.dumps(results,indent=2)+'\n')
        for i in range(len(CASES)):(DATA/f'D_ovito_case_{i}.json').unlink()
