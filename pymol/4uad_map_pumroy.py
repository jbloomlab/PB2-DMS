###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '4uad'

cmd.reinitialize()
cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('antialias','2')
cmd.set('ray_opaque_background','off')
cmd.set('depth_cue', 'off')

### Modify here
url ='https://files.rcsb.org/download/4UAD.pdb'
cmd.load(url, 'orig')

cmd.select('importin', 'chain A')
cmd.select('origPB2', 'chain E')

cmd.load('S009PB2_{0}.pdb'.format(structure, 'S009'))
cmd.super('S009', 'origPB2')

cmd.hide('everything')
cmd.show('cartoon')
# cmd.show('surface')
for c in ['importin', 'S009']:
	cmd.show('surface', '{0}'.format(c))
# cmd.remove('origPB2')
cmd.color('wheat', 'importin')
cmd.color('white', 'S009')

cmd.orient()


### To get PB2 only
cmd.hide('everything')
# cmd.hide('everything', 'origPB2')
# cmd.hide('surface', 'S009')
cmd.show('cartoon', 'importin')
cmd.show('surface', 'S009')

###### END PREAMBLE


metric='Pumroy'
cmd.color('white', 'S009')

# Map sites identified by Pumroy 2015 as interacting with importin-alpha-3

for r in [752, 753, 755]: #major NLS-binding pocket
    cmd.color('green', 'resi {0} and S009'.format(r))
    cmd.show('spheres', 'resi {0} and S009'.format(r))

for r in [738, 739, 737, 718]: #minor NLS-binding site
    cmd.color('cyan', 'resi {0} and S009'.format(r))
    cmd.show('spheres', 'resi {0} and S009'.format(r))

###### SETVIEW

cmd.set_view ('\
    -0.632510662,   -0.737213850,   -0.237579703,\
     0.506035805,   -0.161099315,   -0.847335100,\
     0.586391509,   -0.656170785,    0.474955797,\
    -0.000148103,    0.000024185, -318.995605469,\
   -34.610126495,   13.929591179,   32.111713409,\
   251.499099731,  386.493469238,  -20.000000000 ')


###### POSTAMBLE

cmd.select(None)
cmd.set('specular', "off")
cmd.draw(width=1000, height=1000)
cmd.png('{0}_{1}.png'.format(structure, metric), ray=1)
cmd.rotate('y', angle=180)
cmd.draw(width=1000, height=1000)
cmd.png('{0}_{1}_y.png'.format(structure, metric), ray=1)
cmd.rotate('y', angle=180)
cmd.rotate('x', angle=90)
cmd.draw(width=1000, height=1000)
cmd.png('{0}_{1}_x.png'.format(structure, metric), ray=1)