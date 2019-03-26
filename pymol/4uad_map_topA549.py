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


metric = 'topA549'
cmd.color('white', 'S009')

# knownAdaptive
for r in [256, 9, 526, 271, 660, 534, 535, 158, 286, 684, 701, 189, 702, 192, 199, 74, 714, 588, 591, 598, 740, 489, 627, 636, 253]:
	cmd.color('blue', 'resi {0} and S009'.format(r))
	cmd.show('spheres', 'resi {{0}} and S009'.format(r))

# topA549
for r in [521, 522, 9, 532, 156, 669, 158, 163, 292, 676, 169, 684, 176, 182, 183, 698, 701, 190, 69, 82, 355, 627]:
	cmd.color('red', 'resi {0} and S009'.format(r))
	cmd.show('spheres', 'resi {{0}} and S009'.format(r))

# knownAdaptiveAndTopA549
for r in [9, 684, 627, 701, 158]:
	cmd.color('magenta', 'resi {0} and S009'.format(r))
	cmd.show('spheres', 'resi {{0}} and S009'.format(r))



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