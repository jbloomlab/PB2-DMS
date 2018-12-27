###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '4wsb'

cmd.reinitialize()
cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('antialias','2')
cmd.set('ray_opaque_background','off')
cmd.set('depth_cue', 'off')

### Modify here
url = 'https://files.rcsb.org/download/4WSB.pdb'
cmd.load(url, 'bat')
cmd.select('batPA', 'chain A')
cmd.select('batPB1', 'chain B')
cmd.select('batPB2', 'chain C')

cmd.load('S009PB2_{0}.pdb'.format(structure, 'S009'))
cmd.super('S009', 'bat')

cmd.hide('everything')
cmd.show('cartoon')
# cmd.show('surface')
for c in ['batPA', 'batPB1', 'S009']:
	cmd.show('surface', '{0}'.format(c))
cmd.remove('batPB2')
cmd.color('palegreen', 'batPA')
cmd.color('palecyan', 'batPB1')

cmd.orient()

### To get PB2 only
cmd.hide('everything', 'bat')
cmd.hide('surface', 'S009')
# cmd.show('cartoon', 'bat')
cmd.show('surface', 'S009')

###### END PREAMBLE

metric = 'topA549'
cmd.color('white', 'S009')

# knownAdaptive
for r in [256, 9, 526, 271, 660, 534, 535, 158, 286, 684, 701, 189, 702, 192, 199, 74, 714, 588, 591, 740, 489, 627, 636, 253]:
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
    -0.939411342,    0.224696293,    0.258875728,\
    -0.052751519,    0.651441693,   -0.756858826,\
    -0.338706851,   -0.724658728,   -0.600119412,\
     0.001423955,    0.000171021, -431.928985596,\
  1068.603637695,  -18.134885788,  381.554992676,\
   355.507598877,  508.442962646,  -20.000000000 ')

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