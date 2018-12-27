###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '5d98'

cmd.reinitialize()
cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('antialias','2')
cmd.set('ray_opaque_background','off')
cmd.set('depth_cue', 'off')

### Modify here
url ='https://files.rcsb.org/download/5D98.pdb'
cmd.load(url, 'orig')

cmd.select('origPA', 'chain A')
cmd.select('origPB1', 'chain B')
cmd.select('origPB2', 'chain C')
cmd.select('origPA_2', 'chain D')
cmd.select('origPB1_2', 'chain E')
cmd.select('origPB2_2', 'chain F')
cmd.remove('origPA_2')
cmd.remove('origPB1_2')
cmd.remove('origPB2_2')

cmd.load('S009PB2_{0}.pdb'.format(structure, 'S009'))
cmd.super('S009', 'orig')

cmd.hide('everything')
cmd.show('cartoon')
# cmd.show('surface')
for c in ['origPA', 'origPB1', 'S009']:
	cmd.show('surface', '{0}'.format(c))
cmd.remove('origPB2')
cmd.color('palegreen', 'origPA')
cmd.color('palecyan', 'origPB1')
cmd.color('white', 'S009')

cmd.orient()

### To get PB2 only
cmd.hide('everything', 'orig')
cmd.hide('surface', 'S009')
# cmd.show('cartoon', 'orig')
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
     0.636402488,    0.699312985,   -0.325504273,\
    -0.565354168,    0.135798618,   -0.813591778,\
    -0.524754405,    0.701796234,    0.481783271,\
     0.000000000,    0.000000000, -429.757751465,\
   -69.421409607,    7.238445282,   72.295677185,\
   338.824279785,  520.691223145,  -20.000000000 ')


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