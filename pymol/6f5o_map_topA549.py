###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '6f5o'

cmd.reinitialize()
cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('antialias','2')
cmd.set('ray_opaque_background','off')
cmd.set('depth_cue', 'off')

### Modify here
url ='https://files.rcsb.org/download/6F5O.pdb'
cmd.load(url, 'orig')

cmd.select('origPA', 'chain A')
cmd.select('origPB1', 'chain B')
cmd.select('origPB2', 'chain C')

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

# cmd.set_view ('\
#      0.299299747,   -0.468806446,   -0.831047535,\
#      0.759561479,   -0.410073251,    0.504881859,\
#     -0.577481806,   -0.782342494,    0.233352542,\
#      0.000000000,    0.000000000, -405.080596924,\
#    108.490371704,  107.187545776,  105.390640259,\
#    319.368621826,  490.792572021,  -20.000000000 ')

# Zoomed to region of CTD interaction
cmd.set_view ('\
    -0.324992836,    0.298390329,    0.897407830,\
     0.532062650,   -0.726806164,    0.434348613,\
     0.781849146,    0.618635833,    0.077445090,\
     0.000000000,    0.000000000, -170.138488770,\
   108.915313721,   88.171340942,  145.125610352,\
   134.138488770,  206.138488770,  -20.000000000 ')


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