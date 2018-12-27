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


### Color subunits
metric = 'fullstructure'
# Nter
for r in range(1, 247):
    cmd.color('oxygen', 'resi {0} and S009'.format(r))
# Mid-link
for r in range(247,319) + range(481,538):
    cmd.color('wheat', 'resi {0} and S009'.format(r))    
# Capbinding
for r in range(319, 481):
    cmd.color('selenium', 'resi {0} and S009'.format(r))
# 627
for r in range(538, 680):
    cmd.color('salmon', 'resi {0} and S009'.format(r))
# NLS
for r in range(680, 742):
    cmd.color('firebrick', 'resi {0} and S009'.format(r))




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