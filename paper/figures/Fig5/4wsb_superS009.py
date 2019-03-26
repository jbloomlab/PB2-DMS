###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '4wsb'

cmd.reinitialize()

### Modify here
url = 'https://files.rcsb.org/download/4WSB.pdb'
cmd.load(url, 'bat')
cmd.select('batPA', 'chain A')
cmd.select('batPB1', 'chain B')
cmd.select('batPB2', 'chain C')

cmd.load('S009PB2_{0}.pdb'.format(structure))
cmd.alter('chain ""', 'chain="X"')
cmd.super('S009', 'bat')

cmd.remove('batPB2')

cmd.save('{0}_superS009.pdb'.format(structure))
cmd.save('{0}_superS009_S009only.pdb'.format(structure), 'chain X')

# Get dssp file here: http://www.cmbi.ru.nl/xssp/