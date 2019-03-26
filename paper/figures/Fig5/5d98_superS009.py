###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '5d98'

cmd.reinitialize()

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
cmd.alter('chain ""', 'chain="X"')
cmd.super('S009', 'orig')

cmd.remove('origPB2')

cmd.save('{0}_superS009.pdb'.format(structure))
cmd.save('{0}_superS009_S009only.pdb'.format(structure), 'chain X')

# Get dssp file here: http://www.cmbi.ru.nl/xssp/