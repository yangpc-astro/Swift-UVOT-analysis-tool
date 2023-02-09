#!/usr/bin/env python
# -*- coding utf-8 -*-

'''
Using DS9 to open the image files and check if the source is in the image
'''
import glob
import pyds9
import astroquery
from astroquery.ned import Ned
from astropy.coordinates import SkyCoord
from astropy import units as u
'''
filepaths = glob.glob('/home/pengcheng/Work/data/J1753/00*/uvot/image/*sk.img')

for filepath in filepaths:
	radius = 5.
	zoom = 5
	d = pyds9.DS9()
	#Tell DS9 to open image file at 'filepath'.
	d.set('file %s' %filepath)
	#Tell DS9 to perform ‘zoom in’, 'log and minmax', 'heat'.
	d.set('zoom to fit')
	d.set('zoom %s' %zoom)
	d.set('scale mode minmax')
	d.set('scale log')
	d.set('cmap heat')
	#Tell DS9 to load region ,and center on the source position.
	d.set('regions', 'fk5;circle(268.36787386, -1.45173792,5")')
	d.set('regions', 'fk5;circle(268.3689554,-1.4413813,20")')
	d.set('pan to 268.36787386 -1.45173792 wcs fk5')
	x = input('viewing %s. Hit Enter to continue' %filepath)
'''
filepath = '/home/pengcheng/Work/data/J1753/00030090099/uvot/image/sw00030090099uw1_sk.img'
radius = 5.
zoom = 5
d = pyds9.DS9()
#Tell DS9 to open image file at 'filepath'.
d.set('file %s' %filepath)
#Tell DS9 to perform ‘zoom in’, 'log and minmax', 'heat'.
d.set('zoom to fit')
d.set('zoom %s' %zoom)
d.set('scale mode minmax')
d.set('scale log')
d.set('cmap heat')
#Tell DS9 to load region ,and center on the source position.
d.set('regions', 'fk5;circle(268.36787386, -1.45173792,5")')
d.set('regions', 'fk5;circle(268.3689554,-1.4413813,20")')
d.set('pan to 268.36787386 -1.45173792 wcs fk5')






















