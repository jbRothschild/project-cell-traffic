# growthRate    : in bulk growth rate. Forces aught to slow it down in sim.
# maxLength     : maximum length from center to center (p1 to p2)
# radius        : essentially half width, or radius of spherical tip
# all lengths in micrometers.

DATA_DIR = 'data'

BACT_PARAM = { 0 : { 'growthRate' : 0.013, 'maxLength' : 4.0, 'radius' : 0.5 },
                1 : { 'growthRate' : 0.013, 'maxLength' : 4.0 , 'radius' : 0.5},
                2 : { 'growthRate' : 0.013, 'maxLength' : 4.0 , 'radius' : 0.5},
                3 : { 'growthRate' : 0.006, 'maxLength' : 8.0 , 'radius' : 0.5}
                }

BACT_PLOT = { 0 : 'r', 1 : 'g', 2: 'b'}
