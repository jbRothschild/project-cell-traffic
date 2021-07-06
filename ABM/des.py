from src import agent, environment, default

largeChannel = {'width' : 45.0, 'height' : 12.0}
smallChannel = {'width' : 45.0, 'height' : 6.0}

sim = environment.BiDirMM(**smallChannel)

sim.setup_initial_nbr_cells([2,3])

for bact in sim.bacteriaLst:
    print( bact.__dict__ )
