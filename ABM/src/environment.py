import sys
import numpy as np
from src.agent import Bacteria
from src.default import BACT_PARAM, BACT_PLOT
import matplotlib.pyplot as plt

def ccw(A,B,C):
    # check if these points are listed in counterclockwise order
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

def pnt2line(pnt, start, end):
    # Given a line with coordinates 'start' and 'end' and the
    # coordinates of a point 'pnt' the proc returns the shortest
    # distance from pnt to the line and the coordinates of the
    # nearest point on the line.
    #
    # 1  Convert the line segment to a vector ('line_vec').
    # 2  Create a vector connecting start to pnt ('pnt_vec').
    # 3  Find the length of the line vector ('line_len').
    # 4  Convert line_vec to a unit vector ('line_unitvec').
    # 5  Scale pnt_vec by line_len ('pnt_vec_scaled').
    # 6  Get the dot product of line_unitvec and pnt_vec_scaled ('t').
    # 7  Ensure t is in the range 0 to 1.
    # 8  Use t to get the nearest location on the line to the end
    #    of vector pnt_vec_scaled ('nearest').
    # 9  Calculate the distance from nearest to pnt_vec_scaled.
    # 10 Translate nearest back to the start/end line.
    # Malcolm Kesson 16 Dec 2012
    line_vec = start - end
    pnt_vec = start - pnt
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len
    pnt_vec_scaled = pnt_vec / line_len
    t = np.dot(line_unitvec, pnt_vec_scaled)
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    nearest = t * line_vec
    dist = np.linalg.norm(nearest - pnt_vec)
    nearest = nearest + start
    return dist, pnt, nearest

class BiDirMM():
    def __init__( self, height, width ):
        """
        Biomechanical ordering of dense cell populations, Tsimring et al.
        """
        self.height = height # height of bidirectional mother machine
        self.width = width # width of bidirectional mother machine
        self.bacteriaLst = [] # all the bacteria in the list
        self.beta = 0.8
        self.mass = 1.0
        #self.kn = 2.0 * 10**6 * 98.1 # Something like what they have
        self.kn = 2.2 * 10**2  # 10**7 N m**(-3/2) -> kg /( um min ** 2 )
        self.kt = 0.0 # weird
        self.gamman = 2.2 * 10**2
        self.gammat = 2.2 * 10**2

    def setup_initial_nbr_cells( self, speciesLst ):
        """
        Change to a function factory!

        Input
            speciesLst (list) : List with number of bacteria from each species
        """
        id = 0
        for i, nbrSpecie in enumerate(speciesLst):
            for j in np.arange(0,nbrSpecie):
                maxLen = BACT_PARAM[i]['maxLength']
                # set 1st extremity in the box
                p1 = np.array([np.random.uniform(0.0+maxLen, self.width-maxLen)
                        ,np.random.uniform(0.0+maxLen, self.height-maxLen)])
                # direction of the bacteria
                vec = np.random.rand(2); vec /= np.linalg.norm(vec)
                # set 2nd extermity
                p2 = p1 + ( maxLen * np.random.uniform(0.5,1.0) * vec )
                # create bacteria
                bactDict = {'type' : i, 'id' : str(id), 'p1' : p1
                                , 'p2' : p2, }
                bacteria = Bacteria(**bactDict)
                # add bacteria to simulation
                (self.bacteriaLst).append(bacteria)
                id += 1

    def place_init_cells( self ):
        """
        Place the initial cells
        """
        # cell 1
        # set the extremities in the box
        p1 = np.array([20., 6.0]); p2 = p1 + np.array([3.0,0.0]); id=0
        # create bacteria
        bactDict = {'type' : 0, 'id' : str(id), 'p1' : p1
                        , 'p2' : p2 }
        bacteria = Bacteria(**bactDict)
        # add bacteria to simulation
        (self.bacteriaLst).append(bacteria)

        # cell 2
        # set the extremities in the box
        p2 = np.array([24.25, 6.0]); p1 = p2 + np.array([3.0,0.0]); id=1
        # create bacteria
        bactDict = {'type' : 1, 'id' : str(id), 'p1' : p1
                        , 'p2' : p2 }
        bacteria = Bacteria(**bactDict)
        # add bacteria to simulation
        (self.bacteriaLst).append(bacteria)



    def add_cell(self, cellDict):
        (self.bacteriaLst).append(Bacteria(**cellDict))

    def closest_points_linesegments( self, cell1, cell2 ):
        """
        Finds the closest points between 2 line segments are in 2D.
        In 2D, one of the points has to be the endpoints

        Input
            cell1 : 1st bacteria
            cell2 : 2nd bacteria

        Output
            mindist : distances between points
            pntB    : point on cell1
            pntO    : point on cell2
        """
        # checks if they intersect, if so problem!
        if ( ccw(cell1.p1, cell2.p1, cell2.p2) != ccw(cell1.p2, cell2.p1, cell2.p2)
            and ccw(cell1.p1, cell1.p2, cell2.p1) != ccw(cell1.p1, cell1.p2, cell1.p2)):
            sys.exit("Warning: Cells " + cell1.id + " and " + cell2.id + " intersect!")

        # find shortest distance between the two cells
        dist = np.zeros(4); pnt = np.zeros((4,2)); linepnt = np.zeros((4,2))

        dist[0], pnt[0,:], linepnt[0,:] = pnt2line(cell1.p1, cell2.p1, cell2.p2)
        dist[1], pnt[1,:], linepnt[1,:] = pnt2line(cell1.p2, cell2.p1, cell2.p2)
        dist[2], pnt[2,:], linepnt[2,:] = pnt2line(cell2.p1, cell1.p1, cell1.p2)
        dist[3], pnt[3,:], linepnt[3,:] = pnt2line(cell2.p2, cell1.p1, cell1.p2)
        minidx = np.argmin(dist)

        # whether the point cell1 or cell2
        if minidx < 2:
            cell1Pnt = pnt[minidx,:]; cell2Pnt = linepnt[minidx,:]
        else:
            cell1Pnt = linepnt[minidx,:]; cell2Pnt = pnt[minidx,:]

        return dist[minidx], cell1Pnt, cell2Pnt

    def cell2cell_force( self, cell, othercell, cellPnt, othercellPnt, dist):
        # TODO : Need to figure out how to do each of the forces. Need to have
        # forces between cells
        M_e = self.mass / 2 # m/2
        k_n =  self.kn # (mg/d)
        gamma_n = self.gamman #(g/d)^{1/2}
        gamma_t = self.gammat
        force = np.array([0.0,0.0])
        mu_cc = 0.1
        delta = dist - (cell.radius + othercell.radius)/2.
        n_ij = (cellPnt-othercellPnt)/np.linalg.norm(cellPnt-othercellPnt)
        v_ij = ( ( cell.comVel + cell.angVel*np.array([cellPnt[1],-cellPnt[0]] )
                    ) - ( othercell.comVel + othercell.angVel *
                            np.array([cellPnt[1],-cellPnt[0]]) ) )
        v_n = np.dot(v_ij, n_ij)
        # calculate force
        F_n = k_n * delta**(3./2.) - gamma_n * M_e * delta * v_n
        v_t = v_ij - v_n * n_ij
        if np.linalg.norm(v_t) != 0.0:
            t_ij = v_t/np.linalg.norm(v_t)
        else:
            t_ij=np.array([0.0,0.0])
        F_t = - np.min([gamma_t * M_e * delta**(1./2.) * np.linalg.norm(v_t)
                            , mu_cc * F_n] )

        force = F_n * n_ij + F_t * t_ij

        cell.add_force(force, cellPnt)
        othercell.add_force(-force, othercellPnt)

    def cell2wall_force( self, cell, pval, wallPnt):
        if pval == 'p1':
            cellPnt = cell.p1
        elif pval =='p2':
            cellPnt = cell.p2
        else:
            print('Error! No cell point chosen.')
        M_e = self.mass # m/2
        k_n =  self.kn # (mg/d)
        gamma_n = self.gamman #(g/d)^{1/2}
        gamma_t = self.gammat
        force = np.array([0.0,0.0])
        mu_cw = 0.8

        delta = np.abs(cellPnt[1] - wallPnt) - (cell.radius)
        n_ij = (cellPnt - np.array([cellPnt[0],wallPnt]))/(cellPnt[1] - wallPnt)
        v_ij = ( cell.comVel + cell.angVel*np.array([cellPnt[1],-cellPnt[0]] ))
        v_n = np.dot(v_ij, n_ij)
        # calculate force
        F_n = k_n * delta**(3./2.) - gamma_n * M_e * delta * v_n
        v_t = v_ij - v_n * n_ij
        if np.linalg.norm(v_t) != 0.0:
            t_ij = v_t/np.linalg.norm(v_t)
        else:
            t_ij=np.array([0.0,0.0])
        F_t = - np.min([gamma_t * M_e * delta**(1./2.) * np.linalg.norm(v_t)
                            , mu_cw * F_n ])

        force = F_n * n_ij + F_t * t_ij

        cell.add_force(force, cellPnt)

        # set to true that it touches a wall
        cell.touches_wall()

    def step( self, dt ):

        # each cell needs to check interaction between it and all other cells
        # TODO : Must be a better way to do this.
        for i, cell in enumerate(self.bacteriaLst):
            # Check if cells touch other cell, if true do cell2cell_force()
            for othercell in self.bacteriaLst[i+1:]:
                dist, cellPnt, othercellPnt =\
                            self.closest_points_linesegments( cell, othercell )
                if dist <= cell.radius + othercell.radius:
                    self.cell2cell_force(cell, othercell, cellPnt, othercellPnt, dist)
            if cell.p1[1]-cell.radius < 0.0:
                self.cell2wall_force(cell, 'p1', 0.0)
            elif cell.p1[1]+cell.radius > self.height:
                self.cell2wall_force(cell, 'p1', self.height)
            else:
                pass
            if cell.p2[1]-cell.radius < 0.0:
                self.cell2wall_force(cell, 'p2', 0.0)
            elif cell.p2[1]+cell.radius > self.height:
                self.cell2wall_force(cell, 'p2', self.height)
            else:
                pass

        for i, cell in enumerate(self.bacteriaLst):
            #cell.add_force(np.array([1000,0]),cell.p1)
            cell.integrate_forces(dt, self)

        self.bacteriaLst = [x for x in self.bacteriaLst if not x.out(self)]

        return self.height, self.width, self.bacteriaLst

    def state( self ):
        """
        returns some state or quantitative structure of the whole
        simulation at that time. Might be healthier than saving all
        of the cells in the system.
        """

        return 0

    def plot_bacteria( self, filename=None, save=False ):

        fig = plt.figure(figsize=(self.width/2.,self.height/2.))
        for bacteria in self.bacteriaLst:
            plt.plot([bacteria.p1[0],bacteria.p2[0]],[bacteria.p1[1],bacteria.p2[1]]
                                                ,lw=30
                                                , solid_capstyle='round'
                                                , color=BACT_PLOT[bacteria.type])
        plt.ylim([0,self.height])
        plt.xlim([0,self.width])
        plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False) # labels along the bottom edge are off
        if filename != None:
            # save frame
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()

        return 0
