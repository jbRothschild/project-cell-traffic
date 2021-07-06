import numpy as np
from src.agent import Bacteria
from src.default import BACT_PARAM

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
    pnt_vec_scaled = pnt_vec, 1.0/line_len
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
        self.height = height # height of bidirectional mother machine
        self.width = width # width of bidirectional mother machine
        self.bacteriaLst = [] # all the bacteria in the list
        self.beta = 0.8
        self.kn = 2.0 * 10**6 * 98.1 # Something like what they have
        self.kt = 0.0 # weird
        self.mun = 2.2 * 10**2 * 98.1
        self.mut = 2.2 * 10**2 * 98.1

    def setup_initial_nbr_cells( self, speciesLst ):
        """
        Change to a function factory!
        """
        id = 0
        for i, nbrSpecie in enumerate(speciesLst):
            for j in np.arange(0,nbrSpecie):
                p1 = np.array([np.random.uniform(0.0, self.width)
                        ,np.random.uniform(0.0, self.height)])
                vec = np.random.rand(2); vec /= np.linalg.norm(vec)
                p2 = p1 + ( BACT_PARAM[i]['maxLength']
                            * np.random.uniform(0.5,1.0) * vec )
                bactDict = {'type' : i, 'id' : str(id), 'p1' : p1
                                , 'p2' : p2, }
                bacteria = Bacteria(**bactDict)
                (self.bacteriaLst).append(bacteria)
                id += 1

    def closest_points_linesegments( self, cell1, cell2 ):
        """
        Finds the closest points between 2 line segments are in 2D.
        In 2D, one of the points has to be the

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

        dist[0], pnt[0,:], linepnt[0,:] =  pnt2line(cell1.p1, cell2.p1, cell2.p2)
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

    def forces( self, cell1, cell2 ):
        # TODO : Need to figure out how to do each of the forces. Need to have
        # forces between cells and

        cell1.add_force(force, cell1Pnt)
        cell2.add_force(-force, cell2Pnt)

    def step( self, dt ):

        for cell in self.bacteriaLst:
            s

        for i, cell in enumerate(self.bacteriaLst):
            cell.integrate_forces(dt, self)

        self.bacteriaLst = [x for x in self.bacteriaLst if not x.out(self)]

        return height, width, self.bacteriaLst

    def state( self ):
        """
        returns some state or quantitative structure of the whole
        simulation at that time. Might be healthier than saving all
        of the cells in the system.
        """

        return 0

    def plot_bacteria( self ):

        return 0
