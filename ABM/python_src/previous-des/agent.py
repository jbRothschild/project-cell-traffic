import numpy as np
from src.default import BACT_PARAM


class Bacteria():
    def __init__(self, type, id, p1, p2, comVel=np.array([0.0, 0.0]), angVel=0.0):
        self.id = id  # some way to identify it and it's progeny '1'
        self.p1 = p1  # start of cell
        self.p2 = p2  # end of cell
        self.type = type  # species type of cell
        # distributed randomly on creation of cell
        self.radius = BACT_PARAM[self.type]['radius']
        self.splitLength = np.random.normal(BACT_PARAM[self.type]['maxLength'], 0.05 * BACT_PARAM[self.type]['maxLength'])
        self.comVel = comVel  # center of mass velocity
        self.angVel = angVel  # angular velocity
        self.forces = []  # list of forces on cell
        self.pnts = []  # list of points at which the forces are applied
        self.center_length()

    def center_length(self):
        # automatically sets length and center of cell
        self.center = (self.p1 + self.p2) / 2.0
        self.length = np.linalg.norm(self.p1 - self.p2)

    def out(self, env):
        # return True if it's out of the simualtion
        if (self.center[0] < 0.0 or self.center[0] > env.width):
            return True
        else:
            return False

    def splitting_velocities(self, newp1, newp2):
        # new CoM from newp1 and newp2
        newCenter = (newp1 + newp2) / 2.0

        # angVel3d = w = (angVel, 0, 0)
        angVel3d = np.pad(np.array([self.angVel]), [(0, 2)], mode='constant',
                          constant_values=0)
        # 3d vector from old COM to new COM: (0,x,y)
        vec2newCom3d = np.pad(newCenter - self.center, [(1, 0)], mode='constant', constant_values=0)

        # 3d vector from old CoM to new p1 (0,x1,y1)
        vec2newp13d = np.pad(newp1 - self.center, [(1, 0)], mode='constant', constant_values=0)

        # 3d vector from new CoM to new p1 (0,x2,y2)
        vecNewCom2newp13d = vec2newp13d - vec2newCom3d

        # V_newcom = V_com + w x r_newcom = newComVel3d (0,vx,vy)
        newComVel3d = np.pad(self.comVel, [(1, 0)], mode='constant', constant_values=0) + np.cross(angVel3d, vec2newCom3d)

        # V_newp1 = V_com + w x r = newp1Vel3d (0,vx1,vy1)
        newp1Vel3d = np.pad(self.comVel, [(1, 0)], mode='constant', constant_values=0) + np.cross(angVel3d, vec2newp13d)

        # velocity newp1 wrt velocy of newCoM
        newp1relVelnewCom3d = newp1Vel3d - newComVel3d

        # w = (vxr)/r^2
        newAngVel3d = np.cross(vecNewCom2newp13d, newp1relVelnewCom3d)\
            / (np.linalg.norm(vecNewCom2newp13d)**2)

        return newComVel3d[1:], newAngVel3d[0]

    def split(self, variability=0.02):
        """
        Over a certain length, the cell will split into 2 cells. This splitting
        will be almost in half, however added a bit of noise so that it's not
        completely synchronous growth

        Input
            variability : variability in where the center can be, distributed
        """
        # Where the split happens
        splitPnt = self.p1 + (self.p2 - self.p1)*np.random.normal(1./2., variability)

        # daughter parameters. p2 from parent becomes p2 for daughter
        idDaughter = self.id + '1'
        p2toSplit = splitPnt - self.p2
        p2toSplitNorm = np.linalg.norm(p2toSplit)
        p1Daughter = self.p2 + ((p2toSplitNorm - self.radius) * p2toSplit
                                / p2toSplitNorm)
        comVelDaughter, angVelDaughter = self.splitting_velocities(p1Daughter, np.copy(self.p2))
        p2Daughter = np.copy(self.p2)

        # p1 from parent stay, have to calculate new p2
        p1toSplit = splitPnt - self.p1
        p1toSplitNorm = np.linalg.norm(p1toSplit)
        self.p2 = self.p1 + ((p1toSplitNorm - self.radius) * p1toSplit
                             / p1toSplitNorm)
        self.id += '0'
        self.comVel, self.angVel = self.splitting_velocities(
                                                            np.copy(self.p1), np.copy(self.p2))

        # set new center and length of mother
        self.center_length()

        # daughter parameters
        return {'type': int(np.copy(self.type)), 'id': idDaughter, 'p1': p1Daughter, 'p2': p2Daughter, 'comVel': comVelDaughter, 'angVel': angVelDaughter
                }

    def grow(self, dt, env):
        """
        exponential growth, can be changed to depend on multiple additional
        factors, such as pressure, concentration of cells etc.

        Input
            dt      : time interval
            env     : environment the cell is in
        """

        newLength = (self.length
                     * np.exp(BACT_PARAM[self.type]['growthRate']*dt))
        # if length larger than split length, splits into mother/daughter pair
        if newLength > self.splitLength:
            dictDaughter = self.split()
            env.add_cell(dictDaughter)
        # else grows from the center outwards
        else:
            self.p1 = self.center + (self.p1 - self.center) * newLength / self.length
            self.p2 = self.center + (self.p2 - self.center) * newLength / self.length
            self.length = newLength

    def add_force(self, force, pnt):
        """
        Adding a force to the list of forces

        Input
            force   : force to be added
            pnt     : where the force is applied
        """
        # TODO : Create forces
        (self.forces).append(force)
        (self.pnts).append(pnt)

    def integrate_forces(self, dt, env):
        """
        Integrate equations of motion and grow, split etc.

        Input
            dt      : time interval
            env     : environment the cell is in
        """
        acceleration = 0.0
        torqueAcc = 0.0
        if self.forces != []:
            # Forces on CoM
            # mass of spheres M = A * m where A = massReduced
            massReduced = ((2*self.length) / (np.pi * self.radius) + 1)
            # sum of forces =
            acceleration = sum(self.forces) / massReduced

            # Torque CoM
            inertia = 5  # m = 1
            radialVec = self.pnts - self.center
            torqueAcc = np.sum(np.cross(radialVec, self.forces)) / inertia

        self.comVel += dt * (acceleration - env.damping * self.comVel)
        self.angVel += dt * (torqueAcc - env.damping * self.angVel)

        # rotating the cell by an angle
        dtheta = self.angVel * dt
        rotationMatrix = np.array([[np.cos(dtheta), -np.sin(dtheta)], [np.sin(dtheta), np.cos(dtheta)]])

        # displacements
        self.center += dt * self.comVel
        self.p1 = self.center + rotationMatrix.dot(self.p1 - self.center)
        self.p2 = self.center + rotationMatrix.dot(self.p2 - self.center)
        self.center_length()
        self.forces = []
        self.pnts = []
