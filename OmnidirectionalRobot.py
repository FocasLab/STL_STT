#!/opt/homebrew/bin/python3.11
'''script for `Omnidirectional Robot` example'''
import z3
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class STT_Solver():
    '''class for generating STT based on constraints on trajectory'''
    def __init__(self, degree, dimension, time_step, min_tube_thickness, max_tube_thickness):
        self.setpoints = []
        self.obstacles = []
        self.goal = []
        self._start = 0
        self._finish = 0
        self._step = time_step
        self._range = 0
        self._x_start = 0
        self._x_finish = 0
        self._y_start = 0
        self._y_finish = 0
        self._z_start = 0
        self._z_finish = 0
        self.min_tube_thickness = min_tube_thickness
        self.max_tube_thickness = max_tube_thickness
        self.lambda_values = np.arange(0, 1.1, 0.1)
        self.degree = degree
        self.dimension = dimension
        self.solver = z3.Solver()
        z3.set_param("parallel.enable", True)
        self.C = [z3.Real(f'C{i}') for i in range((2 * self.dimension) * (self.degree + 1))]

    def gammas(self, t):
        '''method to calculate tube equations'''
        tubes = [z3.Real(f'g_{i}') for i in range(2 * self.dimension)]

        for i in range(2 * self.dimension):
            tubes[i] = 0
            power = 0
            for j in range(self.degree + 1):
                tubes[i] += ((self.C[j + i * (self.degree + 1)]) * (t ** power))
                power += 1
        return tubes

    def real_gammas(self, t, C_fin):
        '''method to calculate tube equations'''
        real_tubes = np.zeros(2 * self.dimension)

        for i in range(2 * self.dimension):
            power = 0
            for j in range(self.degree + 1): #each tube eq has {degree+1} terms
                real_tubes[i] += ((C_fin[j + i * (self.degree + 1)]) * (t ** power))
                power += 1
        return real_tubes

    def gamma_dot(self, t):
        '''method to calculate tube equations'''
        tubes = [z3.Real(f'gd_{i}') for i in range(2 * self.dimension)]

        for i in range(2 * self.dimension):
            tubes[i] = 0
            power = 0
            for j in range(self.degree + 1):
                if power < 1:
                    tubes[i] += 0
                    power += 1
                else:
                    tubes[i] += power * ((self.C[j + i * (self.degree + 1)]) * (t ** (power - 1)))
                    power += 1
        return tubes

    def real_gamma_dot(self, t, C_fin):
        '''method to calculate tube equations'''
        real_tubes = np.zeros(2 * self.dimension)

        for i in range(2 * self.dimension):
            power = 0
            for j in range(self.degree + 1):
                if power < 1:
                    real_tubes[i] += 0
                    power += 1
                else:
                    real_tubes[i] += power * ((C_fin[j + i * (self.degree + 1)]) * (t ** (power - 1)))
                    power += 1
        return real_tubes

    def general(self):
        '''method for general specifications'''
        temp_t_values = np.arange(self.getStart(), self.getFinish(), self._step)
        for t in temp_t_values:
            gamma1_L = self.gammas(t)[0]
            gamma2_L = self.gammas(t)[1]
            gamma1_U = self.gammas(t)[2]
            gamma2_U = self.gammas(t)[3]
            constraint_x = z3.And((gamma1_U - gamma1_L) > self.min_tube_thickness, (gamma1_U - gamma1_L) < self.max_tube_thickness)
            constraint_y = z3.And((gamma2_U - gamma2_L) > self.min_tube_thickness, (gamma2_U - gamma2_L) < self.max_tube_thickness)
            self.solver.add(constraint_x)
            self.solver.add(constraint_y)

    def plot_for_2D(self, C_fin):
        x_u = np.zeros(self.getRange())
        x_l = np.zeros(self.getRange())
        y_u = np.zeros(self.getRange())
        y_l = np.zeros(self.getRange())

        gd_xu = np.zeros(self.getRange())
        gd_xl = np.zeros(self.getRange())
        gd_yu = np.zeros(self.getRange())
        gd_yl = np.zeros(self.getRange())

        for i in range(self.getRange()):
            x_u[i] = self.real_gammas(i * self._step, C_fin)[2]
            x_l[i] = self.real_gammas(i * self._step, C_fin)[0]
            y_u[i] = self.real_gammas(i * self._step, C_fin)[3]
            y_l[i] = self.real_gammas(i * self._step, C_fin)[1]

            gd_xu[i] = self.real_gamma_dot(i * self._step, C_fin)[2]
            gd_xl[i] = self.real_gamma_dot(i * self._step, C_fin)[0]
            gd_yu[i] = self.real_gamma_dot(i * self._step, C_fin)[3]
            gd_yl[i] = self.real_gamma_dot(i * self._step, C_fin)[1]

        print("gamma_dot for x_upper max = ", gd_xu, max(gd_xu))
        print("gamma_dot for x_lower max = ", gd_xl, max(gd_xl))
        print("gamma_dot for y_upper max = ", gd_yu, max(gd_yu))
        print("gamma_dot for y_lower max = ", gd_yl, max(gd_yl))

        fig, axs = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True)
        ax, bx = axs
        for i in self.setpoints:        # t1    x1/y1   t2     t1   x2/y2  x1/y1
            square_x = patches.Rectangle((i[4], i[0]), i[5] - i[4], i[1] - i[0], edgecolor='green', facecolor='none')
            square_y = patches.Rectangle((i[4], i[2]), i[5] - i[4], i[3] - i[2], edgecolor='green', facecolor='none')
            ax.add_patch(square_x)
            bx.add_patch(square_y)

        for i in self.obstacles:        # t1    x1/y1   t2     t1   x2/y2  x1/y1
            square_x = patches.Rectangle((i[4], i[0]), i[5] - i[4], i[1] - i[0], edgecolor='red', facecolor='none')
            square_y = patches.Rectangle((i[4], i[2]), i[5] - i[4], i[3] - i[2], edgecolor='red', facecolor='none')
            ax.add_patch(square_x)
            bx.add_patch(square_y)

        t = np.linspace(self.getStart(), self.getFinish(), self.getRange())
        print("range: ", self.getRange(), "\nstart: ", self.getStart(), "\nfinish: ", self.getFinish(), "\nstep: ", self._step)

        ax.plot(t, x_u)
        ax.plot(t, x_l)
        bx.plot(t, y_u)
        bx.plot(t, y_l)

        fig2 = plt.figure(2)
        dx = fig2.add_subplot(111, projection='3d')
        dx.set_xlim(0, 20) ## dx.set_xlim(self.get_x_start(), self.get_x_finish())
        dx.set_ylim(0, 20) ## dx.set_ylim(self.get_y_start(), self.get_y_finish())
        dx.set_zlim(0, 20) ## dx.set_zlim(self.getStart(), self.getFinish())
        dx.set_xlabel('X Axis')
        dx.set_ylabel('Y Axis')
        dx.set_zlabel('Time Axis')

        for i in range(self.getRange()):
            vertices = [[x_u[i], y_u[i], i * self._step], [x_l[i], y_u[i], i * self._step], [x_l[i], y_l[i], i * self._step], [x_u[i], y_l[i], i * self._step],
                        [x_u[i], y_u[i], i * self._step], [x_l[i], y_u[i], i * self._step], [x_l[i], y_l[i], i * self._step], [x_u[i], y_l[i], i * self._step]]

            faces = [   [vertices[0], vertices[1], vertices[2], vertices[3]],  # Bottom face
                [vertices[4], vertices[5], vertices[6], vertices[7]],  # Top face
                [vertices[0], vertices[1], vertices[5], vertices[4]],  # Front face
                [vertices[2], vertices[3], vertices[7], vertices[6]],  # Back face
                [vertices[1], vertices[2], vertices[6], vertices[5]],  # Right face
                [vertices[0], vertices[3], vertices[7], vertices[4]]]  # Left face

            dx.add_collection3d(Poly3DCollection(faces, facecolors='blue', edgecolors='blue', alpha=0.25))

        for i in self.obstacles:
            dx.add_collection3d(Poly3DCollection(self.faces(i), facecolors='red', edgecolors='r', alpha=0.25))

        for i in self.setpoints:
            dx.add_collection3d(Poly3DCollection(self.faces(i), facecolors='green', edgecolors='green', alpha=0.25))

        fig3 = plt.figure(3)
        plt.plot(x_u, y_u, marker='o', linestyle='-')
        plt.plot(x_l, y_l, marker='o', linestyle='-')
        plt.title("Plot of array1 vs array2")
        plt.xlabel("X-Axis")
        plt.ylabel("Y-Axis")

    def find_solution(self):
        '''method to plot the tubes'''
        start = time.time()
        print("Solving...")

        self.setAll()
        self.general()

        if self.solver.check() == z3.sat:
            model = self.solver.model()
            xi = np.zeros((2 * self.dimension) * (self.degree + 1))
            Coeffs = []
            C_fin = np.zeros((2 * self.dimension) * (self.degree + 1))
            for i in range(len(self.C)):
                xi[i] = (np.float64(model[self.C[i]].numerator().as_long()))/(np.float64(model[self.C[i]].denominator().as_long()))
                print("{} = {}".format(self.C[i], xi[i]))
                Coeffs.append(xi[i])

            for i in range(len(Coeffs)):
                C_fin[i] = Coeffs[i]

            self.plot_for_2D(C_fin)
            self.print_equation(C_fin)

            end = time.time()
            self.displayTime(start, end)
            plt.show(block=True)

        else:
            print("No solution found.")
            print("range: ", self.getRange(), "\nstart: ", self.getStart(), "\nfinish: ", self.getFinish(), "\nstep: ", self._step)
            end = time.time()
            self.displayTime(start, end)

    def print_equation(self, C):
        for i in range(2 * self.dimension):
            print("gamma", i, "= ", end = "")
            power = 0
            for j in range(self.degree + 1):
                print("C", j + i * (self.degree + 1), "* t.^", power, "+ ", end = "")
                power += 1
            print("\n")

    def faces(self, i):
        vertices = [[i[0], i[2], i[4]], [i[1], i[2], i[4]], [i[1], i[3], i[4]], [i[0], i[3], i[4]],  # Bottom face
                    [i[0], i[2], i[5]], [i[1], i[2], i[5]], [i[1], i[3], i[5]], [i[0], i[3], i[5]]]   # Top face

        # Define the 6 faces of the cube using the vertices
        faces = [   [vertices[0], vertices[1], vertices[2], vertices[3]],  # Bottom face
                    [vertices[4], vertices[5], vertices[6], vertices[7]],  # Top face
                    [vertices[0], vertices[1], vertices[5], vertices[4]],  # Front face
                    [vertices[2], vertices[3], vertices[7], vertices[6]],  # Back face
                    [vertices[1], vertices[2], vertices[6], vertices[5]],  # Right face
                    [vertices[0], vertices[3], vertices[7], vertices[4]]]  # Left face
        return faces

    def setAll(self):
        all_points = self.setpoints + self.obstacles
        x1, x2, y1, y2, t1, t2 = [], [], [], [], [], []
        for i in all_points:
            x1.append(i[0])
            x2.append(i[1])
            y1.append(i[2])
            y2.append(i[3])
            t1.append(i[4])
            t2.append(i[5])

        self.setStart(min(t1))
        self.setFinish(max(t2))
        self.set_x_start(min(x1))
        self.set_x_finish(max(x2))
        self.set_y_start(min(y1))
        self.set_y_finish(max(y2))
        self.setRange(int((self.getFinish() - self.getStart() + self._step) / self._step))

    def displayTime(self, start, end):
        k = int(end - start)
        days = k // (3600 * 24)
        hrs = (k // 3600) - (days * 24)
        mins = (k // 60) - (hrs * 60)
        if end - start < 1:
            secs = (((end - start) * 10000) // 100) / 100
        else:
            secs = k - (mins * 60) - (hrs * 3600) - (days * 24 * 3600)
        print("Time taken: ", days, "days, ", hrs , "hours, ", mins, "minutes, ", secs, "seconds")

    def getStart(self):
        return self._start

    def setStart(self, start_value):
        self._start = start_value

    def getFinish(self):
        return self._finish

    def setFinish(self, finish_value):
        self._finish = finish_value

    def getRange(self):
        return self._range

    def setRange(self, value):
        self._range = value

    def get_x_start(self):
        return self._x_start

    def set_x_start(self, value):
        self._x_start = value

    def get_x_finish(self):
        return self._x_finish

    def set_x_finish(self, value):
        self._x_finish = value

    def get_y_start(self):
        return self._y_start

    def set_y_start(self, value):
        self._y_start = value

    def get_y_finish(self):
        return self._y_finish

    def set_y_finish(self, value):
        self._y_finish = value

    def get_z_start(self):
        return self._z_start

    def set_z_start(self, value):
        self._z_start = value

    def get_z_finish(self):
        return self._z_finish

    def set_z_finish(self, value):
        self._z_finish = value

solver = STT_Solver(degree=9, dimension=2, time_step=0.1, min_tube_thickness=0.2, max_tube_thickness=0.5)

def reach(x1, x2, y1, y2, t1, t2):
    solver.setpoints.append([x1, x2, y1, y2, t1, t2])
    all_constraints = []
    t_values = np.arange(t1, t2, solver._step)
    lambda_low = 0
    lambda_high = 1

    for t in t_values:
        gamma = solver.gammas(t)
        gamma1_L = gamma[0]
        gamma2_L = gamma[1]
        gamma1_U = gamma[2]
        gamma2_U = gamma[3]

        x_low = (lambda_low * gamma1_L + (1 - lambda_low) * gamma1_U)
        y_low = (lambda_low * gamma2_L + (1 - lambda_low) * gamma2_U)
        constraint_low = z3.And(x_low<x2, x_low>x1, y_low<y2, y_low>y1)
        all_constraints.append(constraint_low)

        x_high = (lambda_high * gamma1_L + (1 - lambda_high) * gamma1_U)
        y_high = (lambda_high * gamma2_L + (1 - lambda_high) * gamma2_U)
        constraint_high = z3.And(x_high<x2, x_high>x1, y_high<y2, y_high>y1)
        all_constraints.append(constraint_high)

    print("Added Reach Constraints: ", solver.setpoints)
    end = time.time()
    solver.displayTime(start, end)
    return all_constraints

def avoid(x1, x2, y1, y2, t1, t2):
    solver.obstacles.append([x1, x2, y1, y2, t1, t2])
    all_constraints = []
    t_values = np.arange(t1, t2, solver._step)
    lambda_low = 0
    lambda_high = 1

    for t in t_values:
        gamma = solver.gammas(t)
        gamma1_L = gamma[0]
        gamma2_L = gamma[1]
        gamma1_U = gamma[2]
        gamma2_U = gamma[3]

        x_low = (lambda_low * gamma1_L + (1 - lambda_low) * gamma1_U)
        y_low = (lambda_low * gamma2_L + (1 - lambda_low) * gamma2_U)
        constraint_low = z3.Or(z3.Or(x_low>x2, x_low<x1), z3.Or(y_low>y2, y_low<y1))
        all_constraints.append(constraint_low)

        x_high = (lambda_high * gamma1_L + (1 - lambda_high) * gamma1_U)
        y_high = (lambda_high * gamma2_L + (1 - lambda_high) * gamma2_U)
        constraint_high = z3.Or(z3.Or(x_high>x2, x_high<x1), z3.Or(y_high>y2, y_high<y1))
        all_constraints.append(constraint_high)

    print("Added Avoid Constraints: ", solver.obstacles)
    end = time.time()
    solver.displayTime(start, end)
    return all_constraints

#------------------ SPECIFICATION -------------------#
#         x1  x2    y1  y2        time
# reach {(4 , 5)  ,(4 , 5) }    ] 0-1s
# reach {(8 , 9)  ,(14, 15)}    ] 3-4s
# reach {(14, 15) ,(8 , 9) }    ] 6-7s
# reach {(8 , 9)  ,(14, 15)}    ] 9-10s
# reach {(14, 15) ,(8 , 9) }    ] 12-13s
# reach {(18, 19) ,(18, 19)}    ] 17-19s
# avoid {(5 , 7)  ,(8 , 10)}    ] 0-20s
# avoid {(15 , 17),(12 ,14)}    ] 0-20s
#----------------------------------------------------#

start = time.time()

S_constraints_list = reach(4, 5, 4, 5, 0, 1)
T1_constraints_list = reach(8, 9, 14, 15, 3, 4)
T2_constraints_list = reach(14, 15, 8, 9, 6, 7)
T3_constraints_list = reach(8, 9, 14, 15, 9, 10)
T4_constraints_list = reach(14, 15, 8, 9, 12, 13)
G_constraints_list = reach(18, 19, 18, 19, 17, 18)
O1_constraints_list = avoid(5, 7, 8, 10, 0, 20)
O2_constraints_list = avoid(15, 17, 12, 14, 0, 20)

for S in S_constraints_list:
    solver.solver.add(S)

for T1 in T1_constraints_list:
    solver.solver.add(T1)

for T2 in T2_constraints_list:
    solver.solver.add(T2)

for T3 in T3_constraints_list:
    solver.solver.add(T3)

for T4 in T4_constraints_list:
    solver.solver.add(T4)

for G in G_constraints_list:
    solver.solver.add(G)

for O1 in O1_constraints_list:
    solver.solver.add(O1)

for O2 in O2_constraints_list:
    solver.solver.add(O2)

solver.find_solution()

# around 20 mins