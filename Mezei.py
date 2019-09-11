##############################################################
#                                                            #
#       Mezei's algorithm for dense fluids from              #
#       Molecular Physics: An International Journal          #
#       at the Interface Between Chemistry and               #
#       Physics: A cavity-biased (T, V, μ) Monte Carlo       #
#       method for the computer simulation of                #
#       fluids. (1980)                                       #
#                                                            #
##############################################################

import random as rn
import math as m
import numpy as np
    
def Mezei(ChemPot, V, T, R_Cut = 3.0):
    """     AMOUNT OF STEPS     """
    MC_Relaxation_Steps = 100000
    MC_Equilibrium_Steps = 250000
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps
    MC_Measurement = 100

    x, y, z = [], [], []
    L = m.pow(V, 1. / 3.)
    Beta = 1. / T
    Displacement, N_Displacement, N_Accepted_Displacement = L / 10., 0, 0
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0
    Energy, N_Measurements = 0., 0
    Energy_Sum, N_Sum = 0., 0.
    for i in range(MC_Steps):
        
        """ PRINTS TO SCREEN SUMMARY """
        if i > 0 and i < MC_Relaxation_Steps and i % int(.01 * MC_Relaxation_Steps) == 0:
            print(int(100*i / MC_Relaxation_Steps), "% Relaxation")
            print("U / N = ", Energy/len(x))
            print("N = ", len(x))
            print("Density = ", len(x)/V)
            print("Max Displacement = ", Displacement)
            print("")

        if i > 0 and i > MC_Relaxation_Steps and i % int(0.01 * (MC_Steps - MC_Relaxation_Steps)) == 0:
            print(100*(i - MC_Relaxation_Steps) / (MC_Equilibrium_Steps), "% Equilibrium")
            print("U / N = ", Energy/len(x))
            print("N = ", len(x))
            print("Density = ", len(x)/V)
            print("Max Displacement = ", Displacement)
            print("")

        if i == MC_Relaxation_Steps:
            print("~~~  STARTING MEASUREMENT STEPS  ~~~")

        #print(i, "/", MC_Steps)

        if len(x) > 1:
            """            MOVEMENT           """
            N_Movement += 1
            N_Displacement += 1
            j = rn.randrange(0, len(x), 1)
            Energy_Old = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
            x_Old, y_Old, z_Old = x[j], y[j], z[j]

            """     DISPLACEMENT AND PERIODIC BOUNDARY CONDITIONS """
            x[j] += Displacement * (rn.random() - 0.5)
            if x[j] <  -L / 2.:
                x[j] += L
            elif x[j] > L / 2.:
                x[j] -= L
                
            y[j] += Displacement * (rn.random() - 0.5)
            if y[j] <  -L / 2.:
                y[j] += L
            elif y[j] > L / 2.:
                y[j] -= L

            z[j] += Displacement * (rn.random() - 0.5)
            if z[j] <  -L / 2.:
                z[j] += L
            elif z[j] > L / 2.:
                z[j] -= L
            
            Energy_New = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
            Delta_E = Energy_New - Energy_Old
            if rn.random() < m.exp(-Beta * Delta_E):
                N_Movement_Accepted += 1
                N_Accepted_Displacement += 1
                Energy += Delta_E
            else:
                N_Movement_Rejected
                x[j] = x_Old
                y[j] = y_Old
                z[j] = z_Old

        RN = rn.random()
        if RN < 0.5:
            """            INSERTION          """
            N_Insertion += 1
            x_Insertion = L*(rn.random() - 0.5)
            y_Insertion = L*(rn.random() - 0.5)
            z_Insertion = L*(rn.random() - 0.5)
            Energy_Insertion = Energy_Virial(L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
            if rn.random() < m.exp(Beta * (ChemPot - Energy_Insertion) + m.log(V / (len(x) + 1))):
                N_Insertion_Accepted += 1
                x.append(x_Insertion)
                y.append(y_Insertion)
                z.append(z_Insertion)
                Energy += Energy_Insertion
            else:
                N_Insertion_Rejected += 1

        elif len(x) > 1:
            """            REMOVAL            """
            N_Removal += 1
            j = rn.randrange(0, len(x), 1)
            Energy_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
            if rn.random() < m.exp(-Beta * (ChemPot - Energy_Removal) + m.log(len(x) / V)):
                N_Removal_Accepted += 1
                x.pop(j)
                y.pop(j)
                z.pop(j)
                Energy -= Energy_Removal
            else:
                N_Removal_Rejected += 1
        if i > 0:
            if m.fmod(i, MC_Measurement) == 0:
                """     COMPUTES AVERAGES    """
                N_Measurements += 1
                Energy_Sum += Energy/len(x)
                N_Sum += len(x)

                """ PARTICLE DISPLACEMENT """
                if 1. * N_Accepted_Displacement / N_Displacement > 0.55:
                    Displacement *= 1.05
                else:
                    Displacement *= 0.95
                N_Accepted_Displacement, N_Displacement = 0, 0
    
    print("< E / N > = ", Energy_Sum / N_Measurements)
    print("< Density > = ", N_Sum / (N_Measurements * V))
    print("< N > = ", N_Sum/N_Measurements)
    print("Movements: ", N_Movement)
    print("     Accepted: ", N_Movement_Accepted)
    print("     Rejected: ", N_Movement_Rejected)
                

def u(r2):
    return 4.0*(m.pow(1. / r2, 6.0) - m.pow(1. / r2, 3.))


def Energy_Virial(L, R_Cut, rx, ry, rz, x, y, z):
    Energy, Virial = 0., 0.
    for i in range(len(x)):
        Delta_x = rx - x[i]
        Delta_y = ry - y[i]
        Delta_z = rz - z[i]        
        """    PERIODIC BOUNDARY CONDITIONS    """
        if Delta_x > L / 2.0:
            Delta_x -= L
        elif Delta_x < -L / 2.0:
            Delta_x += L

        if Delta_y > L / 2.0:
                Delta_y -= L
        elif Delta_y < -L / 2.0:
            Delta_y += L

        if Delta_z > L / 2.0:
            Delta_z -= L
        elif Delta_z < -L / 2.0:
            Delta_z += L
        r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
        if r2 != 0.0:
            if r2 < m.pow(R_Cut, 2):
                Energy += u(r2)
    
    return Energy

Mezei(-3.893, 1500, 2.0)
