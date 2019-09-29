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

import math as m
import numpy as np
import os

from random import random, randrange
from statistics import mean, pstdev
from timeit import default_timer as timer
           
def Mezei(ChemPot, V, T, R_Cut = 3.0):
    start = timer()
    """     AMOUNT OF STEPS     """
    MC_Relaxation_Steps = 5000
    MC_Equilibrium_Steps = 25000
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps
    MC_Measurement = 100

    x, y, z = [], [], []
    L = m.pow(V, 1. / 3.)
    Beta = 1. / T
    Pc, Pc_Sum, Pc_N = dict(), dict(), dict()
    Displacement, N_Displacement, N_Displacement_Accepted = L / 8, 0, 0
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0
    Energy, Virial = 0., 0.
    Energy_Array, Pressure_Array, Density_Array = [], [], []

    print("#" * 40 + "\n#" + " " * 38 + "#")
    print("#" + " " * 4 + "Monte Carlo + Excluded Volume" + " " * 5 + "#" + "\n#" + " " * 38 + "#")
    print("#" * 40)

    Output_Route = "Output/ChPot_%.3f_T_%.2f" % (ChemPot, T)
    if not os.path.exists(Output_Route):
        os.makedirs(Output_Route)
    Average_Energy_File = open("%s/Average_Energy.dat" % Output_Route, "w+")
    Average_Energy_File.write("#\t< E / N >\n")
    Average_Pressure_File = open("%s/Average_Pressure.dat" % Output_Route, "w+")
    Average_Pressure_File.write("#\t< P >\n")
    Average_Density_File = open("%s/Average_Density.dat" % Output_Route, "w+")
    Average_Density_File.write("#\t< Density >\n")
    
    for i in range(MC_Steps):    
        """ PRINTS TO SCREEN SUMMARY """
        if i > 0 and i < MC_Relaxation_Steps and i % int(.01 * MC_Relaxation_Steps) == 0:
            print("%d%% Relaxation" % (100*i / MC_Relaxation_Steps))
            print("U / N = %.6f" % (Energy/len(x)) )
            print("P = %.6f" % ((len(x) * T  - Virial / 3.0) / V) )
            print("N = %d" % len(x))
            print("Density = %.6f" % (len(x)/V) )
            print("Max Displacement = %.6f" % Displacement)
            print("Movements: %d" % N_Movement)
            print("     Accepted: %d" % N_Movement_Accepted)
            print("     Rejected: %d" % N_Movement_Rejected)
            print("Insertions: %d" % N_Insertion)
            print("     Accepted: %d" % N_Insertion_Accepted)
            print("     Rejected: %d" % N_Insertion_Rejected)
            print("Removal: %d" % N_Removal)
            print("     Accepted: %d" % N_Removal_Accepted)
            print("     Rejected: %d" % N_Removal_Rejected)
            print("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0

        if i > 0 and i > MC_Relaxation_Steps and i % int(0.01 * (MC_Steps - MC_Relaxation_Steps)) == 0:
            print("%d%% Equilibrium" % (100*(i - MC_Relaxation_Steps) / (MC_Equilibrium_Steps)))
            print("U / N = %.6f" % (Energy/len(x)) )
            print("P = %.6f" % ((len(x) * T  - Virial / 3.0) / V) )
            print("N = %d" % len(x))
            print("Density = %.6f" % (len(x)/V) )
            print("Max Displacement = %.6f" % Displacement)
            print("Movements: %d" % N_Movement)
            print("     Accepted: %d" % N_Movement_Accepted)
            print("     Rejected: %d" % N_Movement_Rejected)
            print("Insertions: %d" % N_Insertion)
            print("     Accepted: %d" % N_Insertion_Accepted)
            print("     Rejected: %d" % N_Insertion_Rejected)
            print("Removal: %d" % N_Removal)
            print("     Accepted: %d" % N_Removal_Accepted)
            print("     Rejected: %d" % N_Removal_Rejected)
            print("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0

        if i == MC_Relaxation_Steps:
            print("~~~  STARTING MEASUREMENT STEPS  ~~~")

        RN = randrange(1, 4)
        if RN == 1 and len(x) > 1:
            """            MOVEMENT           """
            N_Movement += 1
            N_Displacement += 1
            Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted = Movement(L, Beta, Displacement, Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)

        
        if RN == 2:
            """            INSERTION          """
            N_Insertion += 1
            Grid = 10
            Delta = L / (Grid + 1)
            EVMPS = ExcludedVolume(L, Grid, Delta, x, y, z)
            if np.sum(EVMPS) < m.pow(Grid, 3):
                Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected = Insertion_Mezei(EVMPS, Beta, ChemPot, Grid, Delta, L, V, R_Cut, Pc, Pc_Sum, Pc_N, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z)                
            else:        
                Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected = Insertion(L, V, ChemPot, Beta, R_Cut, x, y, z, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected)

        if RN == 3 and len(x) > 1:
            """            REMOVAL            """
            N_Removal += 1
            j = randrange(0, len(x), 1)
            """   IN CASE THE VALUE REQUIRED DOESN'T EXIST YET,
                    A LINEAR INTERPOLATION IS CARRIED ON USING THE TWO CLOSEST VALUES.  """
            if (len(x) - 1) not in Pc:
                Pc_Interpolation = Interpolation(Pc, len(x))
            else:
                Pc_Interpolation = Pc[len(x) - 1]
            if random() > m.pow(1 - m.pow(Pc_Interpolation, len(x) - 1), m.pow(Grid, 3)):
                Energy_Removal, Virial_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
                if random() < (len(x) / (V * Pc_Interpolation)) * m.exp(Beta * (Energy_Removal - ChemPot)):
                    N_Removal_Accepted += 1
                    x.pop(j)
                    y.pop(j)
                    z.pop(j)
                    Energy -= Energy_Removal
                    Virial -= Virial_Removal
                else:
                    N_Removal_Rejected += 1
            else: 
                Energy_Removal, Virial_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
                if random() < m.exp(Beta * (Energy_Removal - ChemPot) + m.log(len(x) / V)):
                    N_Removal_Accepted += 1
                    x.pop(j)
                    y.pop(j)
                    z.pop(j)
                    Energy -= Energy_Removal
                    Virial -= Virial_Removal
                else:
                    N_Removal_Rejected += 1

        if i > 0:
            if m.fmod(i, MC_Measurement) == 0:
                if i >= MC_Relaxation_Steps:
                    """     COMPUTES AVERAGES    """
                    Energy_Array.append( Energy/len(x) )
                    Average_Energy_File.write("%d\t%f\n" % (len(Energy_Array), mean(Energy_Array)))
                    Pressure_Array.append( (len(x) * T - Virial / 3.0) / V )
                    Average_Pressure_File.write("%d\t%f\n" % (len(Pressure_Array), mean(Pressure_Array)))
                    Density_Array.append( len(x) / V ) 
                    Average_Density_File.write("%d\t%f\n" % (len(Density_Array), mean(Density_Array) ))

                """ PARTICLE DISPLACEMENT """
                if 1. * N_Displacement_Accepted / N_Displacement > 0.55:
                    Displacement *= 1.05
                else:
                    Displacement *= 0.95
                if Displacement < 0.05:
                    Displacement = 0.05
                if Displacement > L / 4.:
                    Displacement = L /4. 
                N_Displacement, N_Displacement_Accepted = 0, 0
    
    Average_Energy_File.close()
    Average_Pressure_File.close()
    Average_Density_File.close()

    Pc_File = open("%s/Pc.dat" % Output_Route, "w+")
    for i in Pc:
        Pc_File.write("%d\t%f\n" % (i, Pc[i]))
    Pc_File.close()

    print("< E / N > = %.6f + %.6f" % (mean(Energy_Array), pstdev(Energy_Array)) )
    print("< P > = %.6f + %.6f" % (mean(Pressure_Array), pstdev(Pressure_Array)) )
    print("< Density > = %.6f + %.6f " % (mean(Density_Array), pstdev(Density_Array)) )
    print("< N > = %.6f + %.6f " % (V*mean(Density_Array), V*pstdev(Density_Array)) )
    dt = timer() - start
    print("Elapsed time: %f s" % dt)

    Summary_File = open("%s/Summary.txt" % Output_Route, "w+")
    Summary_File.write("Mezei algorithm for the Grand Canonical Monte Carlo Simulation\n\n")
    Summary_File.write("~~~~~   INPUT   ~~~~~\n\n")
    Summary_File.write("Chemical Potential = %.6f\n" % ChemPot)
    Summary_File.write("Volume = %.3f (L = %.6f)\n" % (V, L))
    Summary_File.write("Temperature = %.3f\n\n" % T)
    Summary_File.write("Relaxation Steps: %d.   Equilibrium Steps: %d\n\n" % (MC_Relaxation_Steps, MC_Equilibrium_Steps))
    Summary_File.write("~~~~~   OUTPUT  ~~~~~\n\n")
    Summary_File.write("< E / N > = %.6f + %.6f\n" % (mean(Energy_Array), pstdev(Energy_Array)) )
    Summary_File.write("< P > = %.6f + %.6f\n" % (mean(Pressure_Array), pstdev(Pressure_Array)) )
    Summary_File.write("< Density > = %.6f + %.6f\n" % (mean(Density_Array), pstdev(Density_Array)) )
    Summary_File.write("< N > = %.6f + %.6f\n" % (V*mean(Density_Array), V*pstdev(Density_Array)) )

#def u(r2):
#    if r2 <= 1:
#        return m.inf
#    if r2 > 1 and r2 <= 1.5:
#        return -1
#    if r2 > 1.5:
#        return 0

def Movement(L, Beta, Displacement, Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z):
    j = randrange(0, len(x), 1)
    Energy_Old, Virial_Old = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    x_Old, y_Old, z_Old = x[j], y[j], z[j]
    x[j] += Displacement * (random() - 0.5)                
    y[j] += Displacement * (random() - 0.5)
    z[j] += Displacement * (random() - 0.5)
    x[j] = PeriodicBoundaryConditions(L, x[j])
    y[j] = PeriodicBoundaryConditions(L, y[j])
    z[j] = PeriodicBoundaryConditions(L, z[j])
    Energy_New, Virial_New = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    Delta_E = Energy_New - Energy_Old
    Delta_Virial = Virial_New - Virial_Old
    if random() < m.exp(-Beta * Delta_E):
        N_Movement_Accepted += 1
        N_Displacement_Accepted += 1
        Energy += Delta_E
        Virial += Delta_Virial
    else:
        N_Movement_Rejected += 1
        x[j] = x_Old
        y[j] = y_Old
        z[j] = z_Old
    return Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted

def Insertion(L, V, ChemPot, Beta, R_Cut, x, y, z, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected):
    x_Insertion = L*(random() - 0.5)
    y_Insertion = L*(random() - 0.5)
    z_Insertion = L*(random() - 0.5)
    Energy_Insertion, Virial_Insertion = Energy_Virial(L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
    if random() < m.exp(Beta * (ChemPot - Energy_Insertion) + m.log(V / (len(x) + 1))):
        N_Insertion_Accepted += 1
        x.append(x_Insertion)
        y.append(y_Insertion)
        z.append(z_Insertion)
        Energy += Energy_Insertion
        Virial += Virial_Insertion
    else:
        N_Insertion_Rejected += 1
    return Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected

def Insertion_Mezei(EVMPS, Beta, ChemPot, Grid, Delta, L, V, R_Cut, Pc, Pc_Sum, Pc_N, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z):
    aux1 = np.where(EVMPS == 0)
    j = randrange(0, len(aux1[0]), 1)
    i_x = aux1[0][j]
    i_y = aux1[1][j]
    i_z = aux1[2][j]
    x_Insertion = (i_x + 1) * Delta - L / 2.
    y_Insertion = (i_y + 1) * Delta - L / 2.
    z_Insertion = (i_z + 1) * Delta - L / 2.
    Energy_Insertion, Virial_Insertion = Energy_Virial(L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
    """     PROBABILITY Pc      """
    if len(x) not in Pc_N:
        Pc_N[len(x)] = 1
        Pc_Sum[len(x)] = 1 - (np.sum(EVMPS) / m.pow(Grid, 3))
        Pc[len(x)] = 1 - (np.sum(EVMPS) / m.pow(Grid, 3))
    else:
        Pc_N[len(x)] += 1
        Pc_Sum[len(x)] += 1 - (np.sum(EVMPS) / m.pow(Grid, 3)) #N_i/N_t
        Pc[len(x)] = Pc_Sum[len(x)] / Pc_N[len(x)]
    if random() < (V * Pc[len(x)] / (len(x) + 1)) * m.exp(Beta * (ChemPot - Energy_Insertion)):
        N_Insertion_Accepted += 1
        x.append(x_Insertion)
        y.append(y_Insertion)
        z.append(z_Insertion)
        Energy += Energy_Insertion
        Virial += Virial_Insertion
    else:
        N_Insertion_Rejected += 1
    return Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected

def u(r2):
    return 4.0*(m.pow(1. / r2, 6.0) - m.pow(1. / r2, 3.))

def rdu(r2):
    return 4.0*(6*m.pow(1. / r2, 3.0) - 12.0 * m.pow(1.0 / r2, 6.0))

def Energy_Virial(L, R_Cut, rx, ry, rz, x, y, z):
    Energy, Virial = 0., 0.
    for i in range(len(x)):
        Delta_x = rx - x[i]
        Delta_y = ry - y[i]
        Delta_z = rz - z[i]        
        Delta_x = PeriodicBoundaryConditions(L, Delta_x)
        Delta_y = PeriodicBoundaryConditions(L, Delta_y)
        Delta_z = PeriodicBoundaryConditions(L, Delta_z)
        r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
        if r2 != 0.0:
            if r2 < m.pow(R_Cut, 2):
                Energy += u(r2)
                Virial += rdu(r2)
    return Energy, Virial

def PeriodicBoundaryConditions(L, x):
    if x <  -L / 2.:
        x += L
    elif x > L / 2.:
        x -= L
    
    return x

def ExcludedVolume(L, Grid, Delta, x, y, z):
    EVMPS = np.zeros((Grid, Grid, Grid))
    for i_x in range(0, Grid, 1):
        for i_y in range(0, Grid, 1):
            for i_z in range(0, Grid, 1):
                x_Grid = (i_x + 1) * Delta - L / 2.
                y_Grid = (i_y + 1) * Delta - L / 2.
                z_Grid = (i_z + 1) * Delta - L / 2.
                for k in range(0, len(x), 1):
                    Delta_x = x_Grid - x[k]
                    Delta_y = y_Grid - y[k]
                    Delta_z = z_Grid - z[k]
                    r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
                    if r2 < m.pow(0.6, 2):
                        EVMPS[i_x, i_y, i_z] = 1
                        break
    return EVMPS

                #"""     EXCLUDED VOLUME     """
                #V_Excluded = len(x) * (4. / 3.) * m.pi
                #""" CORRECTION DUE TO SPHERE OVERLAPNESS """
                #V_Excluded_Correction = 0
                #for j in range(0, len(x), 1):
                #    for k in range(0, len(x), 1):
                #        if j != k:
                #            Delta_x = x[j] - x[k]
                #            Delta_y = y[j] - y[k]
                #            Delta_z = z[j] - z[k]
                #            r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
                #            if r2 < 1.0: # r2 < 1 implies overlapness 
                #                d = m.sqrt(r2)
                #                V_Excluded_Correction += (m.pi / 12.) * (4 + d) * m.pow((2 - d), 2)
                #            else:
                #                if Delta_x > L - 1.: """  PERIODIC BOUNDARY CONDITION FOR THE OVERLAPNESS CONTRIBUTION  """
                #                    Delta_x -= L
                #                if Delta_y > L - 1.:
                #                    Delta_y -= L
                #                if Delta_z > L - 1.:
                #                    Delta_z -= L
                #                r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
                #                if r2 < 1.0: # r2 < 1 implies overlapness 
                #                    d = m.sqrt(r2)
                #                    V_Excluded_Correction += (m.pi / 12.) * (4 + d) * m.pow((2 - d), 2)
                #Volume_Ratio = (V_Excluded - V_Excluded_Correction) / V
                #if Volume_Ratio > 1 or Volume_Ratio < 0:
                #    raise ValueError("Volume Ratio (", Volume_Ratio, ") can't be negative nor greater than one.")
                #if random() < (V_Excluded - V_Excluded_Correction)/(len(x) + 1) * m.exp(Beta * (ChemPot - Energy_Insertion))
                #    N_Insertion_Accepted += 1
                #    x.append(x_Insertion)
                #    y.append(y_Insertion)
                #    z.append(z_Insertion)
                #    Energy += Energy_Insertion
                #else:
                #    N_Insertion_Rejected += 1

def Interpolaiton(Pc, l):
    if len(Pc) == 1:
        Pc_Interpolation =  Pc[ list( Pc.keys() )[0] ]
    elif len(Pc) > 1:
        lim_inf = l - 1
        while True:
            if lim_inf in Pc:
                break
            lim_inf -= 1

        lim_sup = l + 1
        while True:
            if lim_sup in Pc:
                break
            lim_sup += 1
    elif len(Pc) == 0:
        raise ValueError("Pc has no values still.")

    Pc_Interpolation = Pc[lim_inf] * ((1 - l - lim_inf) / (lim_sup - lim_inf)) + Pc[lim_sup] * ((l - lim_inf) / (lim_sup - lim_inf))
    return Pc_Interpolation

Mezei(17.9183589140279, 250.058699225, 4.0)
