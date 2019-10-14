using Statistics;

function Mezei(ChemPot, V, T, R_Cut = 3.)
    """     CONFIGURATIONAL STEPS       """
    MC_Relaxation_Steps = 15000;
    MC_Equilibrium_Steps = 50000;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;
    MC_Measurement = 10;
    """     VARIABLE INITIALIZATION     """
    x, y, z = Float64[], Float64[], Float64[];
    L = V^(1. / 3.)
    Beta = 1. / T
    Pc_Random = Dict{Int64, Float64}();
    Pc_Random_Sum = Dict{Int64, Float64}();
    Pc_Random_N = Dict{Int64, Int64}();
    Displacement, N_Displacement, N_Displacement_Accepted = L / 8., 0, 0;
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
    Energy, Virial, N_Measurements = 0., 0., 0;
    Energy_Array, Pressure_Array, Density_Array = Float64[], Float64[], Float64[];
    N_Bins = 100;
    g_r = zeros(Float64, N_Bins);
    """     OUTPUT FILES        """
    Output_Route = pwd() * "/Output_Julia/ChemPot_$(round(ChemPot, digits = 2))_T_$(round(T, digits = 2))"
    mkpath("$Output_Route")
    Average_Energy_File = open("$Output_Route/Average_Energy.dat", "w");
    println(Average_Energy_File, "#\t< E / N >\n")
    Average_Pressure_File = open("$Output_Route/Average_Pressure.dat", "w");
    println(Average_Pressure_File, "#\t< P >\n")
    Average_Density_File = open("$Output_Route/Average_Density.dat", "w");
    println(Average_Density_File, "#\t< Density >\n")
    """     SIMULATIONS CYCLES      """
    for i = 1:MC_Steps
        """     PRINTS PROGRESS TO SCREEN   """
        if i < MC_Relaxation_Steps && i % .01MC_Relaxation_Steps == 0
            println("$(convert(Int64, 100i / MC_Relaxation_Steps))% Relaxation")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("P = $(round(length(x) * T - Virial / 3. / V))")
            println("N = $length(x)")
            println("Density = $(len(x) / V)")
            println("Movements: $N_Movement")
            println("   Accepted: $N_Movement_Accepted")
            println("   Rejected: $N_Movement_Rejected")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted")
            println("   Rejected: $N_Insertion_Rejected")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted")
            println("   Rejected: $N_Removal_Rejected")
            println("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
        end

        if i > MC_Relaxation_Steps && i % .01MC_Equilibrium_Steps == 0
            
        end
    end
end

@time Mezei(-3, 250, 2.)