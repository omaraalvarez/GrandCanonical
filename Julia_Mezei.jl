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
    Energy_Array, Pressure_Array, Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
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
            println("P = $(round(length(x) * T - Virial / 3. / V, digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
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
            println("$(convert(Int64, 100(i - MC_Relaxation_Steps) / MC_Equilibrium_Steps))% Equilibrium ($N_Measurements Measurements).")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("P = $(round(length(x) * T - Virial / 3. / V, digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
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

        i == MC_Relaxation_Steps ? println("~~~    STARTING MEASUREMENT STEPS    ~~~") : nothing
        RN = rand(1:3);
        if RN == 1 && length(x) > 1
            N_Movement += 1;
            N_Displacement += 1;
            Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted = Movement(L, Beta, Displacement, Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)
        end
        if RN == 2
            N_Insertion += 1;
            Volume_Ratio, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(L, x, y, z)
            if haskey(length(x), Pc_Random)
                Pc_Random_N[length(x)] = 1
                Pc_Random_Sum[length(x)] = Volume_Ratio;
                Pc_Random[length(x)] = Volume_Ratio;
            else
                Pc_Random_N[length(x)] += 1;
                Pc_Random_Sum[length(x)] += Volume_Ratio;
                Pc_Random[length(x)] = Pc_Random_Sum[length(x)] / Pc_Random_N[length(x)]
            end
            if Volume_Ratio > 0
                Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected = Insertion_Mezei(Beta, ChemPot, L, V, R_Cut, Pc_Random, Pc_Random_Sum, Pc_Random_N, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z, x_Insertion, y_Insertion, z_Insertion)
            else
                Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected = Insertion(L, V, ChemPot, Beta, R_Cut, x, y, z, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected)
            end
        end
        if RN == 3 && length(x) > 1
            N_Removal += 1;
            if length(Pc_Random) == 1:
                Pc_Interpolation = Pc_Random(collect(keys(Pc_Random))[1])
            else
                if haskey(length(x) - 1, Pc_Random)
                    Pc_Interpolation = Pc_Random[length(x) - 1]
                else
                    Pc_Interpolation = Interpolation(Pc_Random, length(x))
                end
            end
            if rand() > (1 - Pc_Interpolation)^1000
                Energy, Virial, N_Removal_Accepted, N_Removal_Rejected = Removal_Mezei(Pc_Interpolation, L, V, Beta, ChemPot, R_Cut, Energy, Virial, N_Removal_Accepted, N_Removal_Rejected, x, y, z);
            else
                Energy, Virial, N_Removal_Accepted, N_Removal_Rejected = Removal(L, V, Beta, ChemPot, R_Cut, Energy, Virial, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
            end
        end
        if i % MC_Measurement == 0
            if i >= MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                Pressure_Array[N_Measurements] = (length(x) * T - Virial / 3.) / V            end
                Density_Array[N_Measurements] length(x) / V;
                g_r += Distribution(N_Bins, L, length(x) / V, x, y, z)
            end
            if i % 10MC_Measurement == 0
                N_Displacement_Accepted / N_Displacement > 0.55 ? Displacement *= 1.05 : Displacement *= 0.95
                Displacement < 0.05 ? Displacement = 0.05 : nothing
                Displacement > L / 4. ? Displacement = L / 4. : nothing
                N_Displacement, N_Displacement_Accepted = 0, 0;
            end
        end
    end
end

function Movement(L, Beta, Displacement, Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)
    j = rand(1:length(x))
    Energy_Old, Virial_Old = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    x_Old, y_Old, z_Old = x[j], y[j], z[j];
    x[j] += Displacement * (random() - 0.5);
    x[j] = PeriodicBoundaryConditions(L, x[j]);
    y[j] += Displacement * (random() - 0.5);
    y[j] = PeriodicBoundaryConditions(L, x[j]);
    z[j] += Displacement * (random() - 0.5);
    z[j] = PeriodicBoundaryConditions(L, x[j]);
    Energy_New, Virial_New = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z);
    Delta_E = Energy_New - Energy_Old;
    Delta_Virial = Virial_New - Virial_Old
    if rand() < exp(-Beta * Delta_E)
        N_Movement_Accepted += 1;
        N_Displacement_Accepted += 1;
        Energy += Delta_E;
        Virial += Delta_Virial;
    else
        N_Movement_Rejected += 1;
        x[j] = x_Old;
        y[j] = y_Old;
        z[j] = z_Old;
    end
    return Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted
end

function Insertion(L, V, ChemPot, Beta, R_Cut, x, y, z, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected)
    x_Insertion = L * (rand() - 0.5)
    y_Insertion = L * (rand() - 0.5)
    z_Insertion = L * (rand() - 0.5)
    Energy_Insertion, Virial_Insertion = Energy_Virial(L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
    if rand() < exp( Beta * (ChemPot - Energy_Insertion) + log(V / (length(x) + 1)) )
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion)
        append!(y, y_Insertion)
        append!(z, z_Insertion)
        Energy += Energy_Insertion;
        Virial += Virial_Insertion;
    else
        N_Insertion_Rejected += 1;
    end
    return Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected
end

function Insertion_Mezei(Beta, ChemPot, L, V, R_Cut, Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z, x_Insertion, y_Insertion, z_Insertion)
    j = rand(1:length(x_Insertion))
    if rand() < (V * Pc_Random[length(x)] / (length(x) + ) ) * exp(Beta * (ChemPot - Energy_Insertion))
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion[j])
        append!(y, y_Insertion[j])
        append!(z, z_Insertion[j])
        Energy += Energy_Insertion;
        Virial += Virial_Insertion
    else
        N_Removal_Rejected += 1;
    end
    return Energy, Virial, N_Insertion_Accepted, N_Insertion_Rejected
end

function Removal(L, V, Beta, ChemPot, R_Cut, Energy, Virial, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
    j = rand(1:length(x));
    Energy_Removal, Virial_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    if rand() < exp( Beta * (Energy_Removal - ChemPot) + log(length(x) / V) )
        N_Removal_Accepted += 1;
        deleteat!(x, j)
        deleteat!(y, j)
        deleteat!(z, j)
        Energy -= Energy_Removal;
        Virial -= Virial_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, Virial, N_Removal_Accepted, N_Removal_Rejected 
end

function Removal_Mezei(Pc_Interpolation, L, V, Beta, ChemPot, R_Cut, Energy, Virial, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
    j = rand(1:length(x))
    Energy_Removal, Virial_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    if rand() < ( length(x) / (V * Pc_Interpolation) ) * exp(Beta * (Energy_Removal - ChemPot))
        N_Removal_Accepted += 1;
        deleteat!(x, j)
        deleteat!(y, j)
        deleteat!(z, j)
        Energy -= Energy_Removal;
        Virial -= Virial_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, Virial, N_Removal_Accepted, N_Removal_Rejected
end

function u(r2)
    if r2 <= 1
        return Inf
    else if r2 <= 1.5^2
        return -1
    else
        return 0
end

function rdu(r2)

function Energy_Virial(L, R_Cut, rx, ry, rz, x, y, z)
    Energy, Virial = 0., 0.;
    for i = 1:length(x)
        Delta_x = rx - x[i];
        Delta_x = PeriodicBoundaryConditions(L, Delta_x);
        Delta_y = ry - y[i];
        Delta_y = PeriodicBoundaryConditions(L, Delta_y);
        Delta_z = rz - z[i];
        Delta_z = PeriodicBoundaryConditions(L, Delta_z);
    end
end

@time Mezei(-3, 250, 2.)