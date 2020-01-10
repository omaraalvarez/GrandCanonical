using Statistics;
using Plots;
using CSV;

function Mezei(ChemPot::Float64, L::Float64, T::Float64, Number_Run::Int64, Total_Run::Int64, R_Cut::Float64 = 3.)
    """     CONFIGURATIONAL STEPS       """
    MC_Relaxation_Steps = 1_000_000;
    MC_Equilibrium_Steps = 5_000_000;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;
    """     VARIABLE INITIALIZATION     """
    Overlap = 1.0;
    σ_p, λ_p = 0.5, 1.5;
    x, y, z = Float64[], Float64[], Float64[];
    V = L ^ 3.;
    MC_Measurement = convert(Int64, floor(V / (10 * (4. / 3.) * π * σ_p^3 )));
    Beta = 1. / T
    Pc, Pc_Sum, Pc_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    #Pc_Grid, Pc_Grid_Sum, Pc_Grid_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    #Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    Displacement, N_Displacement, N_Displacement_Accepted = L / 8., 0, 0;
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
    Energy, N_Measurements = 0., 0;
    Energy_Sum, Density_Sum = 0., 0.;
    Energy_Array, Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Average_Energy_Array, Average_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Error_Energy_Array, Error_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    #Energy, Virial, N_Measurements = 0., 0., 0;
    #Energy_Sum, Pressure_Sum, Density_Sum = 0., 0., 0.;
    #Energy_Array, Pressure_Array, Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    #Average_Energy_Array, Average_Pressure_Array, Average_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    N_Bins = 150;
    N_Image = 1;
    g_r = zeros(Float64, N_Bins);
    """     OUTPUT FILES        """
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))"
    mkpath("$Output_Route/Positions")
    Average_Energy_File = open("$Output_Route/Average_Energy.dat", "w");
    println(Average_Energy_File, "#\t< E / N >")
    Energy_File = open("$Output_Route/Energy.dat", "w");
    println(Energy_File, "#\tE / N ")
    #Average_Pressure_File = open("$Output_Route/Average_Pressure.dat", "w");
    #println(Average_Pressure_File, "#\t< P >")
    #Pressure_File = open("$Output_Route/Pressure.dat", "w");
    #println(Pressure_File, "#\tP")
    Average_Density_File = open("$Output_Route/Average_Density.dat", "w");
    println(Average_Density_File, "#\t< Density >")
    Density_File = open("$Output_Route/Density.dat", "w");
    println(Density_File, "#\tDensity")
    """     SIMULATIONS CYCLES      """
    @inbounds for i = 1:MC_Steps
        """     PRINTS PROGRESS TO SCREEN   """
        if i < MC_Relaxation_Steps && i % .01MC_Relaxation_Steps == 0
            println("$(convert(Int64, 100i / MC_Relaxation_Steps))% Relaxation [$Number_Run / $Total_Run]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            #println("P = $(round((length(x) * T - Virial / 3.) / V, digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Movement")
            println("   Accepted: $N_Movement_Accepted ($(round(100N_Movement_Accepted / N_Movement, digits = 2))%")
            println("   Rejected: $N_Movement_Rejected ($(round(100N_Movement_Rejected / N_Movement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)")
            println("   Rejected: $N_Insertion_Rejected ($(round(100N_Insertion_Rejected / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)")
            println("   Rejected: $N_Removal_Rejected ($(round(100N_Removal_Rejected / N_Removal, digits = 2))%)")
            println("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
        end

        if i > MC_Relaxation_Steps && i % .01MC_Equilibrium_Steps == 0
            println("$(convert(Int64, 100(i - MC_Relaxation_Steps) / MC_Equilibrium_Steps))% Equilibrium ($N_Measurements Measurements). [$Number_Run / $Total_Run]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            #println("P = $(round((length(x) * T - Virial / 3.) / V, digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Movement")
            println("   Accepted: $N_Movement_Accepted ($(round(100N_Movement_Accepted / N_Movement, digits = 2))%")
            println("   Rejected: $N_Movement_Rejected ($(round(100N_Movement_Rejected / N_Movement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)")
            println("   Rejected: $N_Insertion_Rejected ($(round(100N_Insertion_Rejected / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)")
            println("   Rejected: $N_Removal_Rejected ($(round(100N_Removal_Rejected / N_Removal, digits = 2))%)")
            println("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
        end

        if i > MC_Relaxation_Steps && i % .05MC_Equilibrium_Steps == 0
            Positions_File = open("$Output_Route/Positions/Pos_$N_Image.xyz", "w");
            for i = 1:length(x)
                i != length(x) ? println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i]),") : println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i])")
            end
            close(Positions_File)
            N_Image += 1
        end

        i == MC_Relaxation_Steps ? println("~~~    STARTING MEASUREMENT STEPS    ~~~") : nothing
        RN = rand(1:3);
        if RN == 1 && length(x) > 1
            N_Movement += 1;
            N_Displacement += 1;
            Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted = Movement(L, Beta, Displacement, Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)
            #Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted = Movement(L, Beta, Displacement, Energy, Virial, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)
        end
        if RN == 2
            N_Insertion += 1;
            #Pc_Grid, Pc_Grid_Sum, Pc_Grid_N = Grid_Excluded_Volume(Overlap, L, Pc_Grid, Pc_Grid_Sum, Pc_Grid_N, x, y, z)
            Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(Overlap, L, Pc, Pc_Sum, Pc_N, x, y, z)
            #Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N = Analytic_Excluded_Volume(Overlap, L, V, Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N, x, y, z)
            if length(x_Insertion) > 0
                Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion_Mezei(Beta, ChemPot, L, V, R_Cut, Energy, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z, x_Insertion, y_Insertion, z_Insertion, Pc)
            else
                Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion(L, V, ChemPot, Beta, R_Cut, x, y, z, Energy, N_Insertion_Accepted, N_Insertion_Rejected)
            end
        end
        if RN == 3 && length(x) > 1
            N_Removal += 1;
            if length(Pc) == 1
                Pc_Interpolation = Pc(collect(keys(Pc))[1])
            else
                if haskey(Pc, length(x) - 1)
                    Pc_Interpolation = Pc[length(x) - 1]
                else
                    Pc_Interpolation = Interpolation(Pc, length(x))
                end
            end
            if rand() > (1 - Pc_Interpolation)^1000
                Energy, N_Removal_Accepted, N_Removal_Rejected = Removal_Mezei(Pc_Interpolation, L, V, Beta, ChemPot, R_Cut, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z);
            else
                Energy, N_Removal_Accepted, N_Removal_Rejected = Removal(L, V, Beta, ChemPot, R_Cut, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
            end
        end
        if i % MC_Measurement == 0
            if i > MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                println(Energy_File, "$N_Measurements\t$(round(Energy_Array[N_Measurements], digits = 6))")
                Energy_Sum += Energy / length(x);
                Average_Energy_Array[N_Measurements] = Energy_Sum / N_Measurements;
                Error_Energy_Array[N_Measurements] = std(Energy_Array[1 : N_Measurements])
                println(Average_Energy_File, "$N_Measurements\t$(round(Average_Energy_Array[N_Measurements], digits = 6))\t$(round(Error_Energy_Array[N_Measurements], digits = 6))")
                #Pressure_Array[N_Measurements] = (length(x) * T - Virial / 3.) / V
                #println(Pressure_File, "$N_Measurements\t$(round(Pressure_Array[N_Measurements], digits = 6))")
                #Pressure_Sum += (length(x) * T - Virial / 3.) / V;
                #Average_Pressure_Array[N_Measurements] = Pressure_Sum / N_Measurements;
                #println(Average_Pressure_File, "$N_Measurements\t$(round(Average_Pressure_Array[N_Measurements], digits = 6))")
                Density_Array[N_Measurements] = length(x) / V;
                println(Density_File, "$N_Measurements\t$(round(Density_Array[N_Measurements], digits = 6))")
                Density_Sum += length(x) / V;
                Average_Density_Array[N_Measurements] = Density_Sum / N_Measurements;
                Error_Density_Array[N_Measurements] = std(Density_Array[1 : N_Measurements])
                println(Average_Density_File, "$N_Measurements\t$(round(Average_Density_Array[N_Measurements], digits = 6))\t$(round(Error_Density_Array[N_Measurements], digits = 6))")
                g_r += Distribution(N_Bins, L, length(x) / V, x, y, z, R_Cut)
            end
            N_Displacement_Accepted / N_Displacement > 0.55 ? Displacement *= 1.05 : Displacement *= 0.95
            Displacement < 0.05 ? Displacement = 0.05 : nothing
            Displacement > L / 4. ? Displacement = L / 4. : nothing
            N_Displacement, N_Displacement_Accepted = 0, 0;
        end
    end
    close(Average_Energy_File)
    #close(Average_Pressure_File)
    close(Average_Density_File)
    Error_Density_Array[1], Error_Energy_Array[1] = 0., 0.;

    g_r /= N_Measurements;
    #Delta = L / (2N_Bins);
    Delta = R_Cut / N_Bins;
    r = zeros(Float64, N_Bins - 1)
    g_r_File = open("$Output_Route/Radial_Distribution.dat", "w");
    println(Average_Energy_File, "#r\t#Normalized Density\n")
    @inbounds for i = 1:N_Bins - 1
        r[i] = round((i + 0.5)*Delta, digits = 6)
        println(g_r_File, "$(r[i])\t$(round(g_r[i], digits = 6))")
    end
    g_r = g_r[1 : N_Bins - 1]
    close(g_r_File)
    Radial_Distribution_Plot = plot(r, g_r, title = "Chemical Potential = $ChemPot", guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Distance [r]", ylabel = "Normalized Density", width = 3, size = [1920, 1080])
    hline!([1.], style = :dash, width = 3, color = :black)
    savefig(Radial_Distribution_Plot, "$Output_Route/RadialDistribution")

    Pc_File = open("$Output_Route/Pc.dat", "w");
    println(Pc_File, "#N\t#Pc_Random\t#Pc_Grid\t#Pc_Analytic")
    Pc_Array = zeros(Float64, length(Pc) - 1)
    #Pc_Grid_Array = zeros(Float64, length(Pc) - 1)
    #Pc_Analytic_Array = zeros(Float64, length(Pc) - 1)
    @inbounds for i in keys(Pc)
        if i != 0
            Pc_Array[i] = Pc[i];
            #Pc_Grid_Array[i] = Pc_Grid[i];
            #Pc_Analytic_Array[i] = Pc_Analytic[i];
        end
    end
    for i = 1:length(Pc_Array)
        println(Pc_File, "$i\t$(round(Pc_Array[i], digits = 6))")
        #println(Pc_File, "$i\t$(round(Pc_Array[i], digits = 6))\t$(round(Pc_Grid_Array[i], digits = 6))\t$(round(Pc_Analytic_Array[i], digits = 6))")
    end
    close(Pc_File)
    Pc_Plot = plot(Pc_Array, label = "Random", xlabel = "N", ylabel = "Cavity Probability", ylims = (0, 1), width = 3, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, size = [1920, 1080])
    #plot!(Pc_Grid_Array, label = "Grid", width = 3)
    #plot!(Pc_Analytic_Array, label = "Analytic", width = 3)
    savefig(Pc_Plot, "$Output_Route/Pc")

    println("< E / N > = $(round(mean(Energy_Array[1:end - 1]), digits = 6)) ± $(round(std(Energy_Array[1:end - 1]), digits = 6))")
    Energy_Plot = plot(Energy_Array[1:end - 1], legend = false, xlabel = "Measurements", ylabel = "Energy [Unitless]", width = 2, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, size = [1920, 1080])
    hline!([mean(Energy_Array[1:end - 1])], title = "Chemical Potential = $ChemPot", color = :black, width = 2, linestyle = :dash)
    savefig(Energy_Plot, "$Output_Route/Energy")
    Energy_Histogram = histogram(Energy_Array[convert(Int64, floor(MC_Relaxation_Steps/MC_Measurement)):end - 1], title = "Chemical Potential = $ChemPot", bins = 20, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Energy [Unitless]", ylabel = "Frequency", size = [1920, 1080])
    vline!([mean(Energy_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Energy_Histogram, "$Output_Route/Energy_Histogram")
    Average_Energy_Plot = plot(Average_Energy_Array[1:end - 1], ribbon = Error_Energy_Array, fillalpha = 0.2, xlims = (1, N_Measurements), title = "Chemical Potential = $ChemPot", guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Measurements", ylabel = "< Energy > [Unitless]", width = 3, size = [1920, 1080])
    hline!([mean(Energy_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Average_Energy_Plot, "$Output_Route/Average_Energy")

    #println("< P > = $(round(mean(Pressure_Array), digits = 6)) ± $(round(std(Pressure_Array), digits = 6))")
    #Pressure_Plot = plot(Pressure_Array, legend = false, xlabel = "Measurements", ylabel = "Pressure [Unitless]", width = 2, size = [1200, 800])
    #hline!([mean(Pressure_Array)], color = :black, width = 2, linestyle = :dash)
    #savefig(Pressure_Plot, "$Output_Route/Pressure")
    #Pressure_Histogram = histogram(Pressure_Array[convert(Int64, floor(MC_Relaxation_Steps/MC_Measurement)):end], bins = 20, legend = false, xlabel = "Pressure [Unitless]", ylabel = "Frequency", size = [1200, 800])
    #vline!([mean(Pressure_Array)], color = :black, width = 2, linestyle = :dash)
    #savefig(Pressure_Histogram, "$Output_Route/Pressure_Histogram")
    #Average_Pressure_Plot = plot(Average_Pressure_Array, legend = false, xlabel = "Measurements", ylabel = "< Pressure > [Unitless]", width = 3, size = [1200, 800])
    #hline!([mean(Pressure_Array)], color = :black, width = 2, linestyle = :dash)
    #savefig(Average_Pressure_Plot, "$Output_Route/Average_Pressure")

    println("< N > = $(round(V*mean(Density_Array[1:end - 1]), digits = 6)) ± $(round(V*std(Density_Array[1:end - 1]), digits = 6))")
    println("< Density > = $(round(mean(Density_Array[1:end - 1]), digits = 6)) ± $(round(std(Density_Array[1:end - 1]), digits = 6))")
    Density_Plot = plot(Density_Array[1:end - 1], legend = false, xlabel = "Measurements", ylabel = "Density [Unitless]", width = 2, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, size = [1920, 1080])
    hline!([mean(Density_Array[1:end - 1])], title = "Chemical Potential = $ChemPot", color = :black, width = 2, linestyle = :dash)
    savefig(Density_Plot, "$Output_Route/Density")
    Density_Histogram = histogram(Density_Array[convert(Int64, floor(MC_Relaxation_Steps/MC_Measurement)):end - 1], title = "Chemical Potential = $ChemPot", bins = 20, legend = false, xlabel = "Density [Unitless]", ylabel = "Frequency", guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, size = [1920, 1080])
    vline!([mean(Density_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Density_Histogram, "$Output_Route/Density_Histogram")
    Average_Density_Plot = plot(Average_Density_Array[1:end - 1], ribbon = Error_Density_Array, fillalpha = 0.2, xlims = (1, N_Measurements), title = "Chemical Potential = $ChemPot", legend = false, xlabel = "Measurements", ylabel= "< Density > [Unitless]", width = 3, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, size = [1920, 1080])
    hline!([mean(Density_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Average_Density_Plot, "$Output_Route/Average_Density")

    Summary_File = open("$Output_Route/Summary.dat", "w")
    println(Summary_File, "< E / N > = $(round(mean(Energy_Array), digits = 6)) ± $(round(std(Energy_Array), digits = 6))")
    #println(Summary_File, "< P > = $(round(mean(Pressure_Array), digits = 6)) ± $(round(std(Pressure_Array), digits = 6))")
    println(Summary_File, "< N > = $(round(V*mean(Density_Array), digits = 6)) ± $(round(V*std(Density_Array), digits = 6))")
    println(Summary_File, "< Density > = $(round(mean(Density_Array), digits = 6)) ± $(round(std(Density_Array), digits = 6))")
    close(Summary_File)

    Povray_ini(ChemPot, T, N_Image - 1)
    Povray_Pov(L, ChemPot, T)
    #run(`povray $Output_Route/Positions/MC_Animation.ini`)
    
    return mean(Density_Array[1:end - 1]), std(Density_Array[1:end - 1]), Energy_Plot, Energy_Histogram, Average_Energy_Plot, Density_Plot, Density_Histogram, Average_Density_Plot, Radial_Distribution_Plot
end

function Movement(L::Float64, Beta::Float64, Displacement::Float64, Energy::Float64, N_Movement_Accepted::Int64, N_Movement_Rejected::Int64, N_Displacement_Accepted::Int64, R_Cut::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x))
    Energy_Old = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    x_Old, y_Old, z_Old = x[j], y[j], z[j];
    x[j] += Displacement * (rand() - 0.5);
    x[j] = PeriodicBoundaryConditions(L, x[j]);
    y[j] += Displacement * (rand() - 0.5);
    y[j] = PeriodicBoundaryConditions(L, y[j]);
    z[j] += Displacement * (rand() - 0.5);
    z[j] = PeriodicBoundaryConditions(L, z[j]);
    Energy_New = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z);
    Delta_E = Energy_New - Energy_Old;
    #Delta_Virial = Virial_New - Virial_Old
    if rand() < exp(-Beta * Delta_E)
        N_Movement_Accepted += 1;
        N_Displacement_Accepted += 1;
        Energy += Delta_E;
        #Virial += Delta_Virial;
    else
        N_Movement_Rejected += 1;
        x[j] = x_Old;
        y[j] = y_Old;
        z[j] = z_Old;
    end
    return Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted
end

function Insertion(L::Float64, V::Float64, ChemPot::Float64, Beta::Float64, R_Cut::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, Energy::Float64, N_Insertion_Accepted::Int64, N_Insertion_Rejected::Int64)
    x_Insertion = L * (rand() - 0.5)
    y_Insertion = L * (rand() - 0.5)
    z_Insertion = L * (rand() - 0.5)
    Energy_Insertion = Energy_Virial(L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
    if rand() < exp( Beta * (ChemPot - Energy_Insertion) + log(V / (length(x) + 1)) )
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion)
        append!(y, y_Insertion)
        append!(z, z_Insertion)
        Energy += Energy_Insertion;
        #Virial += Virial_Insertion;
    else
        N_Insertion_Rejected += 1;
    end
    return Energy, N_Insertion_Accepted, N_Insertion_Rejected
end

function Insertion_Mezei(Beta::Float64, ChemPot::Float64, L::Float64, V::Float64, R_Cut::Float64, Energy::Float64, N_Insertion_Accepted::Int64, N_Insertion_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, x_Insertion::Array{Float64, 1}, y_Insertion::Array{Float64, 1}, z_Insertion::Array{Float64, 1}, Pc::Dict{Int64, Float64})
    j = rand(1:length(x_Insertion))
    Energy_Insertion = Energy_Virial(L, R_Cut, x_Insertion[j], y_Insertion[j], z_Insertion[j], x, y, z);
    if rand() < (V * Pc[length(x)] / (length(x) + 1) ) * exp(Beta * (ChemPot - Energy_Insertion))
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion[j])
        append!(y, y_Insertion[j])
        append!(z, z_Insertion[j])
        Energy += Energy_Insertion;
        #Virial += Virial_Insertion
    else
        N_Insertion_Rejected += 1;
    end
    return Energy, N_Insertion_Accepted, N_Insertion_Rejected
end

function Removal(L::Float64, V::Float64, Beta::Float64, ChemPot::Float64, R_Cut::Float64, Energy::Float64, N_Removal_Accepted::Int64, N_Removal_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x));
    Energy_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    if rand() < exp( Beta * (Energy_Removal - ChemPot) + log(length(x) / V) )
        N_Removal_Accepted += 1;
        deleteat!(x, j)
        deleteat!(y, j)
        deleteat!(z, j)
        Energy -= Energy_Removal;
        #Virial -= Virial_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, N_Removal_Accepted, N_Removal_Rejected 
end

function Removal_Mezei(Pc_Interpolation::Float64, L::Float64, V::Float64, Beta::Float64, ChemPot::Float64, R_Cut::Float64, Energy::Float64, N_Removal_Accepted::Int64, N_Removal_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x))
    Energy_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
    if rand() < ( length(x) / (V * Pc_Interpolation) ) * exp(Beta * (Energy_Removal - ChemPot))
        N_Removal_Accepted += 1;
        deleteat!(x, j)
        deleteat!(y, j)
        deleteat!(z, j)
        Energy -= Energy_Removal;
        #Virial -= Virial_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, N_Removal_Accepted, N_Removal_Rejected
end

function u_SquareWell(r2::Float64, σ::Float64, λ::Float64, e::Float64 = 1.)
    if r2 <= (2σ)^2
        return Inf
    elseif r2 <= λ^2
        return -e
    else
        return 0
    end
end

#function rdu(r2::Float64)
#    return 4.0(6.0(1.0 / r2)^3.0 - 12.0(1.0 / r2)^6.0)
#end

function Energy_Virial(L::Float64, R_Cut::Float64, rx::Float64, ry::Float64, rz::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    #Energy, Virial = 0., 0.;
    Energy = 0.;
    @inbounds for i = 1:length(x)
        Delta_x = rx - x[i];
        Delta_x = PeriodicBoundaryConditions(L, Delta_x);
        Delta_y = ry - y[i];
        Delta_y = PeriodicBoundaryConditions(L, Delta_y);
        Delta_z = rz - z[i];
        Delta_z = PeriodicBoundaryConditions(L, Delta_z);
        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
        if r2 != 0.
            if r2 < R_Cut^2
                Energy += u_SquareWell(r2, 0.5, 1.5);
                #Virial += rdu(r2);
            end
        end
    end
    return Energy#, Virial
end

function PeriodicBoundaryConditions(L::Float64, x::Float64)
    if x < - L / 2.
        x += L;
    elseif x > L / 2.
        x -= L;
    end
    return x
end

function Random_Excluded_Volume(Overlap::Float64, L::Float64, Pc::Dict{Int64, Float64}, Pc_Sum::Dict{Int64, Float64}, Pc_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    N_Random = 1000;
    N_in = 0;
    x_Insertion, y_Insertion, z_Insertion = Float64[], Float64[], Float64[];
    if length(x) == 1
        for i = 1:N_Random
            x_V = L * (rand() - 0.5);
            y_V = L * (rand() - 0.5);
            z_V = L * (rand() - 0.5);
            N_in += 1;
            append!(x_Insertion, x_V)
            append!(y_Insertion, y_V)
            append!(z_Insertion, z_V)
        end
    else
        @inbounds for i = 1:N_Random
            x_V = L * (rand() - 0.5);
            y_V = L * (rand() - 0.5);
            z_V = L * (rand() - 0.5);
            for j = 1:length(x)
                Delta_x = x_V - x[j];
                Delta_x = PeriodicBoundaryConditions(L, Delta_x);
                Delta_y = y_V - y[j];
                Delta_y = PeriodicBoundaryConditions(L, Delta_y);
                Delta_z = z_V - z[j];
                Delta_z = PeriodicBoundaryConditions(L, Delta_z);
                r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
                if r2 < Overlap^2
                    break
                end
                if j == length(x)
                    N_in += 1;
                    append!(x_Insertion, x_V)
                    append!(y_Insertion, y_V)
                    append!(z_Insertion, z_V)
                end
            end
        end
    end
    Volume_Ratio = N_in /N_Random;
    if !haskey(Pc, length(x))
        Pc_N[length(x)] = 1
        Pc_Sum[length(x)] = Volume_Ratio;
        Pc[length(x)] = Volume_Ratio;
    else
        Pc_N[length(x)] += 1;
        Pc_Sum[length(x)] += Volume_Ratio;
        Pc[length(x)] = Pc_Sum[length(x)] / Pc_N[length(x)]
    end
    return Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion
end

function Grid_Excluded_Volume(Overlap::Float64, L::Float64, Pc_Grid::Dict{Int64, Float64}, Pc_Grid_Sum::Dict{Int64, Float64}, Pc_Grid_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    Grid = 10;
    Delta = L / (Grid + 1);
    EVMPS = zeros(Grid, Grid, Grid);
    @inbounds for i_x = 1:Grid
        @inbounds for i_y = 1:Grid
            @inbounds for i_z = 1:Grid
                x_Grid = i_x * Delta - L / 2.;
                y_Grid = i_y * Delta - L / 2.;
                z_Grid = i_z * Delta - L / 2.;
                @inbounds for i = 1:length(x)
                    Delta_x = x_Grid - x[i];
                    Delta_y = y_Grid - y[i];
                    Delta_z = z_Grid - z[i];
                    r2 = Delta_x^2 + Delta_y^2 + Delta_z^2
                    if r2 < Overlap^2
                        EVMPS[i_x, i_y, i_z] = 1;
                        break
                    end
                end
            end
        end
    end
    Volume_Ratio = 1 - sum(EVMPS) / Grid ^ 3;
    if !haskey(Pc_Grid, length(x))
        Pc_Grid_N[length(x)] = 1;
        Pc_Grid_Sum[length(x)] = Volume_Ratio;
        Pc_Grid[length(x)] = Volume_Ratio;
    else
        Pc_Grid_N[length(x)] += 1;
        Pc_Grid_Sum[length(x)] += Volume_Ratio;
        Pc_Grid[length(x)] = Pc_Grid_Sum[length(x)] / Pc_Grid_N[length(x)];
    end
    return Pc_Grid, Pc_Grid_Sum, Pc_Grid_N
end

function Analytic_Excluded_Volume(Overlap::Float64, L::Float64, V::Float64, Pc_Analytic::Dict{Int64, Float64}, Pc_Analytic_Sum::Dict{Int64, Float64}, Pc_Analytic_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    V_Excluded = length(x) * (4. / 3.) * pi * Overlap^3;
    V_Excluded_Correction = 0;
    @inbounds for i = 1:length(x) - 1
        @inbounds for j = i + 1:length(x)
            Delta_x = x[i] - x[j];
            Delta_y = y[i] - y[j];
            Delta_z = z[i] - z[j];
            r2 = Delta_x^2. + Delta_y^2. + Delta_z^2.;
            if r2 < (2Overlap)^2.
                d = sqrt(r2);
                V_Excluded_Correction += (pi / 12.) * (4 * Overlap + d) * (2 * Overlap - d)^2.;
            else
                Delta_x > L - 1 ? Delta_x -= L : nothing
                Delta_y > L - 1 ? Delta_y -= L : nothing
                Delta_z > L - 1 ? Delta_z -= L : nothing
                r2 = Delta_x^2. + Delta_y^2. + Delta_z^2.
                if r2 < (2Overlap)^2.
                    d = sqrt(r2);
                    V_Excluded_Correction += (pi / 12.) * (4 * Overlap + d) * (2 * Overlap - d)^2.;
                end
            end
        end
    end
    Volume_Ratio = 1 - (V_Excluded - V_Excluded_Correction) / V;
    Volume_Ratio > 1 || Volume_Ratio < 0 ? error("Volume Ratio ($Volume_Ratio) can't be negative or greater than one.") : nothing
    if !haskey(Pc_Analytic, length(x))
        Pc_Analytic_N[length(x)] = 1;
        Pc_Analytic_Sum[length(x)] = Volume_Ratio;
        Pc_Analytic[length(x)] = Volume_Ratio;
    else
        Pc_Analytic_N[length(x)] += 1;
        Pc_Analytic_Sum[length(x)] += Volume_Ratio;
        Pc_Analytic[length(x)] = Pc_Analytic_Sum[length(x)] / Pc_Analytic_N[length(x)];
    end
    return Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N
end

function Interpolation(Pc::Dict{Int64, Float64}, l::Float64)
    if length(Pc) == 1
        Pc_Interpolation = collect(keys(Pc))[1];
    elseif length(Pc) > 1
        lim_inf = l -1
        while true
            if haskey(Pc, lim_inf)
                break
            end
            lim_inf -= 1;
        end

        lim_sup = l + 1;
        while true
            if haskey(Pc, lim_sup)
                break
            end
            lim_sup += 1;
        end
    elseif length(Pc) == 0
        error("Pc has no values still.")
    end
    Pc_Interpolation = Pc[lim_inf] * ((1 - l - lim_inf) / (lim_sup - lim_inf)) + Pc[lim_sup] * ((l - lim_inf) / (lim_sup - lim_inf))
    return Pc_Interpolation
end

function Distribution(N_Bins::Int64, L::Float64, Density::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Float64)
    Delta = L / (2N_Bins);
    Delta = R_Cut / N_Bins;
    g_r = zeros(Float64, N_Bins);
    @inbounds for i = 1:length(x) - 1
        @inbounds for j = i + 1:length(x)
            Delta_x = x[i] - x[j];
            Delta_x = PeriodicBoundaryConditions(L, Delta_x);
            Delta_y = y[i] - y[j];
            Delta_y = PeriodicBoundaryConditions(L, Delta_y);
            Delta_z = z[i] - z[j];
            Delta_z = PeriodicBoundaryConditions(L, Delta_z);
            r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
            #if r2 < (L / 2.)^2
            if r2 < R_Cut^2
                l = convert(Int64, floor(sqrt(r2) / Delta));
                g_r[l] += 2;
            end
        end
    end
    @inbounds for l = 1:N_Bins
        g_r[l] /= (length(x) * 4 * pi * (l * Delta)^2 * Delta * Density)
    end
    return g_r
end

function Povray_Pov(L::Float64, ChemPot::Float64, T::Float64)
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Positions"
    Pov_File = open("$Output_Route/MC_Animation.pov", "w");
    println(Pov_File, "global_settings {\n\tambient_light rgb <0.2, 0.2, 0.2>\tmax_trace_level 15\n}\n")
    println(Pov_File, "background { color rgb <1, 1, 1> }\n")
    println(Pov_File, "#default { finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.7} }\n")
    println(Pov_File, "camera {\n\tperspective\n\tlocation <0, $(-1.55L), 0>\n\tlook_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(-5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(+5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(-5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(+5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "#macro Particle(rx, ry, rz)\n\tintersection {\n\t\tsphere {\n\t\t\t<rx, ry, rz>, 0.5\n\t\t\tpigment {rgbt <107/255, 173/255, 197/255, 0> }\n\t\t}\n\t\tbox {\n\t\t\t<-L/2, -L/2, L/2>,\t<L/2, L/2, -L/2>\n\t\t\tpigment {rgbt <107/255, 173/255, 197/255, 0> }\n\t\t}\n\t}\n#end\n")
    println(Pov_File, "#macro Walls(L, h)\n\tunion {\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, h / 2>, <-L / 2, L / 2, h / 2>, <L / 2, -L / 2, h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t\ttriangle {\n\t\t\t<-L / 2, -L / 2, h / 2>, <-L / 2, L / 2, h / 2>, <L / 2, -L / 2, h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t}\n
        union{\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h / 2>, <-L / 2, L / 2, -h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t\ttriangle {\n\t\t\t<-L / 2, -L / 2, -h / 2>, <-L / 2, L / 2, -h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t}\n#end\n")
    println(Pov_File, "#macro Box(L, h)\n\tunion{\n\t\ttriangle {\n\t\t\t<-L / 2, -L / 2, h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t\ttriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t}
        union {\n\t\ttriangle {\n\t\t\t<L / 2, -L / 2, h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t}\n\n\t
        union {\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <-L / 2, L / 2, -h /2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\ntriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <L / 2, L / 2, h / 2>\n\t\t\tpigment { rgbt <0.75, 0.75, 0.75, 0.5> }\n\t\t}\n\t}\n#end")
    println(Pov_File, "#declare L = $L;\n#declare h = $L;")
    println(Pov_File, """#fopen File_Positions concat("$Output_Route/Pos_", str(clock, 1, 0), ".xyz") read""")
    println(Pov_File, "\t#while (defined( File_Positions ))\n\t\t#read (File_Positions, rx, ry, rz)\n\t\tParticle(rx, ry, -rz)\n\t\t#declare PBC = false;\n\t\t#if (rx > (L - 1) / 2)\n\t\t\t#declare rx = rx - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (rx < -(L - 1) / 2)\n\t\t\t#declare rx = rx + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t
        #if (ry > (L - 1) / 2)\n\t\t\t#declare ry = ry - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (ry < -(L - 1) / 2)\n\t\t\t#declare ry = ry + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (rz > (h - 1) / 2)\n\t\t\t#declare rz = rx - h;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (rz < -(h - 1) / 2)\n\t\t\t#declare rz = rz + h;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (PBC)\n\t\t\tParticle(rx, ry, -rz)\n\t\t#end\n\t#end\n#fclose File")
    println(Pov_File, "Walls(L, h)\nBox(L, h)")
    close(Pov_File)
end

function Povray_ini(ChemPot::Float64, T::Float64, Frames::Int64)
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Positions"
    mkpath("$Output_Route")
    Ini_File = open("$Output_Route/MC_Animation.ini", "w");
    println(Ini_File, "Input_File_Name = $Output_Route/MC_Animation.pov")
    println(Ini_File, "Output_File_Name = $Output_Route/")
    println(Ini_File, "+W800 +H800\n")
    println(Ini_File, "Initial_Frame = 1")
    println(Ini_File, "Final_Frame = $Frames")
    println(Ini_File, "Initial_Clock = 1")
    println(Ini_File, "Final_Clock = $Frames\n")
    println(Ini_File, "Cyclic_Animation = off")
    close(Ini_File)
end

function Cycled_Mezei()
    ChemPot = [5.698524];
    T = 3.;
    L = 10.
    Mean_Density = zeros(Float64, length(ChemPot));
    Std_Density = zeros(Float64, length(ChemPot));
    j = 1;

    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))"
    mkpath("$Output_Route")
    #Density_File = open("$Output_Route/Density_T_$T.dat", "w+");
    #println(Density_File, "ChemPot\tDensity\tErrorDensity")
    p_Radial = Array{Any, 1}(undef, length(ChemPot))
    p_Energy = Array{Any, 1}(undef, length(ChemPot))
    p_Histogram_Energy = Array{Any, 1}(undef, length(ChemPot))
    p_Average_Energy = Array{Any, 1}(undef, length(ChemPot))
    p_Density = Array{Any, 1}(undef, length(ChemPot))
    p_Histogram_Density = Array{Any, 1}(undef, length(ChemPot))
    p_Average_Density = Array{Any, 1}(undef, length(ChemPot))

    for μ in ChemPot
        Mean_Density[j], Std_Density[j], p_Energy[j], p_Histogram_Energy[j], p_Average_Energy[j], p_Density[j], p_Histogram_Density[j], p_Average_Density[j], p_Radial[j] = Mezei(μ, L, T, j, length(ChemPot));
        #println(Density_File, "$μ\t$(Mean_Density[j])\t$(Std_Density[j])")
        j += 1;
    end

    #y = ones(3)
    #title = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, y[2], Plots.text("T = $T", :black, 25)), axis=false, grid = false, leg=false,size=(1920,50))

    #Radial_Plot = plot(p_Radial[1], p_Radial[2], p_Radial[3], p_Radial[4], p_Radial[5], p_Radial[6], p_Radial[7], p_Radial[8], p_Radial[9], p_Radial[10], p_Radial[11], p_Radial[12], p_Radial[13], p_Radial[14], p_Radial[15], p_Radial[16], layout = (4, 4), size = [1920, 1080])
    #Complete_Radial_Plot = plot(title, Radial_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
    #savefig(Complete_Radial_Plot, "$Output_Route/Radial_Distribution_$T")

    #Energy_Plot = plot(p_Energy[1], p_Energy[2], p_Energy[3], p_Energy[4], p_Energy[5], p_Energy[6], p_Energy[7], p_Energy[8], p_Energy[9], p_Energy[10], p_Energy[11], p_Energy[12], p_Energy[13], p_Energy[14], p_Energy[15], p_Energy[16], layout = (4, 4), size = [1920, 1080])
    #Complete_Energy_Plot = plot(title, Energy_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
    #savefig(Complete_Energy_Plot, "$Output_Route/Energy_$T")

    #Energy_Histogram_Plot = plot(p_Histogram_Energy[1], p_Histogram_Energy[2], p_Histogram_Energy[3], p_Histogram_Energy[4], p_Histogram_Energy[5], p_Histogram_Energy[6], p_Histogram_Energy[7], p_Histogram_Energy[8], p_Histogram_Energy[9], p_Histogram_Energy[10], p_Histogram_Energy[11], p_Histogram_Energy[12], p_Histogram_Energy[13], p_Histogram_Energy[14], p_Histogram_Energy[15], p_Histogram_Energy[16], layout = (4, 4), size = [1920, 1080])
    #Complete_Energy_Histogram_Plot = plot(title, Energy_Histogram_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
    #savefig(Complete_Energy_Histogram_Plot, "$Output_Route/Energy_Histogram_$T")

    #Energy_Average_Plot = plot(p_Average_Energy[1], p_Average_Energy[2], p_Average_Energy[3], p_Average_Energy[4], p_Average_Energy[5], p_Average_Energy[6], p_Average_Energy[7], p_Average_Energy[8], p_Average_Energy[9], p_Average_Energy[10], p_Average_Energy[11], p_Average_Energy[12], p_Average_Energy[13], p_Average_Energy[14], p_Average_Energy[15], p_Average_Energy[16], layout = (4, 4), size = [1920, 1080])
    #Complete_Average_Energy_Plot = plot(title, Energy_Average_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
    #savefig(Complete_Average_Energy_Plot, "$Output_Route/Average_Energy_$T")

    #Density_Plot = plot(p_Density[1], p_Density[2], p_Density[3], p_Density[4], p_Density[5], p_Density[6], p_Density[7], p_Density[8], p_Density[9], p_Density[10], p_Density[11], p_Density[12], p_Density[13], p_Density[14], p_Density[15], p_Density[16], layout = (4, 4), size = [1920, 1080])
    #Complete_Density_Plot = plot(title, Density_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
    #savefig(Complete_Density_Plot, "$Output_Route/Density_$T")

    #Density_Histogram_Plot = plot(p_Histogram_Density[1], p_Histogram_Density[2], p_Histogram_Density[3], p_Histogram_Density[4], p_Histogram_Density[5], p_Histogram_Density[6], p_Histogram_Density[7], p_Histogram_Density[8], p_Histogram_Density[9], p_Histogram_Density[10], p_Histogram_Density[11], p_Histogram_Density[12], p_Histogram_Density[13], p_Histogram_Density[14], p_Histogram_Density[15], p_Histogram_Density[16], layout = (4, 4), size = [1920, 1080])
    #Complete_Density_Histogram_Plot = plot(title, Density_Histogram_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
    #savefig(Complete_Density_Histogram_Plot, "$Output_Route/Density_Histogram_$T")

    #Density_Average_Plot = plot(p_Average_Density[1], p_Average_Density[2], p_Average_Density[3], p_Average_Density[4], p_Average_Density[5], p_Average_Density[6], p_Average_Density[7], p_Average_Density[8], p_Average_Density[9], p_Average_Density[10], p_Average_Density[11], p_Average_Density[12], p_Average_Density[13], p_Average_Density[14], p_Average_Density[15], p_Average_Density[16], layout = (4, 4), size = [1920, 1080])
    #Complete_Average_Density_Plot = plot(title, Density_Average_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
    #savefig(Complete_Average_Density_Plot, "$Output_Route/Average_Density_$T")

    #Average_Density_Plot = plot(ChemPot, Mean_Density, yerror = Std_Density, background_color_legend = false, foreground_color_legend = false, legend = :topleft, label = "T = $T", xlabel = "Chemical Potential", ylabel= "< Density > [Unitless]", width = 3, ylims = (0, findmax(Mean_Density)[1] + 0.05 ), size = [1200, 800])
    #savefig(Average_Density_Plot, "$Output_Route/Density_T_$T")
    #close(Density_File)
end

@time Cycled_Mezei()