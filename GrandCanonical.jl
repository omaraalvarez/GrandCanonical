using Statistics;
using Plots; gr()
using Plots.PlotMeasures;
using Test;
using LaTeXStrings;
using StatsPlots;
using Distributions;

function GrandCanonical_MonteCarlo(μ::Type, T::Type, R_Cut::Type = 3.) where {Type <: Real}
    ##################################### CONFIGURATIONAL STEPS #############################
    println("\t\tGRAND CANONICAL MONTE CARLO")
    MC_Relaxation_Steps = 20_000;
    MC_Equilibrium_Steps = 250_000;
    MC_Measurement = 10;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;
    ##################################### VARIABLE INITIALIZATION ###########################
    L, σ, λ = 20., 1.0, 1.5;
    V, Beta, N_Bins, Equilibrium = L ^ 3., 1. / T, 150, false;
    N_Id = convert(Int64, round(0.5 * V));
    Pc, Pc_Sum, Pc_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    Displacement, N_Displacement, N_Displacement_Accepted = L / 8., 0, 0;
    N_Insertion, N_Insertion_Accepted = 0, 0;
    N_Removal, N_Removal_Accepted = 0, 0;
    Energy, Density, N_Measurements = 0., 0., 0;
    Energy_Array, Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Mean_Energy_Array, Mean_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    STD_Energy_Array, STD_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    g_r = zeros(Float64, N_Bins);
    ####################################### OUTPUT ROUTE, INITIAL POSITIONS AND INITIAL ENERGY ###########################
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(μ, digits = 3))"
    mkpath("$Output_Route")
    x, y, z = Float64[], Float64[], Float64[];
    ################################################# SIMULATION CYCLES ###########################################
    @inbounds for k = 1:MC_Steps
        @test all(Array(x) .<= L / 2.) && all(Array(x) .>= -L / 2.)
        @test all(Array(y) .<= L / 2.) && all(Array(y) .>= -L / 2.)
        @test all(Array(z) .<= L / 2.) && all(Array(z) .>= -L / 2.)
        ####################################### PRINTS SIMULATION PROGRESS TO SCREEN ########################################
        if k < MC_Relaxation_Steps && k % .01MC_Relaxation_Steps == 0
            println("$(convert(Int64, 100k / MC_Relaxation_Steps))% Relaxation: [μ = $μ, L = $L, T = $T]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x)) Particles")
            println("ρ = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 2))%)\tRejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100(N_Displacement - N_Displacement_Accepted)/N_Displacement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)\tRejected: $(N_Insertion - N_Insertion_Accepted) ($(round(100(N_Insertion - N_Insertion_Accepted) / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)\tRejected: $(N_Removal - N_Removal_Accepted) ($(round(100(N_Removal - N_Removal_Accepted) / N_Removal, digits = 2))%)\n")
            N_Displacement, N_Displacement_Accepted = 0, 0;
            N_Insertion, N_Insertion_Accepted = 0, 0;
            N_Removal, N_Removal_Accepted = 0, 0;
        end

        if k > MC_Relaxation_Steps && k % .01MC_Equilibrium_Steps == 0
            println("$(convert(Int64, ceil(100(k - MC_Relaxation_Steps) / MC_Equilibrium_Steps)))% Equilibrium: [μ = $μ, L = $L, T = $T, $(N_Measurements + 1) Measurements]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x)) Particles")
            println("ρ = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 2))%)\tRejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100(N_Displacement - N_Displacement_Accepted)/N_Displacement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)\tRejected: $(N_Insertion - N_Insertion_Accepted) ($(round(100(N_Insertion - N_Insertion_Accepted) / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)\tRejected: $(N_Removal - N_Removal_Accepted) ($(round(100(N_Removal - N_Removal_Accepted) / N_Removal, digits = 2))%)\n")
            N_Displacement, N_Displacement_Accepted = 0, 0;
            N_Insertion, N_Insertion_Accepted = 0, 0;
            N_Removal, N_Removal_Accepted = 0, 0;
        end

        k == MC_Relaxation_Steps ? println("- FINISHED RELAXATION STEPS\n") : nothing
        k == MC_Relaxation_Steps ? Equilibrium = true : nothing
        @inbounds for i = 1:N_Id
            RN = rand(1:3);
            if RN == 1 && length(x) > 1
                N_Displacement += 1;
                Energy, N_Displacement_Accepted = Movement(L, Beta, Displacement, Energy, N_Displacement_Accepted, x, y, z)
            end
            if RN == 2
                N_Insertion += 1;
                Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(Equilibrium, σ, L, Pc, Pc_Sum, Pc_N, x, y, z)
                if length(x_Insertion) > 0
                    Energy, N_Insertion_Accepted = Insertion_Mezei(L, μ, Beta, Energy, N_Insertion_Accepted, x, y, z, x_Insertion, y_Insertion, z_Insertion, Pc)
                else
                    Energy, N_Insertion_Accepted = Insertion(L, μ, Beta, Energy, N_Insertion_Accepted, x, y, z)
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
                if rand() > (1 - Pc_Interpolation)^(Equilibrium ? 200 : 1000)
                    Energy, N_Removal_Accepted = Removal_Mezei(L, μ, Beta, Energy, N_Removal_Accepted, x, y, z, Pc_Interpolation);
                else
                    Energy, N_Removal_Accepted = Removal(L, μ, Beta, Energy, N_Removal_Accepted, x, y, z)
                end
            end
        end

        ##################################################### MEASUREMENTS ##################################################
        if k % MC_Measurement == 0
            if k > MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                Density_Array[N_Measurements] = length(x) / V;
                Mean_Energy_Array[N_Measurements] = mean(Energy_Array[1:N_Measurements]);
                Mean_Density_Array[N_Measurements] = mean(Density_Array[1:N_Measurements])
                if N_Measurements > 1
                    STD_Energy_Array[N_Measurements] = std(Energy_Array[1:N_Measurements]);
                    STD_Density_Array[N_Measurements] = std(Density_Array[1:N_Measurements]);
                end
                g_r += RadialDistributionFunction(N_Bins, L, length(x) / V, x, y, z)
            end
            1. * N_Displacement_Accepted / N_Displacement > 0.55 ? Displacement *= 1.05 : Displacement *= 0.95
            Displacement < 0.05 ? Displacement = 0.05 : nothing
            Displacement > L / 4. ? Displacement = L / 4. : nothing

        end  
    end
    ####################################################### END OF SIMULATION CYCLES #############################################
    ########################################################### SUMMARY FILE ##############################################
    println("< U* / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 6)) ± $(round(STD_Energy_Array[N_Measurements], digits = 6))")
    println("< ρ* > = $(round(Mean_Density_Array[N_Measurements], digits = 6)) ± $(round(STD_Density_Array[N_Measurements], digits = 6))")
    Summary_File = open("$Output_Route/Summary.dat", "w+")
    println(Summary_File, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   INPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, " μ = $μ\nL = $L\tV = $V\nT = $T\n$MC_Relaxation_Steps Relaxation Steps.\n$MC_Equilibrium_Steps Equilibrium Steps.\tMeasurements every $MC_Measurement steps.")
    println(Summary_File, "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   OUTPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, "< U* / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 6)) ± $(round(STD_Energy_Array[N_Measurements], digits = 6))")
    println(Summary_File, "< ρ* > = $(round(Mean_Density_Array[N_Measurements], digits = 6)) ± $(round(STD_Density_Array[N_Measurements], digits = 6))")
    close(Summary_File)

    ################################################################ OUTPUT ##############################################################
    Energy_File = open("$Output_Route/Energy.dat", "w+");
    Density_File = open("$Output_Route/Density.dat", "w+");
    println(Energy_File, "Step\tEnergy\tMean_Energy\tStd_Energy")
    println(Density_File, "Step\tDensity\tMean_Density\tStd_Density")
    for i = 1:N_Measurements
        println(Energy_File, "$i\t$(round(Energy_Array[i], digits = 6))\t$(round(Mean_Energy_Array[i], digits = 6))\t$(round(STD_Energy_Array[i], digits = 6))")
        println(Density_File, "$i\t$(round(Density_Array[i], digits = 6))\t$(round(Mean_Density_Array[i], digits = 6))\t$(round(STD_Density_Array[i], digits = 6))")
    end
    close(Energy_File)
    close(Density_File)

    ##################################################### RADIAL DISTRIBUTION FUNCTION ########################################################
    g_r ./= N_Measurements;
    Delta, r = R_Cut / N_Bins, zeros(Float64, N_Bins);
    g_r_File = open("$Output_Route/RadialDistribution.dat", "w+");
    println(g_r_File, "r\tg_r\n")
    @inbounds for i = 1:N_Bins
        r[i] = round((i + 0.5)*Delta, digits = 6);
        println(g_r_File, "$(r[i])\t$(round(g_r[i], digits = 6))")
    end
    close(g_r_File)
    ############################################################ CAVITY PROBABILITY #####################################
    Pc_File = open("$Output_Route/CavityProbability.dat", "w");
    println(Pc_File, "Density\tCavityProbability")
    Pc_Array = zeros(Float64, length(Pc) - 1)
    @inbounds for i in keys(Pc)
        if i != 0
            Pc_Array[i] = Pc[i];
        end
    end
    D_Array = zeros(Float64, length(Pc_Array))
    for i = 1:length(Pc_Array)
        D_Array[i] = i / V;
        println(Pc_File, "$(round(D_Array[i], digits = 6))\t$(round(Pc_Array[i], digits = 6))")
    end
    close(Pc_File)
    Pc_Plot = plot(D_Array, Pc_Array, xlabel = L"\rho^*", ylabel = "Cavity Probability", ylims = (0, 1), yticks = 0:0.1:1, width = 5, label = "Random", guidefontsize = 18, tickfontsize = 14, legendfontsize = 18, foreground_color_legend = false, background_color_legend = false, bottom_margin = 10mm, left_margin = 8mm, widen = true, dpi = 300, size = [1920, 1080])
    savefig(Pc_Plot, "$Output_Route/CavityProbability")

    ###################################################### PLOTS ##############################################
    Radial_Distribution_Plot = plot(r, g_r, xlabel = L"r^*", ylabel = L"g(r^*)", xlim = (0, 3), xticks = 0:0.5:3, width = 3, guidefontsize = 20, tickfontsize = 18, bottom_margin = 7mm, left_margin = 10mm, right_margin = 3mm, legend = false, size = [1920, 1080], dpi = 300)
    hline!([1.0], color = :black, width = 2, linestyle = :dash)
    savefig(Radial_Distribution_Plot, "$Output_Route/RadialDistribution")

    Energy_Plot = plot(Energy_Array, xlabel = "N. Measurements", ylabel = L"U^* / N", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!(Energy_Plot, [Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    Energy_Histogram =  histogram(Energy_Array, normalize = true, xlabel = L"U^* / N", ylabel = "Normalized Frequency", legend = false, framestyle = :box, bins = 20, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(Mean_Energy_Array[N_Measurements], STD_Energy_Array[N_Measurements]), width = 3, linecolor = :black)

    Mean_Energy_Plot = plot(Mean_Energy_Array, ribbon = STD_Energy_Array, xlabel = "N. Measurements", ylabel = L"\langle U^* / N \rangle", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    Energy_Plots = plot(Energy_Plot, Energy_Histogram, Mean_Energy_Plot, layout = (@layout [a{0.3h} ; b c]))
    savefig(Energy_Plots, "$Output_Route/Energy_Plots")

    Density_Plot = plot(Density_Array, xlabel = "N. Measurements", ylabel = L"\rho^*", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!(Density_Plot, [Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    Density_Histogram =  histogram(Density_Array, normalize = true, xlabel = L"\rho^*", ylabel = " Normalized Frequency", legend = false, framestyle = :box, bins = 20, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(Mean_Density_Array[N_Measurements], STD_Density_Array[N_Measurements]), width = 3, linecolor = :black)

    Mean_Density_Plot = plot(Mean_Density_Array, ribbon = STD_Density_Array, xlabel = "N. Measurements", ylabel = L"\langle \rho^* \rangle", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    Density_Plots = plot(Density_Plot, Density_Histogram, Mean_Density_Plot, layout = (@layout [a{0.3h} ; b c]))
    savefig(Density_Plots, "$Output_Route/Density_Plots")
end

function Movement(L::Type, Beta::Type, Displacement::Float64, Energy::Float64, N_Displacement_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x));
    Energy_Old = Energy_Calculation(L, x[i], y[i], z[i], x, y, z);
    x_Old, y_Old, z_Old = x[i], y[i], z[i];
    x[i] += Displacement * (rand() - 0.5);
    y[i] += Displacement * (rand() - 0.5);
    z[i] += Displacement * (rand() - 0.5);
    x[i] = PeriodicBoundaryConditions!(L, x[i]);
    y[i] = PeriodicBoundaryConditions!(L, y[i]);
    z[i] = PeriodicBoundaryConditions!(L, z[i]);
    Energy_New = Energy_Calculation(L, x[i], y[i], z[i], x, y, z);
    Delta_E = Energy_New - Energy_Old;
    if rand() < exp(-Beta * Delta_E)
        N_Displacement_Accepted += 1;
        Energy += Delta_E;
    else
        x[i], y[i], z[i] = x_Old, y_Old, z_Old;
    end
    return Energy, N_Displacement_Accepted
end

function Insertion(L::Type, μ::Type, Beta::Type, Energy::Float64, N_Insertion_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    x_Insertion, y_Insertion, z_Insertion = L * (rand() - 0.5), L * (rand() - 0.5), L * (rand() - 0.5);
    Energy_Insertion = Energy_Calculation(L, x_Insertion, y_Insertion, z_Insertion, x, y, z);
    if rand() < exp( Beta * (μ - Energy_Insertion) + log(L^3. / (length(x) + 1)) )
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion)
        append!(y, y_Insertion)
        append!(z, z_Insertion)
        Energy += Energy_Insertion;
    end
    return Energy, N_Insertion_Accepted
end

function Insertion_Mezei(L::Type, μ::Type, Beta::Type, Energy::Float64, N_Insertion_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, x_Insertion::Array{Float64, 1}, y_Insertion::Array{Float64, 1}, z_Insertion::Array{Float64, 1}, Pc::Dict{Int64, Float64}, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x_Insertion))
    Energy_Insertion = Energy_Calculation(L, x_Insertion[i], y_Insertion[i], z_Insertion[i], x, y, z);
    if rand() < (L^3. * Pc[length(x)] / (length(x) + 1) ) * exp(Beta * (μ - Energy_Insertion))
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion[i])
        append!(y, y_Insertion[i])
        append!(z, z_Insertion[i])
        Energy += Energy_Insertion;
    end
    return Energy, N_Insertion_Accepted
end

function Removal(L::Type, μ::Type, Beta::Type, Energy::Float64, N_Removal_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x));
    Energy_Removal = Energy_Calculation(L, x[i], y[i], z[i], x, y, z);
    if rand() < exp( Beta * (Energy_Removal - μ) + log(length(x) / L^3.) )
        N_Removal_Accepted += 1;
        deleteat!(x, i)
        deleteat!(y, i)
        deleteat!(z, i)
        Energy -= Energy_Removal;
    end
    return Energy, N_Removal_Accepted
end

function Removal_Mezei(L::Type, μ::Type, Beta::Type, Energy::Float64, N_Removal_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, Pc_Interpolation::Float64, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x))
    Energy_Removal = Energy_Calculation(L, x[i], y[i], z[i], x, y, z);
    if rand() < ( length(x) / (L^3. * Pc_Interpolation) ) * exp(Beta * (Energy_Removal - μ))
        N_Removal_Accepted += 1;
        deleteat!(x, i)
        deleteat!(y, i)
        deleteat!(z, i)
        Energy -= Energy_Removal;
    end
    return Energy, N_Removal_Accepted
end

function Energy_Calculation(L::Type, rx::Float64, ry::Float64, rz::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    Energy = 0.;
    @inbounds for i = 1:length(x)
        Delta_x, Delta_y, Delta_z = rx - x[i], ry - y[i], rz - z[i];
        Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
        Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
        Delta_z = PeriodicBoundaryConditions!(L, Delta_z);
        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
        r != 0. ? Energy += U(r) : nothing
    end
    return Energy
end

function PeriodicBoundaryConditions!(L::Type, x::Float64) where {Type <: Real}
    return x - L * round(x / L)
end

function U(r::Float64, σ::Type = 1., λ::Type = 1.5, e::Type = 1.) where {Type <: Real}
    r <= σ ? (return Inf) : r <= λ ? (return -e) : (return 0)
end

function Random_Excluded_Volume(Equilibrium::Bool, σ::Type, L::Type, Pc::Dict{Int64, Float64}, Pc_Sum::Dict{Int64, Float64}, Pc_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}) where {Type <: Real}
    Equilibrium ? N_Random = 200 : N_Random = 1000
    if !haskey(Pc, length(x))
        Equilibrium = false;
        N_Random = 1000;
    end
    N_in = 0;
    x_Insertion, y_Insertion, z_Insertion = Float64[], Float64[], Float64[];
    if length(x) == 1
        for i = 1:N_Random
            x_V, y_V, z_V = L * (rand() - 0.5), L * (rand() - 0.5), L * (rand() - 0.5);
            N_in += 1;
            append!(x_Insertion, x_V)
            append!(y_Insertion, y_V)
            append!(z_Insertion, z_V)
        end
    else
        @inbounds for i = 1:N_Random
            x_V, y_V, z_V = L * (rand() - 0.5), L * (rand() - 0.5), L * (rand() - 0.5);
            for j = 1:length(x)
                Delta_x, Delta_y, Delta_z = x_V - x[j], y_V - y[j], z_V - z[j];
                Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
                Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
                Delta_z = PeriodicBoundaryConditions!(L, Delta_z);
                r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                r <= σ ? break : nothing
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
    if !Equilibrium
        if !haskey(Pc, length(x))
            Pc_N[length(x)] = 1;
            Pc_Sum[length(x)] = Volume_Ratio;
            Pc[length(x)] = Volume_Ratio;
        else
            Pc_N[length(x)] += 1;
            Pc_Sum[length(x)] += Volume_Ratio;
            Pc[length(x)] = Pc_Sum[length(x)] / Pc_N[length(x)];
        end
    end
    return Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion
end

function Interpolation(Pc::Dict{Int64, Float64}, l::Type) where {Type <: Real}
    if length(Pc) == 1
        Pc_Interpolation = collect(keys(Pc))[1];
    elseif length(Pc) > 1
        lim_inf = l - 1;
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

function RadialDistributionFunction(N_Bins::Int64, L::Type, Density::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    Delta = R_Cut / N_Bins;
    g_r = zeros(Float64, N_Bins);
    @inbounds for i = 1:length(x) - 1, j = i + 1:length(x)
        Delta_x, Delta_y, Delta_z = x[i] - x[j], y[i] - y[j], z[i] - z[j];
        Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
        Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
        Delta_z = PeriodicBoundaryConditions!(L, Delta_z);
        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
        if r <= R_Cut
            l = convert(Int64, ceil(r / Delta));
            g_r[l] += 2;
        end
    end
    @inbounds for l = 1:N_Bins
        g_r[l] /= (length(x) * 4 * pi * (l * Delta)^2 * Delta * Density)
    end
    return g_r
end

@time GrandCanonical_MonteCarlo(parse(Float64, ARGS[1]), parse(Float64, ARGS[2]))