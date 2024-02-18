@time using Distributions, NearestNeighbors, StatsBase, DelimitedFiles, ProgressBars, CSV


function MGE(state, District, Location, ε, ID, Num_Attendants, Event_Zone_Size, Micro_Time_Step, brownian, β, INCUBATION, RECOVERY, latent_period, 
    infectious_period, MGE_location, I_District, Num_Attendant_I)
    
        Event_Location = maximum(Location[:, 1]) + 2;
        
        I_MGE = min(Num_Attendant_I, sum(I_District));

        EVENT_ID_I = rand(ID[I_District], I_MGE);

        if MGE_location=="L"
            EVENT_ID_OTHERS = rand(ID[(state .!= 'I').&(District.=='L')], Num_Attendants-I_MGE);
        elseif MGE_location=="H"
            EVENT_ID_OTHERS = rand(ID[(state .!= 'I').&(District.=='H')], Num_Attendants-I_MGE);
        end
        
        EVENT_ID_S = intersect(EVENT_ID_OTHERS, ID[state.=='S']);
        Attendants_ID = cat(EVENT_ID_I, EVENT_ID_OTHERS, dims=1);
        Original_Location = copy(Location[Attendants_ID, :]);
        
        Location[Attendants_ID, :] = [Event_Location Event_Location].+rand(Uniform(0, Event_Zone_Size), Num_Attendants, 2);
        
        EVENT_EXPOSED_TOTAL = zeros(Int64, length(EVENT_ID_S));

        for i=1:Micro_Time_Step
            θ = rand(Uniform(0, 2π), Num_Attendants);
            Location[Attendants_ID, :] = [Event_Location Event_Location].+mod.(Location[Attendants_ID, :]+(rand(brownian, Num_Attendants)').*hcat(θ, θ), Event_Zone_Size);
            Idx = inrange(KDTree(transpose(Location[EVENT_ID_I, :])), transpose(Location[EVENT_ID_S, :]), ε, false);
            contact = length.(Idx);
            bit_infected = rand(length(EVENT_ID_S)) .< (1 .- (1 - β).^contact);
            EVENT_EXPOSED_TOTAL = EVENT_EXPOSED_TOTAL+bit_infected;            
        end

        EVENT_EXPOSED_ID = EVENT_ID_S[EVENT_EXPOSED_TOTAL.>0];
        if !isempty(EVENT_EXPOSED_ID)
            state[EVENT_EXPOSED_ID] .='E';
            INCUBATION[EVENT_EXPOSED_ID] .= round.(rand(latent_period, length(EVENT_EXPOSED_ID)));    
            RECOVERY[EVENT_EXPOSED_ID] .= INCUBATION[EVENT_EXPOSED_ID] + round.(rand(infectious_period, length(EVENT_EXPOSED_ID)));
        end              

        Location[Attendants_ID, :] = Original_Location;
    return Location, state, INCUBATION, RECOVERY;
end


function ABM(Pop_Cluster1, Pop_Cluster2, β, ε, L, End_date, seed_number, tₑ₁, dₑ₁, Lₑ₁, Iₑ₁, Rp, SD_level, MGE_level, Initial_Host, case, MGE_location)
   
    Nₑ₁ = MGE_level*Int(round(Pop_Cluster1*0.005));

    println("SD : $SD_level, MGE : $MGE_level, Movement rate : $Rp");

    n = Pop_Cluster1 + Pop_Cluster2 ; 

    ID = 1:n;
    
    brownian = MvNormal(2, 0.01);
    latent_period = Gamma(4.059, 1.3500); 
    infectious_period = Gamma(5, 1.4); 

    filename = "SD-$SD_level, MGE-$MGE_level, MR-$Rp, $case.txt"

    open(filename, "w") do file

        @time Threads.@threads for seed in ProgressBar(1:seed_number)
            S = [Pop_Cluster1+Pop_Cluster2-1];
            E = [0];
            I = [1];
            R = [0];
            Daily = [1];

            E1 = [0];
            I1 = [1];
            R1 = [0];
            D1 = [1];
            E2 = [0];
            I2 = [0];
            R2 = [0];
            D2 = [0];

            INCUBATION = zeros(Int64, n) .- 1;
            RECOVERY = zeros(Int64, n) .- 1;
            
            t = 0;
        
            state = Array{Char, 1}(undef, n); 
            state .= 'S'; 
        
            state[Initial_Host] = 'I';
            RECOVERY[Initial_Host] = round(rand(infectious_period, 1)[1]) + 1;

            Location = rand(Uniform(0, L), n, 2);
            
            District = Array{Char, 1}(undef, n); 
            District[1:Pop_Cluster1].='L';
            District[Pop_Cluster1+1:end].='H';

            if SD_level==0
                βₒ = β;
            elseif SD_level==1
                βₒ = 0.9β;
            elseif SD_level==2
                βₒ = 0.75β;
            elseif SD_level==3
                βₒ = 0.5β;
            end
            Rp1 = Rp;
            Rp2 = Rp/4;

            thd = 0.02;
            category = -999;

            while t<End_date
                t+=1;
                if sum(state .=='E') + sum(state .=='I') >0 
        
                    To_Cluster2 = rand(ID[(state.!='R').&(District.=='L')], Int(round(Pop_Cluster1*Rp1)));
                            
                    To_Cluster1 = rand(ID[(state.!='R').&(District.=='H')], Int(round(Pop_Cluster2*Rp2)));
                    
                    District[To_Cluster2].='H';
                    District[To_Cluster1].='L';
                    Location[To_Cluster2, :] = [2L 2L].+rand(Uniform(0, L), length(To_Cluster2), 2);
                    Location[To_Cluster1, :] = rand(Uniform(0, L), length(To_Cluster1), 2);        
                    
                    INCUBATION .-= 1;
                    RECOVERY .-= 1;
                    
                    state[(INCUBATION .<= 0).&(state .=='E')] .= 'I';
                    state[(RECOVERY .<= 0).&(state .=='I')] .= 'R';
                        
                    bit_S = (state .== 'S');
                    bit_E = (state .== 'E');
                    bit_I = (state .== 'I');
                    bit_R = (state .== 'R');
        
                    if t!=1
                        push!(S, sum(bit_S));
                        push!(E, sum(bit_E));
                        push!(I, sum(bit_I));
                        push!(Daily, sum(INCUBATION.==0));
                        push!(E1, sum(bit_E.&(District.=='L')));
                        push!(I1, sum(bit_I.&(District.=='L')));
                        push!(R1, sum(bit_R.&(District.=='L')));
                        push!(D1, sum((INCUBATION .== 0).&(District.=='L').&(state.=='I')));
                        push!(E2, sum(bit_E.&(District.=='H')));
                        push!(I2, sum(bit_I.&(District.=='H')));
                        push!(R2, sum(bit_R.&(District.=='H')));
                        push!(D2, sum((INCUBATION .== 0).&(District.=='H').&(state.=='I')));
                        push!(R, sum(bit_R.&(District.=='L'))+sum(bit_R.&(District.=='H')));
                    end
                
                    Movement_Idx = state.!='R';
                    r_L = rand(brownian, sum(Movement_Idx.&(District.=='L')))';
                    r_H = rand(brownian, sum(Movement_Idx.&(District.=='H')))';
                    θ_L = rand(Uniform(0, 2π), sum(Movement_Idx.&(District.=='L')));
                    θ_H = rand(Uniform(0, 2π), sum(Movement_Idx.&(District.=='H')));
                    
                    if MGE_location=="L"
                        I_in_Cluster = (state .=='I').&(District.=='L');
                    elseif MGE_location=="H"
                        I_in_Cluster = (state .=='I').&(District.=='H');
                    end

                    if t in tₑ₁ && sum(I_in_Cluster)>0 && Nₑ₁>0
                        Location, state, INCUBATION, RECOVERY = MGE(
                            state, District, Location, ε, ID, Nₑ₁, Lₑ₁, dₑ₁, brownian, β, 
                            INCUBATION, RECOVERY, latent_period, infectious_period, MGE_location, I_in_Cluster, Iₑ₁);
                    end

                    Location[Movement_Idx.&(District.=='L'), :] = mod.(Location[Movement_Idx.&(District.=='L'), :] + r_L.*hcat(θ_L, θ_L), L);
                    Location[Movement_Idx.&(District.=='H'), :] = [2L 2L].+mod.(Location[Movement_Idx.&(District.=='H'), :] + r_H.*hcat(θ_H, θ_H), L);            
                        
                    Contact1 = length.(inrange(KDTree(transpose(Location[ID[bit_I.&(District.=='L')],:])), transpose(Location[ID[bit_S.&(District.=='L')],:]), ε, false));
                    Contact2 = length.(inrange(KDTree(transpose(Location[ID[bit_I.&(District.=='H')],:])), transpose(Location[ID[bit_S.&(District.=='H')],:]), ε, false));
        
                    bit_infected_Cluster1 = rand(sum(bit_S.&(District.=='L'))) .< (1 .- (1 - βₒ).^Contact1);    
                    bit_infected_Cluster2 = rand(sum(bit_S.&(District.=='H'))) .< (1 .- (1 - βₒ).^Contact2);
                    
                    ID_S_Cluster1 = ID[bit_S.&(District.=='L')];
                    ID_S_Cluster2 = ID[bit_S.&(District.=='H')];
            
                    ID_infected_Cluster1 = ID_S_Cluster1[bit_infected_Cluster1];
                    ID_infected_Cluster2 = ID_S_Cluster2[bit_infected_Cluster2];

                    state[union(ID_infected_Cluster1, ID_infected_Cluster2)].='E';                
                    ID_infected = union(ID_infected_Cluster1, ID_infected_Cluster2);
                    INCUBATION[ID_infected] .= round.(rand(latent_period, sum(bit_infected_Cluster1)+sum(bit_infected_Cluster2)));
                    RECOVERY[ID_infected] .= INCUBATION[ID_infected] + round.(rand(infectious_period, sum(bit_infected_Cluster1)+sum(bit_infected_Cluster2)));
                end
            end

            if R1[end]≤thd*Pop_Cluster1 && R2[end]≤thd*Pop_Cluster2 
                category = 1;
            elseif (R1[end]>thd*Pop_Cluster1 && R2[end]≤thd*Pop_Cluster2)
                category = 2;
            elseif (R1[end]≤thd*Pop_Cluster1 && R2[end]>thd*Pop_Cluster2)
                category = 3;
            elseif (R1[end]>thd*Pop_Cluster1) && (R2[end]>thd*Pop_Cluster2)
                category = 4;
            end
            
            writedlm(file, Int64.([SD_level MGE_level category R1[end] R2[end]]))
        end        
    end 
end 
