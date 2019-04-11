% Models a PV-panel without bypass diodes under various shading conditions
% and all possible SP connection patterns when all cells are connected and
% when at least all shaded (non-ideal) cells are disconnected.
% Results are stored in tab separated textfiles (.txt) and as IV- and 
% PV-graphs in image form (.png)
% For a full explanation of the working of this algorithm please refer to
% my master thesis available at ....
%
% Author: Alexander Eick
% e-mail: alexander.eick@outlook.de
% Release: 2
% Release date: 17/11/2017

% This code is publish under the GNU GPL Version 3
% When reusing the code please mention the source


tic

addpath(genpath('export_fig'))
addpath(genpath('partitions'))

clc
clear

%% Parameters, Kyocera KD140SX-UFBS (single cell)%%
Nc = 36; %number of PV cells
Ns = 36; %number of PV cells in series
Np = 1; %number of PV cells in parallel
PanelHori = 4;  %number of cells horizontally in PV panel grid
PanelVert = 9;  %number of cells vertically in PV panel grid
                %horizontal * vertical must be Nc!!!
Iscr = 8.68; %ref short-circuit current (A)
Vocr = 0.614; %ref open-circuit voltage (V)
Imr = 7.91; %ref MPP current (A)
Vmr = 0.492; %ref MPP voltage (V)
coef_Iscr = 5.2e-3; %temp coeficent of Iscr (A/deg C)
coef_Vocr = -2.211e-3; %temp coeficient of Vocr (V/deg C)
NOCT = 45; %Nominal Operating Cell Tempeature (deg C)
Tr = 25; %ref cell temperaure (deg C)
Gr = 1000; %ref Irradiance (W/m2)
%% Variables %%
simulation = [  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 1  1  1  1  1  1  1  1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5;
                .5 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                .5 .5 1  1  1  1  1  1  1  .5 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                .5 .5 .5 1  1  1  1  1  1  .5 .5 1  1  1  1  1  1  1  .5 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                .5 .5 .5 .5 1  1  1  1  1  .5 .5 .5 1  1  1  1  1  1  .5 .5 1  1  1  1  1  1  1  .5 1  1  1  1  1  1  1  1;
                .5 .5 .5 .5 .5 1  1  1  1  .5 .5 .5 .5 1  1  1  1  1  .5 .5 .5 1  1  1  1  1  1  .5 .5 1  1  1  1  1  1  1;
                .5 .5 .5 .5 .5 .5 1  1  1  .5 .5 .5 .5 .5 1  1  1  1  .5 .5 .5 .5 1  1  1  1  1  .5 .5 .5 1  1  1  1  1  1;
                .5 .5 .5 .5 .5 .5 .5 1  1  .5 .5 .5 .5 .5 .5 1  1  1  .5 .5 .5 .5 .5 1  1  1  1  .5 .5 .5 .5 1  1  1  1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 .5 .5 .5 .5 .5 1  1  .5 .5 .5 .5 .5 .5 1  1  1  .5 .5 .5 .5 .5 1  1  1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 .5 .5 .5 .5 .5 1  1  .5 .5 .5 .5 .5 .5 1  1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 .5 .5 .5 .5 .5 1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5;
                .5 1  1  1  1  1  1  1  1  .5 1  1  1  1  1  1  1  1  .5 1  1  1  1  1  1  1  1  .5 1  1  1  1  1  1  1  1;
                .5 .5 1  1  1  1  1  1  1  .5 .5 1  1  1  1  1  1  1  .5 .5 1  1  1  1  1  1  1  .5 .5 1  1  1  1  1  1  1;
                .5 .5 .5 1  1  1  1  1  1  .5 .5 .5 1  1  1  1  1  1  .5 .5 .5 1  1  1  1  1  1  .5 .5 .5 1  1  1  1  1  1;
                .5 .5 .5 .5 1  1  1  1  1  .5 .5 .5 .5 1  1  1  1  1  .5 .5 .5 .5 1  1  1  1  1  .5 .5 .5 .5 1  1  1  1  1;
                .5 .5 .5 .5 .5 1  1  1  1  .5 .5 .5 .5 .5 1  1  1  1  .5 .5 .5 .5 .5 1  1  1  1  .5 .5 .5 .5 .5 1  1  1  1;
                .5 .5 .5 .5 .5 .5 1  1  1  .5 .5 .5 .5 .5 .5 1  1  1  .5 .5 .5 .5 .5 .5 1  1  1  .5 .5 .5 .5 .5 .5 1  1  1;
                .5 .5 .5 .5 .5 .5 .5 1  1  .5 .5 .5 .5 .5 .5 .5 1  1  .5 .5 .5 .5 .5 .5 .5 1  1  .5 .5 .5 .5 .5 .5 .5 1  1;
                .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 .5 .5 .5 .5 .5 .5 1;
                .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5];
simulation = simulation.';
ovName = 'results/results_v11_0_overall.txt'; 
ovFile = fopen(ovName,'a');
for si = 1:size(simulation,2)%si --> simulation counter
    si
    clearvars -except Nc Np Ns PanelHori PanelVert Iscr Vocr Imr Vmr coef_Iscr coef_Vocr NOCT Tr Gr simulation si ovFile
    margin = 0.2;%margin of Vm to determine best configuration after shading
    shade = 1;%500 - half shaded; 1 - fully shaded
    shade2 = 500;
    fName = sprintf('results/results_v11_%d',si);
    Parameters = zeros(Nc,14);  %initialize with 0 values %% [G, Ta, Tc, Vt,
                                %Isc, Voc, Im, Vm, I0, Rs, Neighbor north,
                                %Neighbor east, Neighbor south, Neighbor west]
    for i = 1:Nc
        Parameters(i,1) = simulation(i,si)*1000; %set G to 1000 (full sun)
        Parameters(i,2) = 25; %set Ta to 25
    end

    %% Newton Raphson %%
    array_max = 0.7; %maximum V until which each cell is evaluated
    array_inc = 1e-3; %I-V curve evaluation steps
    array_size = int16(array_max/array_inc+1);
    V = zeros(Nc,array_size); %initialize with 0 values %% voltage
    I = zeros(Nc,array_size); %initialize with 0 values %% current


    for i = 1:Nc
        V(i,:) = 0:array_inc:array_max;
        [   Parameters(i,3), Parameters(i,4), Parameters(i,5),...
            Parameters(i,6), Parameters(i,7), Parameters(i,8),...
            Parameters(i,9), Parameters(i,10)]...
            =   calc(1, 1, Iscr, Vocr, Imr, Vmr, coef_Iscr, coef_Vocr,...
            NOCT, Tr, Gr, Parameters(i,2), Parameters(i,1));
        sizeV = size(V);
        I(i,1) =    newton_raphson_method_v2(Parameters(i,5), V(i,1), Parameters(i,5),...
            Parameters(i,9), Parameters(i,10), 1, Parameters(i,4));
        for j = 2:sizeV(2)
            I(i,j) =    newton_raphson_method_v2(I(i,j-1), V(i,j),  Parameters(i,5),...
                Parameters(i,9), Parameters(i,10), 1, Parameters(i,4));
        end

        % fill Neighbors matrix %
        Nnorth = i - 1;
        Neast = i + PanelVert;
        Nsouth = i + 1;
        Nwest = i - PanelVert;

        if Nnorth/PanelVert == ceil(Nnorth/PanelVert)
            %set north neighbor to 0 if cell is in the uppermost row
            Parameters(i,11) = 0;
        else if Nnorth > 0
                %set north neighbor
                Parameters(i,11) = Nnorth;
            end
        end

        if Neast <= Nc
            %set east neighbor
            Parameters(i,12) = Neast;
        else
            %set east neighbor to 0 if cell is out of bounds
            Parameters(i,12) = 0;
        end

        if i/PanelVert == ceil(i/PanelVert)
            %set south neighbor to 0 if cell is in the lowermost row
            Parameters(i,13) = 0;
        else if Nsouth <= Nc
                %set south neighbor
                Parameters(i,13) = Nsouth;
            end
        end

        if Nwest > 0
            %set west neighbor
            Parameters(i,14) = Nwest;
        else
            %set west neighbor to 0 if cell is out of bounds
            Parameters(i,14) = 0;
        end
    end    
    %% calculate same number partitions of 1:Nc %%
    %  https://en.wikipedia.org/wiki/Partition_(number_theory)

    %  using code written by John D'Errico 
    %  https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer

    clist = cell(1,Nc);
    maxpos = 0;

    for h=1:Nc
        plist = partitions(h);%find all partitions for 1 --> Nc
        conlist = zeros(size(plist));
        con_i = 0;

        for i=1:size(plist,1)%search for lines in matrix which only have 1 
                             %non-zero value, hence compose the number of 
                             %same number partitions
            k = 0;
            for j=1:size(plist,2)
                if plist(i,j) ~= 0
                    k = k + 1;
                end
                if k > 1%break if more than 1 non-zero entry is found (save time)
                    break
                end
            end
            if k == 1%if 1 non-zero entry is found save this line
                con_i = con_i + 1;
                conlist(con_i,:) = plist(i,:);
            end
        end
        conlist = conlist(1:con_i,:);%crop matrix to same number partitions posibilities
        maxpos = maxpos + con_i;%count posibilities of all numbers
        clist{h} = conlist;%store matrix in list for later use
    end

    clear con_i conlist h i j k plist %clear unused variables to liberate memory

    %% generate I-V curves w/o excluding shaded cell %%

    i_gen = 0;
    OutputIV = zeros(3*maxpos,array_size);%output matrix for I-V curves to write to file
    OutputKey = zeros(3*maxpos,8);%output matrix for key parameters to write to file [Nc,Ns,Np,Isc,Voc,Pm,Im,Vm]
    IVFileWrite = '';%initialize string to define writing in file

    for i = 1:size(clist{Nc},1)
        i_gen = i_gen + 1;
        for j=1:size(clist{Nc},2)
            if clist{Nc}(i,j) ~= 0
                Np = clist{Nc}(i,j);
                Ns = j;
                break
            end
        end
        Vtemp = zeros(Np,array_size);%initialize with 0 values %% temporary voltage matrix while assembling PV panel
        Itemp = Inf(Np,array_size,'like',Vtemp);%initialize with infinite values %% temporary current matrix while assembling PV panel
        Ifinal = zeros(1,array_size);%initialize with 0 values %% final voltage array of assembled PV panel in current Np-Ns configuration
        Vfinal = Inf(1,array_size,'like',Ifinal);%initialize with infinite values %% final current array of assembled PV panel in current Np-Ns configuration

        OutputKey(2*i_gen-1,1) = Nc;%store Nc in output matrix
        OutputKey(2*i_gen-1,2) = Ns;%store Ns in output matrix
        OutputKey(2*i_gen-1,3) = Np;%store Np in output matrix

        for j = 1:Np
            for k = 1:Ns
                Vtemp(j,:) = Vtemp(j,:) + V(((j-1)*Ns)+k,:);%add voltages of series connected cells
                if I(((j-1)*Ns)+k,1) < Itemp(j,1)%find smallest current of series connected cells
                    Itemp(j,:) = I(((j-1)*Ns)+k,:);
                end
            end
            Ifinal = Ifinal + Itemp(j,:);%add currents of parallel connected strings
            if Vtemp(j,1) < Vfinal(1)%find smallest voltage of parallel connected strings
                Vfinal = Vtemp(j,:);
            end
        end

        OutputIV(2*i_gen-1,:) = Vfinal;%store final V array in output matrix
        OutputIV(2*i_gen,:) = Ifinal;%store final I array in output matrix

        OutputKey(2*i_gen-1,4) = Ifinal(1);%find Isc of PV-panel

        I_root = [0,Inf];
        for j = 1:size(Ifinal,2)
            if abs(Ifinal(j)) < I_root(2)%find where I-V curve crosses x-axis
                I_root(1) = j;
                I_root(2) = Ifinal(j);
            end
        end
        OutputKey(2*i_gen-1,5) = Vfinal(I_root(1));%find Voc

        Pfinal = Ifinal .* Vfinal;
        OutputKey(2*i_gen-1,6) = max(Pfinal);%find Pm

        Pm_pos = find(Pfinal==OutputKey(2*i_gen-1,6));%array position of Pm
        OutputKey(2*i_gen-1,7) = Ifinal(Pm_pos);%find Im
        OutputKey(2*i_gen-1,8) = Vfinal(Pm_pos);%find Vm

        IVFileWrite = strcat(IVFileWrite,'%.3f\t%.5f\t');%concat string to write in file with extra variables
    end

    %% generating I-V curves w/ excluding shaded cell %%
    ShadedPos = find(Parameters(:,1)<1000);%find indexes of shaded cells
    Vex = V; %voltage matrix w/o shaded cells
    Iex = I; %current matrix w/o shaded cells
    for j = size(ShadedPos,1):-1:1%exclude shaded cells
        Vex = [Vex(1:ShadedPos(j)-1,:); Vex(ShadedPos(j)+1:Nc-size(ShadedPos,1)+j,:)];
        Iex = [Iex(1:ShadedPos(j)-1,:); Iex(ShadedPos(j)+1:Nc-size(ShadedPos,1)+j,:)];
    end

    for h = Nc-size(ShadedPos,1):-1:1
        for i = 1:size(clist{h},1)
            i_gen = i_gen + 1;
            for j=1:size(clist{h},2)
                if clist{h}(i,j) ~= 0
                    Np = clist{h}(i,j);
                    Ns = j;
                    break
                end
            end
            Vtemp = zeros(Np,array_size);
            Itemp = Inf(Np,array_size,'like',Vtemp);
            Ifinal = zeros(1,array_size);
            Vfinal = Inf(1,array_size,'like',Ifinal);

            OutputKey(2*i_gen-1,1) = h;%store Nc in output matrix
            OutputKey(2*i_gen-1,2) = Ns;%store Ns in output matrix
            OutputKey(2*i_gen-1,3) = Np;%store Np in output matrix

            for j = 1:Np
                for k = 1:Ns
                    Vtemp(j,:) = Vtemp(j,:) + Vex(((j-1)*Ns)+k,:);%add voltages of series connected cells
                    if Iex(((j-1)*Ns)+k,1) < Itemp(j,1)%find smallest current of series connected cells
                        Itemp(j,:) = Iex(((j-1)*Ns)+k,:);
                    end
                end
                Ifinal = Ifinal + Itemp(j,:);%add currents of parallel connected strings
                if Vtemp(j,1) < Vfinal(1)%find smallest volatge of parallel connected stirngs
                    Vfinal = Vtemp(j,:);
                end
            end

            OutputIV(2*i_gen-1,:) = Vfinal;%store final V array in output matrix
            OutputIV(2*i_gen,:) = Ifinal;%store final I array in output matrix

            OutputKey(2*i_gen-1,4) = Ifinal(1);%find Isc of PV-panel

            size_Ifinal = size(Ifinal);
            I_root = [0,Inf];
            for j = 1:size_Ifinal(2)
                if abs(Ifinal(j)) < I_root(2)%find where I-V curve crosses x-axis
                    I_root(1) = j;
                    I_root(2) = Ifinal(j);
                end
            end
            OutputKey(2*i_gen-1,5) = Vfinal(I_root(1));%find Voc

            Pfinal = Ifinal .* Vfinal;
            OutputKey(2*i_gen-1,6) = max(Pfinal);%find Pm

            Pm_pos = find(Pfinal==OutputKey(2*i_gen-1,6));%array position of Pm
            OutputKey(2*i_gen-1,7) = Ifinal(Pm_pos);%find Im
            OutputKey(2*i_gen-1,8) = Vfinal(Pm_pos);%find Vm

            IVFileWrite = strcat(IVFileWrite,'%.3f\t%.5f\t');%concat string to write in file with extra variables
        end
    end

    %% crop final matrix and store IV curve data in file %%
    OutputIV = OutputIV(1:2*i_gen,:);%crop matrix
    OutputKey = OutputKey(1:2*i_gen,:);
    Output = horzcat(OutputKey, OutputIV);
    IVFileWrite = strcat(IVFileWrite,'\n');

    fileID = fopen(sprintf('%s_all.txt',fName),'w');
    fprintf(fileID,IVFileWrite,Output);
    fclose(fileID);

    %% analyze IV curves %%

    %% import standard IV curve %%
    sIV = importdata('results/standard IV curve.txt');%import standard IV data for fully iluminated panel with Np=3, Ns=12
    sNs = sIV(2,1);
    sNp = sIV(3,1);
    sIsc = sIV(4,1);
    sVoc = sIV(5,1);
    sPm = sIV(6,1);
    sIm = sIV(7,1);
    sVm = sIV(8,1);
    sIV = sIV(10:end,:);

    lower_Vm_limit = sVm - sVm * margin;%set lower limit for Vm
    upper_Vm_limit = sVm + sVm * margin;%set upper limit for Vm
    lower_Pm_limit = sPm - sPm * margin;%set lower limit for Pm
    upper_Pm_limit = sPm + sPm * margin;%set upper limit for Pm

    %% create arrays to draw limits within graph %%
    lower_lim = zeros(701,2);
    lower_lim(:,1) = lower_Vm_limit;
    lower_lim(:,2) = 0:1:700;
    upper_lim = zeros(701,2);
    upper_lim(:,1) = upper_Vm_limit;
    upper_lim(:,2) = 0:1:700;
    
    lower_P_lim = zeros(701,2);
    lower_P_lim(:,2) = lower_Pm_limit;
    lower_P_lim(:,1) = 0:1:700;
    upper_P_lim = zeros(701,2);
    upper_P_lim(:,2) = upper_Pm_limit;
    upper_P_lim(:,1) = 0:1:700;

    %% extract IV curves within the Vm limits %%
    OutputWithin = zeros(size(Output));
    IVwithinFileWrite = '';
    i_oW = 1;%counter for IV within limits (i_outputWithin)
    for i_ext = 1:2:size(Output,1)
        if Output(i_ext,8) > lower_Vm_limit && Output(i_ext,8) < upper_Vm_limit
            OutputWithin(i_oW:i_oW+1,:) = Output(i_ext:i_ext+1,:);
            IVwithinFileWrite = strcat(IVwithinFileWrite,'%.3f\t%.5f\t');
            i_oW = i_oW + 2;
        end
    end
    OutputWithin = OutputWithin(1:i_oW-1,:);%crop matrix
    IVwithinFileWrite = strcat(IVwithinFileWrite,'\n');
    
    %% extract IV curves within the Pm limits %%
    OutputWithinP = zeros(size(OutputWithin));
    IVwithinPFileWrite = '';
    i_oWP = 1;
    for i_ext = 1:2:size(OutputWithin,1)
        if OutputWithin(i_ext,6) > lower_Pm_limit && OutputWithin(i_ext,6) < upper_Pm_limit
            OutputWithinP(i_oWP:i_oWP+1,:) = OutputWithin(i_ext:i_ext+1,:);
            IVwithinPFileWrite = strcat(IVwithinPFileWrite,'%.3f\t%.5f\t');
            i_oWP = i_oWP + 2;
        end
    end
    if i_oWP ~= 1
        OutputWithinP = OutputWithinP(1:i_oWP-1,:);%crop matrix to actual size
    else
        OutputWithinP = OutputWithinP(1:i_oWP+1,:);%crop matrix to only 2 lines --> graph drawing will not fail
    end
    IVwithinPFileWrite = strcat(IVwithinPFileWrite,'\n');

    %% export figure and data %%
    figure(1)%plot I-V curves
    clf
    p2 = plot(Output(1,9:end),Output(2,9:end),'-b');
    hold on
    for i = 3:2:size(Output,1)
        plot(Output(i,9:end),Output(i+1,9:end),'-b');
    end
    p3 = plot(OutputWithin(1,9:end),OutputWithin(2,9:end),'-g');
    plot(OutputWithin(1,8),OutputWithin(1,7),'*g');
    for i = 3:2:i_oW-1
        plot(OutputWithin(i,9:end),OutputWithin(i+1,9:end),'-g');
        plot(OutputWithin(i,8),OutputWithin(i,7),'*g');
    end
    p6 = plot(OutputWithin(1,9:end),OutputWithin(2,9:end));
    p7 = plot(OutputWithin(1,8),OutputWithin(1,7));
    set(p6,'Color',[1 0.65 0]);
    set(p7,'MarkerEdgeColor',[1 0.65 0],'MarkerFaceColor',[1 0.65 0],'Marker','*');
    p1 = plot(sIV(:,1),sIV(:,2),'-r');
    plot(sVm,sIm,'*r');
    p4 = plot(lower_lim(:,1),lower_lim(:,2),'-m');
    p5 = plot(upper_lim(:,1),upper_lim(:,2),'-m');

    legend( [p1 p4 p3 p2 p6] , ...
        'standard configuration of Ns=12 & Np=3' , ...
        'lower & upper margin to standard Vm' , ...
        'configurations with Vm within the margin' , ...
        'configurations with Vm outside of the margin' , ...
        'shaded standard configuration', ...
        'location','northeast');
    ylabel('Current (A)')
    xlabel('Voltage (V)')
    axis([0, 20, 0, 350])
    set(gcf, 'Color', 'w');
    export_fig(sprintf('%s_IV--0-20.png',fName),'-r350');%export figure
    axis([0, 8, 0, 40])
    set(gcf, 'Color', 'w');
    export_fig(sprintf('%s_IV--0-8.png',fName),'-r350');%export figure
    
    figure(2)%plot P-V curves
    clf
    p2 = plot(Output(1,9:end),Output(2,9:end).*Output(1,9:end),'-b');
    hold on
    for i_fig = 3:2:size(Output,1)
        plot(Output(i_fig,9:end),Output(i_fig+1,9:end).*Output(i_fig,9:end),'-b');
    end
    p3 = plot(OutputWithin(1,9:end),OutputWithin(2,9:end).*OutputWithin(1,9:end),'-g');
    plot(OutputWithin(1,8),OutputWithin(1,7).*OutputWithin(1,8),'*g');
    for i_fig = 3:2:size(OutputWithin,1)
        plot(OutputWithin(i_fig,9:end),OutputWithin(i_fig+1,9:end).*OutputWithin(i_fig,9:end),'-g');
        plot(OutputWithin(i_fig,8),OutputWithin(i_fig,7).*OutputWithin(i_fig,8),'*g');
    end
    p8 = plot(OutputWithinP(1,9:end),OutputWithinP(2,9:end).*OutputWithinP(1,9:end),'-c');
    if OutputWithinP(1,1) ~= 0
        plot(OutputWithinP(1,8),OutputWithinP(1,6),'*c');
        for i_fig = 3:2:size(OutputWithinP,1)
            plot(OutputWithinP(i_fig,9:end),OutputWithinP(i_fig+1,9:end).*OutputWithinP(i_fig,9:end),'-c');
            plot(OutputWithinP(i_fig,8),OutputWithinP(i_fig,6),'*c');
        end
    end
    p9 = plot(OutputWithin(1,9:end),OutputWithin(2,9:end).*OutputWithin(1,9:end));
    p10 = plot(OutputWithin(1,8),OutputWithin(1,7).*OutputWithin(1,8));
    set(p9,'Color',[1 0.65 0]);
    set(p10,'MarkerEdgeColor',[1 0.65 0],'MarkerFaceColor',[1 0.65 0],'Marker','*');
    p1 = plot(sIV(:,1),sIV(:,2).*sIV(:,1),'-r');
    plot(sVm,sIm*sVm,'*r');
    p4 = plot(lower_lim(:,1),lower_lim(:,2),'-m');
    p5 = plot(upper_lim(:,1),upper_lim(:,2),'-m');
    p6 = plot(lower_P_lim(:,1),lower_P_lim(:,2),'-m');
    p7 = plot(upper_P_lim(:,1),upper_P_lim(:,2),'-m');
    
    legend( [p1 p4 p3 p8 p2 p9] , ...
        'standard configuration of Ns=12 & Np=3' , ...
        'lower & upper margins to standard Vm & Pm' , ...
        'configurations with Vm within the margin' , ...
        'configurations with Pm within the margin' , ...
        'configurations with Vm outside of the margin' , ...
        'shaded standard configuration', ...
        'location','northeast');
    ylabel('Power (W)')
    xlabel('Voltage (V)')
    axis([0, 20, 0, 180])
    set(gcf, 'Color', 'w');
    export_fig(sprintf('%s_PV--0-20.png',fName),'-r350');%export figure
    axis([0, 8, 0, 180])
    set(gcf, 'Color', 'w');
    export_fig(sprintf('%s_PV--0-8.png',fName),'-r350');%export figure
    
    %print data to file
    fileID = fopen(sprintf('%s_within.txt',fName),'w');
    fprintf(fileID,IVwithinFileWrite,OutputWithin);
    fclose(fileID);
    
    fileID = fopen(sprintf('%s_withinP.txt',fName),'w');
    fprintf(fileID,IVwithinPFileWrite,OutputWithinP);
    fclose(fileID);

    %% determine new size %%
    newNs = 1;
    newNp = 1;
    connections = Inf(newNp,newNs);
    prev_P = Inf;
    while (newNs ~= sNs || newNp ~= sNp) && ismember(Inf,connections)
        best_P = 0;
        best_Pindex = 0;
        for i = 1:2:i_oW-1
            if OutputWithin(i,6) > best_P && OutputWithin(i,6) < prev_P
                best_P = OutputWithin(i,6);
                best_Pindex = i;
            end
        end
        prev_P = best_P;
        newNs = OutputWithin(best_Pindex,2);
        newNp = OutputWithin(best_Pindex,3);
        fprintf(ovFile,'%d\t%d\t%d\n',si,newNp,newNs);

        %% find connection pattern %%
        if newNs ~= sNs || newNp ~= sNp
            connections = Inf(newNp,newNs);
            blacklist = cell(newNp,newNs);
            blacklist(:) = {zeros(1,Nc)};
            i_fc = 1;
            i_nc = 2;
            i_nei = 11;%positions of neighbor information in parameter array (11-14)
            i_Np = 1;

            while i_Np <= newNp
                i_fc = 1;
                while i_fc <= Nc%find first cell without shade
                    if      Parameters(i_fc,1) == 1000 && ...%check if cell is fully iluminated
                            ~ismember(i_fc,connections) && ...%check if it is already allocated
                            ~ismember(i_fc,blacklist{i_Np,1})%check if it is blacklisted for this position
                        connections(i_Np,1) = i_fc;
                    else %if cell is not suited increase counter
                        i_fc = i_fc + 1;
                    end
                    if connections(i_Np, 1) ~= Inf %execute only if suitable first cell is found
                        i_nc = 2;
                        while i_nc <= newNs %find other cell connections until number of series connected cells is reached
                            i_nei = 11;
                            while i_nei <= 14 %cycle through all neighboring cells (north, east, south, west)
                                if      Parameters(connections(i_Np, i_nc-1), i_nei) ~= 0 && ...%check if neighbor cell is 0
                                        Parameters(Parameters(connections(i_Np, i_nc-1), i_nei), 1) == 1000 && ...%check if neighbor cell is fully iluminated
                                        Parameters(connections(i_Np, i_nc-1), i_nei) < connections(i_Np, i_nc) && ...%check if neighbor cell smaller than current pick
                                        ~ismember(Parameters(connections(i_Np, i_nc-1), i_nei), connections) && ...%check if neighbor cell is already allocated
                                        ~ismember(Parameters(connections(i_Np, i_nc-1), i_nei),blacklist{i_Np, i_nc})%check if neighbor is on the blacklist for this position
                                    connections(i_Np, i_nc) = Parameters(connections(i_Np, i_nc-1), i_nei);
                                end
                                i_nei = i_nei + 1;
                            end
                            if connections(i_Np, i_nc) == Inf %if suitable cell for this position is not found
                                for i = 1:Nc
                                    if blacklist{i_Np, i_nc-1}(i) == 0
                                        blacklist{i_Np, i_nc-1}(i) = connections(i_Np, i_nc-1); %blacklist cell of previous position
                                        break
                                    end
                                end
                                connections(i_Np, i_nc-1) = Inf; %delete previous position cell
                                blacklist{i_Np, i_nc}(:) = 0; %delete current position blacklist
                                i_nc = i_nc - 1; %return to previous position
                                if i_nc == 1
                                    break %if i_nc = 1 stop this while loop
                                end
                            else
                                i_nc = i_nc + 1;
                            end
                        end
                    end
                    if i_nc == 1 && i_Np ~= 1 && i_fc == Nc 
                        %if not in the loop to determine the first string, delete last
                        %cell assignment of previous string if i_fc is at last cell
                        i_Np = i_Np - 1;
                        for i = 1:Nc
                            if blacklist{i_Np, newNs}(i) == 0
                                blacklist{i_Np, newNs}(i) = connections(i_Np, newNs);
                                break
                            end
                        end
                        connections(i_Np, newNs) = Inf; %delete last cell assignment of previous string
                        blacklist{i_Np+1, i_nc}(:) = 0; %delete current position's blacklist
                        i_fc = 1; %restart 'i_fc' while loop
                    end
                    if ~ismember(Inf, connections(i_Np, :)) %if all connection in this string are found
                        i_Np = i_Np + 1; %start filling the next parallel string
                        break %stop the 'i_fc' while loop
                    end
                end
                if nnz(blacklist{1, 1}) == newNp*newNs
                    break %finish loop if all non-shaded cells are blacklisted on the first position
                end
            end
        end
    end

    %% write pattern to file %%
    ConnFileWrite = '';
    for i = 1:newNp
        ConnFileWrite = strcat(ConnFileWrite, '%d\t');
    end
    ConnFileWrite = strcat(ConnFileWrite, '\n');

    fileID = fopen(sprintf('%s_connection.txt', fName),'w');
    fprintf(fileID,ConnFileWrite, connections);
    fclose(fileID);
end
fclose(ovFile);
toc;