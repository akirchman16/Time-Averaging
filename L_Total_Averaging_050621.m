clearvars;
close all;

% This code creates "average" curves for multiple values of the
% total concentration of RAD51, L_Total. Input different values for L_Total and
% the model will run 10 iterations of the simulation for each L_Total value and
% then create averages for each value. Included in this model is the ratio
% of monomers to dimers.

N = 8660;   %length of DNA lattice
n = 3;  %length of a monomer
w = 1;  %cooperativity parameter
Conc_Values = [0.1,1,10];    %total concentration of RAD51
k_on = 1;   %kinetic rate constants
k_off = 1;
Ratio = 1;   %Percentage of solution which is monomers (0 to 1)
AverageIterations = 10;    %number of iterations at each ratio value

UncoveredLength = 0.34; %length of a DNA nt without RAD51 bound to it (according to van der Heijden paper) - nm
CoveredLength = 0.51;   %length of a DNA nt where RAD51 is bound - nm

L_Prep = ones(1,AverageIterations);  %sorts all values for cooperativity from smallest to largest
Concentrations = [];
for f = 1:length(Conc_Values)
    Concentrations = sort([Concentrations,Conc_Values(f)*L_Prep]);
end

Colors = zeros(numel(Conc_Values),3);   %selects color to plot average curves in for each cooperativity value
for z = 1:numel(Conc_Values)
    Colors(z,:) = [rand,rand,rand]; %random color for each value of L_Total
end

minIterations = 1000;

%Memory Allocation
EventFractions = zeros(numel(Concentrations),7);
FracCover = zeros(numel(Concentrations),minIterations);
DNA_Lengths = zeros(numel(Concentrations),minIterations);
t = zeros(numel(Concentrations),minIterations);
Max_Time = zeros(1,numel(Concentrations));
Equilibrium_Coverage = zeros(1,numel(Concentrations));
L_Monomer = zeros(1,numel(Concentrations));
L_Dimer = zeros(1,numel(Concentrations));
Variable_List = zeros(1,AverageIterations*numel(Concentrations));

Loops = 0;
for L_Total = Concentrations
    Loops = Loops+1;
    Variable_List(Loops) = L_Total;    %lists the respective variable for each loop (helpful when averaging)
    
    L_Monomer(Loops) = Ratio*L_Total;        %Concentration of monomer RAD51
    L_Dimer(Loops) = (1-Ratio)*L_Total;    %Concentration of dimer RAD51
    
    DNA = zeros(1,N);
    BoundAtSpot = zeros(1,N);   %records where monomers are bound on lattice

    %Memroy Allocations
    Populations = zeros(minIterations,7);
    a = zeros(minIterations,7);
    Probabilities = zeros(minIterations,7);
    FiringAmounts = zeros(minIterations,7);
    dt = zeros(1,minIterations);
    j = zeros(1,minIterations);
    Location_History = zeros(7,minIterations);

    Equilibrium = 0;
    Events = 0;
    while max(t(Loops,:)) < 1.25
        Events = Events+1;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        searcH_vars = {'Left','Right','Gap_Size','Left_Available_M','Right_Available_M','Left_Available_D','Right_Available_D','Right_BindingSite_M','Right_BindingSite_D','Gaps_L2_M','Gaps_R2_M','Gaps_L2_D','Gaps_R2_D','Gap_Size2_M','Gap_Size2_D','Gap_SizeI_M','Gap_SizeI_D','Left_I_M','Left_I_D','Isolated_M','SinglyContiguous_M','Doubly_Contiguous_M','Isolated_D','SinglyContiguous_D','DoublyContiguous_D'};
        clear searcH_vars;

        Left = find(diff([1 DNA 1]) == -1);
        Right = find(diff([1 DNA 1]) == 1)-1;
        Gap_Size = Right-Left+1;

        Left_Available_M = Left(Gap_Size >= n);
        Right_Available_M = Right(Gap_Size >= n);
        Left_Available_D = Left(Gap_Size >= 2*n);
        Right_Available_D = Right(Gap_Size >= 2*n);
        Right_BindingSite_M = Right_Available_M-(n-1);
        Right_BindingSite_D = Right_Available_D-(2*n-1);

        %Doubly Contiguous Searches
        Doubly_Contiguous_M = Left_Available_M(Left_Available_M == Right_BindingSite_M & 1 < Left_Available_M & Left_Available_M < N-(n-1));
        Doubly_Contiguous_D = Left_Available_D(Left_Available_D == Right_BindingSite_D & 1 < Left_Available_D & Left_Available_D < N-(2*n-1));

        %Singly Contiguous Searches
        Singly_Contiguous_M = unique([Left_Available_M(~ismember(Left_Available_M,Doubly_Contiguous_M)),Right_Available_M(~ismember(Right_Available_M-(n-1),Doubly_Contiguous_M))-(n-1)]);
        Singly_Contiguous_M(Singly_Contiguous_M == N-(n-1) | Singly_Contiguous_M == 1) = [];
        Singly_Contiguous_D = unique([Left_Available_D(~ismember(Left_Available_D,Doubly_Contiguous_D)),Right_Available_D(~ismember(Right_Available_D-(2*n-1),Doubly_Contiguous_D))-(2*n-1)]);
        Singly_Contiguous_D(Singly_Contiguous_D == N-(2*n-1) | Singly_Contiguous_D == 1) = [];

        %Isolated Searches
        Gaps_L2_M = Left(~ismember(Left,Doubly_Contiguous_M));
        Gaps_L2_D = Left(~ismember(Left,Doubly_Contiguous_M) & ~ismember(Left,Doubly_Contiguous_D));
        Gaps_R2_M = Right(~ismember(Right,Doubly_Contiguous_M+(n-1)));
        Gaps_R2_D = Right(~ismember(Right,Doubly_Contiguous_M+(n-1)) & ~ismember(Right,Doubly_Contiguous_D+(2*n-1)));
        Gap_Size2_M = Gaps_R2_M-Gaps_L2_M+1;
        Gap_Size2_D = Gaps_R2_D-Gaps_L2_D+1;

        Gap_SizeI_M = Gap_Size2_M(Gap_Size2_M > n); %size of gaps where isolated binding can happen
        Gap_SizeI_D = Gap_Size2_D(Gap_Size2_D > 2*n);
        Left_I_M = Gaps_L2_M(ismember(Gap_Size2_M,Gap_SizeI_M));    %left position of gaps where isolated binding can happen
        Left_I_D = Gaps_L2_D(ismember(Gap_Size2_D,Gap_SizeI_D));

        Isolated_M = zeros(1,sum(Gap_SizeI_M)-((n+1)*numel(Gap_SizeI_M))+logical(ismember(N,Left_I_M+Gap_SizeI_M-1))+logical(sum(DNA) == 0));   %memory allocation
        Isolated_D = zeros(1,sum(Gap_SizeI_D)-((2*n+1)*numel(Gap_SizeI_D))+logical(ismember(N,Left_I_D+Gap_SizeI_D-1))+logical(sum(DNA) == 0));
        for i = 1:numel(Left_I_M)
            Isolated_M(find(Isolated_M == 0,1):find(Isolated_M == 0,1)+Gap_SizeI_M(i)-(n+1)-1+logical(N == Left_I_M(i)+Gap_SizeI_M(i)-1)+logical(Left_I_M(i) == 1)) = Left_I_M(i)+(1-logical(Left_I_M(i) == 1):Gap_SizeI_M(i)-(n+1)+logical(N == Left_I_M(i)+Gap_SizeI_M(i)-1));
        end
        for k = 1:numel(Left_I_D)
            Isolated_D(find(Isolated_D == 0,1):find(Isolated_D == 0,1)+Gap_SizeI_D(k)-(2*n+1)-1+logical(N == Left_I_D(k)+Gap_SizeI_D(k)-1)+logical(Left_I_D(k) == 1)) = Left_I_D(k)+(1-logical(Left_I_D(k) == 1):Gap_SizeI_D(k)-(2*n+1)+logical(N == Left_I_D(k)+Gap_SizeI_D(k)-1));
        end

        %Population Numbers
        xB_IM = length(Isolated_M);
        xB_SCM = length(Singly_Contiguous_M);
        xB_DCM = length(Doubly_Contiguous_M);
        xB_ID = length(Isolated_D);
        xB_SCD = length(Singly_Contiguous_D);
        xB_DCD = length(Doubly_Contiguous_D);
        xAB = sum(DNA)/n;

        Populations(Events,:) = [xB_IM,xB_SCM,xB_DCM,xB_ID,xB_SCD,xB_DCD,xAB];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        a(Events,:) = Populations(Events,:).*[k_on*L_Monomer(Loops),k_on*L_Monomer(Loops)*w,k_on*L_Monomer(Loops)*(w^2),k_on*L_Dimer(Loops),k_on*L_Dimer(Loops)*w,k_on*L_Dimer(Loops)*(w^2),k_off];
        Probabilities(Events,:) = a(Events,:)./sum(a(Events,:));

        r = [rand,rand,rand,rand,rand,rand,rand];
        tau = (1./a(Events,:)).*log(1./r);
        FiringAmounts(Events,:) = a(Events,:).*tau;
        dt(Events) = min(tau);
        j(Events) = find(tau == min(tau));

        if j(Events) == 1        %Isolated monomer binding
            Bind_Spot_I_M = Isolated_M(randi(xB_IM));
            DNA(Bind_Spot_I_M:Bind_Spot_I_M+(n-1)) = 1;
            Location_History(1,Events) = Bind_Spot_I_M;
            BoundAtSpot(Bind_Spot_I_M) = 1;
        elseif j(Events) == 2    %Singly-Contiguous monomer binding
            Bind_Spot_SC_M = Singly_Contiguous_M(randi(xB_SCM));
            DNA(Bind_Spot_SC_M:Bind_Spot_SC_M+(n-1)) = 1;
            Location_History(2,Events) = Bind_Spot_SC_M;
            BoundAtSpot(Bind_Spot_SC_M) = 1;
        elseif j(Events) == 3    %Doubly-Contiguous monomer binding
            Bind_Spot_DC_M = Doubly_Contiguous_M(randi(xB_DCM));
            DNA(Bind_Spot_DC_M:Bind_Spot_DC_M+(n-1)) = 1;
            Location_History(3,Events) = Bind_Spot_DC_M;
            BoundAtSpot(Bind_Spot_DC_M) = 1;
        elseif j(Events) == 4    %Isolated dimer binding
            Bind_Spot_I_D = Isolated_D(randi(xB_ID));
            DNA(Bind_Spot_I_D:Bind_Spot_I_D+(2*n-1)) = 1;
            Location_History(4,Events) = Bind_Spot_I_D;
            BoundAtSpot(Bind_Spot_I_D) = 1;
            BoundAtSpot(Bind_Spot_I_D+n) = 1;
        elseif j(Events) == 5    %Singly-Contiguous dimer binding
            Bind_Spot_SC_D = Singly_Contiguous_D(randi(xB_SCD));
            DNA(Bind_Spot_SC_D:Bind_Spot_SC_D+(2*n-1)) = 1;
            Location_History(5,Events) = Bind_Spot_SC_D;
            BoundAtSpot(Bind_Spot_SC_D) = 1;
            BoundAtSpot(Bind_Spot_SC_D+n) = 1;
        elseif j(Events) == 6    %Doubly-Contiguous dimer binding
            Bind_Spot_DC_D = Doubly_Contiguous_D(randi(xB_DCD));
            DNA(Bind_Spot_DC_D:Bind_Spot_DC_D+(2*n-1)) = 1;
            Location_History(6,Events) = Bind_Spot_DC_D;
            BoundAtSpot(Bind_Spot_DC_D) = 1;
            BoundAtSpot(Bind_Spot_DC_D+n) = 1;
        elseif j(Events) == 7    %Monomer unbinding
            Bound_Locations = find(BoundAtSpot == 1);
            Unbind_Spot = Bound_Locations(randi(length(Bound_Locations)));
            DNA(Unbind_Spot:Unbind_Spot+(n-1)) = 0;
            Location_History(7,Events) = Unbind_Spot;
            BoundAtSpot(Unbind_Spot) = 0;
        end

        FracCover(Loops,Events+1) = sum(DNA)/N;
        t(Loops,Events+1) = t(Loops,Events)+dt(Events);
        
        DNA_Length(Loops,Events+1) = (numel(find(DNA==1))*CoveredLength)+(numel(find(DNA==0))*UncoveredLength); %calculated length of DNA strand based on lengths given previously

    %     Testing for Equilibrium
        if Events > minIterations && Equilibrium == 0
            FracCoverStates = (FracCover(Loops,(Events-floor(0.25*Events)):1:Events));
            FracCoverChange = abs(diff(FracCoverStates));   %Difference between each state for the last 1/4 of the simulation
            if (mean(FracCoverChange) <= 2*n/N) && (abs(FracCover(Loops,Events-floor(0.25*Events))-FracCover(Loops,Events)) <= 5*n/N)
                Equilibrium = 1;
%                 Irrelevant_t = find(t(Loops,2:end) == 0);
%                 t(Loops,Irrelevant_t) = NaN;
%                 Irrelevant_FracCover = find(FracCover(Loops,2:end) == 0);
%                 FracCover(Loops,Irrelevant_FracCover) = NaN;
            end
        end
        
        if ~isempty(find(BoundAtSpot == 1 & DNA ~= 1, 1))
            disp(['STOP - ERROR AT ', num2str(find(BoundAtSpot == 1 & DNA ~= 1, 1))]);
            flag = 1;
            break
        end
    end

    EventFractions(Loops,:) = [numel(find(j==1)),numel(find(j==2)),numel(find(j==3)),numel(find(j==4)),numel(find(j==5)),numel(find(j==6)),numel(find(j==7))]./Events;
    Max_Time(Loops) = max(t(Loops,:));
   
    figure(1);
%     subplot(2,1,1);
    hold on;
    scatter(t(Loops,:),FracCover(Loops,:),1,Colors(Conc_Values == L_Total,:),'filled','HandleVisibility','off');
    ylabel('Fractional Coverage');
    xlim([0 1.25]);
    ylim([0 1]);
    title('Saturation of DNA Lattice');
    box on;
%     subplot(2,1,2);
%     hold on;
%     scatter(t(Loops,:),DNA_Length(Loops,:)/1000,1,Colors(Conc_Values == w,:),'filled','HandleVisibility','off');
%     ylabel('Length (\mum)');
    xlabel('Time, t');
%     xlim([0 1.25]);
%     ylim([N*UncoveredLength/1000 N*CoveredLength/1000]);
%     title('Length of DNA Strand');
%     box on;
    
    if flag == 1
        break
    end
end

AllTime_MaxTime = max(Max_Time);

Concentrations = unique(Conc_Values);
for b = Concentrations
    clear TimeMatrix;
    clear FracCoverMatrix;
    clear TimeArray;
    clear FracCoverArray;
    clear MaximumTime;
    clear TimeBinLength;
    clear TimeBins;
    
    LoopNumbers = find(Variable_List == b);   %loop numbers where ratio is equal to given b value
    Zeros = zeros(1,length(LoopNumbers));   %zeros to be added for the initial conditions of each loop
    TimeMatrix(1:length(LoopNumbers),:) = t(LoopNumbers,:); %time profiles for this given ratio
    FracCoverMatrix(1:length(LoopNumbers),:) = FracCover(LoopNumbers,:);    %saturation profiles for given ratio
    DNA_LengthMatrix(1:length(LoopNumbers),:) = DNA_Length(LoopNumbers,:);  %growth profiles for the given ratio
    TimeArray = reshape(TimeMatrix,[1,numel(TimeMatrix)]);  %reshape profiles as 1-row arrays
    FracCoverArray = reshape(FracCoverMatrix,[1,numel(FracCoverMatrix)]);
    StrandLengthArray = reshape(DNA_LengthMatrix,[1,numel(DNA_LengthMatrix)]);
%     There should only be 3 zeros at the beginning of the Time and
%     FracCover Arrays
    TimeArray(TimeArray == 0) = []; %removes all zeros
    TimeArray = [Zeros,TimeArray];  %re-adds corresponding zeros for first step of each loop
    FracCoverArray(FracCoverArray == 0) = [];   %removes all zeros
    FracCoverArray = [Zeros,FracCoverArray];    %re-adds corresponding zeros to zero time value
    StrandLengthArray(StrandLengthArray == 0) = []; %removes all zeros
    StrandLengthArray = [Zeros,StrandLengthArray];  %re-adds corresponding zeros to to zero time value
%     Create Time Bins for each value
    MaximumTime = max(TimeArray);  %maximum time
    TimeBinLength = MaximumTime/500;    %width of time bins
    TimeBins = 0:TimeBinLength:MaximumTime; %actual time bins

%     Memory Allocation for next for loop
    PosInBin = 0;
    AvgTimeInBin = zeros(1,length(TimeBins)-1);
    AvgFracCoverInBin = zeros(1,length(TimeBins)-1);
    AvgStrandLengthInBin = zeros(1,length(TimeBins)-1);
    
    for d = 1:length(TimeBins)-1
        PosInBin = find(TimeArray <= TimeBins(d+1) & TimeArray > TimeBins(d));    %positions of values that are within a given time bin
        AvgTimeInBin(d) = sum(TimeArray(PosInBin))/numel(TimeArray(PosInBin));
        AvgFracCoverInBin(d) = sum(FracCoverArray(PosInBin))/numel(FracCoverArray(PosInBin));
        AvgStrandLengthInBin(d) = sum(StrandLengthArray(PosInBin)/numel(PosInBin));
    end

    figure(2);
%     subplot(2,1,1);
    hold on;
    scatter(AvgTimeInBin,AvgFracCoverInBin,10,Colors(Concentrations == b,:),'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.2);
    ylabel('Fractional Coverage');
    xlim([0 AllTime_MaxTime]);
    ylim([0 1]);
    title(['Avg. Saturation of DNA Lattice (n = ', num2str(AverageIterations),')']);
    box on;
%     subplot(2,1,2);
%     hold on;
%     scatter(AvgTimeInBin,AvgStrandLengthInBin/1000,10,Colors(Cooperativities == b,:),'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.5);
    xlabel('Time, t');
%     ylabel('Length (\mum)');
%     xlim([0 AllTime_MaxTime]);
%     ylim([N*UncoveredLength/1000 N*CoveredLength/1000]);
%     title('Length of DNA Strand');
%     box on;
end

Total_Events = zeros(1,length(Concentrations));
Legend = cell(length(Concentrations),1);
for c = 1:length(Concentrations)
    Legend{c} = ['L_T_o_t = ', num2str(Concentrations(c)),' (L_M = ', num2str(Concentrations(c)*Ratio), ', L_D = ', num2str(Concentrations(c)-(Concentrations(c)*Ratio)), ')'];
end
figure(2);
% subplot(2,1,1);
hleg = legend(Legend,'location','southeast');
htitle = get(hleg,'Title');
set(htitle,'String','Total RAD51 Concentration, L_T_o_t');