%==========================================================================
%% Multi-pulse Time-Slip simulation
% 
%==========================================================================
%
%   10/2023 - VM (vmalis@ucsd.edu)   

clear all

%% input parameters:

% number of T-SLIP setups: tag pulses timing in [ms], e.g.
% single tag pulse => one timepoint
% five tag pulses  => five timepoints
tags(1).time = 25;
tags(2).time = [25,45];
tags(3).time = [25,45,65];
tags(4).time = [25,45,65,85];

% tissue relacation times: muscle, bone, blood in [ms]
tissues=3;
% composite tissue fractions, for example: muscle, bone, blood
frac(1) = 1;
frac(2) = 0;
frac(3) = 0;
%T1=[1350,365,1650];
T1=[1500,365,1650];
T2=[32,133,250];
velocity=[0.1, 0, 0]; %   mm/ms

%T1=[1420,2500,3000];
%T2=[32,1000,1000];


% geometry:
y_FOV    = 400;     %   mm
Tslip_h  = 200;     %   mm
Tslip_y1 = 0;       %   mm
delta    = 50;      %   mm
roi_h    = 200;     %   mm


% Time step [ms]
dt = 0.1;     
% Total simulation time in milliseconds [ms]
t_total = 5000;

% noise level in %
sigma=5;


%..........................................................................
% other misc.

% coordinates of ROI
roi_y1   = Tslip_y1+Tslip_h+delta;
roi_ymin = roi_y1;
roi_ymax = roi_y1+roi_h;

% random noise have to be generated multiple times
% for each acquisition: (control + first tag = 2) + each tag
noise=zeros(t_total/dt,1+size(tags,2));
for a=1:2+size(tags,2)
    r = (sigma/50)*randi(100,1,t_total/dt)/100;
    noise(:,a)=r'-mean(r,'all');
end

no_noise_index=3:2:a;
noise(:,no_noise_index)=0;

GKMd=zeros(t_total/dt,size(tags,2));        % General Kinetic Model dark
GKMb=zeros(t_total/dt,size(tags,2));        % General Kinetic Model bright

%% Simulation
for simNumber=1:size(tags,2)

    formatSpec      = 'Simulation: %d of %d is in progress...';
    progressString  = sprintf(formatSpec,simNumber,size(tags,2));
    disp(progressString)
    
    % RF pulse timings (first pulse is Non-selective, rest are selective)
    non_select_time = 5;
    tag_pulse_times = tags(simNumber).time;
    pulse_times=[non_select_time,tag_pulse_times];
    
    % prelocate coordinates of tslip
    tslip_ymin = nan(t_total/dt,length(tag_pulse_times));
    tslip_ymax = nan(t_total/dt,length(tag_pulse_times));
    
    % Bloch
    [mz,~]=blochMz(T1,T2,dt,t_total,pulse_times);
    control=ones(size(mz,1),1);
    
    % GKM
    m = control - mz(:,1,end-1);
    GKMd(:,simNumber)= ASL_gkm(1,T1(1),500,2000,dt,t_total,m,0.5*1e-3);
    m = abs(mz(:,1,1) - mz(:,1,end));
    GKMb(:,simNumber)= ASL_gkm(1,T1(1),500,2000,dt,t_total,m,0.5*1e-3);
    
    
    %% start calculating from NS pulse applied
    for tissue=1:tissues
    
        tslip_ymin(tag_pulse_times(1)/dt,1)=Tslip_y1;
        tslip_ymax(tag_pulse_times(1)/dt,1)=Tslip_y1+Tslip_h;
        n=2;
         
            for t=tag_pulse_times(1)/dt+1:t_total/dt
                
                phase=0;
                Flow(t) = velocity(tissue);
                
                switch tissue
                    case 1
        
                        Flow(t) = velocity(tissue);
        
                        %bpm=75;
                        %Flow(t) = velocity(tissue) * sin(2*pi*bpm/60/1000*t*dt+phase);
                        %if Flow(t)    < 0
                        %    Flow(t)=Flow(t)*0.1;
                        %end
        
                    case 2
        
                        Flow(t) = velocity(tissue);
                        
                        %bpm=75;
                        %Flow(t) = velocity(tissue) * sin(2*pi*bpm/60/1000*t*dt+phase);
                        %if Flow(t)    < 0
                        %   Flow(t)=Flow(t)*0.1;
                        %end
        
                    case 3
        
                        Flow(t) = velocity(tissue);
                        
                        %bpm=75;
                        %Flow(t) = velocity(tissue) * sin(2*pi*bpm/60/1000*t*dt+phase);
                        %if Flow(t)    < 0
                        %    Flow(t)=Flow(t)*0.1;
                        %end
        
                    
                    %case N
                end
        
                % example of analytical function:
                   
            
                for pulse=1:n-1
                tslip_ymin(t,pulse)=tslip_ymin(t-1,pulse)+Flow(t)*dt;
                tslip_ymax(t,pulse)=tslip_ymax(t-1,pulse)+Flow(t)*dt;
                end
            
                if n<=length(tag_pulse_times) && t*dt==tag_pulse_times(n)
                    tslip_ymin(t,n)=Tslip_y1;
                    tslip_ymax(t,n)=Tslip_y1+Tslip_h;
                n=n+1;
                end
            end
            
            % create grid
            dy=dt*velocity(tissue);
            fractionPulses=zeros(t_total/dt,length(tag_pulse_times));
            if dy~=0
                gridY=0:dy:y_FOV;
                % update grid for every time step and add how many pulses that portion
                % experienced
                roi_idxMax = findClosestElementIndex(gridY,roi_ymax);
                roi_idxMin = findClosestElementIndex(gridY,roi_ymin);
                Signal_total=roi_idxMax-roi_idxMin+1;
                pulses_applied=0;
                
                for t=tag_pulse_times(1)/dt:t_total/dt
                    Tslip=zeros(1,length(gridY));
                    % check if new pulse is applied
                    if pulses_applied<length(tag_pulse_times) &&...
                    t*dt==tag_pulse_times(pulses_applied+1)
                        pulses_applied=pulses_applied+1;
                    end
                    
                    for pulse=1:pulses_applied
                        tslip_idxMax = findClosestElementIndex(gridY,tslip_ymax(t,pulse));
                        tslip_idxMin = findClosestElementIndex(gridY,tslip_ymin(t,pulse));
                        % add pulse
                        Tslip(tslip_idxMin:tslip_idxMax)=Tslip(tslip_idxMin:tslip_idxMax)+1;
                    end
                    
                    ROI = Tslip(roi_idxMin:roi_idxMax);
                    % calculate fraction inside roi
                
                    for pulse=1:pulses_applied
                        logicalArray = (ROI == pulse);
                        pulse_signal = sum(logicalArray);
                        fractionPulses(t,pulse) = double(pulse_signal)/double(Signal_total);
                    end
                end
            end
        
        
        fractionNS=ones(t_total/dt,1)-sum(fractionPulses,2);
        fractionNS(1:non_select_time/dt)=0;
        
        
        frac_pulses(tissue,simNumber).fractionPulses=fractionPulses;
        frac_NS(tissue,simNumber).fractionNS=fractionNS;
        
    end

end

%% signal calculation:

BrightON        =   zeros(size(mz,1),tissues,simNumber);
BrightOFF       =   zeros(size(mz,1),tissues);
DarkON          =   BrightON;
DarkOFF         =   BrightOFF;

BrightON_noise  =   BrightON;
BrightOFF_noise =   BrightOFF;
DarkON_noise    =   BrightON;
DarkOFF_noise   =   BrightOFF;

SIRBright       =   zeros(size(mz,1),simNumber);
SIRBright_noise =   SIRBright;

SIRDark         =   SIRBright;
SIRDark_noise   =   SIRBright;

%%

for simulation=1:simNumber
    
    for tissue=1:tissues



        fractionNS=frac_NS(tissue,simulation).fractionNS;
        
        SignalSelective=zeros(size(mz,1),1);
        SignalSelective_noise=zeros(size(mz,1),1);

        for i=1:simulation
            fractionS=frac_pulses(tissue,simulation).fractionPulses(:,i);
            SignalSelective=SignalSelective+mz(:,tissue,i+1).*fractionS;
            SignalSelective_noise=SignalSelective_noise+(mz(:,tissue,i+1)+noise(:,simulation+2)).*fractionS;
        end

        BrightON(:,tissue,simulation)         = frac(tissue)*(squeeze((mz(:,tissue,1))).*fractionNS + SignalSelective);
        BrightON_noise(:,tissue,simulation)   = frac(tissue)*(squeeze((mz(:,tissue,1)+noise(:,2))).*fractionNS + SignalSelective_noise);

        BrightOFF(:,tissue)             = frac(tissue)*squeeze((mz(:,tissue,1)));
        BrightOFF_noise(:,tissue)       = frac(tissue)*squeeze((mz(:,tissue,1)+noise(:,2)));
    
        SignalSelective=zeros(size(mz,1),1);
        SignalSelective_noise=zeros(size(mz,1),1);

        for i=1:simulation
            fractionS=frac_pulses(tissue,simulation).fractionPulses(:,i);
            SignalSelective=SignalSelective+(mz(:,tissue,i)).*fractionS;
            SignalSelective_noise=SignalSelective_noise+(mz(:,tissue,i)+noise(:,simulation+1)).*fractionS;
        end

        
        DarkON(:,tissue,simulation)         = frac(tissue)*(control.*fractionNS + SignalSelective);
        DarkON_noise(:,tissue,simulation)   = frac(tissue)*((control+noise(:,1)).*fractionNS + SignalSelective_noise);
        
        DarkOFF(:,tissue)                   = frac(tissue)*control;
        DarkOFF_noise(:,tissue)             = frac(tissue)*(control+noise(:,1));


    end

end

ON_bright        =   abs(squeeze(sum(BrightON,2)));
ON_bright_noise  =   abs(squeeze(sum(BrightON_noise,2)));
OFF_bright       =   abs(squeeze(sum(BrightOFF,2)));
OFF_bright_noise =   abs(squeeze(sum(BrightOFF_noise,2)));    

ON_dark          =   abs(squeeze(sum(DarkON,2)));
ON_dark_noise    =   abs(squeeze(sum(DarkON_noise,2)));
OFF_dark         =   abs(squeeze(sum(DarkOFF,2)));
OFF_dark_noise   =   abs(squeeze(sum(DarkOFF_noise,2)));



for simulation=1:simNumber

    SIRBright(:,simulation)          = abs(OFF_bright-ON_bright(:,simulation))./control;
    SIRBright_noise(:,simulation)    = abs(OFF_bright_noise-ON_bright_noise(:,i))./(control+noise(:,1));
    
    SIRDark(:,simulation)            = abs(OFF_dark-ON_dark(:,simulation))./control;
    SIRDark_noise(:,simulation)      = abs(OFF_dark_noise-ON_dark_noise(:,simulation))./(control+noise(:,1));

end

SIRBright(1:non_select_time/dt,:)=0;
SIRBright_noise(1:non_select_time/dt,:)=0;

SIRDark(1:non_select_time/dt,:)=0;
SIRDark_noise(1:non_select_time/dt,:)=0;


%% Plots

%% all bright
figure('Color', 'white'); % Set background to white
plot([0:dt:t_total-dt]/1000,(SIRBright(:,1)),'LineStyle','-','LineWidth',2);
legStr{1}='1-tag';
grid on
set(gca, 'TickLabelInterpreter', 'latex'); % Using LaTeX interpreter for tick labels
hold on
for i=2:size(tags,2)
    plot([0:dt:t_total-dt]/1000,(SIRBright(:,i)),'LineStyle','-','LineWidth',2);
    legStr{i}=strcat(num2str(i),'-tag');
end
legend(legStr,'Interpreter','latex','FontSize',18)
ylim([0,1])
xlim([0,t_total/1000])
ylabel('Signal Increase','FontSize',14,'Interpreter','latex');xlabel('time [ms]','FontSize',14,'Interpreter','latex')
title('Bright','FontSize',16,'Interpreter','latex')

export_fig bright.pdf
close


figure('Color', 'white'); % Set background to white
plot([0:dt:t_total-dt]/1000,(GKMb(:,1)),'LineStyle','-','LineWidth',2);
legStr{1}='1-tag';
grid on
set(gca, 'TickLabelInterpreter', 'latex'); % Using LaTeX interpreter for tick labels
hold on
for i=2:size(tags,2)
    plot([0:dt:t_total-dt]/1000,(GKMb(:,i)),'LineStyle','-','LineWidth',2);
    legStr{i}=strcat(num2str(i),'-tag');
end
legend(legStr,'Interpreter','latex','FontSize',18)
ylim([0,1])
xlim([0,t_total/1000])
ylabel('Signal Increase','FontSize',14,'Interpreter','latex');xlabel('time [ms]','FontSize',14,'Interpreter','latex')
title('Bright GKM','FontSize',16,'Interpreter','latex')

export_fig brightGKM.pdf
close


%% all dark
figure('Color', 'white'); % Set background to white
plot([0:dt:t_total-dt]/1000,(SIRDark(:,1)),'LineStyle','-','LineWidth',2);
legStr{1}='1-tag';
grid on
set(gca, 'TickLabelInterpreter', 'latex'); % Using LaTeX interpreter for tick labels
hold on
for i=2:size(tags,2)
    plot([0:dt:t_total-dt]/1000,(SIRDark(:,i)),'LineStyle','-','LineWidth',2);
    legStr{i}=strcat(num2str(i),'-tag');
end
legend(legStr,'Interpreter','latex','FontSize',18)
ylim([0,1])
xlim([0,t_total/1000])
ylabel('Signal Increase','FontSize',14,'Interpreter','latex');xlabel('time [ms]','FontSize',14,'Interpreter','latex')
title('Dark','FontSize',16,'Interpreter','latex')

export_fig dark.pdf
close


figure('Color', 'white'); % Set background to white
plot([0:dt:t_total-dt]/1000,(GKMd(:,1)),'LineStyle','-','LineWidth',2);
legStr{1}='1-tag';
grid on
set(gca, 'TickLabelInterpreter', 'latex'); % Using LaTeX interpreter for tick labels
hold on
for i=2:size(tags,2)
    plot([0:dt:t_total-dt]/1000,(GKMd(:,i)),'LineStyle','-','LineWidth',2);
    legStr{i}=strcat(num2str(i),'-tag');
end
legend(legStr,'Interpreter','latex','FontSize',18)
ylim([0,1])
xlim([0,t_total/1000])
ylabel('Signal Increase','FontSize',14,'Interpreter','latex');xlabel('time [ms]','FontSize',14,'Interpreter','latex')
title('Dark GKM','FontSize',16,'Interpreter','latex')

export_fig darkGKM.pdf
close


%% bright vs dark
for i=1:size(tags,2)

    figure('Color', 'white'); % Set background to white

    plot([0:dt:t_total-dt]/1000,(SIRBright(:,i)),'LineStyle','-','LineWidth',2);
    hold on
    plot([0:dt:t_total-dt]/1000,(SIRDark(:,i)),'LineStyle','-','LineWidth',2);

    grid on
    set(gca, 'TickLabelInterpreter', 'latex'); % Using LaTeX interpreter for tick labels

    titlStr=strcat(num2str(i),'-tag');

    legend({'bright','dark'},'Interpreter','latex','FontSize',18)

    ylim([0,1])
    xlim([0,t_total/1000])

    ylabel('Signal Increase','FontSize',14,'Interpreter','latex');xlabel('time [ms]','FontSize',14,'Interpreter','latex')
    
    title(titlStr,'FontSize',16,'Interpreter','latex')

    filename=strcat(num2str(i),'tags','.pdf');

    export_fig(filename)
    close

end


%% bright vs dark GKM
for i=1:size(tags,2)

    figure('Color', 'white'); % Set background to white

    plot([0:dt:t_total-dt]/1000,(GKMb(:,i)),'LineStyle','-','LineWidth',2);
    hold on
    plot([0:dt:t_total-dt]/1000,(GKMd(:,i)),'LineStyle','-','LineWidth',2);

    grid on
    set(gca, 'TickLabelInterpreter', 'latex'); % Using LaTeX interpreter for tick labels

    titlStr=strcat(num2str(i),'-tag');

    legend({'bright GKM','dark GKM'},'Interpreter','latex','FontSize',18)

    ylim([0,1])
    xlim([0,t_total/1000])

    ylabel('Signal Increase','FontSize',14,'Interpreter','latex');xlabel('time [ms]','FontSize',14,'Interpreter','latex')
    
    title(titlStr,'FontSize',16,'Interpreter','latex')

    filename=strcat(num2str(i),'tagsGKM','.pdf');

    export_fig(filename)
    close

end


%% passThrough vs GKM
