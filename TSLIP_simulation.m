%==========================================================================
%% Multi-pulse Time-Slip simulation
% 
%==========================================================================
%
%   10/2023 - VM (vmalis@ucsd.edu)   


%% input parameters:

% number of T-SLIP setups: tag pulses timing in [ms], e.g.
% single tag pulse => one timepoint
% five tag pulses  => five timepoints
tags(1).time = 25;
tags(2).time = [25,45];
tags(3).time = [25,45,65];
tags(4).time = [25,45,65,85];


% composite tissue fractions, for example: bone, blood and muslce
muscle_frac = 0.25;
bone_frac   = 0.5;
blood_frac  = 0.25;


% tissue relacation times: muscle, bone, blood in [ms]
T1=[1420,365,1550];
T2=[32,133,250];


% geometry:
y_FOV    = 400;     %   mm
Tslip_h  = 200;     %   mm
Tslip_y1 = 0;       %   mm
delta    = 50;      %   mm
roi_h    = 200;     %   mm
velocity = 0.25;    %   mm/ms


% Time step [ms]
dt = 0.1;     
% Total simulation time in milliseconds [ms]
t_total = 10000;

%..........................................................................

% other misc.
% coordinates of ROI
roi_y1   = Tslip_y1+Tslip_h+delta;
roi_ymin = roi_y1;
roi_ymax = roi_y1+roi_h;

% random noise
r = 0.03*randi(100,1,t_total/dt)/100;
noise=r';



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






%% start calculating from NS pulse applied
tslip_ymin(tag_pulse_times(1)/dt,1)=Tslip_y1;
tslip_ymax(tag_pulse_times(1)/dt,1)=Tslip_y1+Tslip_h;
n=2;

for t=tag_pulse_times(1)/dt+1:t_total/dt
    
    phase=0;
    %bloodFlow(t) = velocity;
    
    % example of analytical function:
       bpm=75;
       bloodFlow(t) = velocity * sin(2*pi*bpm/60/1000*t*dt+phase);
       if bloodFlow(t)    < 0
          bloodFlow(t)=bloodFlow(t)*0.1;
       end

    for pulse=1:n-1
    tslip_ymin(t,pulse)=tslip_ymin(t-1,pulse)+bloodFlow(t)*dt;
    tslip_ymax(t,pulse)=tslip_ymax(t-1,pulse)+bloodFlow(t)*dt;
    end

    if n<=length(tag_pulse_times) && t*dt==tag_pulse_times(n)
        tslip_ymin(t,n)=Tslip_y1;
        tslip_ymax(t,n)=Tslip_y1+Tslip_h;
    n=n+1;
    end
end

% create grid
dy=dt*velocity;
gridY=0:dy:y_FOV;
fractionPulses=zeros(t_total/dt,length(tag_pulse_times));

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

fractionNS=ones(t_total/dt,1)-sum(fractionPulses,2);
fractionNS(1:non_select_time/dt)=0;

%% signal calculation:
% Signal Tag ON
S_tagON=fractionNS.*mz(:,3,1);

for pulse=1:length(tag_pulse_times)
    S_tagON=S_tagON+fractionPulses(:,pulse).*mz(:,3,pulse+1);
end

% Signal Tag OFF
fractionNS=ones(t_total/dt,1);
fractionNS(1:non_select_time/dt)=0;
S_tagOFF=fractionNS.*mz(:,3,1);

% Subtraction
X=abs(S_tagON-S_tagOFF);


MultiOff=abs(muscle_frac*mz(:,1,1))+abs(bone_frac*mz(:,2,1))+abs(blood_frac*fractionNS.*mz(:,3,1));


SIR=abs(S_tagON-S_tagOFF)./abs(MultiOff(end));      % predicted
SIRn=(abs(S_tagON-S_tagOFF)+noise)./abs(MultiOff);  % with noise


%write to results
Results(simNumber).SIR=SIR;     
Results(simNumber).SIRn=SIRn;

end




%% example of SIR plot
figure('Color', 'white'); % Set background to white
plot([0:dt:t_total-dt]/1000,(Results(1).SIR),'LineStyle','-','LineWidth',2);
legStr{1}='1-tag';
grid on
set(gca, 'TickLabelInterpreter', 'latex'); % Using LaTeX interpreter for tick labels
hold on
for i=2:size(tags,2)
    plot([0:dt:t_total-dt]/1000,(Results(i).SIR),'LineStyle','-','LineWidth',2);
    legStr{i}=strcat(num2str(i),'-tag');
end
legend(legStr,'Interpreter','latex','FontSize',18)
ylim([0,2])
xlim([0,7])
ylabel('Signal Increase','FontSize',14,'Interpreter','latex');xlabel('time [ms]','FontSize',14,'Interpreter','latex')
