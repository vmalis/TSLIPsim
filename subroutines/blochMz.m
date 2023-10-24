function [mz_results,mxy_results]=blochMz(T1,T2,dt,t_total,pulse_times)

    % Constants
    gamma = 2.675e8;    % Gyromagnetic ratio in rad/(s*T)
    B0 = 3;             % Magnetic field strength in Tesla
    
    % Tissue properties
    NumTissues = 3;
    num_pulse           = length(pulse_times);     % number of tag pulses
    pulse_flip_angle    = pi;                      % 180-degree pulse
    
    % Preallocate arrays to store results
    mz_results  = zeros(t_total / dt,3,num_pulse);
    mxy_results = zeros(t_total / dt,3,num_pulse);
    % Bloch simulation for NonSelective pulse
    time_axis     = 0:dt:t_total-dt;     % Define time_axis for 1 pulse
    Mxy(1:3)      = 0;                % Initialize magnetization
    dMxy_dt (1:3) = 0;                % Initialize deltas
    Mz(1:3)       = 1;                % Initialize magnetization
    dMz_dt (1:3)  = 0;                % Initialize deltas
    
    for pulse=1:length(pulse_times)
        for tissue=1:NumTissues
            Mxy(1:3)      = 0;                % Initialize magnetization
            dMxy_dt (1:3) = 0;                % Initialize deltas
            Mz(1:3)       = 1;                % Initialize magnetization
            dMz_dt (1:3)  = 0;                % Initialize deltas
        
            for t = 1:t_total/dt
                % Calculate relaxation effects
                dMz_dt(tissue) = (1 - Mz(tissue)) / T1(tissue);
                dMxy_dt(tissue)= -Mz(tissue) / T2(tissue);
                % Apply RF pulse for 1 pulse
                if any(abs(t*dt - pulse_times(1:pulse)) < dt/2)
                Mz(tissue) = -Mz(tissue); % Flip magnetization by -1
                end
                % Update spin states using Euler's method
                Mz(tissue) = Mz(tissue) + dMz_dt(tissue) * dt;
                Mxy(tissue) = Mxy(tissue) + dMxy_dt(tissue) * dt;
                % Store results for 1 pulse
                mz_results(t,tissue,pulse) = Mz(tissue);
                mxy_results(t,tissue,pulse) = Mxy(tissue);
            end
        end
    end

end