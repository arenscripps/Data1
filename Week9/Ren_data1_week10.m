%%% Ren_data1_week9.m
%
%   Look at coherence between two records of your choice
%
%   A. Ren, December 5, 2017

%% attempt to use glider data and calculate for wavenumber

%% glider data
%   two choices, Line 80 and Line 66 == choose like 66 for now

%ncdisp('CUGN_line_66.nc')

%% load data

% filename = 'CUGN_line_66.nc';
% 
% L66_time = ncread(filename, 'time_uv');
% L66_v = ncread(filename, 'v_depth_mean');
% L66_u = ncread(filename, 'u_depth_mean');
% L66_lat = ncread(filename, 'lat_uv');
% L66_lon = ncread(filename, 'lon_uv');
% 
% %  adjust time to matlab datetime
% L66_time_mt = datetime(L66_time, 'ConvertFrom', 'posixtime', ...
%     'Timezone', 'America/Los_Angeles');
% %  now remove the timezone which helps with conversions and calculations
% L66_time_mt.TimeZone = '';

%% diagnostics
% time_start = L66_time_mt(1);
% time_end = L66_time_mt(end);
% 
% L66_lat(1:10)
% L66_lon(1:10)
% diff(L66_lat(1:10))
% diff(L66_lon(1:10))
% gdist(L66_lat(1), L66_lon(1), L66_lat(2), L66_lon(2))
% gdist(L66_lat(2), L66_lon(2), L66_lat(3), L66_lon(3))
% gdist(L66_lat(4), L66_lon(4), L66_lat(5), L66_lon(5))
% L66_lon1 = ncread(filename, 'lon');
% diff(L66_lon1)
% gdist(L66_lat1(1), L66_lon1(1), L66_lat1(2), L66_lon1(2))
% gdist(L66_lat1(2), L66_lon1(2), L66_lat1(3), L66_lon1(3))
% L66_time_all = ncread(filename, 'time');
% diff(L66_time_all)
% L66_time_all(1:10)
% L66_time_all = datetime(L66_time_all, 'ConvertFrom', 'posixtime', ...
% 'Timezone', 'America/Los_Angeles');
% L66_time_all(1:10)
% diff(L66_time_all(1:10))

%   Because the glider data are not collected uniformly in space or time,
%   we have to table this for the assignment.  Use climatology/objectively
%   mapped products instead.

%%  Glider Climatology Data
%
%   Investigate geostrophic velocity:
%   Take Line 67 "total" data and investigate coherence between 80 meters
%   depth and 200 meters depth.  I expect to see some coherence on large
%   scales in terms of the size of the major currents, but I expect the
%   spectra should show different wavenumber peaks.

%% load data
filename = 'total_z_66.nc';

L66_dist = ncread(filename, 'distance');
L66_time = ncread(filename, 'time');
L66_time_mt = datetime(L66_time, 'ConvertFrom', 'posixtime', ...
    'TimeZone', 'America/Los_Angeles');
L66_time_mt.TimeZone = '';

L66_dep = ncread(filename, 'depth');
L66_gvel = ncread(filename, 'geostrophic_velocity');

%% checking data
%  we need the distance from shore variable to be evenly spaced, which is
%  probably true since this is a data product

figure
plot(L66_dist)
% yes, straight line means constant distance between
diff(L66_dist(1:10))

%  the distance between each data point is 5 km
samplingdist = 5; %km
%   the nyquist wavenumber is 1/(2*samplingdist) = 1/10; 1 cycle per 10 km

%% basics

L66_vel80 = squeeze(L66_gvel(8, :, :));  % matrix is time versus distance
L66_vel80 = L66_vel80';                  % change to columns of segments
L66_vel200 = squeeze(L66_gvel(20, :, :));
L66_vel200 = L66_vel200';

%% what do data look like?

%   plotting some characteristic segments (first 8 transects in time)
figure('Name', 'Select Transects')
subplot(2, 1, 1)
    plot(L66_dist, L66_vel80(:, 1:8))
    title('Alongshore velocity at 80 m')
    ylabel('Geostrophic Velocity (m/s)')
    xlabel('Distance from Shore (km)')
subplot(2, 1, 2)
    plot(L66_dist, L66_vel200(:, 1:8))
    title('Alongshore velocity at 200 m')
    ylabel('Geostrophic Velocity (m/s)')
    xlabel('Distance from Shore (km)')
    
%   plotting all data
figure('Name', 'All Data')
subplot(2, 1, 1)
    plot(L66_vel80(:))
    title('Alongshore velocity at 80 m')
    ylabel('Geostrophic Velocity (m/s)')
    xlabel('Spatial Scale (km)')
subplot(2, 1, 2)
    plot(L66_vel200(:))
    title('Alongshore velocity at 200 m')
    ylabel('Geostrophic Velocity (m/s)')
    xlabel('Spatial Scale (km)')
    

%% there are NaNs
%   get rid of first 5 data points closest to shore (1-25 km)
vel80 = L66_vel80(6:81, :);
vel200 = L66_vel200(6:81, :);
distfft = L66_dist(6:81);

%   get rid of any remaining segments with NaNs
%   1) identify such segments
indexnan1 = any(isnan(vel80));
indexnan2 = any(isnan(vel200));

vel80_nonan = vel80(:, ~indexnan1);
vel200_nonan = vel200(:, ~indexnan2);

sum(sum(isnan(vel80_nonan)))
sum(sum(isnan(vel200_nonan)))

%% detrend
%  looks like I should just demean - fairly stationary
vel80_dt = detrend(vel80_nonan, 'constant');  % remove mean
vel200_dt = detrend(vel200_nonan, 'constant');

%% not going to window ??
filteron = 'hanning';

N = length(distfft);  % how many data points per segment

switch filteron
    case 'hanning'
vel80_h = vel80_dt.* (hann(N)*ones(1, size(vel80_dt, 2)));
vel200_h = vel200_dt.* (hann(N)*ones(1, size(vel80_dt, 2)));
    case 'none'
end

%%  Take FFT
switch filteron
    case 'hanning'
        A_vel80 = fft(vel80_h);
        A_vel200 = fft(vel200_h);
    case 'none'
        A_vel80 = fft(vel80_dt);
        A_vel200 = fft(vel200_dt);
end

%% Calculate spectra
amp_vel80 = abs(A_vel80(1:(N/2) + 1, :)).^2 ;  % even N
amp_vel80(2:end-1, :) = 2*amp_vel80(2:end-1, :); % mult by 2
alpha = 5; %km %lengthofsample/N = deltaL in units you want
amp_vel80 = amp_vel80 * alpha/N;              % normalize

amp_vel200 = abs(A_vel200(1:(N/2) + 1, :)).^2 ;  % even N
amp_vel200(2:end-1, :) = 2*amp_vel200(2:end-1, :); % mult by 2
alpha = 5; %km %lengthofsample/N = deltaL in units you want
amp_vel200 = amp_vel200 * alpha/N;              % normalize

% average over segments
amp_vel80_m = mean(amp_vel80, 2);
amp_vel200_m = mean(amp_vel200, 2);

% xaxis
scale = 1/5; %km  1 measurement / 5 km * ()= cycles per km
wavenumberaxis = scale* (0:N/2)/N;

% error bar
totalindepseg = size(A_vel80, 2);

switch filteron
    case 'hanning'
        dof = 36/19 * totalindepseg;
    case 'none'
        dof = 2*totalindepseg;
end
err_high = dof/chi2inv(0.05/2, dof);
err_low = dof/chi2inv(1-0.05/2, dof);

%% Figures Spectra

figure('Name', 'Velocity at 80 m')
    loglog(wavenumberaxis, amp_vel80_m, 'LineWidth', 1.2)
    hold on
    loglog([wavenumberaxis(2) wavenumberaxis(2)], ...
        ([err_low err_high]*amp_vel80_m(2) *1), 'k', 'LineWidth', 1.8)
    loglog(wavenumberaxis(2:10), (6e-5)*wavenumberaxis(2:10).^-(2),...
        '--', 'LineWidth', 1)    
    loglog(wavenumberaxis(3:22), (2e-8)*wavenumberaxis(3:22).^-(11/3),...
        '--', 'LineWidth', 1)
    loglog(wavenumberaxis(11:30), (0.3e-9)*wavenumberaxis(11:30).^-5,...
        '--', 'LineWidth', 1)
    grid on
    title('Alongshore Velocity at 80 m')
    ylabel('[m/s]^{2} / cycles per km')
    xlabel('Cycles per km')
    legend('Spectrum', 'error bar', 'K^{-2}', 'K^{-(11/3)}', 'K^{-5}') 
    
figure('Name', 'Velocity at 200 m')
    loglog(wavenumberaxis, amp_vel200_m, 'LineWidth', 1.2)
    hold on
    loglog([wavenumberaxis(2) wavenumberaxis(2)], ...
        ([err_low err_high]*amp_vel200_m(2) *1), 'k', 'LineWidth', 1.8)
    loglog(wavenumberaxis(2:10), (1e-5)*wavenumberaxis(2:10).^-(2),...
        '--', 'LineWidth', 1)    
    loglog(wavenumberaxis(3:22), (1e-8)*wavenumberaxis(3:22).^-(11/3),...
        '--', 'LineWidth', 1)
    loglog(wavenumberaxis(11:30), (1e-10)*wavenumberaxis(11:30).^-5,...
        '--', 'LineWidth', 1)
    grid on
    title('Alongshore Velocity at 200 m')
    ylabel('[m/s]^{2} / cycles per km')
    xlabel('Cycles per km')
    legend('Spectrum', 'error bar', 'K^{-2}', 'K^{-(11/3)}', 'K^{-5}') 

%% Coherence

C_xy = conj(A_vel80(1:(N/2) + 1,:)).* A_vel200(1:(N/2) + 1,:)*alpha/N;
C_xy(2:end-1,:) = 2*C_xy(2:end-1,:);
% 
coherence = abs(mean(C_xy,2))./sqrt(amp_vel80_m .* amp_vel200_m);

alp = 0.05;
segs = size(A_vel80, 2) - 1;
threshold = sqrt( 1-alp^(1/segs) );
figure
    semilogx((0:N/2)/N*scale,coherence,...
        [wavenumberaxis(2) wavenumberaxis(end)],[threshold threshold])
    grid on
    ylabel('Coherence')
    xlabel('Cycles per km')

%   and Phase
seg = segs+1;
phase=atan2(-imag(mean(C_xy,2)),real(mean(C_xy,2)));
delta_phase = sqrt((1-coherence.^2)./(abs(coherence).^2*2*seg));
figure
    semilogx(wavenumberaxis, [phase phase+delta_phase phase-delta_phase])
    axis([wavenumberaxis(1) wavenumberaxis(end) -pi pi])
    grid on
    ylabel('Phase (radians)')
    xlabel('Cycles per km')
    
%  subplot figure
figure
subplot(2, 1, 1)
    semilogx((0:N/2)/N*scale,coherence,...
        [wavenumberaxis(2) wavenumberaxis(end)],[threshold threshold])
    grid on
    ylabel('Coherence')
    xlabel('Cycles per km')
subplot(2, 1, 2)
    semilogx(wavenumberaxis, [phase phase+delta_phase phase-delta_phase])
    axis([wavenumberaxis(1) wavenumberaxis(end) -pi pi])
    grid on
    ylabel('Phase (radians)')
    xlabel('Cycles per km')
    
%%  Testing if coherence is reasonable
%   (from Sarah Gille)

%   segments with random noise 1000 data points per segment, 300 segments
vel80=randn(1000,300);
vel200=randn(1000,300);

%   fft
A_vel80 = fft(vel80);
A_vel200 = fft(vel200);

%   calculate spectra
N = 300;

amp_vel80 = abs(A_vel80(1:(N/2) + 1, :)).^2 * (1/N);   % no normalization
amp_vel200 = abs(A_vel200(1:(N/2) + 1, :)).^2 * (1/N); % no normalization
amp_vel80_m = mean(amp_vel80, 2);
amp_vel200_m = mean(amp_vel200, 2);

wavenumberaxis = (0:N/2)/N;

%   coherence
C_xy = conj(A_vel80(1:(N/2) + 1,:)).* A_vel200(1:(N/2) + 1,:)/N;
coherence = abs(mean(C_xy,2))./sqrt(amp_vel80_m .* amp_vel200_m);

%   threshold
alp = 0.05;
segs = N - 1;
threshold = sqrt( 1-alp^(1/segs) );

%   figure
figure
    semilogx((0:N/2)/N,coherence,...
        [wavenumberaxis(2) wavenumberaxis(end)],[threshold threshold])
    grid on
    ylim([0 1])

