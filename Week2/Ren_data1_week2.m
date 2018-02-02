%%% Ren_data1_week2.m
%
%   HW week 2, Look at temperature from Scripps Pier, manual versus
%   automatic sampling, compare
%
%   A. Ren October 12, 2017

%%  Load Data
%   Manual Scripps Shore Station Data is in an excel file:
load('scrippspier_manual_t.mat', 'scrippsmanualdata', ...
    'header_scrippsmanualdata')           % manually import and create .mat

yearvector = scrippsmanualdata(:, 1);
monthvector = scrippsmanualdata(:, 2);
dayvector = scrippsmanualdata(:, 3);
hourvector = floor(scrippsmanualdata(:, 4)/100);
minutevector = scrippsmanualdata(:, 4) - (hourvector*100);
secondvector = repmat([0], size(yearvector, 1), 1);

scrippsmanual_time = datetime(yearvector, monthvector, dayvector, ...
                              hourvector, minutevector, secondvector, ...
                              'TimeZone', 'America/Los_Angeles');
%  note that this includes naT, "not a time" values, for when time was NaN
%  in the original data

scrippsmanual_stemp = scrippsmanualdata(:, 6);
scrippsmanual_btemp = scrippsmanualdata(:, 8);
%  for now, ignore the flag system, col 5, 7, 9 - for time, stemp, btemp
                          
%   Automatic Scripps Pier data from thredds server:
yearlist = 2005:2015;
SP_temp_2015 = [];
SP_time_2015 = [];
SP_depth_2015 = [];

for n = 1:length(yearlist);
    yearof = yearlist(n);
    fileloca = sprintf('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-%04d.nc', yearof);
    %ncdisp(fileloca)

    newtemp = ncread(fileloca, 'temperature');
    newtime = ncread(fileloca, 'time');
    newdepth = ncread(fileloca, 'depth');
    SP_temp_2015 = [SP_temp_2015; newtemp];
    SP_time_2015 = [SP_time_2015; newtime];
    SP_depth_2015 = [SP_depth_2015; newdepth];
end         
%there is an outlier from visual inspection
% remove
badvalue = find(SP_temp_2015 > 30);
SP_temp_2015(badvalue) = NaN;
SP_time_2015(badvalue) = NaN;
SP_depth_2015(badvalue) = NaN;


SP_temp_2015 = double(SP_temp_2015);
SP_time_2015 = double(SP_time_2015);
SP_depth_2015 = double(SP_depth_2015);

SP_time_2015_mtlabtime = datetime(SP_time_2015, 'ConvertFrom', ...
    'posixtime', 'TimeZone', 'America/Los_Angeles');

save('scrippspier_auto_t.mat', 'SP_temp_2015', 'SP_time_2015', ...
    'SP_depth_2015', 'SP_time_2015_mtlabtime')

%%  Basic Plot

figure('Name', 'Temperature, Raw Comparison')

% temperature: automatic
plot(SP_time_2015_mtlabtime, SP_temp_2015, 'k', 'LineWidth', 0.5)

hold on

% find where 2005 starts
startday = find(scrippsmanual_time > ...
           datetime([2005 01 01], 'TimeZone', 'America/Los_Angeles'), ...
           1, 'first');

% surface temperature: manual
plot(scrippsmanual_time(startday:end), ...
    scrippsmanual_stemp(startday:end), 'r', 'LineWidth', 0.5)

% bottom temperature: manual
plot(scrippsmanual_time(startday:end), ...
    scrippsmanual_btemp(startday:end), 'b', 'LineWidth', 0.5)

xlim([datenum(2004, 12, 31) datenum(2016, 1, 2)])
grid on
xlabel('Time')
ylabel('Temperature (degC)')

%%   Lots of data, pick one year to zoom in, say 2013

% find time interval
startday1 = find(...
    SP_time_2015_mtlabtime > ...
    datetime([2013 01 01], 'TimeZone', 'America/Los_Angeles'), 1, 'first');
endday1 = find(...
    SP_time_2015_mtlabtime > ...
    datetime([2013 12 31], 'TimeZone', 'America/Los_Angeles'), 1, 'first');

startday2 = find(...
    scrippsmanual_time > ...
    datetime([2013 01 01], 'TimeZone', 'America/Los_Angeles'), 1, 'first');
endday2 = find(scrippsmanual_time > ...
    datetime([2013 12 31], 'TimeZone', 'America/Los_Angeles'), 1, 'first');

%%   Plot 2013
figure('Name', 'Temperature, raw, 2013 comparison')

% temperature: automatic
plot(SP_time_2015_mtlabtime(startday1:endday1), ...
    SP_temp_2015(startday1:endday1), 'k', 'LineWidth', 0.5)

hold on

% surface temperature: manual
plot(scrippsmanual_time(startday2:endday2), ...
    scrippsmanual_stemp(startday2:endday2), 'r', 'LineWidth', 0.7)

% bottom temperature: manual
plot(scrippsmanual_time(startday2:endday2), ...
    scrippsmanual_btemp(startday2:endday2), 'Color', [135 206 250]/255, 'LineWidth', 0.5)

grid on
xlabel('Time')
ylabel('Temperature (degC)')

%%   Manual Temperature Record : use average of surface and bottom

scrippsmanual_avtemp = (scrippsmanual_stemp + scrippsmanual_btemp)/2;

%%  Means and Error of Mean
%   since manual record ends on Oct 31, 2015, use automatic record until
%   that day
%   since automatic record starts on Jun 16, 2005, use manual record
%   starting then

endday4 = find(SP_time_2015_mtlabtime > ...
               datetime([2015 10 31 9 25 0], 'TimeZone', 'America/Los_Angeles'), ...
               1, 'first');
startday4 = find(scrippsmanual_time > ...
               datetime([2005 06 16], 'TimeZone', 'America/Los_Angeles'), ...
               1, 'first');

manualdata = scrippsmanual_avtemp(startday4:end);
automadata = SP_temp_2015(1:endday4);

%   calculate mean 
SP_man_av = mean(manualdata, 'omitnan');  disp('mean: '); disp(SP_man_av)
SP_aut_av = mean(automadata, 'omitnan');  disp(SP_aut_av)

%   how many data? N?
SP_man_N = sum(~isnan(manualdata));  disp('N: '); disp(SP_man_N)
SP_aut_N = sum(~isnan(automadata));  disp(SP_aut_N)

%   error of mean
SP_man_av_std = std(manualdata, 'omitnan')/ sqrt(SP_man_N);
SP_aut_av_std = std(automadata, 'omitnan')/ sqrt(SP_aut_N);

disp('error of mean: '); disp(SP_man_av_std); disp(SP_aut_av_std)

%%  Variance / Standard deviation

SP_man_std = std(manualdata, 'omitnan'); disp(SP_man_std)
SP_aut_std = std(automadata, 'omitnan'); disp(SP_aut_std)

%   subsample automatic data to find nearest neighbor to manual data
manualtime = scrippsmanual_time(startday4:end);
SP_aut_subs_i = [];
for n = 1:length(manualtime)
    subtractvalue = abs(SP_time_2015_mtlabtime - manualtime(n));
    subtractvalue = seconds(subtractvalue);
    indexsmall = find(subtractvalue == min(subtractvalue));
    if length(indexsmall) >1 % sometimes there is a duplicate "closest" pt
        indexsmall = indexsmall(1);
    elseif isempty(indexsmall);
        %disp('help'); disp(n);
        indexsmall = NaN;
    end
    SP_aut_subs_i = [SP_aut_subs_i, indexsmall];
end

SP_aut_subs_std = std(automadata(~isnan(SP_aut_subs_i)), 'omitnan');
disp(SP_aut_subs_std)

%% Theoretical PDFs

%  uniform
%  from Bendat and Piersol, mu = (a+b)/2; var = (b-a)^2/12;

%  SP_aut_std^2 = (b-a)^2/12
%  b - a = sqrt(12*(SP_aut_std^2))

%  b + a = 2* SP_aut_av;

b = (...
    sqrt(...
    12*(SP_aut_std^2)...
    )...
    +(2* SP_aut_av)...
    )/2;
a = 2* SP_aut_av - b;
disp('a: '); disp(a)
disp('b: '); disp(b)

%  and the probability density function is 1/b-a for x between a and b.
pd = 1/(b-a);
disp(pd)

% gaussian
x = 1:0.1:35;
y = (1/(SP_aut_std*sqrt(2*pi)))*exp((-(x-SP_aut_av).^2)/(2*SP_aut_std^2));

figure
plot(x, y)
ylabel('Probability Density')
xlabel('Temperature (degC)')

% coin toss
% assume that a and b are sigma away from mu; each at 0.5 pd

%% pdfs
edges = 10:1:26;
figure
h = histogram(manualdata, edges, 'Normalization', 'pdf');
hold on
h2 = histogram(automadata, edges, 'Normalization', 'pdf');
ylabel('Probability Density')
xlabel('Temperature (degC)')

%% crude cdfs

cdf_manual = cumsum(h.Values);
cdf_automa = cumsum(h2.Values);

figure
plot(11:26, cdf_manual)
hold on
plot(11:26, cdf_automa)
ylabel('Cumulative Probability')
xlabel('Temperature (degC)')

diff = abs(cdf_automa - cdf_manual);
maxdiff = max(diff);

%  parameter (from wikipedia)
calpha = 1.63; % alpha = 0.01
testD = calpha * sqrt( (SP_man_N + SP_aut_N)/ (SP_man_N * SP_aut_N));

disp('maximum difference (D) : '); disp(maxdiff);
disp('D parameter test : '); disp(testD);





