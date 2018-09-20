% Licensed under GNU GPL v3
% Copyright 2018 - Flávio Eler de Melo (flavio.eler@gmail.com)
function generate_plots(model,truth,meas,est_phd,est_cphd,est_dgcphd)

% Save figures?
saveFigures = false;

% Set path for figures
FIGURES_PATH = 'figures/';

numFilters = 3;
colors = lines(numFilters);

iptsetpref('ImshowBorder','tight');

%% Plot results
[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);

titleFontSize = 22;
legendFontSize = 24;
labelFontSize = 20;

markerSize = 8;

% Handlers and labels
hnd = zeros(3,1);
lbl = cell(3,1);

% Plot tracks for PHD
% figure('outerposition',[0 0 1080 1080],'Color','white');
figure('position',[140 140 940 940],'outerposition',[140 140 940 940],'Color','white','Resize','off');

% Plot ground truth
limit = [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
for i = 1:truth.total_tracks
    Pt = X_track(:,k_birth(i):1:k_death(i),i); Pt = Pt([1 3],:);
    hnd(1) = plot( Pt(1,:), Pt(2,:), 'k-', 'LineWidth', 2); hold on;
    plot( Pt(1,1), Pt(2,1), 'ko', 'MarkerSize', markerSize, 'LineWidth', 1.5);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'k^', 'MarkerSize', markerSize, 'LineWidth', 1.5);
end
axis equal; axis(limit);
lbl{1} = 'Ground truth';

% Plot measurements
for k = 1:meas.K
    if ~isempty(meas.Z{k})
        hnd(2) = line(meas.Z{k}(1,:),meas.Z{k}(2,:),'LineStyle','none','Marker','x','Markersize',markerSize,'Color',0.7*ones(1,3));
    end
end
lbl{2} = 'Measurements';

% Plot estimates
for k = 1:meas.K
    if ~isempty(est_phd.X{k})
        P = est_phd.X{k}([1 3],:);
        hnd(3) = line(P(1,:),P(2,:),'LineStyle','none','Marker','.','Markersize',markerSize,'Color','blue');
    end
end
lbl{3} = 'Estimates';

h_l = legend(hnd,lbl,'Location', 'Northwest', 'Orientation', 'Vertical');
set(h_l,'FontSize',legendFontSize);
% set(gca,'FontSize',24);
ylabel('Coordinate y [m]','FontSize',labelFontSize);
xlabel('Coordinate x [m]','FontSize',labelFontSize);
title('Tracks - PHD filter', 'FontSize', titleFontSize);

set(gca,'LooseInset',get(gca,'TightInset'));
if saveFigures
    export_fig([FIGURES_PATH, 'Tracks_PHD.pdf'], '-pdf', '-native', '-nocrop');
end


% Plot tracks for CPHD
% figure('outerposition',[0 0 1080 1080],'Color','white');
figure('position',[140 140 940 940],'outerposition',[140 140 940 940],'Color','white','Resize','off');

% Plot ground truth
limit = [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
for i = 1:truth.total_tracks
    Pt = X_track(:,k_birth(i):1:k_death(i),i); Pt = Pt([1 3],:);
    hnd(1) = plot( Pt(1,:), Pt(2,:), 'k-', 'LineWidth', 2); hold on;
    plot( Pt(1,1), Pt(2,1), 'ko', 'MarkerSize', markerSize, 'LineWidth', 1.5);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'k^', 'MarkerSize', markerSize, 'LineWidth', 1.5);
end
axis equal; axis(limit);
lbl{1} = 'Ground truth';

% Plot measurements
for k = 1:meas.K
    if ~isempty(meas.Z{k})
        hnd(2) = line(meas.Z{k}(1,:),meas.Z{k}(2,:),'LineStyle','none','Marker','x','Markersize',markerSize,'Color',0.7*ones(1,3));
    end
end
lbl{2} = 'Measurements';

% Plot estimates
for k = 1:meas.K
    if ~isempty(est_cphd.X{k})
        P = est_cphd.X{k}([1 3],:);
        hnd(3) = line(P(1,:),P(2,:),'LineStyle','none','Marker','.','Markersize',markerSize,'Color','blue');
    end
end
lbl{3} = 'Estimates';

h_l = legend(hnd,lbl,'Location', 'Northwest', 'Orientation', 'Vertical');
set(h_l,'FontSize',legendFontSize);
% set(gca,'FontSize',24);
ylabel('Coordinate y [m]','FontSize',labelFontSize);
xlabel('Coordinate x [m]','FontSize',labelFontSize);
title('Tracks - CPHD filter', 'FontSize', titleFontSize);

set(gca,'LooseInset',get(gca,'TightInset'));
if saveFigures
    export_fig([FIGURES_PATH, 'Tracks_CPHD.pdf'], '-pdf', '-native', '-nocrop');
end

% Plot tracks for DG-CPHD
% figure('outerposition',[0 0 1080 1080],'Color','white');
figure('position',[140 140 940 940],'outerposition',[140 140 940 940],'Color','white','Resize','off');

% Plot ground truth
limit = [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) ];
for i = 1:truth.total_tracks
    Pt = X_track(:,k_birth(i):1:k_death(i),i); Pt = Pt([1 3],:);
    hnd(1) = plot( Pt(1,:), Pt(2,:), 'k-', 'LineWidth', 2); hold on;
    plot( Pt(1,1), Pt(2,1), 'ko', 'MarkerSize', markerSize, 'LineWidth', 1.5);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'k^', 'MarkerSize', markerSize, 'LineWidth', 1.5);
end
axis equal; axis(limit);
lbl{1} = 'Ground truth';

% Plot measurements
for k = 1:meas.K
    if ~isempty(meas.Z{k})
        hnd(2) = line(meas.Z{k}(1,:),meas.Z{k}(2,:),'LineStyle','none','Marker','x','Markersize',markerSize,'Color',0.7*ones(1,3));
    end
end
lbl{2} = 'Measurements';

% Plot estimates
for k = 1:meas.K
    if ~isempty(est_dgcphd.X{k})
        P = est_dgcphd.X{k}([1 3],:);
        hnd(3) = line(P(1,:),P(2,:),'LineStyle','none','Marker','.','Markersize',markerSize,'Color','blue');
    end
end
lbl{3} = 'Estimates';

h_l = legend(hnd,lbl,'Location', 'Northwest', 'Orientation', 'Vertical');
set(h_l,'FontSize',legendFontSize);
% set(gca,'FontSize',24);
ylabel('Coordinate y [m]','FontSize',labelFontSize);
xlabel('Coordinate x [m]','FontSize',labelFontSize);
title('Tracks - DG-CPHD filter', 'FontSize', titleFontSize);

set(gca,'LooseInset',get(gca,'TightInset'));
if saveFigures
    export_fig([FIGURES_PATH, 'Tracks_DG-CPHD.pdf'], '-pdf', '-native', '-nocrop');
end

if saveFigures
    close all;
end
pause(0.01);

%% Overall performance index
% Inserted for the plots used in the paper
filterTypeString = {...
    'PHD filter', ...
    'CPHD filter', ...
    'DG-CPHD filter'};
colors = lines(6);

% Plot Cardinality
N_max = max(truth.N);
limits = [0 N_max+ceil(sqrt(N_max))+2];
% limits = [0 20];
iFilter = 1;

axisFontSize = 14;
labelFontSize = 20;
legendFontSize = 20;
titleFontSize = 20;

figure('units','normalized','outerposition',[0 0 1/2 1]);
set(gcf, 'PaperPositionMode', 'auto', 'Color', 'w');

CardMean = [est_phd.N,est_cphd.N,est_dgcphd.N]';
CardStdDev = sqrt( [est_phd.S,est_cphd.S,est_dgcphd.S]' );

subplot(3,1,1); box on; hold on; grid on;
stairs(1:meas.K,truth.N,'k','LineWidth',2);
plot(1:meas.K, CardMean(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2);
plot(1:meas.K, CardMean(iFilter,:)+CardStdDev(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2, 'LineStyle', '--');
plot(1:meas.K, CardMean(iFilter,:)-CardStdDev(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2, 'LineStyle', '--');
set(gca,'FontSize',axisFontSize);
ylim(limits);
% xlim([0 8]);
xlabel(' ');
ylabel('Cardinality','FontSize',labelFontSize);
title('PHD filter','FontSize',titleFontSize);
hl = legend(gca,'True','Mean','Mean with standard deviation','Location','Southeast'); iFilter = iFilter+1;
set(hl,'FontSize',legendFontSize);
set(gca, 'LooseInset', [0,0,0,0],'OuterPosition',[0 2/3 1 1/3]);

subplot(3,1,2); box on; hold on; grid on;
stairs(1:meas.K,truth.N,'k','LineWidth',2);
plot(1:meas.K, CardMean(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2);
plot(1:meas.K, CardMean(iFilter,:)+CardStdDev(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2, 'LineStyle', '--');
plot(1:meas.K, CardMean(iFilter,:)-CardStdDev(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2, 'LineStyle', '--');
set(gca,'FontSize',axisFontSize);
ylim(limits);
% xlim([0 8]);
xlabel(' ');
ylabel('Cardinality','FontSize',labelFontSize);
title('CPHD filter','FontSize',titleFontSize);
hl = legend(gca,'True','Mean','Mean with standard deviation','Location','Southeast'); iFilter = iFilter+1;
set(hl,'FontSize',legendFontSize);
set(gca, 'LooseInset', [0,0,0,0],'OuterPosition',[0 1/3 1 1/3]);

subplot(3,1,3); box on; hold on; grid on;
stairs(1:meas.K,truth.N,'k','LineWidth',2);
plot(1:meas.K, CardMean(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2);
plot(1:meas.K, CardMean(iFilter,:)+CardStdDev(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2, 'LineStyle', '--');
plot(1:meas.K, CardMean(iFilter,:)-CardStdDev(iFilter,:), 'Color', colors(iFilter,:), 'LineWidth', 2, 'LineStyle', '--');
set(gca,'FontSize',axisFontSize);
ylim(limits);
% xlim([0 8]);
ylabel('Cardinality','FontSize',labelFontSize);
xlabel('Time [s]','FontSize',labelFontSize);
title('Discrete-Gamma CPHD filter','FontSize',titleFontSize);
hl = legend(gca,'True','Mean','Mean with standard deviation','Location','Southeast');
set(hl,'FontSize',legendFontSize);
set(gca, 'LooseInset', [0,0,0,0],'OuterPosition',[0 0 1 1/3]);

% set(gca,'LooseInset',get(gca,'TightInset'));

if saveFigures
    export_fig([FIGURES_PATH, 'Mean_cardinality_versus_time' '.pdf'], '-pdf', '-native', '-nocrop');
end

% Compute MOSPA
MOSPA = zeros(3,truth.K);
%plot error
ospa_vals= zeros(truth.K,3);
ospa_c= 100;
ospa_p= 1;
for k=1:meas.K
    [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_phd.X{k},[1 3]),ospa_c,ospa_p);
end
MOSPA(1,:) = ospa_vals(:,1)';
for k=1:meas.K
    [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_cphd.X{k},[1 3]),ospa_c,ospa_p);
end
MOSPA(2,:) = ospa_vals(:,1)';
for k=1:meas.K
    [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est_dgcphd.X{k},[1 3]),ospa_c,ospa_p);
end
MOSPA(3,:) = ospa_vals(:,1)';

% Plot MOSPA x time
legendFontSize = 40;
labelFontSize = 38;
axisFontSize = 26;
markers = {'+','o','x','*','square','diamond','pentagram','hexagram'};
figure('units','normalized','outerposition',[0 0 1 1],'Color','w');
hnd = zeros(numFilters,1);
lbl = cell(numFilters,1);
for iFilter = 1:numFilters
    lbl{iFilter,1} = filterTypeString{iFilter};
    hnd(iFilter) = plot(1:meas.K,MOSPA(iFilter,:),'Color',colors(iFilter,:),'Marker',markers{iFilter}, 'LineWidth', 3, 'MarkerSize', 12);
    hold on;
end
h_l = legend(hnd,lbl,'Location', 'Northeast', 'Orientation', 'Vertical');
set(h_l,'FontSize',legendFontSize);
set(gca,'FontSize',axisFontSize);
% xlim([0 8]);
ylabel('Mean OSPA metric','FontSize',labelFontSize);
xlabel('Time [s]','FontSize',labelFontSize);

grid on;

set(gcf, 'PaperPositionMode', 'auto', 'Color', 'w');
set(gca,'LooseInset',get(gca,'TightInset'));

if saveFigures
    export_fig([FIGURES_PATH, 'MOSPA_versus_time' '.pdf'], '-pdf', '-native', '-nocrop');
end

if saveFigures
    close all;
end

function [X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks)

K= size(X,1);
x_dim= size(X{K},1);
k=K-1; while x_dim==0, x_dim= size(X{k},1); k= k-1; end
X_track= zeros(x_dim,K,total_tracks);
k_birth= zeros(total_tracks,1);
k_death= zeros(total_tracks,1);

max_idx= 0;
for k=1:K
    if ~isempty(X{k})
        X_track(:,k,track_list{k})= X{k};
    end
    if max(track_list{k})> max_idx %new target born?
        idx= find(track_list{k}> max_idx);
        k_birth(track_list{k}(idx))= k;
    end
    if ~isempty(track_list{k}), max_idx= max(track_list{k}); end
    k_death(track_list{k})= k;
end

function Xc= get_comps(X,c)

if isempty(X)
    Xc= [];
else
    Xc= X(c,:);
end
