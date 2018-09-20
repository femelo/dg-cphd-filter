% The original code is freely available at http://ba-tuong.vo-au.com/codes.html

function truth = gen_truth(model)

%variables
truth.K= model.num_time_steps;        %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
% nbirths= 32;
bounds = [800; 5; 800; 5];
i = 1;

N = model.num_targets;
N0 = round(N/4);
N1 = round(N/2)-N0;
N2 = round(3*N/4)-(N0+N1);
N3 = N-(N0+N1+N2);

% 0
for j = 1:N0
    xstart(:,i)  = -bounds +2*bounds.*rand(4,1);
    % xstart(:,i)  = bounds.*randn(4,1);
    tbirth(i)  = 1;
    if j <= 5
        tdeath(i) = 80;
    else
        tdeath(i)  = truth.K+1;
    end
    i = i+1;
end

% 1
for j = 1:N1
    xstart(:,i)  = -bounds +2*bounds.*rand(4,1);
    tbirth(i)  = 20;
    tdeath(i)  = truth.K+1;
    i = i+1;
end

xstart(:,i)  = -bounds +2*bounds.*rand(4,1); tbirth(i)  = 20; tdeath(i) = truth.K+1; i = i+1;
xstart(:,i)  = -bounds +2*bounds.*rand(4,1); tbirth(i)  = 20; tdeath(i) = truth.K+1; i = i+1;

% 2
for j = 1:N2
    xstart(:,i)  = -bounds +2*bounds.*rand(4,1);
    tbirth(i)  = 40;
    tdeath(i)  = truth.K+1;
    i = i+1;
end

xstart(:,i)  = -bounds +2*bounds.*rand(4,1); tbirth(i)  = 40; tdeath(i) = truth.K+1; i = i+1;

% 3
for j = 1:N3
    xstart(:,i)  = -bounds +2*bounds.*rand(4,1);
    tbirth(i)  = 60;
    tdeath(i)  = truth.K+1;
    i = i+1;
end

xstart(:,i)  = -bounds +2*bounds.*rand(4,1); tbirth(i)  = 60; tdeath(i) = truth.K+1; i = i+1;
xstart(:,i)  = -bounds +2*bounds.*rand(4,1); tbirth(i)  = 60; tdeath(i) = truth.K+1; i = i+1;

nbirths = i-1;

%generate the tracks
for targetnum=1:nbirths
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        targetstate = gen_newstate_fn(model,targetstate,'noiseless');
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
     end
end
truth.total_tracks= nbirths;
