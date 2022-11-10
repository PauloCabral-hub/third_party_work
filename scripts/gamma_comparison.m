% Run ltpb_analysis.m before runnning the script

% Filter Data Matrix
subsample = 1; % choose group: 1 for LTPB and 2 for Control
cut = find(OriData(:,12) == subsample);
data = OriData(cut,:);

Aaux = min(data(:,6));
Baux = 1;
for a = 1:size(data,1)
    if data(a,6) == Aaux
       data(a,6) = Baux;
    else
       Aaux = data(a,6);
       Baux = Baux + 1;
       data(a,6) = Baux;
    end
end

% Matlab Parameters for graphics
set(0,'defaultfigurecolor',[1 1 1])
set(0, 'DefaultFigureRenderer', 'painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathtogit = '/home/roberto/Documents/Pos-Doc/BagOfCodesAndData/Git'; % Insert GIT folder location
ntrials = 1000; % Insert the number of trials 
tau = 12; % Insert tree ID
pinpoint = 2; % Supplementary information in points. 1 - success rate, 2 - sample size, 3 - pairing
[ids, vids, tauofid, gametest, trees, treesizes] = data_setaspects(data, ntrials, pathtogit);

t_id = find(tauofid == tau);
szaux1 = length(t_id); szaux2 = treesizes(find(trees == tau));

srate = zeros(length(t_id),1);
tau_repoT = cell(szaux1,szaux2); counts_repoT = zeros(szaux1,szaux2);
tau_repoS = cell(szaux1,szaux2); counts_repoS = zeros(szaux1,szaux2);
tau_repoF = cell(szaux1,szaux2); counts_repoF = zeros(szaux1,szaux2);
counts_total = zeros(szaux1,szaux2);
from = 1; till = 1000;

norm = 0;
limitingl = 0; limitingr = 1.5;
samp_inc = 10; % number of samples necessary for modelling the response times;
control = zeros(szaux1,1);
gamma_repoT = zeros(szaux1,3 ,szaux2); gamma_repoS = zeros(szaux1,3 ,szaux2); gamma_repoF = zeros(szaux1,3 ,szaux2);
steps = [1 2 2 1 3];

for sit = 1:3
for a = 1:length(t_id)
   % Getting participant data
   tree_file_address = ['/home/roberto/Documents/Pos-Doc/BagOfCodesAndData/Git/files_for_reference/tree_behave' num2str(tau) '.txt' ];
   [ctx_rtime, ctx_er, ctx_resp, contexts, ctxrnds, ct_pos] = rtanderperctx(data, t_id(a), from, till, tree_file_address, 0, tau);
   % Getting deterministic success rate
   control(a,1) = sum(ctx_resp{1,1} == 1,1)+ sum(ctx_resp{2,1} == 2,1) + sum(ctx_resp{3,1} == 2,1) + sum(ctx_resp{4,1} == 1,1);
   control(a,1) = control(a,1)/( length(ctx_resp{1,1})+length(ctx_resp{1,1})+length(ctx_resp{3,1})+length(ctx_resp{4,1}) );
   % Getting the success rate
   [chain,responses, ~] = get_seqandresp(data,tau, ids(a), 1, ntrials);
   srate(a,1) = sum(chain == responses)/length(chain);
   % 
   [chain,responses, times] = get_seqandresp(data,tau, ids(a), from, till);
   for b = 1:szaux2
      [ctx_fer,ct_poscell] = lastwas_error(ct_pos, ctx_er, contexts, chain, responses,steps(b)); 
      if sit == 1 % considering all responses
      aux = [1:length(ctx_er{b,1})];
      tau_repoT{a,b} = ctx_rtime{b,1}(aux,1);
      tau_repoT{a,b} = tau_repoT{a,b}(find( (tau_repoT{a,b} <= limitingr)&(tau_repoT{a,b} >= limitingl) ),1);
      counts_repoT(a,b) = length(tau_repoT{a,b});
      [gamma_repoT(a,1,b), gamma_repoT(a,2,b), x, fest, gamma_repoT(a,3,b)] = gamma_model(tau_repoT{a,b}, [0 5], 0.1, 0);
      counts_total(a,b) = length(ctx_rtime{b,1});
      elseif sit == 2  % considering success in the the last responses
      aux = find(ctx_fer{b,1} == 0); 
      tau_repoS{a,b} = ctx_rtime{b,1}(aux,1);
      tau_repoS{a,b} = tau_repoS{a,b}(find( (tau_repoS{a,b} <= limitingr)&(tau_repoS{a,b} >= limitingl) ),1);
      counts_repoS(a,b) = length(tau_repoS{a,b});
          if length(tau_repoS{a,b})>samp_inc
             [gamma_repoS(a,1,b), gamma_repoS(a,2,b), x, fest, gamma_repoS(a,3,b)] = gamma_model(tau_repoS{a,b}, [0 5], 0.1, 0);
          end
      counts_total(a,b) = length(ctx_rtime{b,1});
      else % considering failure in the the last responses
      aux = find(ctx_fer{b,1} == 1); 
      tau_repoF{a,b} = ctx_rtime{b,1}(aux,1);
      tau_repoF{a,b} = tau_repoF{a,b}(find( (tau_repoF{a,b} <= limitingr)&(tau_repoF{a,b} >= limitingl) ),1);
      counts_repoF(a,b) = length(tau_repoF{a,b});
          if length(tau_repoF{a,b} > samp_inc)
             [gamma_repoF(a,1,b), gamma_repoF(a,2,b), x, fest, gamma_repoF(a,3,b)] = gamma_model(tau_repoF{a,b}, [0 5], 0.1, 0); 
          end
      counts_total(a,b) = length(ctx_rtime{b,1});
      end 
   end
end    
    
end

valid = find(isoutlier(control)==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 1; % choose the context

Sgamma= [];
Fgamma = [];
for a = 1:length(t_id)
   if ~ismember(t_id(1,a),find(isoutlier(control)== 1)) 
      if (gamma_repoS(a,1,w) ~= 0)&&(gamma_repoS(a,2,w) ~= 0)
         Sgamma = [Sgamma; gamma_repoS(a,1,w) gamma_repoS(a,2,w)]; 
      end
      if (gamma_repoS(a,1,w) ~= 0)&&(gamma_repoS(a,2,w) ~= 0) % ( gamma_repoF(a,1,w) ~= 0)&&( gamma_repoF(a,2,w) ~= 0)
         Fgamma = [Fgamma; gamma_repoF(a,1,w) gamma_repoF(a,2,w)];    
      end
   end
end

% Representing Gamma as k and theta

Sgamma_alt = zeros(size(Sgamma,1), size(Sgamma,2));
Smedian = zeros(size(Sgamma,1), 1);
for a = 1:size(Sgamma,1)
   Sgamma_alt(a,1) = Sgamma(a,1)*Sgamma(a,2);
   Sgamma_alt(a,2) = Sgamma(a,1)*(Sgamma(a,2)^2);
   Smedian(a,1) = gaminv(0.5,Sgamma(a,1),Sgamma(a,2));
end


Fgamma_alt = zeros(size(Fgamma,1), size(Fgamma,2));
Fmedian = zeros(size(Fgamma,1), 1);
for a = 1:size(Fgamma,1)
   Fgamma_alt(a,1) = Fgamma(a,1)*Fgamma(a,2);
   Fgamma_alt(a,2) = Fgamma(a,1)*(Fgamma(a,2)^2);
   Fmedian(a,1) = gaminv(0.5,Fgamma(a,1),Fgamma(a,2));
end

% Calculating the mean vector (centroid)

centS = zeros(1,2); centS_alt = zeros(1,2);
for a = 1:size(Sgamma,1)
centS = centS+Sgamma(a,:);
centS_alt = centS_alt+Sgamma_alt(a,:);
end
centS = (1/size(Sgamma,1))*centS;
centS_alt = (1/size(Sgamma_alt,1))*centS_alt;

centF = zeros(1,2);centF_alt = zeros(1,2);
for a = 1:size(Fgamma,1)
centF = centF+Fgamma(a,:);
centF_alt = centF_alt+Fgamma_alt(a,:);
end
centF = (1/size(Fgamma,1))*centF;
centF_alt = (1/size(Fgamma_alt,1))*centF_alt;

Tgamma_alt = [Sgamma_alt; Fgamma_alt];

centT_alt = zeros(1,2);
for a = 1:size(Tgamma_alt,1)
centT_alt = centT_alt+Tgamma_alt(a,:);
end
centT_alt = (1/size(Tgamma_alt,1))*centT_alt;

% Distance from different points of the circle

dofcenter = zeros(size(Tgamma_alt,1),1);
for a = 1:size(Tgamma_alt,1)
dofcenter(a,1)= sqrt((centT_alt(1,1) - Tgamma_alt(a,1))^2  + (centT_alt(1,2) - Tgamma_alt(a,2))^2);
end

posit = find(dofcenter == max(dofcenter));

xmin = min(Tgamma_alt(:,1));xmax = max(Tgamma_alt(:,1));
ymin = min(Tgamma_alt(:,2));ymax = max(Tgamma_alt(:,2));

rdiu = dofcenter(posit,1); 

% Plotting the parameters
figure
subplot(1,2,1)
plot(Sgamma(:,1),Sgamma(:,2),'b.','MarkerSize',12)
hold on
plot(Fgamma(:,1),Fgamma(:,2),'r.','MarkerSize',12)

% Mean Vectors (centroids)
plot(centS(1,1),centS(1,2),'bo','MarkerSize',14)
plot(centF(1,1),centF(1,2),'ro','MarkerSize',14)
xlabel('$\hat{k}$','Interpreter','Latex','FontSize',16); ylabel('$\hat{\theta}$','Interpreter','Latex','FontSize',16);
title(['w = ' num2str(contexts{1,w})])
legend('Z^{w,s}','Z^{w,f}','S mean vector', 'F mean vector')
set(gca,'FontName','Times New Roman','FontSize',12, 'FontWeight','bold')

% Distributions with the mean vectors

xS = linspace(0,2,1000);
festS = ( (xS.^(centS(1,1)-1)).*exp(-xS./centS(1,2)) )/( gamma(centS(1,1))*(centS(1,2)^centS(1,1)) );
xF = linspace(0,2,1000);
festF = ( (xF.^(centF(1,1)-1)).*exp(-xF./centF(1,2)) )/( gamma(centF(1,1))*(centF(1,2)^centF(1,1)) );

subplot(1,2,2)
plot(xS,festS,'b-')
hold on
plot(xF,festF,'r-')
legend('S mean vector', 'F mean vector')
xlabel('z (normalized)', 'FontSize',12);
ylabel('$f(z;k,\theta)$','Interpreter','Latex','FontSize',12);
set(gca,'FontName','Times New Roman','FontSize',12, 'FontWeight','bold')
axis square

% Representing gamma as mean and var

% Plotting the parameters
figure
subplot(1,2,1)
plot(Sgamma_alt(:,1),Sgamma_alt(:,2),'b.','MarkerSize',12)
hold on
plot(Fgamma_alt(:,1),Fgamma_alt(:,2),'r.','MarkerSize',12)

% Adding supplementary information to points

for a = 1:size(Sgamma_alt,1)
if pinpoint == 1      % pinpoint success rate
    text(Sgamma_alt(a,1), Sgamma_alt(a,2), [ num2str(valid(a,1)) '(' num2str(srate(valid(a,1)),'%.2f') ')'], 'FontSize', 8);
    text(Fgamma_alt(a,1), Fgamma_alt(a,2), [ num2str(valid(a,1)) '(' num2str(srate(valid(a,1)),'%.2f') ')'], 'FontSize', 8);    
elseif pinpoint == 2  % pinpoint number of samples
    text(Sgamma_alt(a,1), Sgamma_alt(a,2), [ num2str(valid(a,1)) '(' num2str( num2str(size(tau_repoS{valid(a,1),w},1)) ) ')'], 'FontSize', 8);
    text(Fgamma_alt(a,1), Fgamma_alt(a,2), [ num2str(valid(a,1)) '(' num2str( num2str(size(tau_repoF{valid(a,1),w},1)) ) ')'], 'FontSize', 8);    
else                  % pinpoint pairing
    plot([Sgamma_alt(a,1) Fgamma_alt(a,1)], [Sgamma_alt(a,2) Fgamma_alt(a,2)],'k-')    
end

end

% Adding the ID

xlabel('$\hat{\mu}$','Interpreter','Latex','FontSize',14); ylabel('$\hat{\sigma}^{2}$','Interpreter','Latex','FontSize',14);
title(['w = ' num2str(contexts{1,w})])
legend({'Z^{w,s}';'Z^{w,f}'},'Location', 'northwest','FontSize',14)
legend boxoff
set(gca,'FontName','Times New Roman','FontSize',14, 'FontWeight','bold')

% Distributions with the mean vectors

xS = linspace(0,2,1000);
festS = ( (xS.^(centS(1,1)-1)).*exp(-xS./centS(1,2)) )/( gamma(centS(1,1))*(centS(1,2)^centS(1,1)) );
xF = linspace(0,2,1000);
festF = ( (xF.^(centF(1,1)-1)).*exp(-xF./centF(1,2)) )/( gamma(centF(1,1))*(centF(1,2)^centF(1,1)) );

% Calculating the distances

Tgamma_alt = [Sgamma_alt; Fgamma_alt];
reference = [0 , 0];

Sgamma_dist = zeros(size(Sgamma_alt,1),1);
for a = 1:size(Sgamma_alt,1)
    Sgamma_dist(a,1) = sqrt(( Sgamma_alt(a,1) - reference(1,1) )^2+( Sgamma_alt(a,2) - reference(1,2) )^2 );
end

Fgamma_dist = zeros(size(Fgamma_alt,1),1);
for a = 1:size(Fgamma_alt,1)
    Fgamma_dist(a,1) = sqrt((Fgamma_alt(a,1)^2)+(Fgamma_alt(a,2)^2));
end

subplot(1,2,2)

sbox_varsize([1*ones(size(Sgamma_dist,1),1) ; 2*ones(size(Fgamma_dist,1),1)], ...
    [Sgamma_dist; Fgamma_dist],  '', '$D(\hat{\mu},\hat{\sigma}^{2})$', '', {'Z^{w,s}';'Z^{w,f}'}, 0, 0, [])
set(gca,'FontName','Times New Roman','FontSize',14, 'FontWeight','bold')
axis square

[p,h,stats] = ranksum(Sgamma_dist(find(isoutlier(Sgamma_dist)==0),1), Fgamma_dist(find(isoutlier(Fgamma_dist)==0),1)); %#ok<FNDSB>

% Distributions of number of samples

valid = find(isoutlier(control)==0);

boxdata = [counts_repoS(valid,w); counts_repoF(valid,w); ];
boxgroup = [1*ones(length(counts_repoS(valid,w)),1) ; 2*ones(length(counts_repoF(valid,w)),1)];

figure
sbox_varsize(boxgroup, boxdata,  '', 'num. of samples', '', {'Z^{w,s}';'Z^{w,f}'}, 0, 0, [])
title(['w = ' num2str(contexts{1,w})])
set(gca,'FontName','Times New Roman','FontSize',14, 'FontWeight','bold')
axis square