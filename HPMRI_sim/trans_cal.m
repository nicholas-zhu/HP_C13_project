% transaction area testing
% based on trans_calc.m file
% - generate ratio map for this
% variables :
%   area    : area ratio matrix for mu, f0, A0
%   threshold: 0.9
%   A0      : 0.7~2.1 (G)
%   f0      : 5e3~10e4 (kHz)
%   mu      : 10~100;

close all;
clear;
clc;
addpath(genpath('functions'));

gamma_ba = 1.071e3;
f0 = (3:.5:10)*1e3;
mu = 10:2:60;
p = 3;

threshold = 0.95;
ratio = zeros(length(f0),length(mu));
for j = 1 : length(f0)
    for k = 1 : length(mu)
       A0 = min(p*f0(j)/(2*gamma_ba*sqrt(mu(k))),2.5);
       drawnow;
       ratio(j,k) = trans_calc_inv(A0,f0(j),mu(k),threshold);
    end
end
%%
figure;
fig1 = (1-ratio);
[y,x] = meshgrid(mu,f0/1000);
surf(x,y,fig1);colormap hot;
axis([min(x(:)) max(x(:)) min(y(:)) max(y(:)) 0.3 1])
xlabel('Bandwidth/kHz');
ylabel('adiabatic factor');
zlabel('1 - transition area ratio')
%%
for i = 1 : length(f0)
    temp = cat(1,ratio(i,:),1:size(ratio,2))';
    temp = sortrows(temp,1);
    mu_index(i) = round(mean(temp(2:3,2),1));
    % mu_index(i) = find(ratio(i,:) == min(ratio(i,:)),1);
end

mu_f0(1,:) = mu(mu_index);
mu_f0(2,:) =  min(ratio,[],2)*100;
mu_f0(3,:) = min(p*f0./(2*gamma_ba.*sqrt(mu_f0(1,:))),2.5);

mu_f0(4,:) = f0/1e3;

figure
[AX, H1, H2] = plotyy(mu_f0(4,:),mu_f0(2,:),mu_f0(4,:),mu_f0(3,:));
% axis(AX(1),[3 10 0 20])
% axis(AX(2),[3 10 0 3])
set(H1,'LineStyle','-','Color',[1 0.5 0],'LineWidth',1);
set(H2,'LineStyle','--','Color','b','LineWidth',1);


set(AX(1),'XColor','k','YColor',[1 0.5 0]);
set(AX(2),'XColor','k','YColor','b');

HH1 = get(AX(1),'Ylabel');
HH2 = get(AX(2),'Ylabel');
HH = get(AX(1),'Xlabel');

set(HH1,'String','transition area ratio(%)');
set(HH2,'String','RF peak power(G)');
set(HH,'String','RF bandwidth(kHz)');

