% minimize the transition area
% bwd1 map
% thk 15mm, peak 2G, mu: 30~70

bwd1 = 9e3;
bwd2 = 3e3;

mu1 = linspace(6,80,60);
mu2 = linspace(6,80,60);

wid1 = zeros(size(mu1));
wid2 = zeros(size(mu2));
for i = 1:length(wid1)
    wid1(i) = trans_est(1e-2,1.8,bwd1,mu1(i),1.5,.90);
    wid2(i) = trans_est(1e-2,1.2,bwd2,mu2(i),1.5,.90);
end

%% comparison
figure;
subplot(2,2,2);
trans_1 = trans_est(1e-2,1.8,9e3,30,1.5,.9);
subplot(2,2,4);
trans_2 = trans_est(1e-2,1.2,3e3,7,1.5,.9);
subplot(2,2,1);
pulse_1 = pulse_HS ( pi, 30, 9e3, 1e-2, 1e-5);
plot(linspace(0,1e-2,length(pulse_1)),[abs(pulse_1);angle(pulse_1)]);
xlabel('duration/s');
legend('Magnitude','Phase');
subplot(2,2,3);
pulse_2 = pulse_HS ( pi, 7, 3e3, 1e-2, 1e-5);
plot(linspace(0,1e-2,length(pulse_2)),[abs(pulse_2);angle(pulse_2)]);
xlabel('duration/s');
legend('Magnitude','Phase');