% magnitude correction
% noise level simulation
% parameters setting
ADC = 1e-3;
Signal = flip([5,10,50,100,500]);
N = 10000;
b = [50,1000,500];
FA = [1,cos(pi/6),cos(pi/6).^2];

% simulation process
ADC_m = zeros(length(Signal),N);
ADC_c = zeros(length(Signal),N);
ADC_psnr = ADC_m;
ADC_psnr2 = ADC_c;
for i = 1:length(Signal)
    for n = 1:N
        noise = .7071*(randn(1,3)+1i*randn(1,3));
        s = abs(Signal(i).*FA.*exp(-ADC.*b) + noise)./FA;
        x = [1,0];
        [x,~,~,rms] = expfitw(b,s,x,FA);
        ADC_m(i,n) = x(2);
        ADC_psnr(i,n) = rms/s(1);
        x = [1,0];
        [x,~,~,rms] = expfitw(b,sqrt(abs(s.^2-0.5./FA.^2)),x,FA);% FA correction effect
        ADC_c(i,n) = x(2);
        ADC_psnr2(i,n) = rms/s(1);
    end
end


