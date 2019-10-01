function y = fconv(x,Fs,convWidth)
% convWidth = 50;
% Fs = 1000;
a = 1/convWidth;

convWindow = zeros(Fs, 1);
for i = 1:length(convWindow)
    convWindow(i) = (a^2)*i*exp(-a*i);
%     Window(i) = sqrt(1/a/pi)*exp(-(1/a)*((i-3*(1/a))/Fs)^2);
%     Window(i) = 1/2*(1-cos(2*pi*i/(Fs-1)))*exp(-a*abs(Fs-1-2*i)/(Fs-1));
end
convWindow = convWindow/max(convWindow);

% y = conv2(x,convWindow); %lucas code
y = filter(convWindow,1,x); %sijia code
y = y*10;

%% Plot an example of convolution's result
% figure(1000); clf;
% subplot(2,1,1);
% plot(convWindow,'k');
% xlabel('[ms]');
% title(['convolution filter with convWidth=', num2str(convWidth)]);
%
% subplot(2,1,2);
% hold on;
% plot(x,'--k');
% plot(y,'r');
% hold off;
%
% xlabel('time from onset [ms]');
% legend({'original';'convoluted'});
% title(['convolution applied before & after']);
%
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
% filename = ['eg_conv'];
% saveas(gcf,['./Plots/' filename, '.png'],'png');
