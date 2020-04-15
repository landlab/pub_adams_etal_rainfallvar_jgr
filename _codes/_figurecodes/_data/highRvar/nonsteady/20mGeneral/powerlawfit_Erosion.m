function [b,bint,gridlog10A,gridlog10S,gridlog10I]=powerlawfit(aall,sall,uall,fig)


%[b,bint,cov,r]=powerlawfit(aall,sall,uall,fig,thresh)
%... Input data goes here

%[aall,sall,uall]=getslpaup(basenm,ts);

% ct=0;
% 
% for i=1:length(aall)
%     if aall(i)>thresh
%         ct=ct+1;
%         A(ct)=aall(i);
%         I(ct)=uall(i);
%         S(ct)=sall(i);
%     end
% end

Tr = csvread('/Users/Jordan/Desktop/LA_Results_a1/20mGeneral/durations.txt');
P = csvread('/Users/Jordan/Desktop/LA_Results_a1/20mGeneral/intensities.txt');
%d = csvread('/Users/joad3819/Desktop/ID_Results_a1/20mGeneral/depths.txt');
%disp(storm_array);

%Tr = storm_array(1:50,1);
%disp(Tr);
%P = storm_array(:,2);
%disp(mean(Tr));
%disp(mean(P));
%disp(mean(d));
aall = reshape(Tr.', 1, numel(Tr));
aall = aall(1:401);

sall = reshape(P, 1, numel(P));
sall = sall(1:401);
%disp(mean(sall))
Er = csvread('/Users/Jordan/Desktop/LA_Results_a1/20mGeneral/midstream1_median_erosion_rates_tc0.txt');
%Er = csvread('/Users/joad3819/Desktop/outlet_median_erosion_rate.txt')
Er = Er(:,1);
%Er(find(Er < 0)) = 0.;


uall = reshape(Er.', 1, numel(Er));
uall = uall(1:401);
%disp(uall);
ids = find(uall > 0);
%disp(Er);

%disp(Er);
A=aall(ids);
S=sall(ids);
I=uall(ids);
fig=1;

log10A=log10(A);
log10S=log10(S);
log10I=log10(I);

%... Calculate covariance and correlation
[cov,r]=calc_cov(log10A,log10S,log10I);

disp('=====CORRELATION MATRIX======');
str=sprintf('\tlog10A\tlog10S\tlog10I');
disp(str);
str = sprintf('log10A\t%g\t%g\t%g',r(1,1),r(1,2),r(1,3));
disp(str);
str = sprintf('log10S\t%g\t%g\t%g',r(2,1),r(2,2),r(2,3));
disp(str);
str = sprintf('log10I\t%g\t%g\t%g',r(3,1),r(3,2),r(3,3));
disp(str);

[b,bint,r,rint,stats]=regress(log10I',[ones(size(log10I,2),1),log10A',log10S']);

disp(' ');
disp('=====REGRESSION RESULTS=====');
str=sprintf('K* = %g  (95%% conf interval: %g to %g)', 10^b(1), (10).^bint(1,:));
disp(str)
str=sprintf('M* = %g  (95%% conf interval: %g to %g)', b(2), bint(2,:));
disp(str);
str=sprintf('N* = %g  (95%% conf interval: %g to %g)', b(3), bint(3,:));
disp(str);
str=sprintf('R^2 = %g',stats(1));
disp(str);

%... Organize data for 3D plot
minlog10A=floor(min(log10A));
maxlog10A=ceil(max(log10A));
steplog10A=(maxlog10A-minlog10A)/100;
minlog10S=floor(min(log10S));
maxlog10S=ceil(max(log10S));
steplog10S=(maxlog10S-minlog10S)/100;

[gridlog10A,gridlog10S] = meshgrid(minlog10A:steplog10A:maxlog10A,minlog10S:steplog10S:maxlog10S);
gridlog10I = b(1)+ b(2)*gridlog10A + b(3)*gridlog10S;

%... Slope versus area plot
figure(fig)
%subplot(2,2,1)
%... Organize incision rate contours for area vs slope plot
%[C,h] = contour((10).^gridlog10A, (10).^gridlog10S, (10).^gridlog10I);
set(gca,'FontSize',18)
[C,h] = contour(gridlog10A, gridlog10S, gridlog10I);
%use below to make in non-log space, then manually changes axes to log to 
%get normal slope-area plot
%contour((10).^glA,(10).^glS,glI)
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
colormap cool
hold on
plot(log10(A),log10(S),'x');
ylabel('log(precipitation intensity (mm/hr))')
xlabel('log(storm duration (hr)')
title('P vs. Tr with Contours for peak discharge');
%hold off

%... Incision rate versus area
%     subplot(2,2,2)
%     loglog(A,I,'x');
%     hold on
%     xlabel('Log10(Upstream Area km^2)')
%     ylabel('Log10(Incision Rate mm/yr)')
%     title('Incision Rate vs. Upstream Area')
%     hold off

%... Incision rate versus slope
% subplot(2,2,3)
% plot(log10S,log10I,'x');
% hold on
% xlabel('Log10(Channel Gradient)')
% ylabel('Log10(Incision Rate mm/yr)')
% title('Incision Rate vs. Channel Gradient')
% hold off

% ... 3D plot
figure(fig+1)
%subplot(2,2,4)
figure(fig+1)
%plot3(log10A, log10S, log10I,'x');
plot3(sall,aall,uall,'x');
hold on
plot3((10).^gridlog10S, (10).^gridlog10A,  (10).^gridlog10I,'-r');
%plot3(gridlog10A,gridlog10S, gridlog10I,'-r');
grid on
xlabel('Log10(Storm duration (hr))');
ylabel('Log10(Precipitation (mm/hr))');
zlabel('Log10(peak discharge (m^3/s)');
title('Plane is Best-Fit Powerlaw, with Data Shown by Crosses');
hold off

end

% function [cov,r]=calc_cov(x,y,z)
% % given x, y, and z as row matrices of the same size, this function will 
% % calculate and return the covariance matrix (cov) and correlation matrix (r)
% 
% data=[x',y',z'];
% 
% %Calc covariance Matrix
% len=size(data,1);
% 
% j=ones(len,1);
% %row matrix of the mean of the variables
% %m=1/len*j'*data;
% %matrix that is same size as data but contains means
% %M=j*m;
% %matrix of deviations
% %D=data-M;
% %covariance matrix (square of deviations divided by N-1
% %cov=1/(len-1)*D'*D;
% 
% %quick and dirty method to get covariance matrix
% cov=1/(len-1)*data'*(eye(len)-1/len*j*j')*data;
% 
% %now calculate correlation matrix
% %extract the variances from the diagnol of the covariance matrix
% vars=diag(cov); %column matrix
% varm=diag(vars); %diagonal matrix
% s=sqrt(varm);
% w=inv(s);
% r=w*cov*w; %correlation matrix
% end