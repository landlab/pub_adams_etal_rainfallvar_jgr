function [cov,r]=calc_cov(up,ar,slp)
% [cov,r]=calc_cov(up,ar,slp)
% given row matrices containing data on uplift (up), area (ar), and 
% slope (slp) (they should be the same length) this function will 
% calculate and return the 
% covariance matrix (cov) and correlation matrix (r)
%NIC YOU SHOULD BE FEEDING LOG VALUES OF UP, AR, AND SLP

data=[up',ar',slp'];

%Calc covariance Matrix
len=size(data,1);

j=ones(len,1);
%row matrix of the mean of the variables
m=1/len*j'*data;
%matrix that is same size as data but contains means
M=j*m;
%matrix of deviations
D=data-M;
%covariance matrix (square of deviations divided by N-1
cov=1/(len-1)*D'*D;

%quick and dirty method to get covariance matrix
%cov=1/(len-1)*data'*(eye(len)-1/len*j*j')*data;

%now calculate correlation matrix
%extract the variances from the diagnol of the covariance matrix
vars=diag(cov); %column matrix
varm=diag(vars); %diagonal matrix
s=sqrt(varm);
w=inv(s);
r=w*cov*w; %correlation matrix