function [p f_crit t] = tcirc_2(z1,z2,alpha)
% T^2_circ statistic based on Victor & Mast (1991)

% check inputs
if nargin<3
    alpha = 0.95;
end
if alpha<0 || alpha>1 
    error('alpha must be between 0 and 1');
end
if nargin<2
    z2 = 0;
end
if nargin<1
    error('You must provide some data');
end
if isempty(z1)
    z1 = 0;
end
if ~(isnumeric(z1) && isnumeric(z2))
    error('Inputs must be vectors');
end
if numel(z1)>length(z1) || numel(z2)>length(z2)
    error('Inputs must be vectors');
end

% determine number of samples
M1 = length(z1);
M2 = length(z2);

% compute means
m_z1 = mean(z1);
m_z2 = mean(z2);

% compute variance
v_indiv = 1/(2*(M1 + M2 - 2)) * (sum(abs(z1-m_z1).^2)+sum(abs(z2-m_z2).^2));
v_group = M1*M2/(2*(M1 + M2)) * abs(m_z1-m_z2)^2;

% compute the tcirc-statistic
t = (M1 + M2)/(M1*M2) * v_group/v_indiv;
%t = (M1+M2-2) * (abs(m_z1-m_z2))^2/sum(abs(z1-m_z1)^2 + sum(abs(z2-m_z2))^2);

% perform an F-test
f_crit = finv(1-alpha,2,2*M1+2*M2-4);
p = fpdf(t,2,2*M1+2*M2-4);