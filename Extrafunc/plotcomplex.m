function plotcomplex(Data,MM)

% Data is a complex vector
%%
if ~exist('MM','var') || isempty(MM)
    M1 = max(abs(real(Data)));
    M2 = max(abs(imag(Data)));
    MM = max(M1,M2);
end

hold on;
line([-MM MM],[0 0],'linestyle','--','color',[.7 .7 .7],'linewidth',1.5);
line([0 0],[-MM MM],'linestyle','--','color',[.7 .7 .7],'linewidth',1.5);


for k = 1:numel(Data)
    line([0 real(Data(k))],[0 imag(Data(k))],'Color',[.2 .2 .7],'linewidth',1.5);
end


xlim([-MM MM]);
ylim([-MM MM]);

end