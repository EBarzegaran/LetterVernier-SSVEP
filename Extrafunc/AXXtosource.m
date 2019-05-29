function sourceAxx = AXXtosource(Axx,InvM)
% InvM is a MxS matrix, where M is the number of electrodes and S is
% the number of sources

    sourceAxx = Axx;
    
    for t = 1:size(Axx.Cos,3)
        Cos(:,:,t) = Axx.Cos(:,:,t)*InvM;
        Sin(:,:,t) = Axx.Sin(:,:,t)*InvM;
        Wave(:,:,t) = Axx.Wave(:,:,t)*InvM;
    end
    
    sourceAxx.Cos = Cos;
    sourceAxx.Sin = Sin;
    sourceAxx.Wave = Wave;
    sourceAxx.Amp = abs(Cos+(1i*Sin));
end