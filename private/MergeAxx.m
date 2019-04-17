function out_axx = MergeAxx(axxlist)
out_axx = axxlist{1};
    for C = 2:numel(axxlist)
        out_axx = out_axx.MergeTrials(axxlist{C});
    end
end