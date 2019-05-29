function out_axx = MergeAxx(axxlist)
out_axx = axxlist{1};
if ~isempty(out_axx)
    for C = 2:numel(axxlist)
        if ~isempty(axxlist{C}) 
            out_axx = out_axx.MergeTrials(axxlist{C});
        end
    end
else
    if numel(axxlist)>1
        out_axx = MergeAxx(axxlist(2:end));
    else
        out_axx = axxlist;
    end
end
end