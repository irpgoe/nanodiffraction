function opts = update_defaults(defaults,opts)

    fields = fieldnames(defaults);
    for f = 1:numel(fields)
        if ~isfield(opts,fields{f})
            opts.(fields{f}) = defaults.(fields{f});
        end
    end

end