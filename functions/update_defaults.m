function defaults = update_defaults(defaults,opts)

%     fields = fieldnames(defaults);
%     for f = 1:numel(fields)
%         if ~isfield(opts,fields{f})
%             opts.(fields{f}) = defaults.(fields{f});
%         end
%     end
    
    fields = fieldnames(opts);
    for f = 1:numel(fields)
        if ~isfield(defaults,fields{f})
            warning(sprintf('You provided an illegal option: %s',fields{f}));
        else
            defaults.(fields{f}) = opts.(fields{f});
        end
    end

end