function res = contains(rootstring,querystring)
    res = ~isempty(strfind(rootstring,querystring));
end