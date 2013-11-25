function path = fixpath(path)
    % FIXPATH
    %   fixpath(path) ensures that the path ends with a '/'
    
    % ensure the pathectory ends with a "\"
    if ispc
      if ~strcmp(path(end),'\')
          path = [path '\'];
      end
    else
      if ~strcmp(path(end),'/')
          path = [path '/'];
      end
    end
end
