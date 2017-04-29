function [no_output]=matrix_print(tmp)

if size(tmp,2)>1
    for i=1:size(tmp,1);
        for j=1:size(tmp,2);
            fprintf('%12.6f',tmp(i,j));
        end
        fprintf('\n');
    end
end

if size(tmp,2)==1
    fprintf('%12.4f\n',tmp);
end