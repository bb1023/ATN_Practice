function[result]=compute_threshold(Dead,alpha)
    [m,n]=size(Dead);
    result=alpha(n)*ones(1,n);
    for j=1:n
        for i=1:m
            if Dead(i,j)>1e-3
                if j>1
                    result(j-1)=alpha(i);
                else
                    result(1)=-1;
                end
                break
            end
        end
    end
end