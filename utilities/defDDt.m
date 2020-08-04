    function [D,Dt] = defDDt(w,a)
        % defines finite difference operator D^a
        % and its transpose operator Dt=conj(-1^a)div^a
        
        D = @(U) ForwardD(U,w);
        Dt = @(X,Y) Dive(X,Y,w,a);
        
        function [Dux,Duy] = ForwardD(U,w)
            [m,n] = size(U);
            % backward finite difference operator
            Dux = imfilter(U,w,'circular','full'); % column differences
            Dux = Dux(:,1:n);
            Duy = imfilter(U,w','circular','full'); % row differences
            Duy = Duy(1:m,:);
        end
        
        function DtXY = Dive(X,Y,w,a)
            % Transpose of the backward finite difference operator
            d = length(w);
            w = fliplr(w);
            DtX = imfilter(X,w,'circular','full');
            DtX = DtX(:,d:end);
            DtY = imfilter(Y,w','circular','full');
            DtY = DtY(d:end,:);
            DtXY = conj((-1)^a)*(-1)^a*(DtX+DtY); %%add conj((-1)^a)*(-1)^a
        end
    end