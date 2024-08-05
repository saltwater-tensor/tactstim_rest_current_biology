function B = muliplyInterp_local(W, A, nComponents)
    switch (nComponents)
        case {0, 1}
            B = double(W * A);
        case 2
            B = zeros(2 * size(W,1), size(A,2));
            B(1:2:end,:) = double(W * A(1:2:end,:));
            B(2:2:end,:) = double(W * A(2:2:end,:));
        case 3
            B = zeros(3 * size(W,1), size(A,2));
            B(1:3:end,:) = double(W * A(1:3:end,:));
            B(2:3:end,:) = double(W * A(2:3:end,:));
            B(3:3:end,:) = double(W * A(3:3:end,:));
    end
    B = double(B);
end