function D = distance_points(P)
% distanza_punti calcola la matrice delle distanze euclidee
% Input:
%   P -> matrice Nx3 (N punti nello spazio 3D)
% Output:
%   D -> matrice NxN delle distanze

    N = size(P,1);
    D = zeros(N,N);

    for i = 1:N
        for j = i+1:N
            d = norm(P(i,:) - P(j,:));
            D(i,j) = d;
            D(j,i) = d; % simmetria
        end
    end

end