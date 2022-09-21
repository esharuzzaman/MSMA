function [ alpha ] = alpha_i_r( Happ, i )

% Alpha evolves linearly

j = 2*(i~=2)+1*(i==2);
k = 6-i-j;

[mu0, Msat, rhok1, D, Halpha] = Mat_consts();


if abs(Happ(i))<Halpha
        alpha = ((Happ(i))/(2*Halpha))+0.5;
else
    alpha = 1;
end

end

