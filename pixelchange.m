
function [MSD,s] = pixelchange(IM02,dt)

n_tau = round((length(IM02)/2));

IM03 = zeros([size(IM02{1,1}),size(IM02,2)]);
MSD = zeros([size(IM02{1,1}),n_tau]);

for m = 1:size(IM02,2)
    IM03(:,:,m) = IM02{1,m};
end

for n = 0:1:n_tau
    MSD(:,:,n+1) = mean((IM03(:,:,n+1:end) - IM03(:,:,1:end-n)).^2, 3);
end

s = dt*[0:1:n_tau]; 

end

