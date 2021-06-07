function S_x = chan_size(xmesh,S,x)  % computing S wherever you need it
S_x = interp1(xmesh,S,x)';
end
