function [data] = patches2data(Patches)
%patches2data Converts 4D Patch matrix to a 2D data matrix.

size_orig = size(Patches);
data = reshape(Patches, [prod(size_orig(1:3)), size_orig(4)]);

end