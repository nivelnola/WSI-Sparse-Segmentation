function [Patches] = data2patches(data, sizePatches)
%data2patches Converts 2D data matrix to a 4D Patches matrix.

Patches = reshape(data, [sizePatches, sizePatches, 3, size(data,2)]);

        
end