function visualize_dict_layers(Dictionary)
%visualize_dict_layers Displays the dictionary matrix to visualize the
%individual channels

imagesc(Dictionary)
colormap('jet')
title('Unclustered Dictionary')
ylabel('Features')
xlabel('Atoms')
colorbar

end

