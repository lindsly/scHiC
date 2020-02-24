function erez_imagesc()
    erez_colors = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
    colormap(erez_colors)
    axis square
end
