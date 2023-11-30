jdir = dir('*.mat');
for i = 1:length(jdir)
    data = load(jdir(i).name);
    figure(10);clf; p = imagesc(log10(data.collapse_time_square));
    set(p, 'AlphaData', ~isnan(data.collapse_time_square))
    clim([0,4]); colormap(gca, cmocean('matter', 100)); colorbar;
    title(jdir(i).name)

    %add inf data in another coloru
    ax = gca;
    axnew = axes;
    axnew.Position = ax.Position;
    isinfdata = data.collapse_time_square;
    isinfdata(isinf(data.collapse_time_square)) = 1;
    isinfdata(~isinf(data.collapse_time_square)) = nan;
    p = imagesc(isinfdata);
    set(p, 'AlphaData', ~isnan(isinfdata))
    axnew.Visible = 'off';

    colormap(axnew, [0,1,0]);
    drawnow; shg
    pause
end
