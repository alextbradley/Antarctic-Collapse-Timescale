function ss = get_flowline_data(fname, data, ds)
%generate the flowline data for a given shelf. e.g. fname = 'Ross.mat'
%(must be full filepath). Data is the ice sheet data. ds is the spacing
%along flowlines in metres (e.g. ds = 1e3).
%
%Returns the data in a struct, containing flowline co-ordinates and ice
%thickness, strain rate, melt rate, ice speed, dhdt and surface temp along flowline.

%% Load the shelf data
f = load(fname);

% use the in file to generate a subset of the full data
in = f.IN;
[rId, cId] = find( in ) ; %rows and columns of non-zero entries
xminidx = min(rId);
xmaxidx = max(rId);
yminidx = min(cId);
ymaxidx = max(cId); %indices of min and max in rows and columns

%
% restrict co-ordinates and variables to this subset
xs = data.yy(xminidx:xmaxidx);
ys = data.xx(yminidx:ymaxidx);
hs = data.H(xminidx:xmaxidx, yminidx:ymaxidx);
vs = data.v(xminidx:xmaxidx, yminidx:ymaxidx);
us = data.u(xminidx:xmaxidx, yminidx:ymaxidx);

%quantities needed for the timestepping
speeds   = sqrt(us.^2 + vs.^2);                            %ice velocities
eps_xxs = data.eflow(xminidx:xmaxidx, yminidx:ymaxidx);   %strain rate
%eps_xxs = data.strain(xminidx:xmaxidx, yminidx:ymaxidx);   %strain rate
dhdts   = data.dhdtadj(xminidx:xmaxidx, yminidx:ymaxidx); %dhdt
melts   = data.m(xminidx:xmaxidx, yminidx:ymaxidx);       %melt rate
temps   = data.T(xminidx:xmaxidx, yminidx:ymaxidx);  %ice temperature
%% create interpolations of required data
[xxs, yys] = meshgrid(xs,ys);
yys = yys';
xxs = xxs';

h_si      = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(hs(:)));
speed_si  = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(speeds(:)));
eps_xx_si = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(eps_xxs(:)));
dhdt_si   = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(dhdts(:)));
melt_si   = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(melts(:)));
temp_si   = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(temps(:)));

% make a plot showing the original and reconstruced
if 0
    figure(1); clf;
    subplot(2,4,1); contourf(xs, ys, speeds', 20, 'linestyle', 'none'); title('speed orginal')
    subplot(2,4,2); contourf(xxs, yys, speed_si(yys,xxs), 20, 'linestyle', 'none'); title('speed reconstructed')

    subplot(2,4,3); contourf(xs, ys, eps_xxs', 20, 'linestyle', 'none'); title('strain rate orginal')
    subplot(2,4,4); contourf(xxs, yys, eps_xx_si(yys,xxs), 20, 'linestyle', 'none'); title('strain rate reconstructed')

    subplot(2,4,5); contourf(xs, ys, dhdts', 20, 'linestyle', 'none'); title('dhdt orginal')
    subplot(2,4,6); contourf(xxs, yys, dhdt_si(yys,xxs), 20, 'linestyle', 'none'); title('dhdt reconstructed')

    subplot(2,4,7); contourf(xs, ys, melts', 20, 'linestyle', 'none'); title('melt orginal')
    subplot(2,4,8); contourf(xxs, yys, melt_si(yys,xxs), 20, 'linestyle', 'none'); title('melt reconstructed')

end
%% find edge elements and tag as front or gl based on flow direction
[nxs, nys] = size(hs);
edge_idxs = [];
edge_speed = [];
gl_idxs = [];

for ix = 2:(nxs - 1)
    for iy = 2:(nys-1)
        %take a square of thicknesses around this point
        hsqr = hs(ix-1:ix+1, iy-1:iy+1);
        if any(any(isnan(hsqr))) && (~isnan(speeds(ix,iy))) && (speeds(ix,iy) > 0) %if there is a non-shelf point around and it is a shelf point
            edge_idxs = [edge_idxs; ix, iy];
            edge_speed = [edge_speed; -vs(ix,iy), us(ix,iy) ];

            % determine whether its a gl point by looking
            vv = -vs(ix,iy);
            uu = us(ix,iy); %velocity co-ordinates
            dirx = vv/norm(vv);
            diry = uu/norm(uu);
            if ~isnan(hs(ix+dirx,iy)) && ~isnan(hs(ix,iy+diry)) %if we have thickness data in this direction
                gl_idxs = [gl_idxs; ix, iy];

            end
        end
    end
end

if 0
    figure(2); clf; contourf(xs, ys, speeds', 100, 'linestyle', 'none'); c = colorbar;
    c.Label.String ='speed (m/yr)';
    colormap(parula)
    hold on
    scatter(xs(edge_idxs(:,1)), ys(edge_idxs(:,2)), 10, 'r', 'Filled')
    scatter(xs(gl_idxs(:,1)), ys(gl_idxs(:,2)), 10, 'g', 'Filled')
    quiver(xs(edge_idxs(:,1))', ys(edge_idxs(:,2))', edge_speed(:,1), edge_speed(:,2))
    shg
end

%% make scattered interpolants of velocity
[xxs, yys] = meshgrid(xs,ys);
vvs = vs';
uus = us';
yys = yys';
xxs = xxs';
vsi = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(vs(:)));
usi = scatteredInterpolant(double(yys(:)),double(xxs(:)),double(us(:)));

%% take a point and construct the path from here

maxiters = 1000*1e3/ds;

ss = struct(); %initialize structure
count   =1 ;
for i = 1:length(gl_idxs)
    x0 = xs(gl_idxs(i,1));
    y0 = ys(gl_idxs(i,2));
    flowline = [x0, y0];
    xn = x0;
    yn = y0;
    %plot(xn,yn, 'mo', 'markersize', 6, 'MarkerFaceColor','m', 'MarkerEdgeColor','k')


    inshelf = 1;
    niters = 1;
    dpttol = 5*1e2;
    while inshelf && (niters < maxiters)

        vn = vsi(yn,xn); %v vel at this point
        un = usi(yn,xn);
        vel = [-vn,un];
        direction = vel/norm(vel);
        xn = xn + direction(1)*ds;
        yn = yn + direction(2)*ds;
        %plot(xn,yn, 'ko', 'markersize', 5, 'MarkerFaceColor','k')
        dpt = sqrt((xn - flowline(:,1)).^2 + (yn - flowline(:,2)).^2);
        if niters > 1
            dpt = dpt(2:end); %remove the current point
        else
            dpt = 2*dpttol; %set so it doesn't trigger on the first iteration
        end
        if any(isnan(direction)) || any(dpt < dpttol) %if we're out of the shelf or within a distance dpttol of any other point in the flowline
            inshelf = 0;
        else
            flowline = [flowline; xn yn];
        end

        %drawnow
        niters = niters + 1;
        %pause
    end
    plot(flowline(:,1), flowline(:,2), 'k')
    %  drawnow
    %  pause
    %shg

    %


    % store
    if length(flowline) > 4 %only store if sufficiently many pts
        ss(count).flowline = flowline;
        ss(count).h = h_si(flowline(:,2), flowline(:,1));
        ss(count).eps_xx = eps_xx_si(flowline(:,2), flowline(:,1));
        ss(count).melt = melt_si(flowline(:,2), flowline(:,1));
        ss(count).speed = speed_si(flowline(:,2), flowline(:,1));
        ss(count).dhdt = dhdt_si(flowline(:,2), flowline(:,1));
        ss(count).temp = temp_si(flowline(:,2), flowline(:,1));
        count = count + 1;
    end
    if mod(i,100) == 0
        fprintf('progress on getting flowlines is %.1f percent \n', i/length(gl_idxs)*100);
    end
end %end loop over grounding line pts

axis equal

%% reconstruct quantities from flowlines
if 0
    xr = [];
    yr = [];
    hr = [];
    eps_xxr = [];
    speedr = [];
    dhdtr = [];
    meltr = [];
    for is = 1:length(ss)
        ff = ss(is).flowline;
        xr = [xr; ff(:,1)];
        yr = [yr; ff(:,2)]; %get the x and y co-ords of ALL flowline pts
        hr = [hr; (ss(is).h)];
        eps_xxr = [eps_xxr; (ss(is).eps_xx)];
        meltr = [meltr; ss(is).melt];
        dhdtr = [dhdtr; ss(is).dhdt];
        speedr = [speedr; ss(is).speed];
    end

    % plot comparison
    figure(3); clf;
    %melt
    ax(1) = subplot(2,5,1); contourf(xs, ys, melts', 20, 'linestyle', 'none'); title('melt orginal'); colorbar; %axis equal
    ax(2) = subplot(2,5,2); scatter(xr, yr,10,  meltr, 'filled' );  box on; colorbar; title('melt flowline'); %axis equal

    ax(3) = subplot(2,5,3); contourf(xs, ys, speeds', 20, 'linestyle', 'none');colorbar; title('speed orginal')%axis equal
    ax(4) = subplot(2,5,4); scatter(xr, yr,10,  speedr, 'filled' );  box on; colorbar;title('speed flowline') %axis equal

    ax(5) = subplot(2,5,6); contourf(xs, ys, eps_xxs', 20, 'linestyle', 'none');colorbar;title('strain orginal')% axis equal
    ax(6) = subplot(2,5,7); scatter(xr, yr,10,  eps_xxr, 'filled' );  box on; colorbar; title('strain flowline')%axis equal

    ax(7) = subplot(2,5,8); contourf(xs, ys, hs', 20, 'linestyle', 'none');colorbar;title('thickness orginal') %axis equal
    ax(8) = subplot(2,5,9); scatter(xr, yr,10,  hr, 'filled' );  box on; colorbar;title('thickness flowline')% axis equal

    ax(9) = subplot(2,5,5); contourf(xs, ys, dhdts', 20, 'linestyle', 'none');colorbar;title('dhdt orginal') %axis equal
    ax(10) = subplot(2,5,10); scatter(xr, yr,10,  dhdtr, 'filled' );  box on; colorbar; title('dhdt flowline') %axis equal


    for i = 1:5
        ax(2*i).XLim = ax(2*i - 1).XLim;
        ax(2*i).YLim = ax(2*i - 1).YLim;
        ax(2*i).CLim = ax(2*i - 1).CLim;
        ax(2*i).XTick = [];
        ax(2*i).YTick = [];
        ax(2*i - 1).XTick = [];
        ax(2*i - 1).YTick = [];

    end
end

