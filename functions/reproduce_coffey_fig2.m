Tbk = -2 + 273.15;
colmap = [0,0,0; 1,0,0; 0,0,1];
count = 1;
figure(4); clf; hold on
for Ts = [-2,-17, -32]
Tsk = Ts + 273.15; 
anonT_Niall = @(z) (Tsk + (Tbk - Tsk)*z); %Niall uses z = 0 at surface
anonT = @(z) anonT_Niall(1-z); %mine, with z =0 at the base



[Rxx, ds, db] = spacesearch_ModifiedNye(anonT);

figure(4);hold on;
plot(Rxx,db, 'Color',colmap(count,:), 'linewidth', 1.5);
hold on
plot(Rxx,1-ds, '--','color',  colmap(count,:), 'linewidth', 1.5);
shg
xlim([0,2])
count = count+1;
end