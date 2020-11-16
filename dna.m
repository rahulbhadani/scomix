 x = -3:0.001:3;
y = 1./( 1 + exp(-x));
z = -y + 1;

x1 = x(1:500:end);
y1 = y(1:500:end);
z1 = z(1:500:end);
m = 10;
h = figure;
%plot(x1, y1, 'o');
hold on;

for i = 1:length(x1)
     if mod(i,2) == 1
            h2  = plot(repmat(x1(i), [1, length(linspace(y1(i), z1(i), 6))]),linspace(y1(i), z1(i), 6), 'o', 'MarkerSize', m);
            set(h2, 'markerfacecolor', get(h2, 'color')); % Use same color to fill in markers
     end
end


for i = 1:length(x1)
    h1 =  plot(x1(i), y1(i),  'o', 'MarkerSize', m);
    set(h1, 'markerfacecolor', get(h1, 'color')); % Use same color to fill in markers

    m = m+2;
end


for i = 1:length(x1)
    h1 =  plot(x1(i), z1(i),  'o', 'MarkerSize', m);
    set(h1, 'markerfacecolor', get(h1, 'color')); % Use same color to fill in markers
    m = m-2;
end




set(gca, 'CameraUpVector', [sin(30), cos(10), 0]);
%text(2,0.5,'scOMiX', 'FontSize',172)
set(gca,'visible','off')

 set(gcf, 'color', 'none');
 set(gca, 'color', 'none');
 set(0,'DefaultAxesColor','none');
 set(h, 'InvertHardCopy', 'off');