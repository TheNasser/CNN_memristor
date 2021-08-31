
% plots the desired cell's state and output over the simultaion.
function PlotOutput(t_Vec ,V_xMatHist ,V_yMatHist,N,r,c,t)
    Plotx = zeros(length(t_Vec));
    Ploty = zeros(length(t_Vec));
    for j = 1:N
        Plotx(j) = V_xMatHist(r,c,j);
        Ploty(j) = V_yMatHist(r,c,j);
    end
    figure();
    plot(t_Vec,Plotx(:,1));
    hold on;
    plot(t_Vec,Ploty(:,1));
    grid on;
    xlim([0 t ])
    hold off;
    legend('$V_{x}$','$V_{y}$','Interpreter','latex');
    xlabel('time (sec)')
    Ystr = ['State of cell circuit C(',num2str(r),',',num2str(c),')'];
    ylabel(Ystr);
    Tstr = ['Transient waveform of cell circuit C(',num2str(r),',',num2str(c),')'];
    title(Tstr);

end